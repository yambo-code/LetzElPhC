import numpy as np
import re
from scipy.linalg import eigh


np.set_printoptions(suppress=True)

def parser_doubles_from_string(line):
    return list(map(float, re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)))

def string_start_with(line, prefix):
    return line.strip().startswith(prefix)

def read_dyn_qe_old(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    idx = 0
    idx += 2  # skip two comment lines
    read_fbuf = parser_doubles_from_string(lines[idx])
    if len(read_fbuf) != 9:
        raise ValueError("Error reading line 3 in dyn file")

    ntype = int(round(read_fbuf[0]))
    natom = int(round(read_fbuf[1]))
    nmodes = 3 * natom
    ibrav = int(round(read_fbuf[2]))
    idx += 1



    if ibrav == 0:
        idx += 4  # skip 4 lines (lattice vectors)

    atm_mass = np.zeros(natom + ntype)
    atm_mass_type = atm_mass[natom:]

    # read atomic mass for each type
    for i in range(ntype):
        line = lines[idx]
        match = re.match(r"\s*(\d+)\s*'([^']+)'\s*([+-]?[0-9]*[.]?[0-9]+)", line)
        if not match:
            raise ValueError("Failed to read atomic masses from dyn file")
        atm_mass_type[i] = float(match.group(3))
        if abs(atm_mass_type[i]) < 1e-12:
            raise ValueError("Zero mass in dynamical file")
        idx += 1

    # assign masses for each atom
    for i in range(natom):
        read_fbuf = parser_doubles_from_string(lines[idx])
        if len(read_fbuf) != 5:
            raise ValueError("Failed to read atomic mass line")
        itype = int(round(read_fbuf[1]))
        atm_mass[i] = atm_mass_type[itype - 1]
        idx += 1

    nq_found = 0

    qpts = []
    while idx < len(lines):
        if not string_start_with(lines[idx], "Dynamical"):
            idx += 1
            continue

        idx += 2  # skip "Dynamical Matrix..." and empty line
        qpt_tmp = parser_doubles_from_string(lines[idx])
        if len(qpt_tmp) != 3:
            raise ValueError("Error reading qpoint from file")

        qpts.append(qpt_tmp)
        idx += 2  # skip q-point and empty line

        dyn_mat_tmp = np.zeros((nmodes, nmodes), dtype=np.complex128)

        for ia in range(natom):
            for ib in range(natom):
                itmp, jtmp = map(int, lines[idx].split()[:2])
                if itmp != ia + 1 or jtmp != ib + 1:
                    raise ValueError("Mismatch in atom indices in dynmat block")
                idx += 1
                inv_mass_sqrt = 1.0 / np.sqrt(atm_mass[ia] * atm_mass[ib])
                for ix in range(3):
                    read_fbuf = parser_doubles_from_string(lines[idx])
                    if len(read_fbuf) != 6:
                        raise ValueError("Error reading dyn matrix values")
                    read_fbuf = [x * inv_mass_sqrt for x in read_fbuf]
                    for iy in range(3):
                        re_part = read_fbuf[2 * iy]
                        im_part = read_fbuf[2 * iy + 1]
                        row = 3 * ia + ix
                        col = 3 * ib + iy
                        dyn_mat_tmp[row, col] = re_part + 1j * im_part
                    idx += 1

        # Symmetrize
        for i in range(nmodes):
            for j in range(i):
                avg = 0.5 * (dyn_mat_tmp[i, j] + np.conj(dyn_mat_tmp[j, i]))
                dyn_mat_tmp[i, j] = avg
                dyn_mat_tmp[j, i] = np.conj(avg)

        # Diagonalize
        omega2, eigvecs = eigh(dyn_mat_tmp)
        omega_q = np.sqrt(np.abs(omega2))
        omega_q[omega2 < 0] *= -1
        omega = omega_q

        # Normalize eigenvectors
        eig = np.zeros((nmodes, nmodes), dtype=np.complex128)
        for imode in range(nmodes):
            for jmode in range(nmodes):
                ia = jmode // 3
                eig[jmode, imode] = eigvecs[jmode, imode] / np.sqrt(atm_mass[ia])

        pol_vecs = eig

        nq_found += 1

    if nq_found == 0:
        raise ValueError("No dynamical matrices found")

    return omega

if __name__ == "__main__":
    filename = "/Users/murali/phd/one_phonon_raman/hbn/3D_hBN/ph/ph_interpolated/dyn"
    plot_ref = True
    omega = []
    for i in range(131):
        omega.append(read_dyn_qe_old(filename+str(i+1))*109737.315685)
    omega = np.array(omega)
    if plot_ref:
        from matplotlib import pyplot as plt
        data = np.loadtxt("./hbn.freq.gp")

        x = data[:, 0]      # column 1
        ys = data[:, 1:13]  # columns 2–10

        #print(omega[:,:4])
        # Plot all columns
        for i in range(ys.shape[1]):
            plt.plot(x, ys[:, i], 'r')
            plt.plot(x, omega[:, i], '--b')

        plt.xlabel("Column 1")
        plt.ylabel("Columns 2–10")
        plt.legend()
        plt.tight_layout()
        plt.show()
