/*
This file contains function that rotates eigen vectors
q -> Sq
*/
#include "symmetries.h"

void rotate_eig_vecs(struct symmetry* sym, const struct Lattice* lattice,
                     const ELPH_float* qpt, const ELPH_cmplx* eig_q,
                     ELPH_cmplx* restrict eig_Sq)
{
    // qpt is in crsytal coordinates (un rotated q point)

    const ELPH_float* blat = lattice->blat_vec;
    const ND_int natom = lattice->natom;
    const ND_int nmodes = lattice->nmodes;

    // find the equalivalent positions that the symmetry maps
    ELPH_float* atom_pos_crys = calloc(3 * natom, sizeof(ELPH_float));
    ND_int* rot_map = malloc(sizeof(ND_int) * natom);

    matmul_float('N', 'N', lattice->atomic_pos, blat, atom_pos_crys,
                 1.0 / (2 * ELPH_PI), 0.0, 3, 3, 3, natom, 3, 3);

    // apply symmetry to atomic positions
    for (ND_int ia = 0; ia < natom; ++ia)
    {
        const ELPH_float* pos_tmp = lattice->atomic_pos + 3 * ia;

        ELPH_float rot_atom[3] = { 0, 0, 0 };
        // apply rotation
        MatVec3f(sym->Rmat, pos_tmp, false, rot_atom);
        // add translation
        for (int ix = 0; ix < 3; ++ix)
        {
            rot_atom[ix] += sym->tau[ix];
        }
        ELPH_float rot_atom_crys[3] = { 0, 0, 0 };
        // convert to cystal
        MatVec3f(blat, rot_atom, true, rot_atom_crys);
        // remove 2*pi comming from blat
        for (int ix = 0; ix < 3; ++ix)
        {
            rot_atom_crys[ix] /= (2 * ELPH_PI);
        }
        // find the index
        // note we can use find_kidx_in_list function as the
        // functionality is not limited to reciprocal space
        // (contrary to what the name says so)
        rot_map[ia] = find_kidx_in_list(natom, atom_pos_crys,
                                        rot_atom_crys);

        if (rot_map[ia] < 0)
        {
            error_msg("Unable to find the equalivalent atom for the symmetry operation. "
                      "Wrong symmetry operation detected.");
        }
    }
    free(atom_pos_crys);

    ELPH_float qcart[3];
    ELPH_float Sq_cart[3];
    // Sq in cart units

    // convert qpt to cartisian coordinates
    MatVec3f(blat, qpt, false, qcart); // 2pi is included
    // compute S*q
    MatVec3f(sym->Rmat, qcart, false, Sq_cart); // 2pi is included

    ELPH_cmplx qphase = 0;
    for (int ix = 0; ix < 3; ++ix)
    {
        qphase += Sq_cart[ix] * sym->tau[ix];
    }
    qphase = cexp(-I * qphase);

    ELPH_cmplx* exp_Sqr = malloc(natom * sizeof(ELPH_cmplx));
    // exp(I*Sq*atomic_pos)
    ELPH_cmplx* exp_qr = malloc(natom * sizeof(ELPH_cmplx));
    // exp(-I*q*atomic_pos)

    ELPH_cmplx* rot_eig = calloc(nmodes, sizeof(ELPH_cmplx));

    ELPH_cmplx Rmat_cmplx[9];
    for (int i = 0; i < 9; ++i)
    {
        Rmat_cmplx[i] = sym->Rmat[i];
    }

    for (ND_int ia = 0; ia < natom; ++ia)
    {
        const ELPH_float* pos_tmp = lattice->atomic_pos + 3 * ia;
        exp_Sqr[ia] = 0;
        exp_qr[ia] = 0;

        for (int ix = 0; ix < 3; ++ix)
        {
            exp_Sqr[ia] += Sq_cart[ix] * pos_tmp[ix];
            exp_qr[ia] += qcart[ix] * pos_tmp[ix];
        }
        exp_Sqr[ia] = cexp(I * exp_Sqr[ia]);
        exp_qr[ia] = cexp(-I * exp_qr[ia]);
    }

    for (ND_int imode = 0; imode < nmodes; ++imode)
    {
        ELPH_cmplx* restrict eig_Sq_mode = eig_Sq + imode * nmodes;
        const ELPH_cmplx* eig_q_mode = eig_q + imode * nmodes;

        // S@eig_q_mode.T  = (3,natom) eig_q_mode@S^T =
        matmul_cmplx('N', 'T', eig_q_mode, Rmat_cmplx, rot_eig,
                     1.0, 0.0, 3, 3, 3, natom, 3, 3);

        for (ND_int ia = 0; ia < natom; ++ia)
        { // remove q phase
            ELPH_cmplx* restrict rot_eig_ia = rot_eig + ia * 3;
            for (int ix = 0; ix < 3; ++ix)
            {
                rot_eig_ia[ix] *= exp_qr[ia];
                rot_eig_ia[ix] *= qphase;
            }

            // index of atom that the symmetry maps ia atom
            ND_int ia_rot = rot_map[ia];
            ELPH_cmplx* restrict eig_Sq_mode_ia = eig_Sq_mode + ia_rot * 3;

            // add Sq phase back
            for (int ix = 0; ix < 3; ++ix)
            {
                // in case of time reversal symmetry we need to conjugate
                if (sym->time_rev)
                {
                    rot_eig_ia[ix] = -conj(rot_eig_ia[ix]);
                    // the negative sign is from rotation matrix
                }
                eig_Sq_mode_ia[ix] = exp_Sqr[ia_rot] * rot_eig_ia[ix];
            }
        }
    }

    free(rot_map);
    free(rot_eig);
    free(exp_Sqr);
    free(exp_qr);
}
