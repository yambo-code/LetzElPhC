# Exciton-Phonon Coupling and Luminescence

In this advanced tutorial, we will calculate exciton-phonon interactions from first principles by interfacing **DFPT** (Quantum Espresso) and **BSE** (Yambo). We will compute exciton-phonon coupling matrix elements and apply them to phonon-assisted luminescence in **monolayer MoS$_2$** and **bulk hBN**.

## Requirements

*   Knowledge of Yambo GW-BSE and Quantum Espresso DFPT.
*   **Executables**: `pw.x`, `ph.x` (QE), `yambo`, `yambo_ph` (Yambo), `lelphc` (LetzElPhC).
*   **Python**: `yambopy`.

## Workflow Overview

1.  **QE**: SCF, NSCF, DVSCF (Phonons).
2.  **Yambo**: Initialize SAVE, Run BSE (Finite Momentum).
3.  **LetzElPhC**: Compute El-Ph matrix elements.
4.  **Yambopy**: Compute Exciton-Phonon coupling and observables (Luminescence).

---

## Step 0: Pseudopotentials & Convergence

Ensure pseudopotentials and lattice parameters are compatible. Verify `ecutwfc` convergence for DFPT. We assume these tests are done.

## Step 1: SCF Calculation (MoS$_2$)

Run standard SCF with `pw.x` using non-symmorphic symmetries.

**Input `mos2.scf`**:
```fortran
&control
    calculation = "scf",
    prefix = "mos2",
    ...
/&end
&system
    ibrav = 4, celldm(1) = 5.90, ...
    force_symmorphic = .true.
/&end
K_POINTS { automatic }
6 6 1 0 0 0
```
Run: `mpirun -np 4 pw.x -inp mos2.scf > scf.out`

## Step 2: NSCF Calculation for Yambo

Copy `mos2.save` from SCF. Run NSCF with correct k-grid (matching q-grid).

**Input `mos2.nscf`**:
```fortran
&control
    calculation = "nscf",
    nbnd = 250, ...
/&end
K_POINTS { automatic }
6 6 1 0 0 0
```
Run: `mpirun -np 4 pw.x -inp mos2.nscf > nscf.out`

## Step 3: DVSCF Phonon Calculation

Copy `mos2.save` from SCF. Run `ph.x` to get $\Delta V_{SCF}$.

**Input `mos2.dvscf`**:
```fortran
&inputph
  prefix='mos2',
  fildvscf = 'mos2-dvscf',
  electron_phonon = 'dvscf',
  ldisp=.true.,
  nq1=6, nq2=6, nq3=1
/
```
Run: `mpirun -np 8 ph.x -inp mos2.dvscf > dvscf.out`

## Step 4: Create Yambo SAVE

In the NSCF folder:
```bash
p2y
yambo
```

## Step 5: Run a BSE Calculation

We need a **finite-momentum** BSE calculation.

**Input `bse.in`** (generated via `yambo -X s -o b -k sex -y d -r`):
```bash
...BSEmod= "causal"
BSKmod= "SEX"
Lkind="full"          # Important: Include long-range exchange
% BSEQptR
 1 | 7 |              # Transferred momenta range (all Q points)
% ...
```
Run: `mpirun -np 8 yambo -F bse.in -J bse_Lfull -C bse_Lfull`

## Step 6: Obtain El-Ph Matrix Elements

Use `yambopy` to call `lelphc`.

```bash
yambopy l2y -ph path/to/mos2.dvscf -b 25 28 -par 4 2 -lelphc path/to/lelphc
```
*   `-b 25 28`: Band range (1-based index).
*   `-par 4 2`: 4 q-pools, 2 k-pools.

This generates `ndb.elph` (LetzElPhC format) and `SAVE/ndb.elph_gkkp*` (Yambo format).

## Step 7: Exciton-Phonon Coupling

We compute the coupling:

$$
\mathcal{G}^{\mu}_{\lambda} (0, q) = \sum_{cvk} A^{\lambda}_{cv}(k,q) g^{\mu}_{n m}(k,q)
$$

**Python Script (`calculate_excph.py`)**:
```python
from yambopy import YamboLatticeDB, YamboWFDB, LetzElphElectronPhononDB
from yambopy.exciton_phonon.excph_matrix_elements import exciton_phonon_matelem

# Paths
bsepath    = '1L_MoS2/bse-allq_full'
savepath   = '1L_MoS2/SAVE'
ndb_elph   = '1L_MoS2/ndb.elph'

# Load DBs
lattice = YamboLatticeDB.from_db_file(filename=f'{savepath}/ns.db1')
elph    = LetzElphElectronPhononDB(ndb_elph, read_all=False)
wfcs    = YamboWFDB(filename='ns.wf', save=savepath, latdb=lattice, bands_range=[24,28])

# Compute
exph = exciton_phonon_matelem(lattice, elph, wfcs, BSE_dir=bsepath,
                              neigs=12, dmat_mode='save', exph_file='MoS2_Ex-ph.npy')
```
Run: `python calculate_excph.py`

## Step 8: Phonon-Assisted Luminescence (hBN)

We model the luminescence $I(\omega)$ as a second-order process.

### Formula

$$
I(\omega) \propto \sum_{\nu \lambda} \left| \sum_{\alpha} \frac{ D^{\alpha} \mathcal{G}^{\nu}_{\lambda \alpha} }{ E_{\lambda}(q) - E_{\alpha}(0) - \omega_{\nu q} } \right|^2 \times N_{\lambda}(T_{exc}) (n_{\nu q}(T) + 1) \delta( \omega - (E_{\lambda}(q) - \omega_{\nu q}) )
$$

### HBN Example Runs
1.  **SCF/NSCF/DVSCF**: Standard QE runs for hBN.
2.  **BSE (Finite Q)**: `Lkind="full"`, all Q points (`bse_Lfull`).
3.  **BSE (Optical)**: `Lkind="bar"`, Q=1 (`bse_Lbar`).
4.  **LetzElPhC**: `yambopy l2y ...`

### Luminescence Script (`luminescence.py`)
```python
from yambopy.exciton_phonon.excph_luminescence import exc_ph_luminescence
from yambopy.exciton_phonon.excph_input_data import exc_ph_get_inputs

path = '3D_hBN'
input_data = exc_ph_get_inputs(f'{path}/SAVE', path, f'{path}/bse_Lfull',
                               bse_path2=f'{path}/bse_Lbar',
                               nexc_in=12, nexc_out=12,
                               bands_range=[6,10], phonons_range=[0,12])

ph_en, exc_en, exc_en_in, G, exc_dip = input_data

w, PL = exc_ph_luminescence(T_ph=10, ph_energies=ph_en, exc_energies=exc_en,
                            exc_dipoles=exc_dip, exc_ph_mat=G, exc_energies_in=exc_en_in,
                            exc_temp=10, emin=4.4, emax=4.7, broad=0.005)

# Save/Plot results...
```
This computes the luminescence spectrum showing phonon sidebands.
