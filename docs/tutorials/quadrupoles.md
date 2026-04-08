# Extracting Dynamical Quadrupoles for Electron–Phonon Interpolation

In polar and piezoelectric materials, atomic displacements generate microscopic electric quadrupoles [(Royo et al., PRX 2019)](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.9.021050).
In 2D materials, these quadrupoles are also responsible for coupling out-of-plane phonon modes to in-plane
electric fields [(Royo et al., PRX 2021)](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.041027).
To achieve accurate electron–phonon interpolation in such systems, the long-range analytical
contribution of these quadrupoles must be subtracted from the change in potentials on the coarse
grid. This leaves a perfectly smooth, short-range residual that can be safely interpolated onto a
fine grid.

Unlike Born effective charges, which are a byproduct of standard phonon calculations, the quadrupole
tensors require a dedicated extraction procedure. This tutorial covers how to compute them using
Quantum ESPRESSO (≥ 7.5) and the [`multipole.py`](https://gitlab.com/QEF/q-e/-/blob/develop/test-suite/ph_multipole/multipole.py) script.

---

## Step 1: Prerequisites

Before starting, ensure the following are in place:

- A converged SCF calculation with its `.save` folder present in the working directory.
- The SCF input file is named **`scf.in`**.
- Quantum ESPRESSO version **≥ 7.5** is available.

---


## Step 2: Preparing the `ph.in` Template

`multipole.py`  requires a template `ph.in` file. The template is provided below — save it as `ph.in` in the same directory containing your `.save` folder.
```fortran title="ph.in"
phonons on a grid
 &inputph
  reduce_io   = .true.
  outdir      = '.', ! Set you output directory
  prefix      = 'your_prefix', ! your scf prefix
  tr2_ph      = 1e-20,       ! Set this bit high for good values
  search_sym  = .false.,     ! Symmetry must be off for random q-point clouds.
  lmultipole  = .true.,      ! Writes the charge-density response to disk.
  fildrho     = 'drho.dat.{}',
  fildvscf    = 'drhodv.dat.{}'
 /
```

!!! warning ""
    - **`tr2_ph`** : A loose threshold produces noisy response functions that degrade the fit.

---

## Step 3: Generating the DFPT Input Files

After creating the `ph.in` template file, run `multipole.py` in preparation mode from the directory containing your `.save` folder:
```bash
python3 multipole.py -e --order 3 --epsil_order 4 -p --alat 1B --q_min 0.005 --q_max 0.015 --nq 25
```
This command produces `qpoints.dat` and a set of `ph.in.X` input files, one per q-point.

One of the important things to notice is the number of free components for the multipole tensors printed in the summary:
```text
Summary of multipole tensors:
- Born effective charges
  Size: 27  | Free: 4
  Index 0   | Mo  | Free: 2
  Index 2   | Te  | Free: 2
  Index 3   | Te  | Free: 2
- Quadrupoles
  Size: 81  | Free: 12
  Index 0   | Mo  | Free: 4
  Index 2   | Te  | Free: 4
  Index 3   | Te  | Free: 4
- Total number of free components: 16
```

!!! warning "Choosing `--nq`"
    The number of q-points must always exceed the **Total number of free components**. As a rule of thumb,
    use at least **1.5 ×** this value, and verify the quality of the fit afterwards — see
    [Step 6](#step-6-interpreting-the-diagnostics) for guidance.

!!! warning "Non-polar materials"
    Add the flag --non_polar for non_polar materials.
**Flag reference:**

| Flag | Purpose |
|---|---|
| `-p` | Preparation mode — generates `qpoints.dat` and the numbered `ph.in.X` input files. |
| `--order 3` | Fits up to quadrupole order ($\mathcal{O}(q^3)$). |
| `-e` / `--epsil_order 4` | **3D bulk only.** Simultaneously fits and subtracts the macroscopic electric-field divergence alongside the multipoles. |
| `--alat 1B` | Sets the lattice parameter units to Bohr (a.u.). Here, `1B` represents `1 Bohr`. **Replace** `1` with the actual numerical values. You can easily find the value in the output of your `scf` run (search for `lattice parameter (alat)`).  |
| `--q_min` / `--q_max` | Defines the spherical q-point shell (in units of $2\pi/a$). This needs to checked carefully (see [Step 6](#step-6-interpreting-the-diagnostics)). A good starting range is `q_min = 0.001 * alat` and `q_max = 0.003*alat` |
| `--nq 15` | Number of randomised q-points within the shell. |

!!! note "2D materials"
    For 2D materials, it is recommended to perform a separate `scf` run without the `assume_isolated` flag, solely for computing the quadrupoles.


!!! tip "Tip"
    In most cases, it is better to turn off spin-orbit coupling to reduce computational cost.

---

## Step 4: Running the DFPT Calculations

Submit the calculations for each generated q-point:
```bash
srun -n ${SLURM_NTASKS} ph.x -nk 8 < ph.in.X | tee ph.out
```

!!! tip "Computational cost"
    Each q-point requires only a **single** self-consistent response loop, making these runs
    significantly cheaper than a full phonon calculation (phonons require roughly $3 \times natom$). The investment is well worth the
    improvement in phonon deformation potentials and interpolation accuracy.

---

## Step 5: Fitting the Quadrupole Tensors

Once all `drho.dat.X` files have been generated, run the least-squares fit using the **same physics
flags** as in Step 2:
```bash
python3 multipole.py -f --order 3 --epsil_order 4 --alat 1B --nq 25
```

The script will:

1. Apply your crystal's spatial symmetry constraints to the tensor components.
2. Analytically subtract the macroscopic electric-field contribution (if `-e` was specified).
3. Solve an overdetermined least-squares system for the multipole tensors.

After the fit completes, copy the output file to your `ph_save` directory:
```bash
cp quadrupole.fmt /path/to/ph_save/
```

---

## Step 6: Interpreting the Diagnostics

The script prints a diagnostic block similar to the following:
```text
Summary of multipole tensors fitting:
- Free multipole components:    10
- Rank of coefficient matrix:   10
- RMSE:                         8.57e-04 |e| / bohr³
- Mean relative error:          7.20e+07
- Relative norm error:          0.0626
- R² score:                     0.9961
```

Each metric answers a specific question about the quality of your fit. If the fit quality is poor, you need to rerun the calculations with better inputs, mostly adjusting `q_min`, `q_max`, and `nq`.

---

#### 1. Matrix Rank — The Pass/Fail Check

> **Rule:** `Rank of coefficient matrix` ≥ `Free multipole components`

**Free multipole components** is the number of symmetry-independent tensor elements your crystal
requires (e.g., 10 for a wurtzite structure). **Rank** is the number of linearly independent
equations provided by your q-point cloud. If the rank falls short, the system is underdetermined
and the fit is unreliable — increase `--nq` and rerun.

---

#### 2. R² Score — The Gold Standard

> **Target:** R² > 0.99

This measures how well the fitted analytical surface reproduces the raw DFT response data *after*
the electric-field contribution has been removed. Interpretation:

| R² value | Diagnosis |
|---|---|
| ≥ 0.99 | Excellent fit. |
| 0.90 – 0.99 | Acceptable; consider tightening `tr2_ph` or adjusting the q-bracket. |
| < 0.90 | Likely noise-dominated or q-points are outside the valid analytical regime. Adjust `--q_min`/`--q_max`. |
| < 0.0 | The 3D non-analytic electric field was not properly subtracted. Verify the `-e` flag is being used for 3D bulk. |

---

#### 3. Relative Norm Error — The Magnitude Check

> **Target:** < 0.10

This compares the Euclidean norm of the fitted tensor to the norm of the DFT data vector. A value
of 0.06 means the overall magnitude of the fit deviates from the data by only 6% — an exceptionally
tight result for complex polar materials.

---

#### 4. RMSE — The Absolute Scale

> **Target:** 10⁻³ to 10⁻⁵ \|e\| / bohr³

The root-mean-square error gives the absolute physical residual between the fit and the data. Values
in this range confirm that the discrepancy is microscopic relative to the physical scale of the
response.

---

#### 5. Mean Relative Error — The Illusion to Ignore

> **Rule:** Disregard this metric for symmetric and piezoelectric crystals.

Crystal symmetry often forces specific tensor components to be *exactly* zero. In practice, DFT
outputs floating-point noise on the order of $10^{-15}$ rather than a true zero. The
component-by-component division used to compute this metric then divides by near-zero denominators,
causing the mean to explode to tens of millions of percent. **If your R² and RMSE are excellent,
this number is meaningless and should be ignored entirely.**
