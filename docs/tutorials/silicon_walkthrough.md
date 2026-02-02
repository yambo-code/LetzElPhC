# Silicon Walkthrough (QE + LetzElPhC)

This tutorial provides a step-by-step walkthrough of the example located in `examples/qe/silicon`. We will compute the electron-phonon matrix elements for bulk Silicon.

## Prerequisites

*   Quantum Espresso (`pw.x`, `ph.x`)
*   Yambo (`p2y`, `yambo`)
*   LetzElPhC (`lelphc`)
*   The `Si.upf` pseudopotential (included in the example directory).

## Directory Structure

The example directory is organized as follows:

```
examples/qe/silicon/
├── Si.upf          # Pseudopotential
├── run.sh          # All-in-one script
├── scf/            # Step 1: SCF
├── ph/             # Step 2: Phonons
├── nscf/           # Step 3: NSCF & Yambo Setup
└── elph/           # Step 4: LetzElPhC Calculation
```

You can run the entire workflow using `./run.sh`, but we will go through each step manually to understand the process.

---

## Step 1: SCF Calculation

First, we perform a self-consistent field (SCF) calculation to obtain the ground state density.

**Input:** `scf/scf.in`

```fortran
&CONTROL
  prefix = 'si'
  calculation = 'scf'
  ...
/
&SYSTEM
  ibrav = 0, nat = 2, ntyp = 1
  ecutwfc = 40.0
  ...
/
...
K_POINTS automatic
    4  4  4  0  0  0
```

**Command:**
```bash
cd scf
pw.x < scf.in > scf.out
cd ..
```

---

## Step 2: Phonon Calculation

Next, we calculate the phonons and the deformation potentials ($\delta V_{SCF}$) on a uniform q-grid.

**Input:** `ph/ph.in`

```fortran
&inputph 
  nq1 = 2, nq2 = 2, nq3 = 2,   ! 2x2x2 q-grid
  fildyn = 'si.dyn', 
  tr2_ph = 1e-12, 
  prefix = 'si', 
  fildvscf = 'dvscf',          ! Important: Save potential changes
  ldisp = .true., 
  reduce_io = .true.
/ 
```

**Important Notes:**
*   We use a **2x2x2 q-grid**.
*   `fildvscf` is specified to save the potential variation files.

**Command:**
```bash
cd ph
cp -r ../scf/si.* .  # Copy SCF data
ph.x < ph.in > ph.out
cd ..
```

---

## Step 3: NSCF Calculation & Yambo Setup

We need the electronic wavefunctions on a k-grid that is commensurate with our phonon q-grid.

**Input:** `nscf/nscf.in`

```fortran
&CONTROL
  calculation = 'nscf'
  prefix = 'si'
  ...
/
&SYSTEM
  ...
  nbnd = 10                  ! Number of bands
  force_symmorphic = .true.  ! Crucial for Yambo compatibility
/
...
K_POINTS automatic
    4  4  4  0  0  0         ! 4x4x4 k-grid (commensurate with 2x2x2 q-grid)
```

**Command:**
```bash
cd nscf
cp -r ../scf/si.* .  # Copy SCF data
pw.x < nscf.in > nscf.out
```

### Initialize Yambo Database
Once the NSCF run is complete, we enter the save directory and initialize the Yambo database.

```bash
cd si.save
p2y        # Convert QE data to Yambo format
yambo      # Initialize Yambo SAVE directory
cd ../..
```

---

## Step 4: LetzElPhC Preprocessing

Before running the main calculation, we run the preprocessor in the phonon directory to organize the phonon data into a `ph_save` directory.

**Command:**
```bash
cd ph
lelphc -pp --code=qe -F ph.in
cd ..
```
This creates `ph/ph_save`.

---

## Step 5: LetzElPhC Calculation

Finally, we compute the electron-phonon matrix elements.

**Input:** `elph/elph.in`

```ini
nkpool          = 1 
nqpool          = 1

start_bnd       = 2
end_bnd         = 6

save_dir        = ../nscf/si.save/SAVE   # Path to Yambo SAVE
ph_save_dir     = ../ph/ph_save          # Path to ph_save

kernel          = dfpt                   # Use full DFPT screening
convention      = standard               # <k+q|dV|k>
```

**Command:**
```bash
cd elph
lelphc -F elph.in
```

### Output
The code will generate:
*   Standard output with progress information.
*   **`ndb.elph`**: The database containing the electron-phonon matrix elements (in LetzElPhC format).
*   If running with Yambopy or converting later, this can be transformed into Yambo's `ndb.elph_gkkp` databases.
