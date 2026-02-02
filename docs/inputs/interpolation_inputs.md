These variables control the interpolation calculations.

=== "ph_save_dir"
    * **Type:** `STRING`
    * **Default:** `ph_save`
    * **Description:**
        Directory path containing the phonon save data to be read.

=== "ph_save_interpolation_dir"
    * **Type:** `STRING`
    * **Default:** `ph_save_interpolation`
    * **Description:**
        Directory path where the interpolation results will be written.

=== "interpolate_dvscf"
    * **Type:** `LOGICAL`
    * **Default:** `True`
    * **Description:**
        If True, the code interpolates the deformation potential ($\delta V_{scf}$). 

=== "asr"
    * **Type:** `LOGICAL`
    * **Default:** `True`
    * **Description:**
        If `True`, applies the Acoustic Sum Rule (ASR) to enforce translational invariance.

=== "asr_kind"
    * **Type:** `STRING`
    * **Default:** `simple`
    * **Description:**
        Specifies the type of Acoustic Sum Rule to apply.
        * Example: `simple`

=== "loto"
    * **Type:** `LOGICAL`
    * **Default:** `False`
    * **Description:**
        If `True`, applies non-analytic (LO-TO) splitting corrections to the dynamical matrices and potentials.

=== "loto_dir"
    * **Type:** `REAL (Array of 3)`
    * **Default:** `0.0 0.0 0.0`
    * **Description:**
        A vector specifying the direction approaching the $\Gamma$ point ($q \to 0$) for LO-TO splitting.
        * **Format:** Three space-separated floats.
        * **Example:** `loto_dir = 1.0 0.0 0.0` (along x-axis)

=== "nq1"
    * **Type:** `INTEGER`
    * **Default:** `1`
    * **Description:**
        Dimensions of the fine q-point grid used for interpolation.
        Number of q points along reciprocal lattice vector $b_1$.

=== "nq2"
    * **Type:** `INTEGER`
    * **Default:** `1`
    * **Description:**
        Dimensions of the fine q-point grid used for interpolation.
        Number of q points along reciprocal lattice vector $b_2$.

=== "nq3"
    * **Type:** `INTEGER`
    * **Default:** `1`
    * **Description:**
        Dimensions of the fine q-point grid used for interpolation.
        Number of q points along reciprocal lattice vector $b_3$.

=== "write_dVbare"
    * **Type:** `LOGICAL`
    * **Default:** `False`
    * **Description:**
        If `True`, the code will dump the change in the local part of the bare potential to disk.
        Note: in case `write_dVbare = True` and `interpolate_dvscf = False`, dyn* files are not written.

=== "eta_induced"
    * **Type:** `REAL`
    * **Default:** `1.0`
    * **Description:**
        Ewald parameter used for screened long-range interactions (dipoles, quadrupoles, phonon dynamical matrices).
        * Must be a positive value.

=== "eta_ph"
    * **Type:** `REAL`
    * **Default:** `1.0`
    * **Description:**
        Ewald summation parameter used specifically for the interpolation of dynamical matrices.
        * Must be a positive value.

=== "qlist_file"
    * **Type:** `STRING`
    * **Default:** `(Empty)`
    * **Description:**
        Filename containing a specific list of q-points provided in reduced units to interpolate.
        * If provided, the code interpolates over this list instead of the grid defined by `nq1`, `nq2`, `nq3`.
        * **Format:**
            ```
            Number of qpoints given in the list
            Q1x Q1y Q1z
            Q2x Q2y Q2z
            # is comment line and will be skiped
            ```
