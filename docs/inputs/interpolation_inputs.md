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

=== "nosym"
    * **Type:** `LOGICAL`
    * **Default:** `False`
    * **Description:**
        If `True`, no symmetries are used; i.e., dynamical matrices and dvscf are
        printed on the full Brillouin zone.


=== "asr"
    * **Type:** `STRING`
    * **Default:** `no`
    * **Description:**
        Specifies the type of Acoustic Sum Rule to apply.
        * `no`: No acoustic sum rule applied  
        * `simple`: Simple acoustic sum rule (corrects only the diagonal matrix elements)  
        * `crystal`: Acoustic sum rule imposed via optimized correction of the force constants

=== "loto"
    * **Type:** `LOGICAL`
    * **Default:** `False`
    * **Description:**
        If `True`, applies non-analytic (LO-TO) splitting corrections to the dynamical matrices and potentials.

        If `qlist_file` is provided instead of `nq1`, `nq2`, and `nq3`, then `loto` and `loto_dir` variables are ignored.
        See `qlist_file` for details.



=== "loto_dir"
    * **Type:** `REAL (Array of 3)`
    * **Default:** `0.0 0.0 0.0`
    * **Description:**
        A vector (in cartesian units) specifying the direction approaching the $\Gamma$ point ($q \to 0$) for LO-TO splitting.
        Internally, the code normalizes the vector, so the user does not need to provide a normalized vector.

        If the magnitude of the vector is smaller than $10^{-5}$, it is ignored and treated as `loto = false`.

        * **Format:** Three space-separated floats.
        * **Example:** `loto_dir = 1.0 0.0 0.0` (along x-axis)

        If `qlist_file` is provided instead of `nq1`, `nq2`, and `nq3`, then `loto` and `loto_dir` variables are ignored.
        See `qlist_file` for details.

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
            # Use # to comment. It will be skiped when reading
            ```

        !!! warning "Note on LO-TO splitting"
            When `qlist_file` is provided, the `loto` and `loto_dir` variables are ignored. To obtain LO–TO splitting, include a small direction vector at the $\Gamma$ point. Instead of sampling exactly at $\Gamma$, use a displaced point $\vec{q} = \Gamma + \vec{\epsilon}$, where $\lvert \vec{\epsilon} \rvert$ is small and defines the direction of approach. When $\vec{\epsilon}$ vector is converted to caresian units, its magnitude should be small; a recommended value is $|\vec{\epsilon}_{cart}| = 10^{-4}$.
