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
        Specifies the type of Acoustic Sum Rule to apply for force constant matrix.
        * `no`: No acoustic sum rule applied.
        * `simple`: Simple acoustic sum rule (corrects only the diagonal matrix elements).
        * `crystal`: Translational acoustic sum rule imposed via orthogonal projection.
        * `all`: Translational and rotational acoustic sum rules imposed via orthogonal projection.
        * `all_huang`: (recommended) Translational and rotational sum rules, plus Huang invariances, imposed via orthogonal projection.

=== "zasr"
    * **Type:** `STRING`
    * **Default:** `no`
    * **Description:**
        Specifies the type of Acoustic Sum Rule to apply for Born charges.
        * `no`: No acoustic sum rule applied.
        * `simple`: Simple acoustic sum rule (corrects only the diagonal matrix elements).
        * `crystal`: Acoustic sum rule imposed via orthogonal projection.


=== "loto"
    * **Type:** `LOGICAL`
    * **Default:** `False`
    * **Description:**
        If `True`, applies non-analytic (LO-TO) splitting corrections to the dynamical matrices at $\Gamma$ point.

        It should be noted that this flag does not completely switch off LO-TO treatment. 
        It is only intended for switch on/off LO/TO treatment at the Gamma point.
        If you wish to completely disable LO-TO splitting, please remove or rename the `tensors.xml` file in the `ph_save` directory.

        If `qlist_file` is provided instead of `nq1`, `nq2`, and `nq3`, then `loto` and `loto_dir` variables are ignored.
        See `qlist_file` for details.



=== "loto_dir"
    * **Type:** `REAL (Array of 3)`
    * **Default:** `0.0 0.0 0.0`
    * **Description:**
        A vector (in cartesian units) specifying the direction approaching the $\Gamma$ point ($q \to 0$) for LO-TO splitting.
        Internally, the code normalizes the vector, so the user does not need to provide a normalized vector.
        
        If `loto = false`, this flag is ignored.

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
        Ewald parameter used for seperation of long-range and short range potentials. Used for interpolation of $\partial V_{scf}$
        * Must be a positive value. For 2D systems, must be `eta_induced >= 1`
        * In the case of 3D, a Gaussian decay factor $e^{-\frac{|\mathbf{q}+\mathbf{G}|^2}{4\eta_{\text{induced}}}}$ is applied. For 2D, the macroscopic decay is governed by $1 - \tanh\left(\frac{|\mathbf{q}+\mathbf{G}|_{\parallel} L}{2}\right)$, where the effective slab thickness is bounded by the out-of-plane polarizability condition $L = 4\pi\alpha_{\perp} \times 1.001 \times \eta_{\text{induced}}$.

=== "eta_ph"
    * **Type:** `REAL`
    * **Default:** `1.0`
    * **Description:**
        Ewald parameter used for seperation of long-range and short range force constants. Used specifically for the interpolation of dynamical matrices.
        * Must be a positive value. For 2D systems, must be `eta_ph >= 1`
        * In the case of 3D, a Gaussian decay factor $e^{-\frac{|\mathbf{q}+\mathbf{G}|^2}{4\eta_{\text{ph}}}}$ is applied. For 2D, the macroscopic decay is governed by $1 - \tanh\left(\frac{|\mathbf{q}+\mathbf{G}|_{\parallel} L}{2}\right)$, where the effective slab thickness is bounded by the out-of-plane polarizability condition $L = 4\pi\alpha_{\perp} \times 1.001 \times \eta_{\text{ph}}$.

=== "qlist_file"
    * **Type:** `STRING`
    * **Default:** `""`
    * **Description:**
        Filename containing a specific list of q-points provided in reduced units to interpolate.
        * If a non-zero length string is provided, the code interpolates over this list instead of the grid defined by `nq1`, `nq2`, `nq3`.
        * **Format:**
            ```
            Number of qpoints given in the list
            Q1x Q1y Q1z
            Q2x Q2y Q2z
            # Use # to comment. It will be skiped when reading
            ```

        !!! warning "Note on LO-TO splitting"
            When `qlist_file` is provided, the `loto` and `loto_dir` variables are ignored. To obtain LO–TO splitting, include a small direction vector at the $\Gamma$ point. Instead of sampling exactly at $\Gamma$, use a displaced point $\vec{q} = \Gamma + \vec{\epsilon}$, where $\lvert \vec{\epsilon} \rvert$ is small and defines the direction of approach. When $\vec{\epsilon}$ vector is converted to caresian units, its magnitude should be small; a recommended value is $|\vec{\epsilon}_{cart}| = 10^{-4}$.
