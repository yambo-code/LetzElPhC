#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/constants.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/string_func.h"
#include "elphC.h"
#include "phonon.h"

// Context for Force Constant ASR (Implicit MatVec)
struct ASR_Context
{
    // Dimensions
    ND_int natom;
    ND_int nmodes;  // 3 * natom
    ND_int nconstraints;
    ND_int ngrid;
    ND_int frc_size;  // Size of the flattened real force constant vector
    // Geometry & Grid
    const ND_int* qgrid;
    const ELPH_float* atomic_pos;
    const ELPH_float* lat_vecs;
    const ND_int* ws_vecs;
    ND_int n_ws_vecs;
    const ND_int* ws_degen;
    // Mapping arrays
    const ND_int* grid_map;
    const bool* do_grid;
    // Configuration
    enum asr_kind mode;
    bool huang;
    // Derived Constants
    ND_int Gridyz;
};

// Context for Born Charge ASR (Implicit MatVec)
struct Born_Context
{
    ND_int natom;
    ND_int dim_Z;         // 9 * natom
    ND_int nconstraints;  // 9
};

// static function
static void get_huang_indices(const ND_int idx, ND_int* restrict a,
                              ND_int* restrict b, ND_int* restrict c,
                              ND_int* restrict d);

static void born_matvec(const int mode, const double* restrict x,
                        double* restrict y, void* userdata);

static void asr_matvec(const int mode, const double* restrict x,
                       double* restrict y, void* userdata);

static inline void update_fc_Amat_entry(const int mode, const double* x,
                                        double* restrict y,
                                        const ND_int idx_constr,
                                        const ND_int idx_phi, const double val);

static inline ND_int get_fc_index(const struct ASR_Context* ctx,
                                  const ND_int iig, const bool process_R,
                                  const bool is_gamma, const ND_int ia,
                                  const ND_int alpha, const ND_int ja,
                                  const ND_int beta);

// public functions
enum asr_kind asr_kind_from_string(const char* str, bool print_warning)
{
    if (str == NULL)
    {
        return ASR_NONE;
    }

    // Optional: enforce max length (32 chars)
    if (strlen(str) > 32)
    {
        return ASR_NONE;
    }

    if (my_strcasecmp(str, "no") == 0)
    {
        return ASR_NONE;
    }
    if (my_strcasecmp(str, "simple") == 0)
    {
        return ASR_SIMPLE;
    }
    if (my_strcasecmp(str, "crystal") == 0)
    {
        return ASR_CRYSTAL;
    }
    if (my_strcasecmp(str, "all") == 0)
    {
        return ASR_ALL;
    }
    if (my_strcasecmp(str, "all_huang") == 0)
    {
        return ASR_ALL_HUANG;
    }

    if (print_warning)
    {
        fprintf(stdout,
                "Warning : Invalid asr/zasr string : %s. Switching off "
                "Acoustic sum rule.\n",
                str);
    }

    return ASR_NONE;
}

void asr_kind_to_string(enum asr_kind kind, char* buf)
{
    // buf must be alteast of size 32
    if (buf == NULL)
    {
        return;
    }

    if (kind == ASR_NONE)
    {
        strlcpy_custom(buf, "no", 32);
        return;
    }
    if (kind == ASR_SIMPLE)
    {
        strlcpy_custom(buf, "simple", 32);
        return;
    }
    if (kind == ASR_CRYSTAL)
    {
        strlcpy_custom(buf, "crystal", 32);
        return;
    }
    if (kind == ASR_ALL)
    {
        strlcpy_custom(buf, "all", 32);
        return;
    }
    if (kind == ASR_ALL_HUANG)
    {
        strlcpy_custom(buf, "all_huang", 32);
        return;
    }

    strlcpy_custom(buf, "no", 32);
}

void apply_acoustic_sum_rule_born_charges(enum asr_kind mode, ELPH_float* Zborn,
                                          const ND_int natom)
{
    // asr for born charges
    if (mode == ASR_NONE || !Zborn || !natom)
    {
        return;
    }

    if (mode == ASR_SIMPLE)
    {
        /* ==============================
         * SIMPLE RULE
         * ==============================
         */
        for (ND_int i = 0; i < 3; i++)
        {
            for (ND_int j = 0; j < 3; j++)
            {
                ELPH_float sum = 0.0;
                for (ND_int ia = 0; ia < natom; ia++)
                {
                    sum += Zborn[ia * 9 + i * 3 + j];
                }

                ELPH_float avg = sum / ((ELPH_float)natom);

                for (ND_int ia = 0; ia < natom; ia++)
                {
                    Zborn[ia * 9 + i * 3 + j] -= avg;
                }
            }
        }
    }
    else
    {
        /* ==============================
         * CRYSTAL RULE (USING LSMR)
         * ==============================
         */
        // Problem: minimize ||Z - Z_input|| subject to A * Z = 0
        // Solution: Z = Z_input - A^T * lambda
        // Solve A A^T * lambda = A * Z_input using LSMR.

        ND_int dim_Z = 9 * natom;
        ND_int nconstraints = 9;

        // 1. Prepare Input Z in double precision
        double* Z_real = malloc(dim_Z * sizeof(*Z_real));
        CHECK_ALLOC(Z_real);

        for (ND_int i = 0; i < dim_Z; ++i)
        {
            Z_real[i] = Zborn[i];
        }

        // 2. Setup Context
        struct Born_Context ctx = {
            .natom = natom, .dim_Z = dim_Z, .nconstraints = nconstraints};

        // 3. Solve using LSMR
        double* lambda = calloc(nconstraints, sizeof(*lambda));
        CHECK_ALLOC(lambda);

        // let the compiler remove this loop
        for (ND_int i = 0; i < nconstraints; ++i)
        {
            lambda[i] = 0.0;
        }

        ND_int itn = 0;
        // Solving min || A^T*lambda - Z_input ||
        // Note: implicit matvec A^T is passed as 'M' to LSMR.
        // LSMR solves M*x = b => A^T * lambda = Z_input?
        // No, standard projection is:
        // A * (Z_in - A^T * lambda) = 0  => (A A^T) lambda = A * Z_in.
        // Wait, LSMR minimizes || M x - b ||.
        // If we set M = A^T and b = Z_input,
        // It finds x (lambda) that minimizes || A^T lambda - Z_input ||.
        // The residual r = b - M x = Z_input - A^T lambda.
        // This residual 'r' IS strictly the projected vector Z that satisfies
        // A*Z=0 (if A*Z=0 is the nullspace). Yes, projection onto null(A) is
        // equivalent to removing the range(A^T).

        int status =
            lsmr_solver(dim_Z, nconstraints, born_matvec, &ctx, Z_real, 0.0,
                        1e-10, 1e-10, 1e12, nconstraints * 2, lambda, &itn);

        if (status != 0 && status != 1 && status != 2)
        {
            fprintf(stderr,
                    "Warning: Born LSMR solver status %d after %lld iters\n",
                    status, itn);
        }

        // 4. Compute Correction (A^T * lambda)
        double* correction = calloc(dim_Z, sizeof(*correction));
        CHECK_ALLOC(correction);

        born_matvec(0, lambda, correction, &ctx);  // Mode 0 = A^T * lambda

        // 5. Apply correction
        for (ND_int i = 0; i < dim_Z; ++i)
        {
            // Z_projected = Z_input - Correction
            Zborn[i] = (Z_real[i] - correction[i]);
        }

        free(lambda);
        free(correction);
        free(Z_real);
    }
}

void apply_acoustic_sum_rule_fc(enum asr_kind mode, const ND_int* qgrid,
                                const ND_int natom, ELPH_cmplx* frc,
                                const ELPH_float* atomic_pos,
                                const ELPH_float* lat_vecs,
                                const ND_int* ws_vecs, const ND_int n_ws_vecs,
                                const ND_int* ws_degen)
{
    // asr for force constants
    if (mode == ASR_NONE)
    {
        return;
    }

    bool huang = false;
    if (mode == ASR_ALL_HUANG)
    {
        mode = ASR_ALL;
        huang = true;
    }

    const ND_int ngrid = qgrid[0] * qgrid[1] * qgrid[2];
    const ND_int nmodes = 3 * natom;
    const ND_int Gridyz = qgrid[1] * qgrid[2];

    // Simple ASR
    if (mode == ASR_SIMPLE)
    {
        for (ND_int ia = 0; ia < natom; ++ia)
        {
            for (ND_int alpha = 0; alpha < 3; ++alpha)
            {
                for (ND_int beta = 0; beta < 3; ++beta)
                {
                    ELPH_cmplx sum = 0.0;
                    ND_int base_na_i_j = ia * 9 * natom + alpha * nmodes + beta;
                    for (ND_int ig = 0; ig < ngrid; ++ig)
                    {
                        ND_int grid_off = ig * nmodes * nmodes;
                        for (ND_int ib = 0; ib < natom; ++ib)
                        {
                            sum += frc[grid_off + base_na_i_j + ib * 3];
                        }
                    }
                    ND_int self_idx = base_na_i_j + ia * 3;
                    frc[self_idx] -= sum;
                }
            }
        }
        return;
    }

    /* -----------------------------------------------------------
     * Advanced ASR using LSMR (Matrix-Free Orthogonal Projection)
     * ----------------------------------------------------------- */

    // 1. Setup Grid Mapping (Identify independent R points vs -R points)
    ND_int* grid_map = calloc(ngrid, sizeof(*grid_map));
    bool* do_grid = malloc(ngrid * sizeof(*do_grid));
    CHECK_ALLOC(grid_map);
    CHECK_ALLOC(do_grid);

    for (ND_int ig = 0; ig < ngrid; ++ig)
    {
        do_grid[ig] = true;
    }

    ND_int ngrid_independent = 0;
    for (ND_int ig = 0; ig < ngrid; ++ig)
    {
        if (!do_grid[ig])
        {
            continue;
        }

        ND_int Rx = ig / Gridyz;
        ND_int Ry = (ig % Gridyz) / qgrid[2];
        ND_int Rz = (ig % Gridyz) % qgrid[2];

        ND_int mRx = (qgrid[0] - Rx) % qgrid[0];
        ND_int mRy = (qgrid[1] - Ry) % qgrid[1];
        ND_int mRz = (qgrid[2] - Rz) % qgrid[2];
        ND_int mig = mRz + mRy * qgrid[2] + mRx * Gridyz;

        grid_map[ig] = ngrid_independent;
        if (ig != mig)
        {
            do_grid[mig] = false;
            grid_map[mig] = ngrid_independent;
        }
        ++ngrid_independent;
    }

    // 2. Extract Real Force Constants into a flat double array
    const ND_int frc_size = ngrid_independent * nmodes * nmodes;
    double* frc_real = malloc(frc_size * sizeof(*frc_real));
    CHECK_ALLOC(frc_real);

    ND_int iig_tmp = 0;
    for (ND_int ig = 0; ig < ngrid; ++ig)
    {
        if (!do_grid[ig])
        {
            continue;
        }

        // Copy real parts
        for (ND_int i = 0; i < nmodes * nmodes; ++i)
        {
            frc_real[i + iig_tmp * nmodes * nmodes] =
                creal(frc[i + ig * nmodes * nmodes]);
        }
        ++iig_tmp;
    }
    // Enforce Symmetry at Gamma point (R=0)
    // We zero out the upper triangle so the solver treats it as
    // dependent/zero. We average (lower+upper)/2 into lower triangle.
    // now frc
    for (ND_int i = 0; i < nmodes; ++i)
    {
        for (ND_int j = 0; j < i; ++j)
        {
            frc_real[j * nmodes + i] = 0.0;
        }
    }

    // 3. Setup Implicit Context
    ND_int nconstraints = 9 * natom;
    if (mode == ASR_ALL)
    {
        nconstraints += 9 * natom;
        if (huang)
        {
            nconstraints += 15;
        }
    }

    struct ASR_Context ctx = {.natom = natom,
                              .nmodes = nmodes,
                              .nconstraints = nconstraints,
                              .ngrid = ngrid,
                              .frc_size = frc_size,
                              .qgrid = qgrid,
                              .atomic_pos = atomic_pos,
                              .lat_vecs = lat_vecs,
                              .ws_vecs = ws_vecs,
                              .n_ws_vecs = n_ws_vecs,
                              .ws_degen = ws_degen,
                              .grid_map = grid_map,
                              .do_grid = do_grid,
                              .mode = mode,
                              .huang = huang,
                              .Gridyz = Gridyz};

    // 4. Solve using LSMR
    // Solves min || M*x - b || with M = A^T, b = frc_real.
    // Result x is lambda.
    double* lambda = calloc(nconstraints, sizeof(*lambda));
    CHECK_ALLOC(lambda);

    // let the compiler remove this loop
    for (ND_int i = 0; i < nconstraints; ++i)
    {
        lambda[i] = 0.0;
    }

    ND_int itn = 0;
    int status =
        lsmr_solver(frc_size, nconstraints, asr_matvec, &ctx, frc_real, 0.0,
                    1e-10, 1e-10, 1e12, nconstraints * 2, lambda, &itn);

    if (status != 0 && status != 1 && status != 2)
    {
        fprintf(
            stderr,
            "Warning: LSMR solver stopped with status %d after %lld iters\n",
            status, itn);
    }

    // 5. Compute Final Projected Force Constants
    // r = b - M*x = frc_real - A^T * lambda
    double* correction = calloc(frc_size, sizeof(*correction));
    CHECK_ALLOC(correction);

    // Compute A^T * lambda into 'correction'
    asr_matvec(0, lambda, correction, &ctx);

    // Apply correction to internal real array
    for (ND_int i = 0; i < frc_size; ++i)
    {
        frc_real[i] -= correction[i];
    }

    free(correction);
    free(lambda);

    // Now recontruct the full force constant matrix
    for (ND_int ig = 0; ig < ngrid; ++ig)
    {
        if (!do_grid[ig])
        {
            continue;
        }
        ND_int iig = grid_map[ig];

        ND_int Rx = ig / Gridyz;
        ND_int Ry = (ig % Gridyz) / qgrid[2];
        ND_int Rz = (ig % Gridyz) % qgrid[2];

        // Calculate -R index (mig)
        ND_int mRx = (qgrid[0] - Rx) % qgrid[0];
        ND_int mRy = (qgrid[1] - Ry) % qgrid[1];
        ND_int mRz = (qgrid[2] - Rz) % qgrid[2];
        ND_int mig = mRz + mRy * qgrid[2] + mRx * Gridyz;
        //
        for (ND_int i = 0; i < nmodes; ++i)
        {
            for (ND_int j = 0; j < nmodes; ++j)
            {
                frc[ig * nmodes * nmodes + i * nmodes + j] =
                    frc_real[iig * nmodes * nmodes + i * nmodes + j];
                if (ig != mig)
                {
                    frc[mig * nmodes * nmodes + j * nmodes + i] =
                        frc_real[iig * nmodes * nmodes + i * nmodes + j];
                }
            }
        }
    }

    // Symmetrization at Gamma point
    for (ND_int i = 0; i < nmodes; ++i)
    {
        for (ND_int j = 0; j < i; ++j)
        {
            frc[j * nmodes + i] = frc[i * nmodes + j];
        }
    }

    free(do_grid);
    free(grid_map);
    free(frc_real);
}

// static helpers
static void born_matvec(const int mode, const double* restrict x,
                        double* restrict y, void* userdata)
{
    /* =========================================================================
     * BORN CHARGE IMPLICIT MATVEC
     * =========================================================================
     */
    // Solves A * Z = 0
    // A is 9 x (9*natom) matrix.
    // A_{ic, ia_alpha_beta} = delta_{alpha, ic_alpha} * delta_{beta, ic_beta}
    // Effectively: Sum over atoms for each (alpha, beta) component.
    //
    // if mode = 0, Mx else M^T.x
    // Here M = A^T
    struct Born_Context* ctx = userdata;
    ND_int natom = ctx->natom;

    // Initialize output y to 0
    ND_int out_len = (mode == 0) ? (9 * natom) : 9;
    for (ND_int i = 0; i < out_len; ++i)
    {
        y[i] = 0.0;
    }

    // Iterate over constraints (9 components of the tensor)
    for (ND_int ic = 0; ic < 9; ++ic)
    {
        ND_int alpha = ic / 3;
        ND_int beta = ic % 3;

        // Loop over atoms (columns of A)
        for (ND_int ia = 0; ia < natom; ++ia)
        {
            // The index in Z vector for atom ia, component (alpha, beta)
            ND_int z_idx = ia * 9 + 3 * alpha + beta;

            // Note for mode = 0, we do A^Tx, because out M = A^T.
            if (mode == 0)
            {
                // y = A^T * x
                // A^T maps constraint (ic) back to all atoms
                // A_{ic, z_idx} is 1.0.
                // y[z_idx] += A_{ic, z_idx} * x[ic]
                y[z_idx] += x[ic];
            }
            else
            {
                // y = A * x
                // Sum over all atoms for this component
                // y[ic] += A_{ic, z_idx} * x[z_idx]
                y[ic] += x[z_idx];
            }
        }
    }
}

static void asr_matvec(const int mode, const double* restrict x,
                       double* restrict y, void* userdata)
{
    // matrix vector function for force constants
    // M * x. for mode = 0, M = A^T, ese M = A
    struct ASR_Context* ctx = userdata;

    // Unpack dimensions
    const ND_int natom = ctx->natom;
    const ND_int ngrid = ctx->ngrid;
    const ND_int* qgrid = ctx->qgrid;
    const ND_int Gridyz = ctx->Gridyz;

    // Initialize Output
    ND_int out_size = (mode == 0) ? ctx->frc_size : ctx->nconstraints;
    for (ND_int i = 0; i < out_size; ++i)
    {
        y[i] = 0.0;
    }

    // Constraint offsets
    const ND_int off_rot = 9 * natom;
    const ND_int off_huang = off_rot + 9 * natom;

    // --- Main Loop ---
    ND_int iws_vec = 0;

    for (ND_int ig = 0; ig < ngrid; ++ig)
    {
        ND_int iig = ctx->grid_map[ig];
        bool process_R = ctx->do_grid[ig];
        bool is_gamma = (ig == 0);

        // Miller indices
        ND_int Rx = ig / Gridyz;
        ND_int Ry = (ig % Gridyz) / qgrid[2];
        ND_int Rz = (ig % Gridyz) % qgrid[2];
        ND_int mx = get_miller_idx(Rx, qgrid[0]);
        ND_int my = get_miller_idx(Ry, qgrid[1]);
        ND_int mz = get_miller_idx(Rz, qgrid[2]);

        for (ND_int ia = 0; ia < natom; ++ia)
        {
            const ELPH_float* tau_i = ctx->atomic_pos + 3 * ia;

            for (ND_int ja = 0; ja < natom; ++ja)
            {
                const ELPH_float* tau_j = ctx->atomic_pos + 3 * ja;

                ND_int degen_idx = ig * natom * natom + ia * natom + ja;
                const ND_int iws_vecs_degen = ctx->ws_degen[degen_idx];
                const double idegen_fac = 1.0 / ((double)iws_vecs_degen);

                // ---------------------------------------------------------
                // 1. Translational Constraints
                // ---------------------------------------------------------
                if (iws_vecs_degen > 0)
                {
                    for (ND_int alpha = 0; alpha < 3; ++alpha)
                    {
                        for (ND_int beta = 0; beta < 3; ++beta)
                        {
                            // Get index for Phi(ia, alpha, ja, beta)
                            ND_int idx_phi =
                                get_fc_index(ctx, iig, process_R, is_gamma, ia,
                                             alpha, ja, beta);

                            // Constraint index: T(ia, alpha, beta)
                            ND_int idx_constr = 9 * ia + 3 * alpha + beta;

                            // Apply (+1.0)
                            update_fc_Amat_entry(mode, x, y, idx_constr,
                                                 idx_phi, 1.0);
                        }
                    }
                }

                // Loop over WS vectors (Needed for Rotational and Huang)
                if (ctx->mode == ASR_ALL)
                {
                    for (ND_int ii = 0; ii < iws_vecs_degen; ++ii)
                    {
                        // Geometry Calculation
                        ELPH_float Rpt[3];
                        Rpt[0] = mx + ctx->ws_vecs[3 * iws_vec];
                        Rpt[1] = my + ctx->ws_vecs[3 * iws_vec + 1];
                        Rpt[2] = mz + ctx->ws_vecs[3 * iws_vec + 2];

                        ELPH_float tau_Rj[3], tau_Rij[3];
                        MatVec3f(ctx->lat_vecs, Rpt, false, tau_Rj);
                        tau_Rj[0] += tau_j[0];
                        tau_Rj[1] += tau_j[1];
                        tau_Rj[2] += tau_j[2];

                        tau_Rij[0] = tau_i[0] - tau_Rj[0];
                        tau_Rij[1] = tau_i[1] - tau_Rj[1];
                        tau_Rij[2] = tau_i[2] - tau_Rj[2];

                        // -----------------------------------------------------
                        // 2. Rotational Constraints
                        // -----------------------------------------------------
                        for (ND_int ic = 0; ic < 9; ++ic)
                        {
                            // Due to anti-symmetric property of rotational
                            // invariance, we only have 3 combinations of
                            // (beta,gamma) constraints. beta_gamma ->
                            // (beta,gamma): 0 -> (1,0); 1 -> (2,0); 2-> (2,1)
                            const ND_int alpha = ic / 3;
                            const ND_int beta_gamma = ic % 3;
                            const ND_int beta = MIN((beta_gamma + 1), 2);
                            const ND_int gamma = beta_gamma / 2;

                            ND_int idx_constr = off_rot + ia * 9 + ic;

                            double val1 = idegen_fac * tau_Rj[gamma];
                            double val2 = -1.0 * idegen_fac * tau_Rj[beta];

                            // Term 1: Phi(ia, alpha, ja, beta)
                            ND_int idx_phi1 =
                                get_fc_index(ctx, iig, process_R, is_gamma, ia,
                                             alpha, ja, beta);
                            update_fc_Amat_entry(mode, x, y, idx_constr,
                                                 idx_phi1, val1);

                            // Term 2: Phi(ia, alpha, ja, gamma)
                            ND_int idx_phi2 =
                                get_fc_index(ctx, iig, process_R, is_gamma, ia,
                                             alpha, ja, gamma);
                            update_fc_Amat_entry(mode, x, y, idx_constr,
                                                 idx_phi2, val2);
                        }

                        // -----------------------------------------------------
                        // 3. Huang Constraints
                        // -----------------------------------------------------
                        if (ctx->huang)
                        {
                            for (ND_int ic = 0; ic < 15; ++ic)
                            {
                                // Due to anti-symmtric properties
                                // (alpha ,beta) <-> (gamma, delta),
                                // we reduce contraints from 81 -> 36
                                // Further more, due to symmetry properties
                                // due to alpha <-> beta and gamma <-> delta,
                                // we further reduce from 36 -> 15 independent
                                // contraints
                                //
                                ND_int alpha, beta, gamma, delta;
                                get_huang_indices(ic, &alpha, &beta, &gamma,
                                                  &delta);
                                ND_int idx_constr = off_huang + ic;

                                double v1 = idegen_fac * tau_Rij[gamma] *
                                            tau_Rij[delta];
                                double v2 = -1.0 * idegen_fac * tau_Rij[alpha] *
                                            tau_Rij[beta];

                                // Term 1: Phi(ia, alpha, ja, beta)
                                ND_int idx_phi1 =
                                    get_fc_index(ctx, iig, process_R, is_gamma,
                                                 ia, alpha, ja, beta);
                                update_fc_Amat_entry(mode, x, y, idx_constr,
                                                     idx_phi1, v1);

                                // Term 2: Phi(ia, gamma, ja, delta)
                                ND_int idx_phi2 =
                                    get_fc_index(ctx, iig, process_R, is_gamma,
                                                 ia, gamma, ja, delta);
                                update_fc_Amat_entry(mode, x, y, idx_constr,
                                                     idx_phi2, v2);
                            }
                        }
                        ++iws_vec;
                    }
                }
            }
        }
    }
}
//
//
static inline ND_int get_fc_index(const struct ASR_Context* ctx,
                                  const ND_int iig, const bool process_R,
                                  const bool is_gamma, const ND_int ia,
                                  const ND_int alpha, const ND_int ja,
                                  const ND_int beta)
{
    /**
     * @brief Calculates the flat index in the Force Constant vector.
     * Handles the logic for -R mapping and Gamma point triangular storage.
     */
    ND_int r_im, r_jm;

    if (process_R)
    {
        // Standard Case: Phi(ia, alpha, ja, beta)
        r_im = 3 * ia + alpha;
        r_jm = 3 * ja + beta;
    }
    else
    {
        // -R Case: Corresponds to Transpose(Phi) -> Phi(ja, beta, ia, alpha)
        r_im = 3 * ja + beta;
        r_jm = 3 * ia + alpha;
    }

    // Gamma Point Symmetry Enforcment:
    // We only store the lower triangle (j <= i) at Gamma.
    // If the logical request asks for upper, we swap to lower.
    if (is_gamma && r_jm > r_im)
    {
        ND_int tmp = r_im;
        r_im = r_jm;
        r_jm = tmp;
    }

    return iig * ctx->nmodes * ctx->nmodes + r_im * ctx->nmodes + r_jm;
}

static inline void update_fc_Amat_entry(const int mode, const double* x,
                                        double* restrict y,
                                        const ND_int idx_constr,
                                        const ND_int idx_phi, const double val)
{
    /**
     * @brief Applies the sparse matrix value to the vectors based on mode.
     * Mode 0 (A^T): y[phi] += val * x[constr]
     * Mode 1 (A):   y[constr] += val * x[phi]
     */
    if (mode == 0)
    {
        y[idx_phi] += val * x[idx_constr];
    }
    else
    {
        y[idx_constr] += val * x[idx_phi];
    }
}
//
static void get_huang_indices(const ND_int idx, ND_int* restrict a,
                              ND_int* restrict b, ND_int* restrict c,
                              ND_int* restrict d)
{
    // maps [0,14] indices to 15 independent coordinates rank-4 tensor
    // with has properties of huang invariance matrix
    // clang-format off
    switch (idx)
    {
        /* Pair (0,0) combined with others */
        case 0:  *a=0; *b=0; *c=0; *d=1; break;
        case 1:  *a=0; *b=0; *c=0; *d=2; break;
        case 2:  *a=0; *b=0; *c=1; *d=1; break;
        case 3:  *a=0; *b=0; *c=1; *d=2; break;
        case 4:  *a=0; *b=0; *c=2; *d=2; break;

        /* Pair (0,1) combined with others */
        case 5:  *a=0; *b=1; *c=0; *d=2; break;
        case 6:  *a=0; *b=1; *c=1; *d=1; break;
        case 7:  *a=0; *b=1; *c=1; *d=2; break;
        case 8:  *a=0; *b=1; *c=2; *d=2; break;

        /* Pair (0,2) combined with others */
        case 9:  *a=0; *b=2; *c=1; *d=1; break;
        case 10: *a=0; *b=2; *c=1; *d=2; break;
        case 11: *a=0; *b=2; *c=2; *d=2; break;

        /* Pair (1,1) combined with others */
        case 12: *a=1; *b=1; *c=1; *d=2; break;
        case 13: *a=1; *b=1; *c=2; *d=2; break;

        /* Pair (1,2) combined with (2,2) */
        case 14: *a=1; *b=2; *c=2; *d=2; break;

        // you should never reach this
        default: *a=0; *b=0; *c=0; *d=0; break;
    }
    // clang-format on
}

// clang-format off
// For refecence, below commented is using standard QR
/*
void apply_acoustic_sum_rule_born_charges(enum asr_kind mode, ELPH_float* Zborn,
                                          const ND_int natom)
{
    if (mode == ASR_NONE || !Zborn || !natom)
    {
        return;
    }

    if (mode == ASR_SIMPLE)
    {
        for (ND_int i = 0; i < 3; i++)
        {
            for (ND_int j = 0; j < 3; j++)
            {
                ELPH_float sum = 0.0;

                for (ND_int ia = 0; ia < natom; ia++)
                {
                    sum += Zborn[ia * 9 + i * 3 + j];
                }

                ELPH_float avg = sum / ((ELPH_float)natom);

                for (ND_int ia = 0; ia < natom; ia++)
                {
                    Zborn[ia * 9 + i * 3 + j] -= avg;
                }
            }
        }
    }
    else
    {
        // allocate constraint matrix
        ELPH_float* Amat = calloc(81 * natom, sizeof(*Amat));
        CHECK_ALLOC(Amat);

        // let the compilers remove this.
        for (ND_int i = 0; i < (81 * natom); ++i)
        {
            Amat[i] = 0.0;
        }

        for (ND_int ic = 0; ic < 9; ++ic)
        {
            ND_int alpha = ic / 3;
            ND_int beta = ic % 3;
            for (ND_int ia = 0; ia < natom; ++ia)
            {
                Amat[ic * 9 * natom + ia * 9 + 3 * alpha + beta] = 1.0;
            }
        }

        // Now we need to find |Z_input-Z_ideal|_min in such a way
        // that A*Z_input = 0
        ND_int err_info =
            orthogonal_projection(9, 9 * natom, 9 * natom, Amat, Zborn, 1e-5);
        if (err_info)
        {
            error_msg("Orthogonal projection of born charges failed.");
        }

        free(Amat);
    }
}

void apply_acoustic_sum_rule_fc(enum asr_kind mode, const ND_int* qgrid,
                                const ND_int natom, ELPH_cmplx* frc,
                                const ELPH_float* atomic_pos,
                                const ELPH_float* lat_vecs,
                                const ND_int* ws_vecs, const ND_int n_ws_vecs,
                                const ND_int* ws_degen)

{
    // C. Lin et al.: npj Comput. Mater. 8, 236 (2022)
    if (mode == ASR_NONE)
    {
        return;
    }

    bool huang = false;
    if (mode == ASR_ALL_HUANG)
    {
        mode = ASR_ALL;
        huang = true;
    }

    const ND_int ngrid = qgrid[0] * qgrid[1] * qgrid[2];
    const ND_int nmodes = 3 * natom;
    const ND_int Gridyz = qgrid[1] * qgrid[2];

    // Simple ASR
    if (mode == ASR_SIMPLE)
    {
        // (ig, na, alpha, nb, beta)
        for (ND_int ia = 0; ia < natom; ++ia)
        {
            for (ND_int alpha = 0; alpha < 3; ++alpha)
            {
                for (ND_int beta = 0; beta < 3; ++beta)
                {
                    ELPH_cmplx sum = 0.0;
                    ND_int base_na_i_j = ia * 9 * natom + alpha * nmodes + beta;

                    for (ND_int ig = 0; ig < ngrid; ++ig)
                    {
                        ND_int grid_off = ig * nmodes * nmodes;
                        for (ND_int ib = 0; ib < natom; ++ib)
                        {
                            ND_int idx = grid_off + base_na_i_j + ib * 3;
                            sum += frc[idx];
                        }
                    }
                    // remove self interaction
                    ND_int self_idx = base_na_i_j + ia * 3;
                    frc[self_idx] -= sum;
                }
            }
        }
        return;
    }

    // Crystal or all or all_huang
    // Orthogonal projection method
    //
    // Phi_{ai, bj}(R)= Phi_{bj,ai}(-R)
    // if R = -R, so we skip all -R,
    // for R = 0, the force constant matrix is symmetric
    //
    // In should be noted that R = -R => R =0, but in practice,
    // for boundary FFT points, R \eqvil to -R mod grid_size,
    // but this is an artifact and at this R point, we should not
    // symmtrize the force constant matrix.
    //
    // First find number of grid points
    //
    ND_int* grid_map = calloc(ngrid, sizeof(*grid_map));
    CHECK_ALLOC(grid_map);

    bool* do_grid = malloc(ngrid * sizeof(*do_grid));
    CHECK_ALLOC(do_grid);
    // if do_grid[ig] = false, then ig is -R and is mapped
    // to some other R

    for (ND_int ig = 0; ig < ngrid; ++ig)
    {
        do_grid[ig] = true;
    }
    //
    ND_int ngrid_independent = 0;
    for (ND_int ig = 0; ig < ngrid; ++ig)
    {
        if (!do_grid[ig])
        {
            continue;
        }

        ND_int Rx = ig / Gridyz;
        ND_int Ry = (ig % Gridyz) / qgrid[2];
        ND_int Rz = (ig % Gridyz) % qgrid[2];

        // Calculate -R index (mig)
        ND_int mRx = (qgrid[0] - Rx) % qgrid[0];
        ND_int mRy = (qgrid[1] - Ry) % qgrid[1];
        ND_int mRz = (qgrid[2] - Rz) % qgrid[2];
        ND_int mig = mRz + mRy * qgrid[2] + mRx * Gridyz;
        //
        grid_map[ig] = ngrid_independent;
        if (ig != mig)
        {
            do_grid[mig] = false;
            grid_map[mig] = ngrid_independent;
        }
        ++ngrid_independent;
    }
    //
    const ND_int frc_size = ngrid_independent * nmodes * nmodes;
    ELPH_float* frc_real = malloc(frc_size * sizeof(*frc_real));
    CHECK_ALLOC(frc_real);

    // Make the force constant matrix real
    ND_int iig_tmp = 0;
    for (ND_int ig = 0; ig < ngrid; ++ig)
    {
        if (!do_grid[ig])
        {
            continue;
        }
        for (ND_int i = 0; i < nmodes * nmodes; ++i)
        {
            frc_real[i + iig_tmp * nmodes * nmodes] =
                creal(frc[i + ig * nmodes * nmodes]);
        }
        ++iig_tmp;
    }

    ND_int nconstraints = 9 * natom;  // translational
    if (mode == ASR_ALL)
    {
        nconstraints += 9 * natom;  // rotational
        if (huang)
        {
            nconstraints += 15;  // huang
        }
    }

    // Allocate the contraint matrix
    // A(ngrid_independent, ngrid, na , 3, nb, 3)
    ELPH_float* Amat = calloc(nconstraints * frc_size, sizeof(*Amat));
    CHECK_ALLOC(Amat);

    // most compilers will remove this loop. But let's stick to standard.
    for (ND_int i = 0; i < (nconstraints * frc_size); ++i)
    {
        Amat[i] = 0.0;
    }

    // 1) translational
    for (ND_int ig = 0; ig < ngrid; ++ig)
    {
        ND_int iig = grid_map[ig];
        for (ND_int ic = 0; ic < 9 * natom; ++ic)
        {
            //
            const ND_int ia = ic / 9;
            const ND_int alpha = (ic % 9) / 3;
            const ND_int beta = (ic % 9) % 3;
            // (ic, ngrid, na , 3, nb, 3)
            for (ND_int ja = 0; ja < natom; ++ja)
            {
                // For R, set A[ic,igg,ia,alpha,ja,beta] += 1
                ELPH_float* Amat_tmp = Amat + ic * frc_size + ia * 9 * natom +
                                       alpha * nmodes + beta +
                                       iig * nmodes * nmodes + 3 * ja;
                if (!do_grid[ig])
                {
                    // this is -R, so we need to wrap
                    // For -R, set A[ic,igg,ja,beta,ia,alpha] += 1
                    Amat_tmp = Amat + ic * frc_size + ja * 9 * natom +
                               beta * nmodes + alpha + iig * nmodes * nmodes +
                               3 * ia;
                }
                *Amat_tmp += 1.0;
            }
        }
    }

    // 2) rotational sum rules in case requested
    // 3) Huang invariances
    if (mode == ASR_ALL)
    {
        ND_int iws_vec = 0;
        for (ND_int ig = 0; ig < ngrid; ++ig)
        {
            ND_int iig = grid_map[ig];

            ND_int Rx = ig / Gridyz;
            ND_int Ry = (ig % Gridyz) / qgrid[2];
            ND_int Rz = (ig % Gridyz) % qgrid[2];
            //
            Rx = get_miller_idx(Rx, qgrid[0]);
            Ry = get_miller_idx(Ry, qgrid[1]);
            Rz = get_miller_idx(Rz, qgrid[2]);

            for (ND_int ia = 0; ia < natom; ++ia)
            {
                const ELPH_float* tau_i = atomic_pos + 3 * ia;
                for (ND_int ja = 0; ja < natom; ++ja)
                {
                    const ELPH_float* tau_j = atomic_pos + 3 * ja;
                    const ND_int iws_vecs_degen =
                        ws_degen[ig * natom * natom + ia * natom + ja];
                    const ELPH_float idegen_fac =
                        1.0 / ((ELPH_float)iws_vecs_degen);
                    // (ic, ngrid, na , 3, nb, 3)
                    ELPH_float* Amat_tmp =
                        Amat + iig * nmodes * nmodes + ja * 3 + ia * 3 * nmodes;
                    if (!do_grid[ig])
                    {
                        // this is -R, so we need to wrap
                        Amat_tmp = Amat + iig * nmodes * nmodes + ia * 3 +
                                   ja * 3 * nmodes;
                    }

                    for (ND_int ii = 0; ii < iws_vecs_degen; ++ii)
                    {
                        if (iws_vec >= n_ws_vecs)
                        {
                            error_msg("Wigner seitz vectors Out of bound.");
                        }
                        ELPH_float Rpt[3];
                        Rpt[0] = Rx + ws_vecs[3 * iws_vec];
                        Rpt[1] = Ry + ws_vecs[3 * iws_vec + 1];
                        Rpt[2] = Rz + ws_vecs[3 * iws_vec + 2];
                        // convert to cart and add atomic postion
                        ELPH_float tau_Rj[3], tau_Rij[3];
                        MatVec3f(lat_vecs, Rpt, false, tau_Rj);
                        //
                        tau_Rj[0] += tau_j[0];
                        tau_Rj[1] += tau_j[1];
                        tau_Rj[2] += tau_j[2];  //
                        //
                        tau_Rij[0] = tau_i[0] - tau_Rj[0];
                        tau_Rij[1] = tau_i[1] - tau_Rj[1];
                        tau_Rij[2] = tau_i[2] - tau_Rj[2];
                        // (summ over ja, R)
                        // 2) Rotational sum rule
                        for (ND_int ic = 0; ic < 9; ++ic)
                        {
                            const ND_int alpha = ic / 3;
                            const ND_int beta_gamma = ic % 3;
                            // Due to anti-symmetric property of rotational
                            // invariance, we only have 3 combinations of
                            // (beta,gamma) constraints. beta_gamma ->
                            // (beta,gamma): 0 -> (1,0); 1 -> (2,0); 2-> (2,1)
                            const ND_int beta = MIN((beta_gamma + 1), 2);
                            const ND_int gamma = beta_gamma / 2;
                            //
                            // +9natom is to append the rotation after 9*natom
                            // translational rules
                            ELPH_float* Amat_tmp_ic =
                                Amat_tmp + (ic + ia * 9 + 9 * natom) * frc_size;

                            // (j, beta, gamma)
                            if (!do_grid[ig])
                            {
                                // A[ic,igg,ja, beta ia, alpha]
                                Amat_tmp_ic[alpha + beta * nmodes] +=
                                    (idegen_fac * tau_Rj[gamma]);
                                // A[ic,igg,ja, gamma ia, alpha]
                                Amat_tmp_ic[alpha + gamma * nmodes] -=
                                    (idegen_fac * tau_Rj[beta]);
                            }
                            else
                            {
                                // A[ic,igg,ia, alpha ja, beta ]
                                Amat_tmp_ic[alpha * nmodes + beta] +=
                                    (idegen_fac * tau_Rj[gamma]);
                                // A[ic,igg,ia, alpha ja, gamma]
                                Amat_tmp_ic[alpha * nmodes + gamma] -=
                                    (idegen_fac * tau_Rj[beta]);
                            }
                        }
                        //
                        // 3) Huang invariances
                        if (huang)
                        {
                            for (ND_int ic = 0; ic < 15; ++ic)
                            {
                                ELPH_float* Amat_tmp_ic =
                                    Amat_tmp + (ic + 18 * natom) * frc_size;
                                // Due to anti-symmtric properties
                                // (alpha ,beta) <-> (gamma, delta),
                                // we reduce contraints from 81 -> 36
                                // Further more, due to symmetry properties
                                // due to alpha <-> beta and gamma <-> delta,
                                // we further reduce from 36 -> 15 independent
                                // contraints
                                ND_int alpha, beta, gamma, delta;
                                get_huang_indices(ic, &alpha, &beta, &gamma,
                                                  &delta);
                                //
                                //
                                if (!do_grid[ig])
                                {
                                    // for -R
                                    //  A[ic,igg,ja, beta ia, alpha]
                                    Amat_tmp_ic[alpha + beta * nmodes] +=
                                        (idegen_fac * tau_Rij[gamma] *
                                         tau_Rij[delta]);
                                    // A[ic,igg,ja, delta ia, gamma]
                                    Amat_tmp_ic[gamma + delta * nmodes] -=
                                        (idegen_fac * tau_Rij[alpha] *
                                         tau_Rij[beta]);
                                }
                                else
                                {
                                    // for R
                                    // A[ic,igg,ia, alpha ja, beta]
                                    Amat_tmp_ic[alpha * nmodes + beta] +=
                                        (idegen_fac * tau_Rij[gamma] *
                                         tau_Rij[delta]);
                                    // A[ic,igg,ia, gamma ja, delta]
                                    Amat_tmp_ic[gamma * nmodes + delta] -=
                                        (idegen_fac * tau_Rij[alpha] *
                                         tau_Rij[beta]);
                                }
                            }
                        }
                        ++iws_vec;
                    }
                }
            }
        }
    }

    // To, Symmetrize R = 0 point. we remove upper triangular part
    // in A and frc_real
    // 1st A mat
    for (ND_int ic = 0; ic < nconstraints; ++ic)
    {
        ELPH_float* Amat_tmp = Amat + ic * frc_size;
        for (ND_int i = 0; i < nmodes; ++i)
        {
            for (ND_int j = 0; j < i; ++j)
            {
                // A[ic,0,ia,alpha,ja,beta] += A[ic,0,ja,beta,ia,alpha]
                Amat_tmp[i * nmodes + j] += Amat_tmp[j * nmodes + i];
                Amat_tmp[j * nmodes + i] = 0.0;
            }
        }
    }
    // now frc
    for (ND_int i = 0; i < nmodes; ++i)
    {
        for (ND_int j = 0; j < i; ++j)
        {
            frc_real[j * nmodes + i] = 0.0;
        }
    }

    // Now we need to find |\Phi_input-\Phi_ideal|_min in such a way
    // that A\Phi_input = 0
    ND_int err_info = orthogonal_projection(nconstraints, frc_size, frc_size,
                                            Amat, frc_real, 1e-5);
    if (err_info)
    {
        error_msg(
            "Orthogonal projection of force constant Constraint matrix "
            "failed.");
    }

    // Now recontruct the full force constant matrix
    for (ND_int ig = 0; ig < ngrid; ++ig)
    {
        if (!do_grid[ig])
        {
            continue;
        }
        ND_int iig = grid_map[ig];

        ND_int Rx = ig / Gridyz;
        ND_int Ry = (ig % Gridyz) / qgrid[2];
        ND_int Rz = (ig % Gridyz) % qgrid[2];

        // Calculate -R index (mig)
        ND_int mRx = (qgrid[0] - Rx) % qgrid[0];
        ND_int mRy = (qgrid[1] - Ry) % qgrid[1];
        ND_int mRz = (qgrid[2] - Rz) % qgrid[2];
        ND_int mig = mRz + mRy * qgrid[2] + mRx * Gridyz;
        //
        for (ND_int i = 0; i < nmodes; ++i)
        {
            for (ND_int j = 0; j < nmodes; ++j)
            {
                frc[ig * nmodes * nmodes + i * nmodes + j] =
                    frc_real[iig * nmodes * nmodes + i * nmodes + j];
                if (ig != mig)
                {
                    frc[mig * nmodes * nmodes + j * nmodes + i] =
                        frc_real[iig * nmodes * nmodes + i * nmodes + j];
                }
            }
        }
    }

    // Symmetrization at Gamma point
    for (ND_int i = 0; i < nmodes; ++i)
    {
        for (ND_int j = 0; j < i; ++j)
        {
            frc[j * nmodes + i] = frc[i * nmodes + j];
        }
    }

    free(do_grid);
    free(grid_map);
    free(Amat);
    free(frc_real);
}
*/
// clang-format on
