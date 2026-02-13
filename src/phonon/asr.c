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

// static functions
static void get_huang_indices(const ND_int idx, ND_int* restrict a,
                              ND_int* restrict b, ND_int* restrict c,
                              ND_int* restrict d);

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
         * CRYSTAL RULE
         * ==============================
         */
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

// static helpers
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
