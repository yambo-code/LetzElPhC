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
enum asr_kind asr_kind_from_string(const char* str)
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

    return ASR_NONE;
}

void apply_acoustic_sum_rule_born_charges(enum asr_kind mode, ELPH_float* Zborn,
                                          const ND_int nat)
{
    if (mode == ASR_NONE || !Zborn || !nat)
    {
        return;
    }

    ND_int total_elements = 3 * 3 * nat;

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

                for (ND_int na = 0; na < nat; na++)
                {
                    sum += Zborn[na * 9 + i * 3 + j];
                }

                ELPH_float avg = sum / nat;

                for (ND_int na = 0; na < nat; na++)
                {
                    Zborn[na * 9 + i * 3 + j] -= avg;
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
        ND_int zeu_dim = total_elements;
        ND_int zeu_u_dim = 18 * zeu_dim;

        ELPH_float* zeu_u = calloc(zeu_u_dim, sizeof(ELPH_float));
        CHECK_ALLOC(zeu_u);

        ELPH_float* zeu_new = malloc(zeu_dim * sizeof(ELPH_float));
        CHECK_ALLOC(zeu_new);

        memcpy(zeu_new, Zborn, zeu_dim * sizeof(ELPH_float));

        ND_int p = 0;
        /* Initialize translational sum rules in zeu_u */
        for (ND_int i = 0; i < 3; i++)
        {
            for (ND_int j = 0; j < 3; j++)
            {
                for (ND_int iat = 0; iat < nat; iat++)
                {
                    ND_int idx = p * zeu_dim + (iat * 9 + i * 3 + j);
                    zeu_u[idx] = 1.0;
                }
                p++;
            }
        }

        /* Gram-Schmidt Orthonormalization */
        ELPH_float* zeu_w = calloc(zeu_dim, sizeof(ELPH_float));
        CHECK_ALLOC(zeu_w);
        //
        ELPH_float* zeu_x = calloc(zeu_dim, sizeof(ELPH_float));
        CHECK_ALLOC(zeu_x);
        //
        ELPH_float* tempZeu = calloc(zeu_dim, sizeof(ELPH_float));
        CHECK_ALLOC(tempZeu);

        /* zeu_less stores indices, so it must be ND_int */
        ND_int* zeu_less = calloc(18, sizeof(ND_int));
        CHECK_ALLOC(zeu_less);

        ELPH_float scalar;
        ND_int nzeu_less = 0;
        ND_int r;

        for (ND_int k = 0; k < p; k++)
        {
            memcpy(zeu_w, zeu_u + k * zeu_dim,
                   total_elements * sizeof(ELPH_float));
            memcpy(zeu_x, zeu_u + k * zeu_dim,
                   total_elements * sizeof(ELPH_float));

            for (ND_int q = 0; q < k; q++)
            {
                r = 1;
                for (ND_int iZeu_less = 0; iZeu_less < nzeu_less; iZeu_less++)
                {
                    if (zeu_less[iZeu_less] == q)
                    {
                        r = 0;
                    }
                }

                if (r != 0)
                {
                    memcpy(tempZeu, zeu_u + q * zeu_dim,
                           total_elements * sizeof(ELPH_float));

                    /* Dot product (zeu_x . tempZeu) */
                    scalar = 0.0;
                    for (ND_int z = 0; z < total_elements; z++)
                    {
                        scalar += zeu_x[z] * tempZeu[z];
                    }

                    for (ND_int x = 0; x < total_elements; x++)
                    {
                        zeu_w[x] -= scalar * tempZeu[x];
                    }
                }
            }

            /* Dot product (zeu_w . zeu_w) */
            ELPH_float norm2 = 0.0;
            for (ND_int z = 0; z < total_elements; z++)
            {
                norm2 += zeu_w[z] * zeu_w[z];
            }

            if (sqrt(norm2) > ELPH_EPS)
            {
                ELPH_float inv_sqrt_norm = 1.0 / sqrt(norm2);
                for (ND_int x = 0; x < total_elements; x++)
                {
                    zeu_u[k * zeu_dim + x] = zeu_w[x] * inv_sqrt_norm;
                }
            }
            else
            {
                zeu_less[nzeu_less] = k;
                nzeu_less++;
            }
        }

        /* Projection */
        for (ND_int ii = 0; ii < zeu_dim; ++ii)
        {
            zeu_w[ii] = 0.0;
        }

        for (ND_int k = 0; k < p; k++)
        {
            r = 1;
            for (ND_int izeu_less = 0; izeu_less < nzeu_less; izeu_less++)
            {
                if (zeu_less[izeu_less] == k)
                {
                    r = 0;
                }
            }

            if (r != 0)
            {
                memcpy(zeu_x, zeu_u + k * zeu_dim,
                       total_elements * sizeof(ELPH_float));

                /* Dot product (zeu_x . zeu_new) */
                scalar = 0.0;
                for (ND_int z = 0; z < total_elements; z++)
                {
                    scalar += zeu_x[z] * zeu_new[z];
                }

                for (ND_int x = 0; x < total_elements; x++)
                {
                    zeu_w[x] += scalar * zeu_u[k * zeu_dim + x];
                }
            }
        }

        for (ND_int x = 0; x < total_elements; x++)
        {
            Zborn[x] = zeu_new[x] - zeu_w[x];
        }

        free(zeu_u);
        free(zeu_new);
        free(zeu_w);
        free(zeu_x);
        free(tempZeu);
        free(zeu_less);
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

    ND_int ngrid = qgrid[0] * qgrid[1] * qgrid[2];
    ND_int nmodes = 3 * natom;

    // Orthogonal projection method
    const ND_int frc_size = ngrid * nmodes * nmodes;
    ELPH_float* frc_real = malloc(frc_size * sizeof(*frc_real));
    CHECK_ALLOC(frc_real);

    // Make the force constant matrix real
    for (ND_int i = 0; i < frc_size; ++i)
    {
        frc_real[i] = creal(frc[i]);
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

    ELPH_float* Amat = calloc(nconstraints * frc_size, sizeof(*Amat));
    CHECK_ALLOC(Amat);

    // most compilers will remove this loop. But let's stick to standard.
    for (ND_int i = 0; i < (nconstraints * frc_size); ++i)
    {
        Amat[i] = 0.0;
    }

    // 1) translational
    for (ND_int ic = 0; ic < 9 * natom; ++ic)
    {
        //
        const ND_int ia = ic / 9;
        const ND_int alpha = (ic % 9) / 3;
        const ND_int beta = (ic % 9) % 3;
        // (ic, ngrid, na , 3, nb, 3)
        ELPH_float* Amat_tmp = Amat + ic * ngrid * nmodes * nmodes +
                               ia * 9 * natom + alpha * nmodes + beta;
        // R,j'
        for (ND_int ig = 0; ig < ngrid; ++ig)
        {
            for (ND_int ja = 0; ja < natom; ++ja)
            {
                Amat_tmp[ig * nmodes * nmodes + ja * 3] = 1.0;
            }
        }
    }

    ND_int Gridyz = qgrid[1] * qgrid[2];
    // 2) rotational sum rules in case requested
    // 3) Huang invariances
    if (mode == ASR_ALL)
    {
        ND_int iws_vec = 0;
        for (ND_int ig = 0; ig < ngrid; ++ig)
        {
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
                        Amat + ig * nmodes * nmodes + ja * 3 + ia * 3 * nmodes;

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
                                Amat_tmp + alpha * nmodes +
                                (ic + ia * 9 + 9 * natom) * frc_size;

                            // (j, beta, gamma)
                            Amat_tmp_ic[beta] += (idegen_fac * tau_Rj[gamma]);
                            Amat_tmp_ic[gamma] -= (idegen_fac * tau_Rj[beta]);
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
                                Amat_tmp_ic[alpha * nmodes + beta] +=
                                    (idegen_fac * tau_Rij[gamma] *
                                     tau_Rij[delta]);
                                Amat_tmp_ic[gamma * nmodes + delta] -=
                                    (idegen_fac * tau_Rij[alpha] *
                                     tau_Rij[beta]);
                            }
                        }
                        ++iws_vec;
                    }
                }
            }
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

    // copy back the force constants
    for (ND_int i = 0; i < frc_size; ++i)
    {
        frc[i] = frc_real[i];
    }

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
