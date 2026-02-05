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
                                const ND_int nat, ELPH_cmplx* frc,
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

    const ND_int nr1 = qgrid[0];
    const ND_int nr2 = qgrid[1];
    const ND_int nr3 = qgrid[2];
    ND_int n_grid = nr1 * nr2 * nr3;
    ND_int dim_nb = nat * 3;
    ND_int dim_i = 3 * dim_nb;
    ND_int dim_na = nat * dim_i;
    ND_int dim_total = n_grid * dim_na;

    /*    CASE 1: Simple ASR */
    /*    Sum_nb Phi(na, nb) = 0 */
    if (mode == ASR_SIMPLE)
    {
        for (ND_int na = 0; na < nat; na++)
        {
            for (ND_int i = 0; i < 3; i++)
            {
                for (ND_int j = 0; j < 3; j++)
                {
                    ELPH_cmplx sum = 0.0;
                    ND_int base_na_i_j = na * dim_i + i * dim_nb + j;

                    for (ND_int g = 0; g < n_grid; g++)
                    {
                        ND_int grid_off = g * dim_na;
                        for (ND_int nb = 0; nb < nat; nb++)
                        {
                            ND_int idx = grid_off + base_na_i_j + nb * 3;
                            sum += frc[idx];
                        }
                    }
                    // remove self interaction
                    ND_int self_idx = base_na_i_j + na * 3;
                    frc[self_idx] -= sum;
                }
            }
        }
        return;
    }

    /* CASE 2: Projection Methods, cyrstal and all */
    ND_int n_trans = 9;
    ND_int n_rot = 0;
    ND_int n_huang = 0;

    if (mode == ASR_ALL)
    {
        n_rot = 9;
        if (huang)
        {
            n_huang = 15;
        }
    }

    ND_int total_vecs = (n_trans + n_rot) * nat + n_huang;

    ELPH_cmplx** u = malloc(total_vecs * sizeof(ELPH_cmplx*));
    CHECK_ALLOC(u);

    for (ND_int k = 0; k < total_vecs; k++)
    {
        u[k] = malloc(dim_total * sizeof(ELPH_cmplx));
        CHECK_ALLOC(u[k]);
        for (ND_int z = 0; z < dim_total; z++)
        {
            u[k][z] = 0.0;
        }
    }

    // Stream index for ws_vecs
    ND_int ws_stream_idx = 0;

    ND_int p_rot_start = 9 * nat;
    ND_int p_huang_start = p_rot_start + 9 * nat;

    ND_int grid_yz = nr2 * nr3;

    // Loop order must match build_wigner_seitz_vectors to read stream correctly
    for (ND_int g = 0; g < n_grid; g++)
    {
        ND_int rx_raw = g / grid_yz;
        ND_int ry_raw = (g % grid_yz) / nr3;
        ND_int rz_raw = (g % grid_yz) % nr3;

        ND_int rx = get_miller_idx(rx_raw, nr1);
        ND_int ry = get_miller_idx(ry_raw, nr2);
        ND_int rz = get_miller_idx(rz_raw, nr3);

        // Calculate Grid R vector in Cartesian using lat_vecs
        ELPH_float r_vec_cart[3];
        for (int d = 0; d < 3; d++)
        {
            r_vec_cart[d] = lat_vecs[d * 3 + 0] * rx +
                            lat_vecs[d * 3 + 1] * ry + lat_vecs[d * 3 + 2] * rz;
        }

        for (ND_int na = 0; na < nat; na++)
        {
            for (ND_int nb = 0; nb < nat; nb++)
            {
                ND_int blk_offset = g * dim_na + na * dim_i + nb * 3;

                // 1. Translations (Crystal ASR)
                for (ND_int i = 0; i < 3; i++)
                {
                    for (ND_int j = 0; j < 3; j++)
                    {
                        ND_int p_idx = (na * 3 + i) * 3 + j;
                        u[p_idx][blk_offset + i * dim_nb + j] = 1.0;
                    }
                }

                if (mode == ASR_ALL)
                {
                    // Degeneracy Index
                    ND_int degen_idx = g * nat * nat + na * nat + nb;
                    ND_int neq = ws_degen[degen_idx];

                    ELPH_float sum_r[3] = {0.0, 0.0, 0.0};
                    ELPH_float sum_rr[3][3] = {{0}};

                    // Atom position difference
                    ELPH_float tau_diff[3];
                    for (int d = 0; d < 3; d++)
                    {
                        tau_diff[d] =
                            atomic_pos[nb * 3 + d] - atomic_pos[na * 3 + d];
                    }

                    for (ND_int q = 0; q < neq; q++)
                    {
                        // Read Integer Wigner-Seitz Vector
                        if (ws_stream_idx >= n_ws_vecs)
                        {
                            error_msg("Wigner seitz vectors out of bound.");
                        }
                        ND_int t_int[3];
                        t_int[0] = ws_vecs[3 * ws_stream_idx];
                        t_int[1] = ws_vecs[3 * ws_stream_idx + 1];
                        t_int[2] = ws_vecs[3 * ws_stream_idx + 2];

                        ws_stream_idx++;

                        // Convert Integer WS Vector to Cartesian
                        ELPH_float t_vec_cart[3];
                        for (int d = 0; d < 3; d++)
                        {
                            t_vec_cart[d] = lat_vecs[d * 3 + 0] * t_int[0] +
                                            lat_vecs[d * 3 + 1] * t_int[1] +
                                            lat_vecs[d * 3 + 2] * t_int[2];
                        }

                        // r_eff = R_grid + T_ws + (tau_nb - tau_na)
                        ELPH_float r_eff[3];
                        for (int d = 0; d < 3; d++)
                        {
                            r_eff[d] =
                                r_vec_cart[d] + t_vec_cart[d] + tau_diff[d];

                            sum_r[d] += r_eff[d];
                        }

                        if (huang)
                        {
                            for (int a = 0; a < 3; a++)
                            {
                                for (int b = 0; b < 3; b++)
                                {
                                    sum_rr[a][b] += r_eff[a] * r_eff[b];
                                }
                            }
                        }
                    }

                    if (neq > 0)
                    {
                        ELPH_float inv_neq = 1.0 / (ELPH_float)neq;
                        for (int d = 0; d < 3; d++)
                        {
                            sum_r[d] *= inv_neq;
                        }
                        if (huang)
                        {
                            for (int a = 0; a < 3; a++)
                            {
                                for (int b = 0; b < 3; b++)
                                {
                                    sum_rr[a][b] *= inv_neq;
                                }
                            }
                        }
                    }

                    // 2. Rotational Constraints
                    for (ND_int ax = 0; ax < 3; ax++)
                    {
                        ND_int ax1 = (ax + 1) % 3;
                        ND_int ax2 = (ax + 2) % 3;

                        for (ND_int i = 0; i < 3; i++)
                        {
                            ND_int p_rot =
                                p_rot_start + (ax * 3 + i) * nat + na;
                            u[p_rot][blk_offset + i * dim_nb + ax1] -=
                                (sum_r[ax2]);
                            u[p_rot][blk_offset + i * dim_nb + ax2] +=
                                (sum_r[ax1]);
                        }
                    }

                    // 3. Huang Constraints (Explicit expansion, no macros)
                    if (huang)
                    {
                        // 1. yx
                        u[p_huang_start + 0][blk_offset + 0 * dim_nb + 0] -=
                            (sum_rr[0][1]);

                        u[p_huang_start + 0][blk_offset + 0 * dim_nb + 1] +=
                            (sum_rr[0][0]);

                        // 2. zx
                        u[p_huang_start + 1][blk_offset + 0 * dim_nb + 0] -=
                            (sum_rr[0][2]);
                        u[p_huang_start + 1][blk_offset + 0 * dim_nb + 2] +=
                            (sum_rr[0][0]);

                        // 3. xx-yy
                        u[p_huang_start + 2][blk_offset + 0 * dim_nb + 0] -=
                            (sum_rr[1][1]);
                        u[p_huang_start + 2][blk_offset + 1 * dim_nb + 1] +=
                            (sum_rr[0][0]);

                        // 4. yz (pair 1)
                        u[p_huang_start + 3][blk_offset + 0 * dim_nb + 0] -=
                            (sum_rr[1][2]);
                        u[p_huang_start + 3][blk_offset + 1 * dim_nb + 2] +=
                            (sum_rr[0][0]);

                        // 5. xx-zz
                        u[p_huang_start + 4][blk_offset + 0 * dim_nb + 0] -=
                            (sum_rr[2][2]);
                        u[p_huang_start + 4][blk_offset + 2 * dim_nb + 2] +=
                            (sum_rr[0][0]);

                        // 6. xy
                        u[p_huang_start + 5][blk_offset + 0 * dim_nb + 1] -=
                            (sum_rr[0][2]);
                        u[p_huang_start + 5][blk_offset + 0 * dim_nb + 2] +=
                            (sum_rr[0][1]);

                        // 7. xy (pair 2)
                        u[p_huang_start + 6][blk_offset + 0 * dim_nb + 1] -=
                            (sum_rr[1][1]);
                        u[p_huang_start + 6][blk_offset + 1 * dim_nb + 1] +=
                            (sum_rr[0][1]);

                        // 8. yz (pair 2)
                        u[p_huang_start + 7][blk_offset + 0 * dim_nb + 1] -=
                            (sum_rr[1][2]);
                        u[p_huang_start + 7][blk_offset + 1 * dim_nb + 2] +=
                            (sum_rr[0][1]);

                        // 9. zz
                        u[p_huang_start + 8][blk_offset + 0 * dim_nb + 1] -=
                            (sum_rr[2][2]);
                        u[p_huang_start + 8][blk_offset + 2 * dim_nb + 2] +=
                            (sum_rr[0][1]);

                        // 10. yy
                        u[p_huang_start + 9][blk_offset + 0 * dim_nb + 2] -=
                            (sum_rr[1][1]);
                        u[p_huang_start + 9][blk_offset + 1 * dim_nb + 1] +=
                            (sum_rr[0][2]);

                        // 11. yz (pair 3)
                        u[p_huang_start + 10][blk_offset + 0 * dim_nb + 2] -=
                            (sum_rr[1][2]);
                        u[p_huang_start + 10][blk_offset + 1 * dim_nb + 2] +=
                            (sum_rr[0][2]);

                        // 12. zz (pair 2)
                        u[p_huang_start + 11][blk_offset + 0 * dim_nb + 2] -=
                            (sum_rr[2][2]);
                        u[p_huang_start + 11][blk_offset + 2 * dim_nb + 2] +=
                            (sum_rr[0][2]);

                        // 13. zy
                        u[p_huang_start + 12][blk_offset + 1 * dim_nb + 1] -=
                            (sum_rr[1][2]);
                        u[p_huang_start + 12][blk_offset + 1 * dim_nb + 2] +=
                            (sum_rr[1][1]);

                        // 14. yy-zz (pair 2)
                        u[p_huang_start + 13][blk_offset + 1 * dim_nb + 1] -=
                            (sum_rr[2][2]);
                        u[p_huang_start + 13][blk_offset + 2 * dim_nb + 2] +=
                            (sum_rr[1][1]);

                        // 15. yz (final)
                        u[p_huang_start + 14][blk_offset + 1 * dim_nb + 2] -=
                            (sum_rr[2][2]);
                        u[p_huang_start + 14][blk_offset + 2 * dim_nb + 2] +=
                            (sum_rr[1][2]);
                    }
                }
            }
        }
    }

    /* --------------------------------------------------------------------- */
    /* 4. Symmetrization */
    /* --------------------------------------------------------------------- */

    for (ND_int n1 = 0; n1 < nr1; n1++)
    {
        ND_int m1 = (nr1 - n1) % nr1;
        for (ND_int n2 = 0; n2 < nr2; n2++)
        {
            ND_int m2 = (nr2 - n2) % nr2;
            for (ND_int n3 = 0; n3 < nr3; n3++)
            {
                ND_int m3 = (nr3 - n3) % nr3;

                ND_int g_A = n1 * nr2 * nr3 + n2 * nr3 + n3;
                ND_int g_B = m1 * nr2 * nr3 + m2 * nr3 + m3;

                for (ND_int na = 0; na < nat; na++)
                {
                    for (ND_int nb = 0; nb < nat; nb++)
                    {
                        for (ND_int i = 0; i < 3; i++)
                        {
                            for (ND_int j = 0; j < 3; j++)
                            {
                                ND_int idx_A = g_A * dim_na + na * dim_i +
                                               i * dim_nb + nb * 3 + j;
                                ND_int idx_B = g_B * dim_na + nb * dim_i +
                                               j * dim_nb + na * 3 + i;

                                if (idx_A <= idx_B)
                                {
                                    ELPH_cmplx avg_f =
                                        0.5 * (frc[idx_A] + frc[idx_B]);
                                    frc[idx_A] = avg_f;
                                    frc[idx_B] = avg_f;

                                    for (ND_int k = 0; k < total_vecs; k++)

                                    {
                                        ELPH_cmplx avg_u =
                                            0.5 * (u[k][idx_A] + u[k][idx_B]);
                                        u[k][idx_A] = avg_u;
                                        u[k][idx_B] = avg_u;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /*  --------------------------------------------------------------------- */
    /*    5. Gram-Schmidt Orthogonalization */
    /*  --------------------------------------------------------------------- */

    ELPH_cmplx* w = malloc(dim_total * sizeof(ELPH_cmplx));
    CHECK_ALLOC(w);

    ND_int active_count = (n_trans + n_rot) * nat + n_huang;

    ND_int* keep_mask = calloc(active_count, sizeof(ND_int));
    CHECK_ALLOC(keep_mask);

    for (ND_int k = 0; k < active_count; k++)
    {
        keep_mask[k] = 1;

        for (ND_int z = 0; z < dim_total; z++)
        {
            w[z] = u[k][z];
        }

        for (ND_int q = 0; q < k; q++)
        {
            if (keep_mask[q])
            {
                ELPH_cmplx dot = Cmplxdot(u[q], w, dim_total);
                for (ND_int z = 0; z < dim_total; z++)
                {
                    w[z] -= dot * u[q][z];
                }
            }
        }

        ELPH_float norm = sqrt(creal(Cmplxdot(w, w, dim_total)));

        if (norm > ELPH_EPS)
        {
            ELPH_float inv_norm = 1.0 / norm;
            for (ND_int z = 0; z < dim_total; z++)
            {
                u[k][z] = w[z] * inv_norm;
            }
        }
        else
        {
            keep_mask[k] = 0;
        }
    }

    /* ---------------------------------------------------------------------*/
    /*     6. Projection and Subtraction */
    /*  ---------------------------------------------------------------------*/

    for (ND_int z = 0; z < dim_total; z++)
    {
        w[z] = 0.0;
    }

    for (ND_int k = 0; k < active_count; k++)
    {
        if (keep_mask[k])
        {
            ELPH_cmplx dot = Cmplxdot(u[k], frc, dim_total);
            for (ND_int z = 0; z < dim_total; z++)
            {
                w[z] += dot * u[k][z];
            }
        }
    }

    for (ND_int z = 0; z < dim_total; z++)
    {
        frc[z] -= w[z];
    }

    free(w);
    free(keep_mask);
    for (ND_int k = 0; k < total_vecs; k++)
    {
        free(u[k]);
    }
    free(u);
}
