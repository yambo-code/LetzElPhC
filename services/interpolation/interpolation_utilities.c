// This file contains some helper functions used in interpolation

#include "interpolation_utilities.h"

#include "elphC.h"
#include "fft/fft.h"
// we must always place complex.h before fftw3.
#include <complex.h>
#include <ctype.h>
#include <fftw3.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/ELPH_timers.h"
#include "common/constants.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/omp_pragma_def.h"

static int qpt_sort_cmp(const void* a, const void* b);

void fft_q2R(ELPH_cmplx* data, const ND_int* qgrid, const ND_int nsets)
{
    // given a data (qx,qy,qz,nsets)
    // performs inplace fourier transform q->R
    ELPH_start_clock("q2R FFT");
    //
    ND_int qx = qgrid[0];
    ND_int qy = qgrid[1];
    ND_int qz = qgrid[2];

    fftw_fun(iodim64) dims[3];
    dims[0].n = qx;
    dims[0].is = qy * qz * nsets;
    dims[0].os = qy * qz * nsets;

    dims[1].n = qy;
    dims[1].is = qz * nsets;
    dims[1].os = qz * nsets;

    dims[2].n = qz;
    dims[2].is = nsets;
    dims[2].os = nsets;

    // How the transforms are repeated over nsets
    fftw_fun(iodim64) howmany_dims[1];
    howmany_dims[0].n = nsets;
    howmany_dims[0].is = 1;
    howmany_dims[0].os = 1;

    // DOnot overwrite the data, so always set to
    // FFTW_ESTIMATE.
    fftw_generic_plan plan = fftw_fun(plan_guru64_dft)(
        3, dims, 1, howmany_dims, data, data, FFTW_FORWARD, FFTW_ESTIMATE);

    if (NULL == plan)
    {
        error_msg("FFT q->R plan failed.");
    }

    fftw_fun(execute)(plan);
    fftw_fun(destroy_plan)(plan);
    // Normalize
    ELPH_float norm_fft = 1.0 / (qx * qy * qz);
    //
    ELPH_OMP_PAR_FOR_SIMD
    for (ND_int i = 0; i < (qx * qy * qz * nsets); ++i)
    {
        data[i] *= norm_fft;
    }
    ELPH_stop_clock("q2R FFT");
    return;
}

void fft_R2q_dyn(const ELPH_cmplx* dataR, const ELPH_float* qpt_crys,
                 const ND_int* qgrid, const ND_int natom, const ND_int* ws_vecs,
                 const ND_int nws, const ND_int* ws_degen, ELPH_cmplx* dataq)
{
    // (Rx, Ry, Rz, nmodes, nmodes)-> (nmodes, nmodes)
    // This is (N^2) fouier transform (slow one).
    // G. Pizzi et al 2020 J. Phys.: Condens. Matter 32 165902

    ELPH_start_clock("R2q dyn invFFT");

    const ND_int qx = qgrid[0];
    const ND_int qy = qgrid[1];
    const ND_int qz = qgrid[2];
    const ND_int nRpts = qx * qy * qz;

    const ND_int nmodes = 3 * natom;

    // set the output buffer to 0;
    //
    for (ND_int i = 0; i < nmodes * nmodes; ++i)
    {
        dataq[i] = 0.0;
    }

    // Do fourier interpolation
    // See Eq. 47 of
    // G. Pizzi et al 2020 J. Phys.: Condens. Matter 32 165902
    ND_int Gridyz = qy * qz;
    ND_int iws_vec = 0;
    //
    for (ND_int i = 0; i < nRpts; ++i)
    {
        ND_int Rx = i / Gridyz;
        ND_int Ry = (i % Gridyz) / qgrid[2];
        ND_int Rz = (i % Gridyz) % qgrid[2];
        //
        Rx = get_miller_idx(Rx, qgrid[0]);
        Ry = get_miller_idx(Ry, qgrid[1]);
        Rz = get_miller_idx(Rz, qgrid[2]);

        ELPH_float Rpt[3] = {Rx, Ry, Rz};
        const ELPH_cmplx eiqR =
            cexp(I * 2 * ELPH_PI * dot3_macro(qpt_crys, Rpt));
        //
        for (ND_int ia = 0; ia < natom; ++ia)
        {
            for (ND_int ja = 0; ja < natom; ++ja)
            {
                const ND_int iws_vecs_degen =
                    ws_degen[i * natom * natom + ia * natom + ja];
                // compute \sum e^iqGr
                ELPH_cmplx eiqTR = 0.0;
                for (ND_int ii = 0; ii < iws_vecs_degen; ++ii)
                {
                    if (iws_vec >= nws)
                    {
                        error_msg("Wigner seitz vectors Out of bound.");
                    }
                    const ELPH_float dot_TRq =
                        qpt_crys[0] * ws_vecs[3 * iws_vec] +
                        qpt_crys[1] * ws_vecs[3 * iws_vec + 1] +
                        qpt_crys[2] * ws_vecs[3 * iws_vec + 2];
                    //
                    eiqTR += cexp(I * 2 * ELPH_PI * dot_TRq);
                    //
                    ++iws_vec;
                }
                eiqTR *= eiqR;
                eiqTR *= (1.0 / iws_vecs_degen);
                //
                for (ND_int ixy = 0; ixy < 9; ++ixy)
                {
                    const ND_int ix = ixy / 3;
                    const ND_int iy = ixy % 3;
                    //
                    const ND_int im = ix + 3 * ia;
                    const ND_int jm = iy + 3 * ja;
                    dataq[im * nmodes + jm] +=
                        (eiqTR * dataR[i * nmodes * nmodes + im * nmodes + jm]);
                }
            }
        }
    }
    ELPH_stop_clock("R2q dyn invFFT");
    return;
}

void fft_R2q_dvscf(const ELPH_cmplx* dataR, const ELPH_float* qpt_crys,
                   const ND_int* qgrid, const ND_int natom, const ND_int nsets,
                   const ND_int* ws_vecs, const ND_int nws,
                   const ND_int* ws_degen, ELPH_cmplx* dataq)
{
    // (Rx, Ry, Rz, nmodes, nsets)-> (nmodes, nsets)
    // This is (N^2) fouier transform (slow one).
    // G. Pizzi et al 2020 J. Phys.: Condens. Matter 32 165902
    ELPH_start_clock("R2q dvscf invFFT");
    const ND_int qx = qgrid[0];
    const ND_int qy = qgrid[1];
    const ND_int qz = qgrid[2];
    const ND_int nRpts = qx * qy * qz;

    const ND_int nmodes = 3 * natom;

    // set the output buffer to 0;
    //
    for (ND_int i = 0; i < nmodes * nsets; ++i)
    {
        dataq[i] = 0.0;
    }

    ND_int iws_vec = 0;
    ND_int iws_vecs_degen = 1;
    ELPH_cmplx eiqTR = 0.0;
    // Do fourier interpolation
    // See Eq. 47 of
    // G. Pizzi et al 2020 J. Phys.: Condens. Matter 32 165902
    ND_int Gridyz = qy * qz;
    for (ND_int i = 0; i < nRpts; ++i)
    {
        ND_int Rx = i / Gridyz;
        ND_int Ry = (i % Gridyz) / qgrid[2];
        ND_int Rz = (i % Gridyz) % qgrid[2];
        //
        Rx = get_miller_idx(Rx, qgrid[0]);
        Ry = get_miller_idx(Ry, qgrid[1]);
        Rz = get_miller_idx(Rz, qgrid[2]);

        ELPH_float Rpt[3] = {Rx, Ry, Rz};
        ELPH_cmplx eiqR = cexp(I * 2 * ELPH_PI * dot3_macro(qpt_crys, Rpt));
        //
        for (ND_int ia = 0; ia < natom; ++ia)
        {
            iws_vecs_degen = ws_degen[i * natom + ia];
            // compute \sum e^iqGr
            eiqTR = 0.0;
            for (ND_int ii = 0; ii < iws_vecs_degen; ++ii)
            {
                if (iws_vec >= nws)
                {
                    error_msg("Wigner seitz vectors Out of bound.");
                }
                ELPH_float dot_TRq = qpt_crys[0] * ws_vecs[3 * iws_vec] +
                                     qpt_crys[1] * ws_vecs[3 * iws_vec + 1] +
                                     qpt_crys[2] * ws_vecs[3 * iws_vec + 2];
                eiqTR += cexp(I * 2 * ELPH_PI * dot_TRq);
                ++iws_vec;
            }
            eiqTR *= eiqR;
            eiqTR *= (1.0 / iws_vecs_degen);
            //
            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int isets = 0; isets < 3 * nsets; ++isets)
            {
                dataq[ia * 3 * nsets + isets] +=
                    (eiqTR * dataR[(i * nmodes + ia * 3) * nsets + isets]);
            }
        }
    }
    ELPH_stop_clock("R2q dvscf invFFT");
    return;
}

void Sorted_qpts_idxs(const ND_int nqpts, ELPH_float* qpts, ND_int* indices)
{
    /*
    given list of qpts, sort then in fft_grid order.
    Return : sorted indices

    // qpts in crystal coordinates
    */
    if (NULL == qpts)
    {
        return;
    }

    ELPH_float** qvec_ptrs = malloc(sizeof(ELPH_float*) * nqpts);
    CHECK_ALLOC(qvec_ptrs);

    /* fill the struct */
    for (ND_int i = 0; i < nqpts; ++i)
    {
        qvec_ptrs[i] = qpts + 3 * i;
    }

    qsort(qvec_ptrs, nqpts, sizeof(ELPH_float*), qpt_sort_cmp);

    // store the sorted indices
    for (ND_int i = 0; i < nqpts; ++i)
    {
        indices[i] = (qvec_ptrs[i] - qpts) / 3;
    }

    free(qvec_ptrs);
}

void rearrange_qpt_grid(const ND_int nqpts, const ELPH_cmplx* in_buf,
                        const ND_int* idx, ELPH_cmplx* out_buf)
{
    // out[i] = in[idx[i]]
    for (ND_int i = 0; i < nqpts; ++i)
    {
        out_buf[i] = in_buf[idx[i]];
    }
}

void find_qpt_grid(const ND_int nqpts, const ELPH_float* qpts, ND_int* q_grid)
{
    // given list of qpoints, it finds the phonon grid.
    // q_grid is 3 ints

    bool is_gamma_centred = false;
    for (ND_int i = 0; i < 3; ++i)
    {
        q_grid[i] = 1;
    }

    for (ND_int i = 0; i < nqpts; ++i)
    {
        bool is_gamma = true;
        for (ND_int ix = 0; ix < 3; ++ix)
        {
            ELPH_float qx = qpts[3 * i + ix] - floor(qpts[3 * i + ix]);

            if (fabs(qx) > ELPH_EPS && fabs(qx - 1) > ELPH_EPS)
            {
                is_gamma = false;
                qx = 1.0 / qx;
                ND_int qx_int = rint(qx);
                qx -= qx_int;

                if (fabs(qx) < ELPH_EPS)
                {
                    q_grid[ix] = (qx_int > q_grid[ix]) ? qx_int : q_grid[ix];
                }
            }
        }
        if (is_gamma)
        {
            is_gamma_centred = true;
        }
    }

    if (!is_gamma_centred)
    {
        error_msg("Only Gamma centred grids supported.");
    }
    // finally do a sanity check i.e the product of three dimensions are equal
    // to nqpts
    ND_int nqpt_prod = q_grid[0] * q_grid[1] * q_grid[2];

    if (nqpts != nqpt_prod)
    {
        error_msg("Number of qpts != product of found qpt grid.");
    }
}

ELPH_float* parse_qpt_entries(const char* filename, ND_int* count_out)
{
    FILE* fp = fopen(filename, "r");
    if (!fp)
    {
        fprintf(stderr, "Error : Unable to open file %s.\n", filename);
        return NULL;
    }

    char line[1024];
    long long expected_entries = -1;
    long long current_entry_idx = 0;
    bool header_found = false;

    ELPH_float* out_buffer = NULL;

    while (fgets(line, sizeof(line), fp))
    {
        char* ptr = line;
        // Skip leading whitespace
        while (isspace((unsigned char)*ptr))
        {
            ptr++;
        }
        // skip empty lines or comment lines
        if (*ptr == '\0' || *ptr == '#')
        {
            continue;
        }

        if (!header_found)
        {
            if (sscanf(ptr, "%lld", &expected_entries) == 1)
            {
                if (expected_entries <= 0)
                {
                    fclose(fp);
                    return NULL;
                }

                out_buffer = malloc(expected_entries * 3 * sizeof(*out_buffer));
                if (!out_buffer)
                {
                    fclose(fp);
                    return NULL;
                }
                header_found = true;
            }
            else
            {
                fclose(fp);
                return NULL;
            }
        }
        else
        {
            // Processing data entries
            if (current_entry_idx >= expected_entries)
            {
                free(out_buffer);
                fclose(fp);
                return NULL;
            }

            double tmp_buffer[3];
            if (sscanf(ptr, "%lf %lf %lf", &tmp_buffer[0], &tmp_buffer[1],
                       &tmp_buffer[2]) == 3)
            {
                long long offset = current_entry_idx * 3;
                out_buffer[offset + 0] = tmp_buffer[0];
                out_buffer[offset + 1] = tmp_buffer[1];
                out_buffer[offset + 2] = tmp_buffer[2];
                current_entry_idx++;
            }
            else
            {
                free(out_buffer);
                fclose(fp);
                return NULL;
            }
        }
    }

    fclose(fp);
    //
    if (header_found && current_entry_idx == expected_entries)
    {
        if (count_out)
        {
            *count_out = expected_entries;
        }
        return out_buffer;
    }
    else
    {
        fprintf(stderr,
                "Error : In file %s, Number of qpt entries are not equal to "
                "excepted "
                "number.\n",
                filename);
        if (out_buffer)
        {
            free(out_buffer);
        }
        return NULL;
    }
}

/*  Static functions */
static int qpt_sort_cmp(const void* a, const void* b)
{
    /*
    First sorts along x, then y and finally z
    This will ensure that the fastest will be z, followed by y and then x
    */
    ELPH_float* v1 = *(ELPH_float**)a;
    ELPH_float* v2 = *(ELPH_float**)b;

    ELPH_float kdiff[3];

    // first bring them to [0,1)
    for (int i = 0; i < 3; ++i)
    {
        ELPH_float k1 = v1[i];
        k1 -= floor(k1);
        k1 += ELPH_EPS;
        k1 -= floor(k1);

        ELPH_float k2 = v2[i];
        k2 -= floor(k2);
        k2 += ELPH_EPS;
        k2 -= floor(k2);

        kdiff[i] = k1 - k2;
    }

    if (fabs(kdiff[0]) < ELPH_EPS)
    {
        if (fabs(kdiff[1]) < ELPH_EPS)
        {
            if (fabs(kdiff[2]) < ELPH_EPS)
            {
                return 0;
            }
            else
            {
                return kdiff[2] > 0.0 ? 1 : -1;
            }
        }
        else
        {
            return kdiff[1] > 0.0 ? 1 : -1;
        }
    }
    else
    {
        return kdiff[0] > 0.0 ? 1 : -1;
    }
}
