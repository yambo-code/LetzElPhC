#include <complex.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/constants.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/omp_pragma_def.h"
#include "elphC.h"
#include "symmetries.h"
#include "wfc/wfc.h"
/*
This file contains function to rotate dvscf
*/

void rotate_dvscf(const ELPH_float* qpt, const ELPH_cmplx* dvscf_in,
                  struct symmetry* sym, const struct Lattice* lattice,
                  const bool composite_form, ELPH_cmplx* restrict dvscf_out,
                  MPI_Comm commK)
{
    // dvscf : (nmodes, nmag, Nx, Ny, Nz_loc)
    // composite_form : if true. dvscf is in V*s_0 + Bx*sx + By*Sy + Bz*sz
    // where s_0,x,y,z are pauli matrices.
    // if composite_form is false then [Vx, Bx, By, Bz] are given i.e
    // The output will always be in the form of the input
    // indivividual consitututes are given (like in quantum espresso)
    //
    const ND_int* fft_dims = lattice->fft_dims;
    const ND_int nmodes = lattice->nmodes;
    const ND_int nmag = lattice->nmag;

    // before anything, do basic checks.

    int comm_rank, comm_size;
    //
    int mpi_error = MPI_Comm_rank(commK, &comm_rank);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_size(commK, &comm_size);
    MPI_error_msg(mpi_error);

    ND_int nfft = fft_dims[0] * fft_dims[1] * fft_dims[2];

    ND_int nfft_loc = fft_dims[0] * fft_dims[1] * lattice->nfftz_loc;
    // Check if this overflows int (because we need use this in mpi functions)
    CHECK_OVERFLOW_ERROR(nfft_loc, INT_MAX);
    CHECK_OVERFLOW_ERROR(fft_dims[2], INT_MAX);

    // first setup few things
    ELPH_cmplx* tmp_buf = NULL;
    ND_int* rot_idx = NULL;

    int* recv_counts = NULL;
    int* recv_disp = NULL;

    if (!comm_rank)
    {
        recv_counts = malloc(sizeof(*recv_counts) * comm_size);
        CHECK_ALLOC(recv_counts);
        //
        recv_disp = malloc(sizeof(*recv_disp) * comm_size);
        CHECK_ALLOC(recv_disp);

        rot_idx = malloc(nfft * sizeof(*rot_idx));
        CHECK_ALLOC(rot_idx);

        tmp_buf = malloc(2 * nfft * sizeof(*tmp_buf));
        CHECK_ALLOC(tmp_buf);
    }

    int nz_loc = lattice->nfftz_loc;
    int nz_loc_shift = lattice->nfftz_loc_shift;

    mpi_error =
        MPI_Gather(&nz_loc, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, commK);
    MPI_error_msg(mpi_error);

    mpi_error =
        MPI_Gather(&nz_loc_shift, 1, MPI_INT, recv_disp, 1, MPI_INT, 0, commK);
    MPI_error_msg(mpi_error);

    // Create MPI datatypes
    MPI_Datatype column_type, column_type_loc;
    MPI_Datatype column_type_r, column_type_loc_r;
    //
    // Create a gatter/scatter data-type for root
    mpi_error = MPI_Type_vector(fft_dims[0] * fft_dims[1], 1, fft_dims[2],
                                ELPH_MPI_cmplx, &column_type);
    MPI_error_msg(mpi_error);
    // Create a local mpi data to send or get
    mpi_error = MPI_Type_vector(fft_dims[0] * fft_dims[1], 1, nz_loc,
                                ELPH_MPI_cmplx, &column_type_loc);
    MPI_error_msg(mpi_error);
    //
    mpi_error = MPI_Type_commit(&column_type);
    MPI_error_msg(mpi_error);
    //
    mpi_error = MPI_Type_commit(&column_type_loc);
    MPI_error_msg(mpi_error);

    // change the extend
    mpi_error = MPI_Type_create_resized(column_type, 0, sizeof(*tmp_buf),
                                        &column_type_r);
    MPI_error_msg(mpi_error);
    //
    mpi_error = MPI_Type_create_resized(column_type_loc, 0, sizeof(*tmp_buf),
                                        &column_type_loc_r);
    MPI_error_msg(mpi_error);
    //
    mpi_error = MPI_Type_commit(&column_type_r);
    MPI_error_msg(mpi_error);
    //
    mpi_error = MPI_Type_commit(&column_type_loc_r);
    MPI_error_msg(mpi_error);

    // next do the actual thing.
    //

    if (!comm_rank)
    {
        ELPH_float sym_red[9], tau_red[3];
        const ELPH_float* blat = lattice->blat_vec;

        // convert to reduced coordinates.
        MatVec3f(blat, sym->tau, true, tau_red);
        for (int xi = 0; xi < 3; ++xi)
        {
            tau_red[xi] = tau_red[xi] / (2 * ELPH_PI);
            if (sym->time_rev)
            {
                tau_red[xi] = -tau_red[xi];
            }
        }

        // get the inverse of symmetry matrix in reduced cooordiantes
        {
            // S^{-1}_{red} = blat.T @ sym_{cart}.T @ lat
            ELPH_float tmp_sym[9];
            Gemm3x3f(sym->Rmat, 'T', lattice->alat_vec, 'N', tmp_sym);
            // because blat will have (2*pi) factor, we remove it first
            for (int xi = 0; xi < 9; ++xi)
            {
                tmp_sym[xi] = tmp_sym[xi] / (2 * ELPH_PI);
                if (sym->time_rev)
                {
                    tmp_sym[xi] = -tmp_sym[xi];
                }
            }
            Gemm3x3f(blat, 'T', tmp_sym, 'N', sym_red);
        }

        ND_int Nyz = fft_dims[1] * fft_dims[2];
        // first find the rotation indices i.e S^{r-tau}
        // check if symmetry is compatible with fft grid
        volatile bool symm_fft_compat = true;
        // volatile is important to insure we reload everytime (for example if
        // other treads change the value)
        ELPH_OMP_PAR_FOR
        for (ND_int i = 0; i < nfft; ++i)
        {
            if (!symm_fft_compat)
            {
                continue;
            }
            ND_int ix = i / Nyz;
            ND_int iy = (i % Nyz) / fft_dims[2];
            ND_int iz = (i % Nyz) % fft_dims[2];
            ELPH_float tmp_idxs[3] = {ix, iy, iz};
            for (int xi = 0; xi < 3; ++xi)
            {
                tmp_idxs[xi] = tmp_idxs[xi] / fft_dims[xi];
                tmp_idxs[xi] -= tau_red[xi];
            }
            ELPH_float tmp_idxs_rot[3];
            MatVec3f(sym_red, tmp_idxs, false, tmp_idxs_rot);
            // bring it [0,1) and multiply with Nx, Ny, Nz
            for (int xi = 0; xi < 3; ++xi)
            {
                tmp_idxs_rot[xi] -= floor(tmp_idxs_rot[xi]);
                // if there are any 0.9999999999, add a small tolerence and
                // wrap again to [0,1).
                tmp_idxs_rot[xi] -= floor(tmp_idxs_rot[xi] + ELPH_EPS);
                tmp_idxs_rot[xi] -= ELPH_EPS;
                tmp_idxs_rot[xi] *= fft_dims[xi];

                ELPH_float tmp_diff = tmp_idxs_rot[xi] - rint(tmp_idxs_rot[xi]);
                // check now if it is an integer.
                // if (fabs(tmp_diff) > ELPH_EPS)
                // NM : Going below 1-e3 on single precision is giving an error.
                if (fabs(tmp_diff) > 1e-3)
                {
                    if (symm_fft_compat)
                    {
                        ELPH_OMP_PAR_CRITICAL { symm_fft_compat = false; }
                    }
                    /* error_msg( */
                    /*     "Symmetry operation is incompatible with FFT grid.");
                     */
                }
            }
            if (!symm_fft_compat)
            {
                continue;
            }
            ix = rint(tmp_idxs_rot[0]);
            iy = rint(tmp_idxs_rot[1]);
            iz = rint(tmp_idxs_rot[2]);
            rot_idx[i] = ix * Nyz + iy * fft_dims[2] + iz;
            if (rot_idx[i] >= nfft)
            {
                error_msg("rotated FFT index out of bounds");
            }
        }

        if (!symm_fft_compat)
        {
            error_msg("Symmetry operation is incompatible with FFT grid.");
        }
    }

    // compute phase factor from fractional translation
    ELPH_float qcart[3] = {0.0, 0.0, 0.0};
    ELPH_float Rqcart[3] = {0.0, 0.0, 0.0};
    MatVec3f(lattice->blat_vec, qpt, false, qcart);
    MatVec3f(sym->Rmat, qcart, false, Rqcart);
    const ELPH_cmplx exp_iRqv = cexp(-I * (dot3_macro(Rqcart, sym->tau)));

    // rotate the fft grid
    for (ND_int iset = 0; iset < nmodes * nmag; ++iset)
    {
        // gather the grid to root.
        ELPH_cmplx* tmp_buf1 = tmp_buf;
        ELPH_cmplx* tmp_buf2 = tmp_buf ? tmp_buf + nfft : NULL;

        mpi_error = MPI_Gatherv(dvscf_in + iset * nfft_loc, nz_loc,
                                column_type_loc_r, tmp_buf1, recv_counts,
                                recv_disp, column_type_r, 0, commK);
        MPI_error_msg(mpi_error);

        if (!comm_rank)
        {
            for (ND_int i = 0; i < nfft; ++i)
            {
                tmp_buf2[i] = exp_iRqv * tmp_buf1[rot_idx[i]];
            }
        }

        // scatter back
        mpi_error = MPI_Scatterv(tmp_buf2, recv_counts, recv_disp,
                                 column_type_r, dvscf_out + iset * nfft_loc,
                                 nz_loc, column_type_loc_r, 0, commK);
        MPI_error_msg(mpi_error);
    }

    // spinor rotation
    if (4 == nmag)
    {
        // we need to do a spinor rotation
        ELPH_cmplx sumat[4];
        SU2mat(sym->Rmat, lattice->nspinor, false, sym->time_rev, sumat);

        for (ND_int imode = 0; imode < nmodes; ++imode)
        {
            ELPH_cmplx* dV_tmp_mode = dvscf_out + imode * nmag * nfft_loc;

            ELPH_cmplx* dV0 = dV_tmp_mode + 0 * nfft_loc;
            ELPH_cmplx* dV1 = dV_tmp_mode + 1 * nfft_loc;
            ELPH_cmplx* dV2 = dV_tmp_mode + 2 * nfft_loc;
            ELPH_cmplx* dV3 = dV_tmp_mode + 3 * nfft_loc;

            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int i = 0; i < nfft_loc; ++i)
            {
                // if not in composite form, convert to composite form
                if (!composite_form)
                {
                    ELPH_cmplx Vxc_t, Bx_t, By_t, Bz_t;
                    Vxc_t = dV0[i];
                    Bx_t = dV1[i];
                    By_t = dV2[i];
                    Bz_t = dV3[i];

                    dV0[i] = Vxc_t + Bz_t;
                    dV1[i] = Bx_t - I * By_t;
                    dV2[i] = Bx_t + I * By_t;
                    dV3[i] = Vxc_t - Bz_t;
                }
                // perform spinoral rotation i.e SU dV SU^dagger
                // // Compute sumat @ dV @ sumat^dagger
                ELPH_cmplx M0 =
                    sumat[0] *
                        (dV0[i] * conj(sumat[0]) + dV1[i] * conj(sumat[1])) +
                    sumat[1] *
                        (dV2[i] * conj(sumat[0]) + dV3[i] * conj(sumat[1]));

                ELPH_cmplx M1 =
                    sumat[0] *
                        (dV0[i] * conj(sumat[2]) + dV1[i] * conj(sumat[3])) +
                    sumat[1] *
                        (dV2[i] * conj(sumat[2]) + dV3[i] * conj(sumat[3]));

                ELPH_cmplx M2 =
                    sumat[2] *
                        (dV0[i] * conj(sumat[0]) + dV1[i] * conj(sumat[1])) +
                    sumat[3] *
                        (dV2[i] * conj(sumat[0]) + dV3[i] * conj(sumat[1]));

                ELPH_cmplx M3 =
                    sumat[2] *
                        (dV0[i] * conj(sumat[2]) + dV1[i] * conj(sumat[3])) +
                    sumat[3] *
                        (dV2[i] * conj(sumat[2]) + dV3[i] * conj(sumat[3]));

                if (composite_form)
                {
                    dV0[i] = M0;
                    dV1[i] = M1;
                    dV2[i] = M2;
                    dV3[i] = M3;
                }
                else
                {
                    dV0[i] = (M0 + M3) / 2.0;
                    dV1[i] = (M1 + M2) / 2.0;
                    dV2[i] = -I * (M2 - M1) / 2.0;
                    dV3[i] = (M0 - M3) / 2.0;
                }
            }
        }
    }

    if (sym->time_rev)
    {
        for (ND_int i = 0; i < (nmodes * nmag * nfft_loc); ++i)
        {
            dvscf_out[i] = conj(dvscf_out[i]);
        }
    }

    // free stuff
    // first free the derived types
    //
    mpi_error = MPI_Type_free(&column_type_loc_r);
    MPI_error_msg(mpi_error);
    //
    mpi_error = MPI_Type_free(&column_type_r);
    MPI_error_msg(mpi_error);
    //
    // Next free the main types
    // Note these must be free only after resized types
    mpi_error = MPI_Type_free(&column_type);
    MPI_error_msg(mpi_error);
    //
    mpi_error = MPI_Type_free(&column_type_loc);
    MPI_error_msg(mpi_error);
    //
    free(recv_counts);
    free(recv_disp);
    free(rot_idx);
    free(tmp_buf);
}
