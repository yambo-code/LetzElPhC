/*
 * License-Identifier: GPL
 *
 * Copyright (C) 2025 The Yambo Team
 *
 * Authors (see AUTHORS file for details): RR AM
 *
 * Gather distributed dVscf z-slabs to commK rank 0 and apply a forward 3D
 * FFT, producing dV_q^nu(G) ready for the elph_dvG_fill_fn callback.
 */

#include "dvG_utils.h"
#include "fft/fft.h"
#include "common/error.h"
#include <stdlib.h>
#include <string.h>

ELPH_cmplx* gather_dVscf_and_fft(const ELPH_cmplx* dVscf,
                                   const struct Lattice* lattice,
                                   MPI_Comm commK)
{
    ND_int nmodes  = lattice->nmodes;
    ND_int nmag    = lattice->nmag;
    ND_int Nx      = lattice->fft_dims[0];
    ND_int Ny      = lattice->fft_dims[1];
    ND_int Nz      = lattice->fft_dims[2];
    ND_int Nz_loc  = lattice->nfftz_loc;
    ND_int shift   = lattice->nfftz_loc_shift;
    ND_int Nxy     = Nx * Ny;
    ND_int NmodMag = nmodes * nmag;

    int commK_rank, commK_size;
    MPI_Comm_rank(commK, &commK_rank);
    MPI_Comm_size(commK, &commK_size);

    /* Collect per-rank Nz_loc and shifts. */
    int* all_Nz_loc = malloc(commK_size * sizeof(int));
    int* all_shifts = malloc(commK_size * sizeof(int));
    CHECK_ALLOC(all_Nz_loc);
    CHECK_ALLOC(all_shifts);
    int my_Nz_loc = (int)Nz_loc;
    int my_shift  = (int)shift;
    MPI_Allgather(&my_Nz_loc, 1, MPI_INT, all_Nz_loc, 1, MPI_INT, commK);
    MPI_Allgather(&my_shift,  1, MPI_INT, all_shifts,  1, MPI_INT, commK);

    /*
     * Pack local slab: rearrange from
     *   dVscf[nmodes][nmag][Nx][Ny][Nz_loc]  (C row-major)
     * to
     *   packed[Nz_loc][nmodes][nmag][Nx][Ny]
     * so that the z-extent is the leading dimension for contiguous Gatherv.
     */
    ND_int local_count = Nz_loc * NmodMag * Nxy;
    ELPH_cmplx* packed = malloc(local_count * sizeof(ELPH_cmplx));
    CHECK_ALLOC(packed);

    for (ND_int iv = 0; iv < nmodes; ++iv)
    for (ND_int is = 0; is < nmag;   ++is)
    for (ND_int ix = 0; ix < Nx;     ++ix)
    for (ND_int iy = 0; iy < Ny;     ++iy)
    for (ND_int iz = 0; iz < Nz_loc; ++iz)
    {
        packed[iz * NmodMag*Nxy + iv * nmag*Nxy + is * Nxy + ix * Ny + iy] =
            dVscf[iv * nmag*Nxy*Nz_loc + is * Nxy*Nz_loc
                  + ix * Ny*Nz_loc + iy * Nz_loc + iz];
    }

    /* Build Gatherv send/recv counts and displacements (only on rank 0). */
    int* sendcounts  = NULL;
    int* displs      = NULL;
    ELPH_cmplx* gathered = NULL;
    ND_int NmodMagNxy = NmodMag * Nxy;

    if (commK_rank == 0)
    {
        sendcounts = malloc(commK_size * sizeof(int));
        displs     = malloc(commK_size * sizeof(int));
        CHECK_ALLOC(sendcounts);
        CHECK_ALLOC(displs);
        for (int r = 0; r < commK_size; ++r)
        {
            sendcounts[r] = all_Nz_loc[r] * (int)NmodMagNxy;
            displs[r]     = all_shifts[r]  * (int)NmodMagNxy;
        }
        gathered = malloc((size_t)Nz * NmodMag * Nxy * sizeof(ELPH_cmplx));
        CHECK_ALLOC(gathered);
    }

    MPI_Gatherv(packed, (int)local_count, ELPH_MPI_cmplx,
                gathered, sendcounts, displs, ELPH_MPI_cmplx,
                0, commK);

    free(packed);
    free(all_Nz_loc);
    free(all_shifts);
    if (commK_rank == 0)
    {
        free(sendcounts);
        free(displs);
    }

    if (commK_rank != 0)
        return NULL;

    /*
     * Transpose gathered[Nz][nmodes][nmag][Nx][Ny] back to
     * dVG[nmodes][nmag][Nx][Ny][Nz].
     */
    ELPH_cmplx* dVG = malloc((size_t)NmodMag * Nxy * Nz * sizeof(ELPH_cmplx));
    CHECK_ALLOC(dVG);

    for (ND_int iz = 0; iz < Nz;     ++iz)
    for (ND_int iv = 0; iv < nmodes; ++iv)
    for (ND_int is = 0; is < nmag;   ++is)
    for (ND_int ix = 0; ix < Nx;     ++ix)
    for (ND_int iy = 0; iy < Ny;     ++iy)
    {
        dVG[iv * nmag*Nxy*Nz + is * Nxy*Nz + ix * Ny*Nz + iy * Nz + iz] =
            gathered[iz * NmodMagNxy + iv * nmag*Nxy + is * Nxy + ix * Ny + iy];
    }
    free(gathered);

    /*
     * Forward 3D FFT for each (mode, spin) slice.
     * fftw plan_dft_3d uses row-major [n0][n1][n2] = [Nx][Ny][Nz].
     */
    ND_int slice_sz = Nxy * Nz;
    for (ND_int iv = 0; iv < nmodes; ++iv)
    for (ND_int is = 0; is < nmag;   ++is)
    {
        ELPH_cmplx* sl = dVG + (iv * nmag + is) * slice_sz;
        fftw_generic_plan plan = fftw_fun(plan_dft_3d)(
            (int)Nx, (int)Ny, (int)Nz,
            sl, sl,
            FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_fun(execute)(plan);
        fftw_fun(destroy_plan)(plan);
    }

    return dVG;
}
