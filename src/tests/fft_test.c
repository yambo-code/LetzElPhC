#include "../fft/fft.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "../common/parallel.h"
#include "../elphC.h"

enum
{
    NS_PER_SECOND = 1000000000
};

void sub_timespec(struct timespec t1, struct timespec t2, struct timespec* td)
{
    td->tv_nsec = t2.tv_nsec - t1.tv_nsec;
    td->tv_sec = t2.tv_sec - t1.tv_sec;
    if (td->tv_sec > 0 && td->tv_nsec < 0)
    {
        td->tv_nsec += NS_PER_SECOND;
        td->tv_sec--;
    }
    else if (td->tv_sec < 0 && td->tv_nsec > 0)
    {
        td->tv_nsec -= NS_PER_SECOND;
        td->tv_sec++;
    }
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    if (argc != 4)
    {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    ND_int fft_dims[3];

    fft_dims[0] = atoi(argv[1]);
    fft_dims[1] = atoi(argv[2]);
    fft_dims[2] = atoi(argv[3]);
    ND_int FFTx, FFTy, FFTz;
    FFTx = fft_dims[0];
    FFTy = fft_dims[1];
    FFTz = fft_dims[2];

    ND_int Gxy_shift, G_vecs_xy;
    ND_int zshift, zloc;

    G_vecs_xy = get_mpi_local_size_idx(fft_dims[0] * fft_dims[1], &Gxy_shift,
                                       MPI_COMM_WORLD);
    zloc = get_mpi_local_size_idx(fft_dims[2], &zshift, MPI_COMM_WORLD);

    ND_int size_G_vecs = G_vecs_xy * fft_dims[2];
    ELPH_cmplx* dVlocG =
        calloc(2 * size_G_vecs, sizeof(ELPH_cmplx));  // 3*natom* ix_s*jy_s*kz
    int* gvecs = malloc(3 * size_G_vecs * sizeof(int));
    // set gvecs
    for (ND_int ig = 0; ig < G_vecs_xy; ++ig)
    {
        ND_int fft_glob_idx = Gxy_shift + ig;
        int ix = fft_glob_idx / FFTy;
        int jy = fft_glob_idx % FFTy;

        for (ND_int kz = 0; kz < FFTz; ++kz)
        {
            int* restrict g_temp_set = gvecs + FFTz * 3 * ig + kz * 3;
            g_temp_set[0] = ix;
            g_temp_set[1] = jy;
            g_temp_set[2] = kz;
        }
    }

    ELPH_cmplx* Vlocr = malloc(4 * fft_dims[0] * fft_dims[1] * zloc *
                               sizeof(ELPH_cmplx));  // 3*natom* ix_s*jy_s*kz
    ELPH_cmplx* Psir = malloc(2 * fft_dims[0] * fft_dims[1] * zloc *
                              sizeof(ELPH_cmplx));  // 3*natom* ix_s*jy_s*kz

    // fill data
    for (ND_int i = 0; i < (4 * fft_dims[0] * fft_dims[1] * zloc); ++i)
    {
        Vlocr[i] = sin(i);
    }
    for (ND_int i = 0; i < (2 * fft_dims[0] * fft_dims[1] * zloc); ++i)
    {
        Psir[i] = cos(i);
    }

    struct ELPH_fft_plan fft_plan;

    wfc_plan(&fft_plan, size_G_vecs, zloc, G_vecs_xy, gvecs, fft_dims,
             FFTW_MEASURE, MPI_COMM_WORLD);

    struct timespec start, finish, delta;  // timing vars

    MPI_Barrier(MPI_COMM_WORLD);
    clock_gettime(CLOCK_REALTIME, &start);

    fft_convolution3D(&fft_plan, 2, 4, Vlocr, Psir, dVlocG, false);

    MPI_Barrier(MPI_COMM_WORLD);
    clock_gettime(CLOCK_REALTIME, &finish);  // timing end:
    sub_timespec(start, finish, &delta);     // timing end:
    if (!world_rank)
    {
        printf("Time taken for run : %d.%.9ld\n", (int)delta.tv_sec,
               delta.tv_nsec);
    }

    wfc_destroy_plan(&fft_plan);

    free(dVlocG);
    free(gvecs);
    free(Vlocr);
    free(Psir);

    MPI_Finalize();

    return 0;
}