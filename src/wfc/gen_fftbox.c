#include <math.h>
#include <mpi.h>

#include "common/error.h"
#include "elphC.h"
#include "wfc.h"

void get_fft_box(const ELPH_float EcutRy, const ELPH_float* blat,
                 ND_int* fft_box, MPI_Comm commK)
{
    /**
     * @brief Calculate FFT box (Nx, Ny, Nz). The dimensions are determined by:
     * 1. The energy cutoff (EcutRy) which defines the maximum kinetic energy.
     * 2. Each cpu will atleast have one Z stick i.e Nx * Ny > comm_size.
     *    So this is cpu dependent
     *
     * The initial dimensions are calculated as:
     *    N_i = 2*round(sqrt(Ecut)/|b_i| + 1) + 1
     * If the initial grid would be too small for parallelization (Nx*Ny <
     * num_processes), the grid is scaled up by a factor to ensure each process
     * gets at least one Gz stick.
     *
     * @param[in] EcutRy      Energy cutoff in Rydberg units
     * @param[in] blat        Array of reciprocal lattice vectors (9 elements,
     * 3x3 matrix)
     * @param[out] fft_box    Array to store the computed FFT box dimensions
     * (Nx, Ny, Nz)
     * @param[in] commK       MPI communicator for k-point parallelization
     *
     * @note The function uses MPI Bcast so it must be called by all processes
     * in commK
     */
    const ELPH_float b0_norm =
        sqrt(blat[0] * blat[0] + blat[3] * blat[3] + blat[6] * blat[6]);
    const ELPH_float b1_norm =
        sqrt(blat[1] * blat[1] + blat[4] * blat[4] + blat[7] * blat[7]);
    const ELPH_float b2_norm =
        sqrt(blat[2] * blat[2] + blat[5] * blat[5] + blat[8] * blat[8]);

    const ELPH_float Ecut_sqrt = sqrt(EcutRy);
    //
    ND_int Nx = rint(Ecut_sqrt / b0_norm + 1);
    ND_int Ny = rint(Ecut_sqrt / b1_norm + 1);
    ND_int Nz = rint(Ecut_sqrt / b2_norm + 1);
    // we have -ve index gvecs and 0, so 2*N + 1
    Nx = 2 * Nx + 1;
    Ny = 2 * Ny + 1;
    Nz = 2 * Nz + 1;
    //
    int comm_size, comm_rank;
    int mpi_error = MPI_Comm_size(commK, &comm_size);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_rank(commK, &comm_rank);
    MPI_error_msg(mpi_error);

    // check if all cpus can accomidate atleast one gstick
    if (Nx * Ny < comm_size)
    {
        // the factor 1.01 is for numerical round off.
        double factor = sqrt((1.01 * (double)comm_size) / (Nx * Ny * 1.0));
        Nx = rint(1 + factor * Nx);
        Ny = rint(1 + factor * Ny);
        // make them odd
        if (0 == Nx % 2)
        {
            ++Nx;
        }
        if (0 == Ny % 2)
        {
            ++Ny;
        }
    }
    //
    fft_box[0] = Nx;
    fft_box[1] = Ny;
    fft_box[2] = Nz;

    // Make sure we have same number on all cpus. Does it happen?
    // Maybe who knows if we 4.499999999 and 4.5000000001
    // will yeild different numbers.
    mpi_error = MPI_Allreduce(MPI_IN_PLACE, fft_box, 3, ELPH_MPI_ND_INT,
                              MPI_MAX, commK);
    MPI_error_msg(mpi_error);

    return;
}
