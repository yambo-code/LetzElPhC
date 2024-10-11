#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include "../common/dtypes.h"
#include "../common/error.h"
#include "../common/numerical_func.h"
#include "../common/parallel.h"
#include "../elphC.h"
#include "../wfc/wfc.h"
#include "dvloc.h"

/* create a interpolation table on a coarse grid for short local potential in
 * Gspace*/
void create_vlocg_table(const struct Lattice* lattice, struct Pseudo* pseudo,
                        const struct ELPH_MPI_Comms* Comm)
{
    int mpi_error;
    // pseudo->vloc_table->qmax_abs
    const struct local_pseudo* loc_pseudo = pseudo->loc_pseudo;
    const ND_int ntype = pseudo->ntype;
    const ND_int ngrid_max = pseudo->ngrid_max;
    const ELPH_float volume = lattice->volume;

    ELPH_float gmax, gmin;
    ND_int npts;

    // create a small scope
    {
        const ELPH_float* blat = lattice->blat_vec;  // 2*pi is included

        ELPH_float bi[3], bmin, bmax, Nmax;

        bi[0] = sqrt((blat[0] * blat[0]) + (blat[3] * blat[3]) +
                     (blat[6] * blat[6]));
        bi[1] = sqrt((blat[1] * blat[1]) + (blat[4] * blat[4]) +
                     (blat[7] * blat[7]));
        bi[2] = sqrt((blat[2] * blat[2]) + (blat[5] * blat[5]) +
                     (blat[8] * blat[8]));
        bmin = bi[0];
        bmax = bi[0];
        Nmax = lattice->fft_dims[0];

        for (int i = 0; i < 3; ++i)
        {
            if (bi[i] < bmin)
            {
                bmin = bi[i];
            }
            if (bi[i] > bmax)
            {
                bmax = bi[i];
            }
            if (lattice->fft_dims[i] > Nmax)
            {
                Nmax = lattice->fft_dims[i];
            }
        }

        /*
        We try to do this only once for any q, so we should get the max possible
        value of |G+q| |G+q| <= |G| + |q|. If all three values of q (in crystal
        coordinates) are restricted to [qmin,qmax] then |q| <= 3*b_max*|qmax|
        G_max = |N1/2*b1 + N2/2*b2 + N3/2*b3| < N_max/2 * 3 * b_max
        |G|<= 1.5*N_max*b_max
        => |G+q| <= b_max*(1.5*N_max+3*|qmax|)
        So we choose Interpolation range to be [0,(1.5*N_max+3*|qmax|)*b_max]
        // coarse spacing = b_min, npts ~ (1.5*N_max + 3*|qmax|)*b_max/b_min.
        */

        if (bmin > 0.001)
        {
            bmin = 0.001;  // set some bare minimum
        }

        gmax = (1.7 * (Nmax + 4.0) + 3 * fabs(pseudo->vloc_table->qmax_abs)) *
               bmax;
        // we choose 1.7 instead of 1.5
        npts = ceil(gmax / bmin);
        gmin = 0.0;
    }  // end of scope

    ELPH_float* xins = malloc(sizeof(ELPH_float) * npts);
    CHECK_ALLOC(xins);

    ELPH_float* yins = malloc(sizeof(ELPH_float) * npts * ntype);
    CHECK_ALLOC(yins);

    pseudo->vloc_table->vploc_co = malloc(sizeof(ELPH_float) * npts * ntype);
    CHECK_ALLOC(pseudo->vloc_table->vploc_co);

    pseudo->vloc_table->npts_co = npts;
    pseudo->vloc_table->g_co = xins;
    pseudo->vloc_table->vlocg = yins;

    ELPH_float diff = (gmax - gmin) / (npts - 1.0);
    for (ND_int i = 0; i < npts; ++i)
    {
        xins[i] = gmin + i * diff;
    }
    pseudo->vloc_table->dg = diff;

    ND_int n_shift;
    int npts_loc = get_mpi_local_size_idx(npts, &n_shift, Comm->commW);

    ELPH_float* work_array =
        malloc(sizeof(ELPH_float) * (ngrid_max + npts_loc + 1));
    CHECK_ALLOC(work_array);

    ELPH_float* vlocg_cpu = work_array + ngrid_max;

    int* counts_recv = malloc(sizeof(int) * 2 * Comm->commW_size);
    CHECK_ALLOC(counts_recv);

    int* displacements = counts_recv + Comm->commW_size;

    mpi_error = MPI_Allgather(&npts_loc, 1, MPI_INT, counts_recv, 1, MPI_INT,
                              Comm->commW);
    MPI_error_msg(mpi_error);

    int n_shift_temp = n_shift;
    mpi_error = MPI_Allgather(&n_shift_temp, 1, MPI_INT, displacements, 1,
                              MPI_INT, Comm->commW);
    MPI_error_msg(mpi_error);

    for (ND_int itype = 0; itype < ntype; ++itype)
    {
        ELPH_float* vlocg_atom = yins + itype * npts;
        ELPH_float* dyy = pseudo->vloc_table->vploc_co + itype * npts;

        const struct local_pseudo* loc_pseudo_type = loc_pseudo + itype;
        // loop not thread safe
        for (ND_int i = 0; i < npts_loc; ++i)
        {
            vlocg_cpu[i] = Vloc_Gspace(
                work_array, lattice->dimension, xins[i + n_shift],
                loc_pseudo_type->ngrid, loc_pseudo_type->Vloc_atomic,
                loc_pseudo_type->r_grid, loc_pseudo_type->rab_grid,
                loc_pseudo_type->Zval, 1.0, volume);
        }
        mpi_error = MPI_Allgatherv(vlocg_cpu, npts_loc, ELPH_MPI_float,
                                   vlocg_atom, counts_recv, displacements,
                                   ELPH_MPI_float, Comm->commW);
        MPI_error_msg(mpi_error);

        prepare_spline(npts, xins, vlocg_atom, dyy);
    }

    free(counts_recv);
    free(work_array);
}

void free_vlocg_table(struct Vloc_table* vloc_table)
{
    // free the allocate Vloc_table
    free(vloc_table->g_co);
    free(vloc_table->vlocg);
    free(vloc_table->vploc_co);
    vloc_table->g_co = NULL;
    vloc_table->vlocg = NULL;
    vloc_table->vploc_co = NULL;
}
