#include <ctype.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#include "common/cwalk/cwalk.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/parallel.h"
#include "common/string_func.h"
#include "elphC.h"
#include "io/ezxml/ezxml.h"
#include "io/io.h"
#include "qe_io.h"

void get_interpolation_data_from_qe(struct Lattice* lattice,
                                    struct Phonon* phonon,
                                    struct Pseudo* pseudo,
                                    const char* ph_save_dir, ELPH_float** Zvals,
                                    ELPH_float* alat,
                                    const struct ELPH_MPI_Comms* Comm)
{
    /*
       This functions gets basic ground state data to run interpolation. We
       donot have yambo SAVE, so we need to get few more quantities than in
       get_data_from_qe.c

       allocated memory must be freed outside this function
       */
    //
    get_data_from_qe(lattice, phonon, pseudo, ph_save_dir, alat, Comm);
    // we need to set set nnfftz_loc and nfftz_loc_shift

    *Zvals = NULL;
    ELPH_float* Zvalance = malloc(sizeof(*Zvalance) * lattice->natom);
    CHECK_ALLOC(Zvalance);
    *Zvals = Zvalance;
    //
    for (ND_int ia = 0; ia < lattice->natom; ++ia)
    {
        ND_int itype = lattice->atom_type[ia];
        Zvalance[ia] = pseudo->loc_pseudo[itype].Zval;
    }

    lattice->nfftz_loc = get_mpi_local_size_idx(
        lattice->fft_dims[2], &(lattice->nfftz_loc_shift), Comm->commK);

    if (lattice->nfftz_loc < 1)
    {
        error_msg(
            "Some cpus do not contain plane waves. Over parallelization !.");
    }
}
