#include "elph.h"
#include "../common/progess_bar.h"

/*
 * This function contain the wrapper functions to compute and
 * write or read dmat functions to file)
 */

void compute_and_write_dmats(const char* file_name, const struct WFC* wfcs,
                             const struct Lattice* lattice,
                             const ND_int nph_sym, const struct symmetry* sym_data,
                             const struct ELPH_MPI_Comms* Comm)
{
    ND_int nk_totalBZ = lattice->nkpts_BZ;
    ELPH_cmplx* Dkmn_rep_ptr = NULL;

    int ncid, varid, nc_err;

    size_t startp[6] = { 0, 0, 0, 0, 0, 0 };
    size_t countp[6] = { 1, 1, lattice->nspin, lattice->nbnds, lattice->nbnds, 2 };

    ND_int nk_chunk_size = NC4_DEFAULT_CHUCK_KB * 1024; // now this is in bytes
    // scale with complex number size to get the number of elements
    nk_chunk_size /= (sizeof(ELPH_cmplx) * lattice->nspin * lattice->nbnds * lattice->nbnds);
    // chuck the varaible elph_mat with atmost default size
    if (nk_chunk_size == 0)
    {
        nk_chunk_size = 1;
    }
    else if (nk_chunk_size > nk_totalBZ)
    {
        nk_chunk_size = nk_totalBZ;
    }

    if (Comm->commK_rank == 0)
    {
        // we overwrite any existing file
        if ((nc_err = nc_create_par(file_name, NC_NETCDF4 | NC_CLOBBER, Comm->commR,
                                    MPI_INFO_NULL, &ncid)))
        {
            fprintf(stderr, "Error creating Dmat file");
            ERR(nc_err);
        }
        // set no fill mode (to avoid writting twice)
        if ((nc_err = ncsetfill(ncid, NC_NOFILL)))
        {
            fprintf(stderr, "Error setting nc_fill to dmat file.");
            ERR(nc_err);
        }

        def_ncVar(
            ncid, &varid, 6, ELPH_NC4_IO_FLOAT,
            (ND_int[]) { nph_sym, nk_totalBZ, lattice->nspin, lattice->nbnds,
                         lattice->nbnds, 2 },
            "Dmats",
            (char*[]) { "nsym_ph", "nkpts", "nspin", "Rk_band", "k_band", "re_im" },
            (size_t[]) { 1, nk_chunk_size, lattice->nspin, lattice->nbnds, lattice->nbnds, 2 });

        // Make the access INDEPENDENT as not all can call the put_var function
        // simultaneously
        if ((nc_err = nc_var_par_access(ncid, varid, NC_INDEPENDENT)))
        {
            ERR(nc_err);
        }

        Dkmn_rep_ptr = calloc(lattice->nspin * lattice->nbnds * lattice->nbnds,
                              sizeof(ELPH_cmplx));
        CHECK_ALLOC(Dkmn_rep_ptr);
    }

    // for computation of Dmats, we use all the nodes
    // ("nsym", "nk", "nspin", "nbndb", "nbnda")
    ND_int dmat_shift;
    ND_int ndmats = distribute_to_grps(nph_sym * nk_totalBZ, Comm->nqpools * Comm->nkpools,
                                       Comm->commW_rank / Comm->commK_size, &dmat_shift);

    // start the progress bar for dmats
    struct progress_bar pbar[1];
    start_progressbar(pbar, Comm->commW_rank, ndmats);

    for (ND_int idmat = 0; idmat < ndmats; ++idmat)
    {
        ND_int isym = (idmat + dmat_shift) / nk_totalBZ;
        ND_int ikBZ = (idmat + dmat_shift) % nk_totalBZ;

        startp[0] = isym;
        startp[1] = ikBZ;

        // compute the dmats
        electronic_reps(wfcs, lattice, sym_data[isym].Rmat, sym_data[isym].tau,
                        sym_data[isym].time_rev, ikBZ, Dkmn_rep_ptr, Comm);

        if (Comm->commK_rank == 0)
        {
            // write data to file
            if ((nc_err = nc_put_vara(ncid, varid, startp, countp, Dkmn_rep_ptr)))
            {
                ERR(nc_err);
            }
        }

        // update the progress bar
        print_progressbar(pbar);
    }
    if (Comm->commK_rank == 0)
    {
        free(Dkmn_rep_ptr);

        if ((nc_err = nc_close(ncid)))
        {
            ERR(nc_err);
        }
    }
}
