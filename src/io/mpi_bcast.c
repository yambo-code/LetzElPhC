#include <limits.h>

#include "io.h"

void Bcast_wfc(struct WFC* wfc, const struct Lattice* lattice,
               const struct Pseudo* pseudo, bool alloc_mem, int root,
               MPI_Comm comm)
{
    // lattice and pseudo must be on all process
    int mpi_error;
    int my_rank;

    mpi_error = MPI_Comm_rank(comm, &my_rank);

    // npw_total
    mpi_error = MPI_Bcast(&(wfc->npw_total), 1, ELPH_MPI_ND_INT, root, comm);
    // npw_loc
    mpi_error = MPI_Bcast(&(wfc->npw_loc), 1, ELPH_MPI_ND_INT, root, comm);
    /*
    Bcast wfc
    */
    int nspin = lattice->nspin;
    int nbnds = lattice->nbnds;
    int nspinor = lattice->nspinor;

    ND_int nltimesj = pseudo->nltimesj;
    ND_int ntype = pseudo->ntype;

    ND_int wfc_size = nspin * nbnds * nspinor * wfc->npw_loc;
    ND_int gvec_size = 3 * wfc->npw_loc;
    ND_int Fk_size = nltimesj * ntype * wfc->npw_loc;

    if (alloc_mem && (my_rank != root))
    {
        wfc->wfc = malloc(sizeof(ELPH_cmplx) * wfc_size);
        CHECK_ALLOC(wfc->wfc);

        wfc->gvec = malloc(sizeof(ELPH_float) * gvec_size);
        CHECK_ALLOC(wfc->gvec);

        wfc->Fk = malloc(sizeof(ELPH_float) * Fk_size);
        CHECK_ALLOC(wfc->Fk);
    }

    mpi_error = MPI_Bcast(wfc->wfc, wfc_size, ELPH_MPI_cmplx, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(wfc->gvec, gvec_size, ELPH_MPI_float, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(wfc->Fk, Fk_size, ELPH_MPI_float, root, comm);
    MPI_error_msg(mpi_error);
}

void Bcast_symmetries(ND_int nsyms, struct symmetry* symms, int root,
                      MPI_Comm comm)
{
    // This function will Bcast struct symmetrys
    // Note this will not create any data so pass allocated array

    int mpi_error;
    struct symmetry dummy;

    int lengths[4] = { 9, 3, 1, 1 };
    MPI_Aint displacements[4], base_address;

    // get address of members
    mpi_error = MPI_Get_address(&dummy, &base_address);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Get_address(dummy.Rmat, &displacements[0]); // 9
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Get_address(dummy.tau, &displacements[1]); // 3
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Get_address(&dummy.inv_idx, &displacements[2]); // 1
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Get_address(&dummy.time_rev, &displacements[3]); // 1
    MPI_error_msg(mpi_error);

    for (int i = 0; i < 4; ++i)
    {
        displacements[i] = MPI_Aint_diff(displacements[i], base_address);
    }

    MPI_Datatype types[4] = { ELPH_MPI_float, ELPH_MPI_float, ELPH_MPI_ND_INT,
                              MPI_C_BOOL };

    MPI_Datatype symm_type;
    mpi_error = MPI_Type_create_struct(4, lengths, displacements, types, &symm_type);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Type_commit(&symm_type);
    MPI_error_msg(mpi_error);

    // Bcast
    mpi_error = MPI_Bcast(symms, nsyms, symm_type, root, comm);
    MPI_error_msg(mpi_error);

    // free the type
    mpi_error = MPI_Type_free(&symm_type);
    MPI_error_msg(mpi_error);
}

void Bcast_local_pseudo(struct local_pseudo* loc_pseudo, bool alloc_mem,
                        int root, MPI_Comm comm)
{
    // First get the ngrid
    int mpi_error;
    int my_rank;

    mpi_error = MPI_Comm_rank(comm, &my_rank);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&loc_pseudo->ngrid, 1, ELPH_MPI_ND_INT, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&loc_pseudo->Zval, 1, ELPH_MPI_float, root, comm);
    MPI_error_msg(mpi_error);

    if (alloc_mem && (my_rank != root))
    {
        loc_pseudo->Vloc_atomic = malloc(sizeof(ELPH_float) * loc_pseudo->ngrid);
        CHECK_ALLOC(loc_pseudo->Vloc_atomic);

        loc_pseudo->r_grid = malloc(sizeof(ELPH_float) * loc_pseudo->ngrid);
        CHECK_ALLOC(loc_pseudo->r_grid);

        loc_pseudo->rab_grid = malloc(sizeof(ELPH_float) * loc_pseudo->ngrid);
        CHECK_ALLOC(loc_pseudo->rab_grid);
    }

    mpi_error = MPI_Bcast(loc_pseudo->Vloc_atomic, loc_pseudo->ngrid,
                          ELPH_MPI_float, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(loc_pseudo->r_grid, loc_pseudo->ngrid, ELPH_MPI_float,
                          root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(loc_pseudo->rab_grid, loc_pseudo->ngrid,
                          ELPH_MPI_float, root, comm);
    MPI_error_msg(mpi_error);
}
