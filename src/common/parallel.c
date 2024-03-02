/*
This file contains function which distributes cpus
*/
#include "parallel.h"

/*get block size and starting idx of dimension that is distrbuted amoung cpus*/
ND_int get_mpi_local_size_idx(const ND_int n, ND_int* start_idx, MPI_Comm Comm)
{
    int my_rank, total_size;
    int mpi_error;
    mpi_error = MPI_Comm_rank(Comm, &my_rank);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_size(Comm, &total_size);
    MPI_error_msg(mpi_error);

    ND_int n_q = n / total_size;
    ND_int n_r = n % total_size;

    ND_int n_this_cpu = n_q;
    if (my_rank < n_r)
    {
        ++n_this_cpu;
    }

    if (start_idx != NULL)
    {
        *start_idx = n_q * my_rank;
        if (my_rank < n_r)
        {
            *start_idx += my_rank;
        }
        else
        {
            *start_idx += n_r;
        }
    }
    return n_this_cpu;
}

ND_int distribute_to_grps(const ND_int n, const ND_int ngrp, const ND_int igrp,
                          ND_int* shift)
{
    /*
    This function is used to distribute k/q points between  n groups
    return k/q for each group, and the shift
    */

    ND_int n_q = n / ngrp;
    ND_int n_r = n % ngrp;

    ND_int n_this_grp = n_q;
    if (igrp < n_r)
    {
        ++n_this_grp;
    }

    *shift = n_q * igrp;
    if (igrp < n_r)
    {
        *shift += igrp;
    }
    else
    {
        *shift += n_r;
    }

    return n_this_grp;
}

/* allocate comms */
void create_parallel_comms(const int nqpools, const int nkpools,
                           const MPI_Comm MPI_world_comm,
                           struct ELPH_MPI_Comms* Comm)
{
    /*
    !! Warning : Order of arrangement of cpus is very important. the key given
    to comm_split must be the rank of world comm
    */
    int mpi_error;
    // first set the basic things
    Comm->commW = MPI_world_comm;
    Comm->nqpools = nqpools;
    Comm->nkpools = nkpools;

    mpi_error = MPI_Comm_rank(Comm->commW, &Comm->commW_rank);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_size(Comm->commW, &Comm->commW_size);
    MPI_error_msg(mpi_error);

    // check Comm->commW_size is divisible by nqpools*nkpools
    if (nqpools < 1 || nkpools < 1)
    {
        error_msg(" number of pools must be greater than 0 ");
    }
    if (Comm->commW_size % (nqpools * nkpools) != 0)
    {
        error_msg(" product of kpools and qpools must divide total cpus.");
    }

    Comm->commQ_size = Comm->commW_size / nqpools;
    Comm->commK_size = Comm->commQ_size / nkpools;

    /* create comms */
    // commQ
    mpi_error = MPI_Comm_split(Comm->commW, Comm->commW_rank / Comm->commQ_size,
                               Comm->commW_rank, &Comm->commQ);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_rank(Comm->commQ, &Comm->commQ_rank);
    MPI_error_msg(mpi_error);

    // commK
    mpi_error = MPI_Comm_split(Comm->commQ, Comm->commQ_rank / Comm->commK_size,
                               Comm->commW_rank, &Comm->commK);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_rank(Comm->commK, &Comm->commK_rank);
    MPI_error_msg(mpi_error);
    // commR
    mpi_error = MPI_Comm_split(Comm->commW, Comm->commK_rank, Comm->commW_rank,
                               &Comm->commR);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_rank(Comm->commR, &Comm->commR_rank);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_size(Comm->commR, &Comm->commR_size);
    MPI_error_msg(mpi_error);
    // commRq
    mpi_error = MPI_Comm_split(Comm->commQ, Comm->commK_rank, Comm->commW_rank,
                               &Comm->commRq);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_rank(Comm->commRq, &Comm->commRq_rank);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_size(Comm->commRq, &Comm->commRq_size);
    MPI_error_msg(mpi_error);

    // Sanity check
    if (Comm->commK_rank != Comm->commW_rank % Comm->commK_size)
    {
        error_msg("MPI Comm 1 creation failed.");
    }
    if (Comm->commR_size != (nqpools * nkpools))
    {
        error_msg("MPI Comm 2 creation failed.");
    }
    if (Comm->commRq_size != nkpools)
    {
        error_msg("MPI Comm 3 creation failed.");
    }
}

/* free allocated comms */
void free_parallel_comms(struct ELPH_MPI_Comms* Comm)
{
    int mpi_error;
    mpi_error = MPI_Comm_free(&Comm->commQ);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_free(&Comm->commK);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_free(&Comm->commR);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_free(&Comm->commRq);
    MPI_error_msg(mpi_error);
}
