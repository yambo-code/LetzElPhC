/*
This file contains function which distributes cpus 
*/
#include "parallel.h"


/*get block size and starting idx of dimension that is distrbuted amoung cpus*/
ND_int get_mpi_local_size_idx(const ND_int n, ND_int * start_idx,  MPI_Comm Comm)
{   
    int my_rank, total_size;

    MPI_Comm_rank(Comm, &my_rank);
    MPI_Comm_size(Comm, &total_size);

    ND_int n_q = n/total_size ;
    ND_int n_r = n%total_size ;
    
    ND_int n_this_cpu = n_q;
    if (my_rank < n_r) ++n_this_cpu;

    if (start_idx != NULL)
    {
        *start_idx = n_q*my_rank;
        if (my_rank < n_r) *start_idx += my_rank;
        else *start_idx += n_r;
    }
    return n_this_cpu;
}



/* allocate comms */
void create_parallel_comms(int nqpools, int nkpools, MPI_Comm * commQ, MPI_Comm * commK)
{
    /*
    !! Warning : Order of arrangement of cpus is very important. the key given to 
    comm_split must be the rank of world comm
    */
    int my_rank, total_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_size);

    // check is total_size is divisible by nqpools*nkpools
    if (nqpools  == 0 || nkpools == 0) error_msg(" number of pools must be greate than 0 ");
    if (total_size%(nqpools*nkpools) != 0 ) error_msg(" product of kpools and q pools must divide total cpus ");

    int qsize = total_size/nqpools ;
    int ksize = qsize/nkpools;

    /* first make nqpools groups */
    MPI_Comm_split(MPI_COMM_WORLD, my_rank/qsize, my_rank, commQ);

    int krank;
    MPI_Comm_rank(*commQ, &krank);
    // now create k pools
    MPI_Comm_split(*commQ, krank/ksize, my_rank, commK);

}

/* free allocated comms */
void free_parallel_comms(MPI_Comm * commQ, MPI_Comm * commK)
{
    MPI_Comm_free(commK);
    MPI_Comm_free(commQ);
}
