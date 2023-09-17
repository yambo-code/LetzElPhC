/*
This file contains function which distributes cpus 
*/
#include "parallel.h"

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