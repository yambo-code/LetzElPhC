#pragma once 
#include "../elphC.h"

struct ELPH_MPI_Comms
{   
    /* This struct contains all the information about communications 
    used in the code*/
    MPI_Comm commW; 
    // always set to MPI_COMM_"W"ORLD (do not free!)
    MPI_Comm commQ; 
    // COMM for "Q"pool
    MPI_Comm commK; 
    // COMM for "K"pool
    MPI_Comm commR; 
    // equal "R"ank cpus from all kpools in COMM_WORLD form commR 
    MPI_Comm commRq; 
    // equal "R"ank cpus from all kpools in comm"Q" form commRk 
    int nqpools; // total number of q pools
    int nkpools;  // total number of k pools
    // ranks of comms
    int commW_rank;
    int commQ_rank;
    int commK_rank; 
    int commR_rank;
    int commRq_rank;
    // size of comm
    int commW_size;
    int commQ_size;
    int commK_size;
    int commR_size;
    int commRq_size;
    /*
    for Example if there are total 8 cpus in MPI_COMM_WORLD
    divide as 2 qpools and 2 kpools
    
                            Global {0,1,2,3,4,5,6,7}
            _________________|__________________
            |                                   |
    {0,1,2,3} (commQ)                   {4,5,6,7} (commQ)
    _________|__________                _________|__________  
    |                   |               |                   |
    {0,1}(commK)     {2,3}(commK)     {4,5} (commK)       {6,7} (commK)

    npw_cpus = size of commK
    commR contains  {0,2,4,6} and {1,3,5,7} i.e rank%npw_cpus is used to make the 
    division
    commRq contains {0,2} , {1,3}, in 1st commQ and {4,5} and {5,7} in 2nd commQ
    */
};


ND_int get_mpi_local_size_idx(const ND_int n, ND_int * start_idx,  MPI_Comm Comm);

void create_parallel_comms(const int nqpools, const int nkpools, \
        const MPI_Comm MPI_world_comm, struct ELPH_MPI_Comms * Comm);

void free_parallel_comms(struct ELPH_MPI_Comms * Comm);

ND_int distribute_to_grps(const ND_int n, const ND_int ngrp, const ND_int igrp, \
                        ND_int * shift);



