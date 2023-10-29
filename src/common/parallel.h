#pragma once 
#include "../elphC.h"

ND_int get_mpi_local_size_idx(const ND_int n, ND_int * start_idx,  MPI_Comm Comm);

void create_parallel_comms(int nqpools, int nkpools, MPI_Comm * commQ, MPI_Comm * commK);
void free_parallel_comms(MPI_Comm * commQ, MPI_Comm * commK);



