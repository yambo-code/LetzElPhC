#pragma once 
#include "../elphC.h"



void create_parallel_comms(int nqpools, int nkpools, MPI_Comm * commQ, MPI_Comm * commK);
void free_parallel_comms(MPI_Comm * commQ, MPI_Comm * commK);



