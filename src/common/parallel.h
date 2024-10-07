#pragma once
#include "../elphC.h"
#include "dtypes.h"
#include <mpi.h>

ND_int get_mpi_local_size_idx(const ND_int n, ND_int* start_idx, MPI_Comm Comm);

void create_parallel_comms(const int nqpools, const int nkpools,
                           const MPI_Comm MPI_world_comm,
                           struct ELPH_MPI_Comms* Comm);

void free_parallel_comms(struct ELPH_MPI_Comms* Comm);

ND_int distribute_to_grps(const ND_int n, const ND_int ngrp, const ND_int igrp,
                          ND_int* shift);
