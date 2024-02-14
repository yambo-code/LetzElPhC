#pragma once
#include "../common/dtypes.h"
#include "../common/numerical_func.h"
#include "../common/parallel.h"
#include "../common/string_func.h"
#include "../elphC.h"
#include "../symmetries/symmetries.h"

// mpi_bcast.c
// Bcast functions
void Bcast_wfc(struct WFC* wfc, const struct Lattice* lattice,
               const struct Pseudo* pseudo, bool alloc_mem, int root,
               MPI_Comm comm);

void Bcast_symmetries(ND_int nsyms, struct symmetry* symms, int root,
                      MPI_Comm comm);

void Bcast_local_pseudo(struct local_pseudo* loc_pseudo, bool alloc_mem,
                        int root, MPI_Comm comm);
