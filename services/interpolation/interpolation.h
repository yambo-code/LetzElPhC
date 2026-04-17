#pragma once

#include <mpi.h>

#include "common/dtypes.h"
#include "elphC.h"

void interpolation_driver(const char* ELPH_input_file,
                          enum ELPH_dft_code dft_code, MPI_Comm comm_world);
