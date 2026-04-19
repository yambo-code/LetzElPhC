/*
 * License-Identifier: GPL
 *
 * Copyright (C) 2025 The Yambo Team
 *
 * Authors (see AUTHORS file for details): RR AM
 *
 * Thin bridge: called from Fortran via ISO_C_BINDING.
 * Converts a Fortran MPI communicator handle (MPI_Fint) to a C MPI_Comm
 * and forwards the call to LetzElPhC's elph_driver().
 */

#include "elph.h"   /* already pulls in common/dtypes.h and elphC.h */
#include <mpi.h>

void elph_driver_f2c(const char* input_file, int dft_code, MPI_Fint f_comm)
{
    MPI_Comm c_comm = MPI_Comm_f2c(f_comm);
    elph_driver(input_file, (enum ELPH_dft_code)dft_code, c_comm);
}
