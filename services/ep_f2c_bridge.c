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

/*
 * Callback-enabled bridge: fill_fn_ptr is a Fortran procedure pointer
 * (passed as a C function pointer via ISO_C_BINDING's c_funloc).
 */
void elph_driver_cb_f2c(const char* input_file, int dft_code, MPI_Fint f_comm,
                         void* fill_fn_ptr)
{
    MPI_Comm c_comm = MPI_Comm_f2c(f_comm);
    elph_driver_cb(input_file, (enum ELPH_dft_code)dft_code, c_comm,
                   (elph_fill_fn)fill_fn_ptr);
}

/*
 * Extended callback bridge: gkkp fill_fn + dvG fill callback.
 * Either pointer may be NULL (pass c_null_funptr from Fortran) to skip that output.
 */
void elph_driver_cb2_f2c(const char* input_file, int dft_code, MPI_Fint f_comm,
                          void* fill_fn_ptr, void* dvG_fill_fn_ptr)
{
    MPI_Comm c_comm = MPI_Comm_f2c(f_comm);
    elph_driver_cb2(input_file, (enum ELPH_dft_code)dft_code, c_comm,
                    (elph_fill_fn)fill_fn_ptr,
                    (elph_dvG_fill_fn)dvG_fill_fn_ptr);
}
