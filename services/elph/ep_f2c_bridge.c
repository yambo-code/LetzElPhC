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
 *
 */

#include "elph.h"   /* already pulls in common/dtypes.h and elphC.h */
#include "common/print_info.h"
#include <elph/yambo.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include "parser/parser.h"

/*
 * Redirect stdout to the log file for the duration of the LetzElPhC call.
 */
static int saved_stdout_fd = -1;

static void open_letz_log(MPI_Comm comm, const char* log_path)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0 && log_path != NULL && log_path[0] != '\0')
    {
        FILE* fp = fopen(log_path, "w");
        if (fp)
        {
            elph_set_log_file(fp);
            fflush(stdout);
            saved_stdout_fd = dup(STDOUT_FILENO);
            dup2(fileno(fp), STDOUT_FILENO);
        }
    }
}

static void close_letz_log(void)
{
    FILE* fp = elph_get_log_file();
    if (fp != NULL)
    {
        fflush(fp);
        if (saved_stdout_fd >= 0)
        {
            dup2(saved_stdout_fd, STDOUT_FILENO);
            close(saved_stdout_fd);
            saved_stdout_fd = -1;
        }
        fclose(fp);
        elph_set_log_file(NULL);
    }
}

void elph_driver_f2c(const char* input_file, int dft_code, MPI_Fint f_comm,
                     const char* log_path)
{
    MPI_Comm c_comm = MPI_Comm_f2c(f_comm);
    open_letz_log(c_comm, log_path);
    elph_driver(input_file, (enum ELPH_dft_code)dft_code, c_comm);
    close_letz_log();
}

/*
 * Extended callback bridge: callbacks passed from Fortran to handle Linux linking.
 * Communicators (f_comm_q, f_comm_k) passed from Y6 PAR schemes for MPI distribution.
 */
void elph_driver_cb2_f2c(struct elph_usr_input* input_data, struct Y6_info* y6_data,
                         int dft_code, void* fill_fn_ptr, void* dvG_fill_fn_ptr,
                         const char* log_path, int i_control, struct Y6_parallel* y6_par,
                         MPI_Fint f_comm_world, MPI_Fint f_comm_q, MPI_Fint f_comm_k, int bz_mode_code)
{
    MPI_Comm c_comm_world = MPI_Comm_f2c(f_comm_world);
    MPI_Comm c_comm_q = MPI_Comm_f2c(f_comm_q);
    MPI_Comm c_comm_k = MPI_Comm_f2c(f_comm_k);

    open_letz_log(c_comm_world, log_path);

    /* Y6 mode: use provided communicators and callbacks from Fortran */
    elph_driver_cb2(input_data,y6_data,y6_par,(enum ELPH_dft_code)dft_code,
                    (elph_gkkp_fill_fn)fill_fn_ptr,
                    (elph_dvG_fill_fn)dvG_fill_fn_ptr,i_control,c_comm_world,bz_mode_code);

    close_letz_log();
}
