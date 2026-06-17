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

/* Fortran bind(C) callbacks — linked from the ep packages library. */
extern void elph_coll_fill_gkkp(int, int, const void*, int, int, int, int, int, int, int, const void*);
extern void elph_coll_fill_dvg(int, const void*, const void*, int, int, int, int, int, int);

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
 * Callback-enabled bridge: fill_fn_ptr argument is kept for API compatibility
 * but ignored — elph_coll_fill_gkkp is resolved directly by the linker.
 */
void elph_driver_cb_f2c(const char* input_file, int dft_code, MPI_Fint f_comm,
                         void* fill_fn_ptr, const char* log_path)
{
    MPI_Comm c_comm = MPI_Comm_f2c(f_comm);
    open_letz_log(c_comm, log_path);
    elph_driver_cb(input_file, (enum ELPH_dft_code)dft_code, c_comm,
                   (elph_fill_fn)elph_coll_fill_gkkp);
    close_letz_log();
}

/*
 * Extended callback bridge: callbacks passed from Fortran to handle Linux linking.
 * Communicators (f_comm_q, f_comm_k) passed from Y6 PAR schemes for MPI distribution.
 */
void elph_driver_cb2_f2c(struct elph_usr_input* input_data, struct Y6_info* y6_data,
                         int dft_code, void* fill_fn_ptr, void* dvG_fill_fn_ptr,
                         const char* log_path, int i_control,
                         int NQ_todo, int* Q_todo, int NK_todo , int* K_todo,
                         MPI_Fint f_comm_world, MPI_Fint f_comm_q, MPI_Fint f_comm_k)
{
    MPI_Comm c_comm_world = MPI_Comm_f2c(f_comm_world);
    MPI_Comm c_comm_q = MPI_Comm_f2c(f_comm_q);
    MPI_Comm c_comm_k = MPI_Comm_f2c(f_comm_k);

    open_letz_log(c_comm_world, log_path);

    struct Y6_parallel_work* y6_work = malloc(sizeof(struct Y6_parallel_work));

    y6_work->Q  = malloc(NQ_todo * sizeof(int));
    y6_work->NQ = NQ_todo;
    for (int i = 0; i < y6_work->NQ; i++) {
      y6_work->Q[i]  = Q_todo[i]-1;
    }

    y6_work->K  = malloc(NK_todo * sizeof(int));
    y6_work->NK = NK_todo;
    for (int i = 0; i < y6_work->NK; i++) {
      y6_work->K[i]  = K_todo[i]-1;
    }

    /* Y6 mode: use provided communicators and callbacks from Fortran */
    elph_driver_cb2(input_data,y6_data,y6_work,(enum ELPH_dft_code)dft_code,
                    (elph_fill_fn)fill_fn_ptr,
                    (elph_dvG_fill_fn)dvG_fill_fn_ptr,i_control,c_comm_world);

    close_letz_log();
}
