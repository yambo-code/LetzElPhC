/* This file contains functions common for all types*/
#include "error.h"

#include <mpi.h>

void elph_error_msg(const char* error_msg, const char* file, const char* line,
                    const char* func_name)
{
    int my_rank, size_comm;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    fprintf(stderr, "*************************************\n");
    fprintf(stderr,
            "# [ Error !!!] :  File : %s, in function : %s at line : %s \n"
            "Error msg : %s \n",
            file, func_name, line, error_msg);
    fprintf(stderr, "*************************************\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

void MPI_error_msg(int err_code, const char* file, const char* line,
                   const char* func_name)
{
    char message[MPI_MAX_ERROR_STRING + 1];
    int resultlen;
    MPI_Error_string(err_code, message, &resultlen);
    message[resultlen] = '\0';
    elph_error_msg(message, file, line, func_name);
}
