/* This file contains functions common for all types*/
#include "nd_error.h"

#if defined(COMPILE_ONCE)
#ifndef WITH_MPI_ERROR 
void nd_error_msg(const char * error_msg, const char * func_name){
    fprintf(stderr, "*************************************\n");
    fprintf(stderr, "# [ Error !!!] :  In function : %s , Error msg : %s \n", func_name, error_msg);
    fprintf(stderr, "*************************************\n");
    exit(EXIT_FAILURE);
}
#else // MPI
#include <mpi.h>

void nd_error_msg(const char * error_msg, const char * func_name){
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    fprintf(stderr, "*************************************\n");
    fprintf(stderr, "# [ Error !!!] :  In function : %s , Error msg : %s \n", func_name, error_msg);
    fprintf(stderr, "*************************************\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
}
#endif // MPI
#endif //COMPILE_ONCE
