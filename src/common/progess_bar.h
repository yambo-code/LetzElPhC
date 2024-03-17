#pragma once
#include "../elphC.h"

#define PROGRESS_BAR_LEN 5
#define PROGRESS_BAR_SIGN "#####"
#define PROGRESS_BAR_SPACE "     "

struct progress_bar
{
    int rank; // only pbar->rank ==0 process prints the bar
    ND_int iiter; // current iteration
    ND_int niter; // number of iterations
    ND_int isign; // number of signs printed
    double prev_time; // time when previous iter was finished
    double elap_time; // total time elasped
    double max_iter_time; // max time taken by one iteration
};

void start_progressbar(struct progress_bar* pbar,
                       int mpi_rank, ND_int niter);

void print_progressbar(struct progress_bar* pbar);
