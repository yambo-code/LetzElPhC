#include "progess_bar.h"
#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

void start_progressbar(struct progress_bar* pbar, int mpi_rank, ND_int niter)
{
    // initialize the values
    pbar->rank = mpi_rank;
    pbar->iiter = 0;
    pbar->niter = niter;
    pbar->isign = 0;
    pbar->prev_time = MPI_Wtime();
    pbar->elap_time = 0.0;
    pbar->max_iter_time = 0.0;

    if (pbar->rank)
    {
        return;
    }

    // print out empty progress bar
    time_t result = time(NULL);
    if (result != (time_t)(-1))
    {
        const char* dty_str = asctime(localtime(&result));
        int nnew_line = strcspn(dty_str, "\r\n");
        fprintf(stdout, "%.*s : [", nnew_line, dty_str);
    }

    for (int i = 0; i < PROGRESS_BAR_LEN; ++i)
    {
        fprintf(stdout, PROGRESS_BAR_SPACE);
    }
    fprintf(stdout, "] %6.2f%% \n", 0.0);
    fflush(stdout);
}

void print_progressbar(struct progress_bar* pbar)
{
    // only pbar->rank ==0 process prints the bar
    if (pbar->rank)
    {
        return;
    }
    ++pbar->iiter;

    float progress = ((float)pbar->iiter) / pbar->niter;

    int l_bar = progress * PROGRESS_BAR_LEN;
    // length to the current iteration bar

    if (l_bar == pbar->isign)
    {
        return;
    }
    else
    {
        ++pbar->isign;
    }

    time_t result = time(NULL);
    if (result != (time_t)(-1))
    {
        const char* dty_str = asctime(localtime(&result));
        int nnew_line = strcspn(dty_str, "\r\n");
        fprintf(stdout, "%.*s : [", nnew_line, dty_str);
    }
    else
    {
        fprintf(stdout, "Progress at XXX XXX XX XX:XX:XX XXXX : [");
    }

    for (int i = 0; i < l_bar; ++i)
    {
        fprintf(stdout, PROGRESS_BAR_SIGN);
    }
    for (int i = l_bar; i < PROGRESS_BAR_LEN; ++i)
    {
        fprintf(stdout, PROGRESS_BAR_SPACE);
    }

    double current_time = MPI_Wtime();
    double iter_time = current_time - pbar->prev_time;
    if (iter_time > pbar->max_iter_time)
    {
        pbar->max_iter_time = iter_time;
    }
    pbar->elap_time += iter_time;
    pbar->prev_time = current_time;

    int rem_time = pbar->max_iter_time * PROGRESS_BAR_LEN - pbar->elap_time;

    if (pbar->iiter == pbar->niter)
    {
        rem_time = 0;
    }

    fprintf(stdout, "] %6.2f%%  ETC: %4dh %2dm %2ds\n",
            progress * 100, rem_time / 3600,
            (rem_time % 3600) / 60, (rem_time % 3600) % 60);

    if (pbar->iiter == pbar->niter)
    {
        fprintf(stdout, "\n");
    }
    fflush(stdout);
}
