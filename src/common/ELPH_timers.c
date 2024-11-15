// Funtions related to timings in the code
//
//

#include "ELPH_timers.h"

#include <mpi.h>
#include <stdio.h>

#include "ELPH_hash_map/ELPH_hmap.h"

static map_double_t timer_map;
// Data structure to store hash table
//
//

void init_ELPH_clocks(void)
{
    // must be called to inititiate time clocks
    map_init(&timer_map);
    // Also at the same time initiate Total time
    //
    /* double tic = MPI_Wtime(); */
    /* map_set(&timer_map, "Total_Wtime", tic); */
}

void ELPH_start_clock(const char *str)
{
    if (!str)
    {
        return;
    }
    double tic = MPI_Wtime();
    // check if this tag already exists
    double *etime = map_get(&timer_map, str);

    if (etime)
    {
        tic -= *etime;
    }
    tic = -tic;
    map_set(&timer_map, str, tic);
}

void ELPH_stop_clock(const char *str)
{
    if (!str)
    {
        return;
    }
    double tok = MPI_Wtime();
    // check if this tag already exists
    double *etime = map_get(&timer_map, str);

    if (etime)
    {
        tok += *etime;
    }
    else
    {
        // clock does not exist
        return;
    }
    map_set(&timer_map, str, tok);
}

void cleanup_ELPH_clocks(void)
{
    // cleanup all the clock info.
    // nothing should be called after this.
    map_deinit(&timer_map);
}

void print_ELPH_clock_summary(void)
{
    const char *key;
    map_iter_t iter = map_iter(&timer_map);
    fputs("\n", stdout);
    fputs("============== Wall times ==============\n", stdout);

    while ((key = map_next(&timer_map, &iter)))
    {
        fprintf(stdout, "%-20s  : %.4f s.\n", key, *map_get(&timer_map, key));
    }
    fputs("\n", stdout);
}
