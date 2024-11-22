// Funtions related to timings in the code
//
//

#include "ELPH_timers.h"

#include <mpi.h>
#include <stdio.h>
#include <string.h>

#include "ELPH_hash_map/ELPH_hmap.h"

struct ELPH_timer
{
    double Wtime;
    size_t count;
};

typedef map_t(struct ELPH_timer) map_timer_t;

static map_timer_t timer_map;
// Data structure to store hash table
//
//

void init_ELPH_clocks(void)
{
    // must be called to inititiate time clocks
    map_init(&timer_map);
    // Also at the same time initiate Total time
    ELPH_start_clock("Total time");
}

void ELPH_start_clock(const char *str)
{
    if (!str)
    {
        return;
    }
    double tic = MPI_Wtime();
    // check if this tag already exists
    struct ELPH_timer *etime = map_get(&timer_map, str);

    size_t count_tmp = 0;

    if (etime)
    {
        tic -= etime->Wtime;
        count_tmp = etime->count;
    }
    tic = -tic;

    struct ELPH_timer time_set;
    time_set.count = count_tmp;
    time_set.Wtime = tic;
    map_set(&timer_map, str, time_set);
}

void ELPH_stop_clock(const char *str)
{
    if (!str)
    {
        return;
    }
    double tok = MPI_Wtime();
    // check if this tag already exists
    struct ELPH_timer *etime = map_get(&timer_map, str);

    size_t count_tmp = 0;
    if (etime)
    {
        tok += etime->Wtime;
        count_tmp = etime->count + 1;
    }
    else
    {
        // clock does not exist
        return;
    }

    struct ELPH_timer time_set;
    time_set.count = count_tmp;
    time_set.Wtime = tok;

    map_set(&timer_map, str, time_set);
}

void cleanup_ELPH_clocks(void)
{
    // cleanup all the clock info.
    // nothing should be called after this.
    map_deinit(&timer_map);
}

void print_ELPH_clock_summary(void)
{
    // first end the total time
    ELPH_stop_clock("Total time");

    const char *key;
    map_iter_t iter = map_iter(&timer_map);
    fputs("\n", stdout);
    fputs("================== Wall times ==================\n", stdout);
    fprintf(stdout, "%-20s  : Wtime (s)    (  counts  )\n", "Function");
    fputs("------------------------------------------------\n", stdout);
    while ((key = map_next(&timer_map, &iter)))
    {
        if (0 == strcmp(key, "Total time"))
        {
            continue;
        }
        struct ELPH_timer *etime = map_get(&timer_map, key);
        fprintf(stdout, "%-20s  : %.6f     ( %8zu )\n", key, etime->Wtime,
                etime->count);
    }
    fputs("------------------------------------------------\n", stdout);
    // print total time
    key = "Total time";
    fprintf(stdout, "%-20s  : %.6f s.\n", key, map_get(&timer_map, key)->Wtime);
    fputs("\n", stdout);
}
