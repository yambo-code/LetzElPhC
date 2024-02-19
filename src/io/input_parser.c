#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "inih/ini.h"
#include "io.h"

#define MAX_STR_INPUT 1400

static void Bcast_input_data(struct usr_input* input, int root, MPI_Comm comm);
static int handler(void* user, const char* section, const char* name,
                   const char* value);

// function to alloc, initiate usr_input
void init_usr_input(struct usr_input** input)
{
    // this function also sets defaults for the user input file
    *input = malloc(sizeof(struct usr_input));
    struct usr_input* inp = *input;

    inp->save_dir = calloc(MAX_STR_INPUT, 1);
    inp->ph_save_dir = inp->save_dir + 600;
    inp->kernel = inp->save_dir + 1200;

    // defaults
    inp->nkpool = 1;
    inp->nqpool = 1;
    inp->start_bnd = 0;
    inp->end_bnd = 0;
    strcpy(inp->save_dir, "./SAVE");
    strcpy(inp->ph_save_dir, "./ph_save");
    strcpy(inp->kernel, "dfpt");
}

// function to free usr_input struct data
void free_usr_input(struct usr_input* input)
{
    free(input->save_dir);
    free(input);
}

static void Bcast_input_data(struct usr_input* input, int root, MPI_Comm comm)
{
    int mpi_error;

    // all char * will be bcasted in one single call
    mpi_error = MPI_Bcast(input->save_dir, MAX_STR_INPUT, MPI_CHAR, root, comm);
    // ph_save_dir and kernel are also broadcasted when entire save_dir is
    // Bcasted
    mpi_error = MPI_Bcast(&input->nkpool, 1, MPI_INT, root, comm);
    mpi_error = MPI_Bcast(&input->nqpool, 1, MPI_INT, root, comm);
    mpi_error = MPI_Bcast(&input->start_bnd, 1, MPI_INT, root, comm);
    mpi_error = MPI_Bcast(&input->end_bnd, 1, MPI_INT, root, comm);
}

static int handler(void* user, const char* section, const char* name,
                   const char* value)
{
    struct usr_input* inp = user;

    // check if value is just an empty string
    size_t nospace_len = 0;
    size_t val_len = strlen(value);
    for (size_t i = 0; i < val_len; ++i)
    {
        if (value[i] != ' ')
        {
            ++nospace_len;
        }
    }

    if (nospace_len == 0)
    {
        printf("Invalid input for %s ", name);
        error_msg("Invalid input");
    }

    if (strcmp(name, "nkpool") == 0)
    {
        inp->nkpool = atoi(value);
    }
    else if (strcmp(name, "nqpool") == 0)
    {
        inp->nqpool = atoi(value);
    }
    else if (strcmp(name, "start_bnd") == 0)
    {
        inp->start_bnd = atoi(value);
    }
    else if (strcmp(name, "end_bnd") == 0)
    {
        inp->end_bnd = atoi(value);
    }
    else if (strcmp(name, "save_dir") == 0)
    {
        strcpy(inp->save_dir, value);
    }
    else if (strcmp(name, "ph_save_dir") == 0)
    {
        strcpy(inp->ph_save_dir, value);
    }
    else if (strcmp(name, "kernel") == 0)
    {
        strcpy(inp->kernel, value);
        for (char* p = inp->kernel; *p; ++p)
        {
            *p = tolower(*p);
        }
    }
    else
    {
        error_msg("Invalid variable in input file.");
    }

    return 1;
}

void read_input_file(const char* input_file, struct usr_input** input_data,
                     MPI_Comm MPI_world_comm)
{
    // input_data must be free outside

    init_usr_input(input_data);

    int mpi_world_rank, mpi_error;
    mpi_error = MPI_Comm_rank(MPI_world_comm, &mpi_world_rank);

    if (mpi_world_rank == 0)
    {
        if (ini_parse(input_file, handler, *input_data) < 0)
        {
            error_msg("Cannot open input file.");
        }
    }
    // broad cast;
    Bcast_input_data(*input_data, 0, MPI_world_comm);
}
