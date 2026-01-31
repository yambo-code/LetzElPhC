#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/dtypes.h"
#include "common/error.h"
#include "common/string_func.h"
#include "parser.h"
#include "parser/inih/ini.h"

#define READ_STR_LEN 600

static void Bcast_elph_input_data(struct elph_usr_input* input, int root,
                                  MPI_Comm comm);
static int elph_input_handler(void* user, const char* section, const char* name,
                              const char* value);

// function to alloc, initiate elph_usr_input
void init_elph_usr_input(struct elph_usr_input** input)
{
    // this function also sets defaults for the user input file
    *input = malloc(sizeof(struct elph_usr_input));
    CHECK_ALLOC(*input);

    struct elph_usr_input* inp = *input;

    inp->save_dir = calloc(READ_STR_LEN * 3, 1);
    CHECK_ALLOC(inp->save_dir);

    inp->ph_save_dir = inp->save_dir + READ_STR_LEN;
    inp->kernel_str = inp->save_dir + 2 * READ_STR_LEN;

    // defaults
    inp->nkpool = 1;
    inp->nqpool = 1;
    inp->start_bnd = 0;
    inp->end_bnd = 0;
    strncpy_custom(inp->save_dir, "SAVE", READ_STR_LEN);
    strncpy_custom(inp->ph_save_dir, "ph_save", READ_STR_LEN);
    strncpy_custom(inp->kernel_str, "dfpt", READ_STR_LEN);
    inp->kminusq = false;  // default is standard
}

// function to free elph_usr_input struct data
void free_elph_usr_input(struct elph_usr_input* input)
{
    free(input->save_dir);
    free(input);
}

static void Bcast_elph_input_data(struct elph_usr_input* input, int root,
                                  MPI_Comm comm)
{
    int mpi_error;

    // all char * will be bcasted in one single call
    mpi_error =
        MPI_Bcast(input->save_dir, READ_STR_LEN * 3, MPI_CHAR, root, comm);
    MPI_error_msg(mpi_error);
    // ph_save_dir and kernel_str are also broadcasted when entire save_dir is
    // Bcasted
    mpi_error = MPI_Bcast(&input->nkpool, 1, MPI_INT, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&input->nqpool, 1, MPI_INT, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&input->start_bnd, 1, MPI_INT, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&input->end_bnd, 1, MPI_INT, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&input->kminusq, 1, MPI_C_BOOL, root, comm);
    MPI_error_msg(mpi_error);
}

static int elph_input_handler(void* user, const char* section, const char* name,
                              const char* value)
{
    UNUSED_VAR(section);
    // All the new inputs are added here.
    // Note : Set the defaults in init_elph_usr_input function.
    struct elph_usr_input* inp = user;

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

    // add inputs from here. use else if
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
        strncpy_custom(inp->save_dir, value, READ_STR_LEN);
        strip_quotes_in_string(inp->save_dir);
    }
    else if (strcmp(name, "ph_save_dir") == 0)
    {
        strncpy_custom(inp->ph_save_dir, value, READ_STR_LEN);
        strip_quotes_in_string(inp->ph_save_dir);
    }
    else if (strcmp(name, "kernel") == 0)
    {
        strncpy_custom(inp->kernel_str, value, READ_STR_LEN);
        strip_quotes_in_string(inp->kernel_str);
        lowercase_str(inp->kernel_str);
    }
    else if (strcmp(name, "convention") == 0)
    {
        if (strstr(value, "yambo"))
        {
            inp->kminusq = true;
        }
        else if (strstr(value, "standard"))
        {
            inp->kminusq = false;
        }
        else
        {
            error_msg(
                "Invalid value for convention in the input file."
                " Only yambo or standard are accepted.");
        }
    }
    else
    {
        error_msg("Invalid variable in the elph input file.");
    }

    return 1;
}

void read_elph_input_file(const char* input_file,
                          struct elph_usr_input** input_data,
                          MPI_Comm MPI_world_comm)
{
    // input_data must be free outside

    init_elph_usr_input(input_data);

    int mpi_world_rank, mpi_error;

    mpi_error = MPI_Comm_rank(MPI_world_comm, &mpi_world_rank);
    MPI_error_msg(mpi_error);

    if (mpi_world_rank == 0)
    {
        if (ini_parse(input_file, elph_input_handler, *input_data) < 0)
        {
            error_msg("Cannot open input file.");
        }
    }
    // broad cast;
    Bcast_elph_input_data(*input_data, 0, MPI_world_comm);
}
