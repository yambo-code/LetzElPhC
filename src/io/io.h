#pragma once
#include "../common/dtypes.h"
#include "../common/numerical_func.h"
#include "../common/parallel.h"
#include "../common/string_func.h"
#include "../dvloc/dvloc.h"
#include "../elphC.h"
#include "../nonloc/fcoeff.h"
#include "../symmetries/symmetries.h"
#include "mpi_bcast.h"
#include "qe/qe_io.h"

/* This struct contains all the input file details */
struct usr_input
{
    // system varibles
    int nkpool; // k point parallelization
    int nqpool; // q point parallelization
    int start_bnd; // starting band
    int end_bnd; // last band
    char* save_dir; // save dir
    char* ph_save_dir; // ph_save directory
    char* kernel_str; // level of screening to include
    bool kminusq; // true if convention is "yambo" else false
};

#define ERR(e)                                          \
    {                                                   \
        fprintf(stderr, "Error: %s\n", nc_strerror(e)); \
        error_msg("netcdf_error");                      \
    }

void read_and_alloc_save_data(char* SAVEdir, const struct ELPH_MPI_Comms* Comm,
                              ND_int start_band, ND_int end_band,
                              struct WFC** wfcs, char* ph_save_dir,
                              struct Lattice* lattice, struct Pseudo* pseudo,
                              struct Phonon* phonon, enum ELPH_dft_code dft_code);

void free_save_data(struct WFC* wfcs, struct Lattice* lattice,
                    struct Pseudo* pseudo, struct Phonon* phonon);

// parse_upf.c // get pseudo data from upfs
void parse_upf(const char* filename, struct local_pseudo* loc_pseudo);
void get_upf_element(const char* filename, char* atomic_sym);

void init_usr_input(struct usr_input** input);
void free_usr_input(struct usr_input* input);
void read_input_file(const char* input_file, struct usr_input** input_data,
                     MPI_Comm MPI_world_comm);

//======= nc4 function
void def_ncVar(const int ncid, int* varid, ND_int rank, nc_type xtype,
               ND_int* dims, const char* var_name, char** dim_names,
               size_t* chunksize);

void write_basic_data(const int ncid, struct Lattice* lattice,
        struct Phonon* phonon, const char * kernel_str, const char * convention_str);


