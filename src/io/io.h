#pragma once

#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <stddef.h>

#include "common/dtypes.h"
#include "elphC.h"

#define NC4_DEFAULT_CHUCK_KB 2048
// default chunking for large nc varaibles (in Kilobytes)

void read_and_alloc_save_data(char* SAVEdir, const struct ELPH_MPI_Comms* Comm,
                              ND_int start_band, ND_int end_band,
                              struct WFC** wfcs, char* ph_save_dir,
                              struct Lattice* lattice, struct Pseudo* pseudo,
                              struct Phonon* phonon,
                              enum ELPH_dft_code dft_code);

void free_phonon_data(struct Phonon* phonon);

void free_save_data(struct WFC* wfcs, struct Lattice* lattice,
                    struct Pseudo* pseudo, struct Phonon* phonon);

//======= nc4 function
void def_ncVar(const int ncid, int* varid, ND_int rank, nc_type xtype,
               ND_int* dims, const char* var_name, char** dim_names,
               size_t* chunksize);

void write_basic_data(const int ncid, struct Lattice* lattice,
                      struct Phonon* phonon, const char* kernel_str,
                      const char* convention_str);
