#pragma once
#include "../elphC.h"
#include "../common/dtypes.h"

void print_info_msg(int mpi_rank, const char* fmt, ...);

void print_input_info(const char* save_dir, const char* ph_save_dir,
                      const char* kernel, const bool kminusq, const enum ELPH_dft_code dft_code,
                      const struct ELPH_MPI_Comms* Comm);

void print_lattice_info(const struct ELPH_MPI_Comms* Comm, const struct Lattice* lattice);

void print_phonon_info(const struct ELPH_MPI_Comms* Comm, const struct Phonon* phonon);
