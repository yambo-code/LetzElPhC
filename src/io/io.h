#pragma once
#include "../elphC.h"
#include "../common/numerical_func.h"
#include "../common/dtypes.h"
#include "../common/parallel.h"
#include "../nonloc/fcoeff.h"
#include "../symmetries/symmetries.h"
#include "../common/string_func.h"
#include "../dvloc/dvloc.h"
#include "qe/qe_io.h"

/* This struct contains all the input file details */
struct usr_input
{   
    // system varibles
    int     nkpool              ; // k point parallelization
    int     nqpool              ; // q point parallelization
    int     start_bnd           ; // starting band
    int     end_bnd             ; // last band
    char *  save_dir            ; // save dir
    char *  ph_save_dir         ; // ph_save directory
    char *  kernel              ; // level of screening to include
};




#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); error_msg("netcdf_error");}

void read_and_alloc_save_data(char * SAVEdir, const struct ELPH_MPI_Comms * Comm, \
                ND_int start_band, ND_int end_band, struct WFC ** wfcs, \
                char * ph_save_dir, struct Lattice * lattice, \
                struct Pseudo * pseudo, struct Phonon * phonon, char * dft_code);

void free_save_data(struct WFC * wfcs, struct Lattice * lattice, \
                    struct Pseudo * pseudo, struct Phonon * phonon);

void Bcast_ND_arrayFloat(ND_array(Nd_floatS) * array_in, bool alloc_mem, int root, MPI_Comm comm);
void Bcast_ND_arrayCmplx(ND_array(Nd_cmplxS) * array_in, bool alloc_mem, int root, MPI_Comm comm);
void Bcast_wfc(struct WFC * wfc, bool alloc_mem, int root, MPI_Comm comm);


void parse_upf(const char * filename, ELPH_float * Zval, ND_array(Nd_floatS)* Vloc, \
                ND_array(Nd_floatS)* r_grid, ND_array(Nd_floatS)* rab_grid);
void get_upf_element(const char * filename, char* atomic_sym);


void init_usr_input(struct usr_input ** input);
void free_usr_input(struct usr_input *input);
void read_input_file(const char * input_file, struct usr_input ** input_data, \
                    MPI_Comm MPI_world_comm);





