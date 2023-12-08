#pragma once
#include "../elphC.h"
#include "../common/numerical_func.h"
#include "../common/data_structs.h"
#include "../common/parallel.h"
#include "../nonloc/fcoeff.h"
#include "../symmetries/symmetries.h"


/* This struct contains all the input file details */
struct usr_input
{   
    // system varibles
    int     nkpool       ; // k point parallelization
    int     nqpool       ; // q point parallelization
    int     start_bnd    ; // starting band
    int     end_bnd      ; // last band
    char *  save_dir     ; // save dir
    char *  pseudo_dir   ; // pseudo pot dir
    char *  dvscf_file   ; // dvscf directory
    char ** pseudos      ;  // pseudo separated by comma "a.upf,b.upf"
    char    dimension    ; // system dimension 3/2/1 Currently only 2/3 are supported
};



#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);}

void read_and_alloc_save_data(char * SAVEdir, MPI_Comm commQ, MPI_Comm commK,  \
                ND_int start_band, ND_int end_band, struct WFC ** wfcs,char * pseudo_dir, \
                char ** pseudo_pots, struct Lattice * lattice, struct Pseudo * pseudo, \
                const ND_int * FFT_dims);
void free_save_data(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo);


void Bcast_ND_arrayFloat(ND_array(Nd_floatS) * array_in, bool alloc_mem, int root, MPI_Comm comm);
void Bcast_ND_arrayCmplx(ND_array(Nd_cmplxS) * array_in, bool alloc_mem, int root, MPI_Comm comm);
void Bcast_wfc(struct WFC * wfc, bool alloc_mem, int root, MPI_Comm comm);


void parse_upf2(const char * filename, ELPH_float * Zval, ND_array(Nd_floatS)* Vloc, \
                ND_array(Nd_floatS)* r_grid, ND_array(Nd_floatS)* rab_grid);
void get_upf_element(const char * filename, char* atomic_sym);

void get_FFT_dims(const char * file_name, ND_int * nq, ND_int * fft_dims);

void read_dvscfq(const char * file_name, ND_array(Nd_cmplxS) * eigVec,  \
                struct Lattice * lattice, ND_array(Nd_cmplxS) *dVscf, \
                ND_int iq, MPI_Comm commK);


void init_usr_input(struct usr_input ** input);
void free_usr_input(struct usr_input *input);
void read_input_file(const char * input_file, struct usr_input ** input_data);





