#pragma once
#include "../elphC.h"
#include "../common/numerical_func.h"
#include "../common/data_structs.h"
#include "../common/parallel.h"
#include "../nonloc/fcoeff.h"
#include "../symmetries/symmetries.h"


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





