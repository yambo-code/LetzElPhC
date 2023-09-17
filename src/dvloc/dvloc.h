#pragma once
#include "../elphC.h"
#include "../common/data_structs.h"
#include "../wfc/wfc.h"
#include "../common/numerical_func.h"


/*dvpsi.c*/
void dvpsi(const ND_int nmag, const ND_int nspinor, const ND_int nFFT, 
            const ELPH_cmplx * restrict dV_r, ELPH_cmplx * restrict psi_r);

/* VlocG_atom.c */
ELPH_float Vloc_Gspace(ELPH_float * work_arr, const char cutoff, const ELPH_float Gnorm, \
        const ND_array(Nd_floatS)* Vloc_atomic, const ND_array(Nd_floatS)* r_grid, \
        const ND_array(Nd_floatS)* rab_grid, const ELPH_float Zval,const ELPH_float eta, \
        const ELPH_float cutoff_fac);

void dvpsi(const ND_int nmag, const ND_int nspinor, const ND_int nFFT, 
            const ELPH_cmplx * restrict dV_r, const ELPH_cmplx * restrict psi_r, \
            ELPH_cmplx * restrict dVpsi_out);

void add_dvscf(ND_array(Nd_cmplxS) * dVscf, ND_array(Nd_cmplxS) * dVloc); 


void VlocinVFFT(ND_array(Nd_cmplxS) * VlocPot, struct Lattice * lattice, MPI_Comm mpi_comm);


void dVlocq(const ELPH_float * qpt, struct Lattice * lattice, struct Pseudo * pseudo, \
            ND_array(Nd_cmplxS) * eigVec, ND_array(Nd_cmplxS) * Vlocr, MPI_Comm commK);

