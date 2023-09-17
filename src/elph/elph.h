#pragma once
#include "../dvloc/dvloc.h"
#include "../nonloc/Vnonloc.h"



void compute_elph(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo, \
            ELPH_float * qpt, ND_array(Nd_cmplxS) * eigVec, ND_array(Nd_cmplxS) * dVscfq, 
            ELPH_cmplx * elph_kq, MPI_Comm commQ, MPI_Comm commK);


