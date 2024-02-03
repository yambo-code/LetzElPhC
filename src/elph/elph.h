#pragma once
#include "../dvloc/dvloc.h"
#include "../nonloc/Vnonloc.h"
#include "../io/io.h"


void compute_elphq(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo, \
            ELPH_float * qpt, ND_array(Nd_cmplxS) * eigVec, ND_array(Nd_cmplxS) * dVscfq, 
            ELPH_cmplx * elph_kq, const struct ELPH_MPI_Comms * Comm);


void compute_and_write_dmats(const char * file_name, const struct WFC * wfcs, \
        const struct Lattice * lattice, const ND_int nph_sym, \
        const ELPH_float * Rsym_mat,  const ELPH_float * tauR, \
        const bool *tim_revR, const struct ELPH_MPI_Comms * Comm);


