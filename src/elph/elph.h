#pragma once
#include "../dvloc/dvloc.h"
#include "../nonloc/Vnonloc.h"
#include "../io/io.h"


void compute_and_write_elphq(struct WFC * wfcs, struct Lattice * lattice, \
            struct Pseudo * pseudo, struct Phonon  * phonon, \
            const ND_int iqpt, ND_array(Nd_cmplxS) * eigVec, ND_array(Nd_cmplxS) * dVscfq, 
            const int ncid, const int varid, const bool non_loc, const bool kminusq, \
            const struct ELPH_MPI_Comms * Comm);


void compute_and_write_dmats(const char * file_name, const struct WFC * wfcs, \
        const struct Lattice * lattice, const ND_int nph_sym, \
        const ELPH_float * Rsym_mat,  const ELPH_float * tauR, \
        const bool *tim_revR, const struct ELPH_MPI_Comms * Comm);


