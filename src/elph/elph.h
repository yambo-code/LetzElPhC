#pragma once
#include "../dvloc/dvloc.h"
#include "../io/io.h"
#include "../nonloc/Vnonloc.h"

void compute_and_write_elphq(struct WFC* wfcs, struct Lattice* lattice,
                             struct Pseudo* pseudo, struct Phonon* phonon,
                             const ND_int iqpt, ELPH_cmplx* eigVec,
                             ELPH_cmplx* dVscfq, 
                             const int ncid_elph, const int varid_elph, 
                             const int ncid_dmat, const int varid_dmat, 
                             const bool non_loc,
                             const bool kminusq,
                             const struct ELPH_MPI_Comms* Comm);

void compute_and_write_dmats(const char* file_name, const struct WFC* wfcs,
                             const struct Lattice* lattice,
                             const ND_int nph_sym, const struct symmetry* sym_data,
                             const struct ELPH_MPI_Comms* Comm);
