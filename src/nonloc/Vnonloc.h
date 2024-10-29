#pragma once
#include "../common/dtypes.h"
#include "../elphC.h"

void sum_VNL_KKp(ELPH_cmplx* K_ptr, ELPH_cmplx* Kp_ptr, ELPH_cmplx* fcoeff,
                 ND_int nspin, ND_int nbnd, ND_int nspinor,
                 ELPH_cmplx* restrict out);

void add_elphNonLocal(struct WFC* wfcs, struct Lattice* lattice,
                      struct Pseudo* pseudo, int ikq, int ik, int kqsym,
                      int ksym, ELPH_cmplx* eigVec, ELPH_cmplx* elph_kq_mn,
                      const struct ELPH_MPI_Comms* Comm);
