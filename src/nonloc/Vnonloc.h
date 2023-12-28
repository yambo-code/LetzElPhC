#pragma once
#include "../dvloc/dvloc.h"

void sum_VNL_KKp(ELPH_cmplx * K_ptr, ELPH_cmplx *  Kp_ptr, ELPH_cmplx * fcoeff, \
                        ND_int nspin, ND_int nbnd, ND_int nspinor, ELPH_cmplx * restrict out);

void add_elphNonLocal(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo, \
                    int ikq, int ik, int kqsym, int ksym,  ND_array(Nd_cmplxS) * eigVec, \
                    ELPH_cmplx * elph_kq_mn, MPI_Comm commK);

