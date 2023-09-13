#pragma once


void elphNL_q(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo, \
            ELPH_float * qpt, ND_array(Nd_cmplxS) * eigVec, ELPH_cmplx * elph_kq);


void sum_VNL_KKp(ELPH_cmplx * K_ptr, ELPH_cmplx *  Kp_ptr, ELPH_cmplx * restrict fcoeff, \
                        ND_int nspin, ND_int nbnd, ND_int nspinor, ELPH_cmplx * out);

void elphNonLocal(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo, int ikq, int ik, \
                int kqsym, int ksym, const bool trivial_phase, ND_array(Nd_cmplxS) * eigVec, \
                ELPH_cmplx * elph_kq_mn, MPI_Comm commK);