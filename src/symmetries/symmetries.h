#pragma once
#include "../common/numerical_func.h"
#include "../common/parallel.h"
#include "../elphC.h"

#define Function(FUN_NAME, TYPE_SMALL) Function_HIDDEN(FUN_NAME, TYPE_SMALL)
#define Function_HIDDEN(FUN_NAME, TYPE_SMALL) ELPH_##TYPE_SMALL##FUN_NAME

void bz_expand(ND_array(Nd_floatS) * ibz_kpts, ND_array(Nd_floatS) * sym_mats,
               ND_array(Nd_floatS) * lat_vec, ND_array(Nd_floatS) * kpoints,
               nd_arr_i* kmap);

void get_KplusQ_idxs(ND_array(Nd_floatS) * kpoints, int* KplusQidxs,
                     ELPH_float* Q_pt, ND_array(Nd_floatS) * lat_vec,
                     bool Qincrystal);

bool Function(isVECpresent, Nd_cmplxS)(const ND_array(Nd_floatS) * array,
                                       const ELPH_float* vec, ND_int* idx);

void electronic_reps(const struct WFC* wfcs, const struct Lattice* lattice,
                     const ELPH_float* Rsym_mat, const ELPH_float* tauR,
                     const bool tim_revR, const ND_int ikBZ,
                     ELPH_cmplx* Dkmn_rep, const struct ELPH_MPI_Comms* Comm);

void elph_q_rotate(const ELPH_cmplx* elph_mat_q, const struct Lattice* lattice,
                   const ELPH_cmplx* Dmats, const ELPH_float* symS,
                   const bool tim_revS, const ELPH_cmplx fac,
                   const ELPH_float* qpt, ELPH_cmplx* restrict elph_mat_Sq);
