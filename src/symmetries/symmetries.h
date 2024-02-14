#pragma once
#include "../common/numerical_func.h"
#include "../common/parallel.h"
#include "../elphC.h"

ND_int bz_expand(const ND_int Nibz, const ND_int Nsym,
                 const ELPH_float* ibz_kpts, const struct symmetry* symms,
                 const ELPH_float* lat_vec, ELPH_float* kpoints,
                 ND_int* kstar, int* kmap);

void electronic_reps(const struct WFC* wfcs, const struct Lattice* lattice,
                     const ELPH_float* Rsym_mat, const ELPH_float* tauR,
                     const bool tim_revR, const ND_int ikBZ,
                     ELPH_cmplx* Dkmn_rep, const struct ELPH_MPI_Comms* Comm);

void elph_q_rotate(const ELPH_cmplx* elph_mat_q, const struct Lattice* lattice,
                   const ELPH_cmplx* Dmats, const ELPH_float* symS,
                   const bool tim_revS, const ELPH_cmplx fac,
                   const ELPH_float* qpt, ELPH_cmplx* restrict elph_mat_Sq);

ND_int find_inv_symm_idx(ND_int nsym, const ELPH_float* Smat,
                         const ELPH_float* point_group, const bool trans);
