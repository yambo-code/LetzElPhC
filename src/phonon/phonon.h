#pragma once
#include "common/dtypes.h"
#include "elphC.h"

void mass_normalize_pol_vecs(const ELPH_float* atomic_masses,
                             const ND_int nsets, const ND_int natoms,
                             const ELPH_float power, ELPH_cmplx* pol_vecs);

void pol_vecs_to_dyn(const ELPH_float* omega, const ND_int natom,
                     const ELPH_float* atomic_masses, ELPH_cmplx* pol_vecs);

void add_ph_dyn_long_range(const ELPH_float* qpt, struct Lattice* lattice,
                           struct Phonon* phonon, const ND_int* Ggrid,
                           const ND_int sign, const ELPH_float* atomic_masses,
                           ELPH_cmplx* dyn_mat_asr, const ELPH_float eta,
                           ELPH_cmplx* dyn_mat);

void compute_dyn_lr_asr_correction(struct Lattice* lattice,
                                   struct Phonon* phonon, const ND_int* Ggrid,
                                   const ELPH_float* atomic_masses,
                                   const ELPH_float eta,
                                   ELPH_cmplx* dyn_mat_asr);
