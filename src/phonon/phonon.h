#pragma once
#include "common/dtypes.h"
#include "elphC.h"

void mass_normalize_pol_vecs(const ELPH_float* atomic_masses,
                             const ND_int nsets, const ND_int natoms,
                             const ELPH_float power, ELPH_cmplx* pol_vecs);

void mass_normalize_force_constants(const ELPH_float* atomic_masses,
                                    const ND_int nsets, const ND_int natoms,
                                    const ELPH_float power, ELPH_cmplx* frc);

void pol_vecs_to_dyn(const ELPH_float* omega, const ND_int natom,
                     const ELPH_float* atomic_masses, ELPH_cmplx* pol_vecs);

void add_ph_dyn_long_range(const ELPH_float* qpt, struct Lattice* lattice,
                           struct Phonon* phonon, const ND_int* Ggrid,
                           const ND_int sign, const ELPH_float* atomic_masses,
                           ELPH_cmplx* dyn_mat_asr, const ELPH_float eta,
                           ELPH_cmplx* dyn_mat);

void compute_dyn_lr_asr_correction(struct Lattice* lattice,
                                   struct Phonon* phonon, const ND_int* Ggrid,
                                   const ELPH_float eta,
                                   ELPH_cmplx* dyn_mat_asr);

enum asr_kind asr_kind_from_string(const char* str);

void apply_acoustic_sum_rule_fc(enum asr_kind mode, const ND_int* qgrid,
                                const ND_int nat, ELPH_cmplx* frc,
                                const ELPH_float* atomic_pos,
                                const ELPH_float* lat_vecs,
                                const ND_int* ws_vecs, const ND_int n_ws_vecs,
                                const ND_int* ws_degen);

void apply_acoustic_sum_rule_born_charges(enum asr_kind mode, ELPH_float* Zborn,
                                          const ND_int nat);
