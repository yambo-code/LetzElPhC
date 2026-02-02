#pragma once
#include <mpi.h>
#include <stdbool.h>

#include "common/dtypes.h"
#include "elphC.h"

void parse_qexml(const char* xml_file, ND_int* natoms, ELPH_float* lat_vec,
                 ELPH_float* alat, char* dim, bool* is_soc_present,
                 ND_int* nmag, ND_int* fft_dims, ND_int* nph_sym,
                 ELPH_float** ph_sym_mats, ELPH_float** ph_sym_tau,
                 bool** ph_trevs, bool* ph_mag_symm, bool* ph_tim_rev,
                 char** pseudo_dir, char*** pseudo_pots, int* nspinor,
                 ND_int* ntype, int** atomic_type,
                 ELPH_float** atomic_positons);

void get_interpolation_data_from_qe(struct Lattice* lattice,
                                    struct Phonon* phonon,
                                    struct Pseudo* pseudo,
                                    const char* ph_save_dir, ELPH_float** Zvals,
                                    ELPH_float* alat,
                                    const struct ELPH_MPI_Comms* Comm);

void read_ph_tensors_qe(const char* tensor_xml_file, const ND_int natom,
                        struct Phonon* phonon);

void read_quadrupole_fmt(const char* filename, ELPH_float** Qpole_buf,
                         int natom);

void read_pattern_qe(const char* pat_file, struct Lattice* lattice,
                     ELPH_cmplx* pat_vecs);

void read_dvscf_qe(const char* dvscf_file, struct Lattice* lattice,
                   const ELPH_cmplx* eig, const ELPH_cmplx* pats,
                   ELPH_cmplx* dvscf_out, MPI_Comm commK);

void read_qpts_qe(const char* dyn0_file, ND_int* nqpt_iBZ, ND_int* nqpt_fullBZ,
                  ELPH_float** qpts);

ND_int read_dyn_qe(const char* dyn_file, struct Lattice* lattice,
                   ELPH_float* qpts, ELPH_float* omega, ELPH_cmplx* pol_vecs,
                   ELPH_float* amass);

void get_data_from_qe(struct Lattice* lattice, struct Phonon* phonon,
                      struct Pseudo* pseudo, const char* ph_save_dir,
                      ELPH_float* alat_param,
                      const struct ELPH_MPI_Comms* Comm);

void get_dvscf_dyn_qe(const char* ph_save_dir, struct Lattice* lattice,
                      ND_int iq_BZ, ELPH_cmplx* eig, ELPH_cmplx* dvscf,
                      ELPH_float* omega_ph, const struct ELPH_MPI_Comms* Comm);

// write functions
void write_dvscf_qe(const char* dvscf_file, struct Lattice* lattice,
                    const ELPH_cmplx* dvscf_in, MPI_Comm commK);

void write_dyn_qe(const char* file_name, ND_int natom, const ELPH_float* qpts,
                  const ELPH_cmplx* dyn_mat, const ELPH_float* atomic_masses);

void write_qpts_qe(const char* dyn_file, ND_int nqpt_iBZ,
                   const ELPH_float* qpts, const ND_int* qgrid);
