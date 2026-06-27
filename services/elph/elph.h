#pragma once
#include <mpi.h>
#include <stdbool.h>

#include "parser/parser.h"
#include "common/dtypes.h"
#include "elphC.h"
#include "yambo.h"

/*
 * Per-(iq_BZ, ik_BZ) fill callback for yambo COLL integration.
 * Called once per BZ (q,k) pair on commK_rank==0 processes.
 *   iq_BZ, ik_BZ  : 0-based BZ indices
 *   data           : ELPH_cmplx buffer, C row-major (nmodes, nspin, nbnds, nbnds)
 *   nq..nb_start   : full-BZ dimensions (constant across calls); nb_start is 1-based
 *   iq_iBZ         : 0-based iBZ q-point index for this iq_BZ star
 *   qpt_BZ_crys    : crystal-coordinate q-point for iq_BZ (ELPH_float[3])
 */
typedef void (*elph_fill_fn)(int iq_BZ, int ik_BZ,
                              const void* data, int *KplusQidxs,
                              int nq, int nk, int nmodes, int nspin,
                              int nbnds, int nb_start,
                              int iq_iBZ, const void* qpt_BZ_crys);

/*
 * Light callback (reduced signature): just indices and data pointer.
 * Constants (nmodes, nspin, nbnds) known to Fortran side.
 */
typedef void (*elph_fill_fn_light)(int iq_ibz_letz, int iq_bz_letz,
                                    int ik_ibz_letz, int ik_bz_letz,
                                    const void* data);

/*
 * Callback for dV_q^nu(G) in reciprocal space, called once per iBZ q-point
 * from commK rank 0.
 *   iq_iBZ   : 0-based global iBZ q-point index
 *   dVG      : C row-major (nmodes, nmag, nfft_x, nfft_y, nfft_z) after forward FFT
 *   ph_freqs : phonon frequencies omega_ph[nmodes] in Hartree (ELPH_float*)
 *   nq_iBZ   : total number of iBZ q-points
 */
typedef void (*elph_dvG_fill_fn)(int iq_iBZ,
                                  const void* dVG,
                                  const void* ph_freqs,
                                  int nq_iBZ, int nmodes, int nmag,
                                  int nfft_x, int nfft_y, int nfft_z);

void elph_driver(const char* ELPH_input_file, enum ELPH_dft_code dft_code,
                 MPI_Comm comm_world);

/* Callback-enabled variant: skips ndb.elph; calls fill_fn per (q,k) instead. */
void elph_driver_cb(const char* ELPH_input_file, enum ELPH_dft_code dft_code,
                    MPI_Comm comm_world, elph_fill_fn fill_fn);

/*
 * Extended callback variant: like elph_driver_cb plus calls dvG_fill_fn once per
 * iBZ q-point with the full dV_q^nu(G) potential in reciprocal space.
 * Either callback may be NULL to skip that output.
 * comm_q, comm_k: Y6 PAR communicators for q,k distribution.
 */
void elph_driver_cb2(struct elph_usr_input* input_data, struct Y6_info* y6_data, struct Y6_parallel_work* y6_work, enum ELPH_dft_code dft_code,
                     elph_fill_fn_light fill_fn,
                     elph_dvG_fill_fn dvG_fill_fn,int i_control,
                     MPI_Comm comm_world, int bz_mode_code);

void compute_and_write_elphq(struct WFC* wfcs, struct Lattice* lattice,
                             struct Pseudo* pseudo, struct Phonon* phonon,
                             const ND_int iqpt, ELPH_cmplx* eigVec,
                             ELPH_cmplx* dVscfq, const int ncid_elph,
                             const int varid_elph, const int ncid_dmat,
                             const int varid_dmat, const bool non_loc,
                             const bool kminusq,
                             const struct ELPH_MPI_Comms* Comm,
                             elph_fill_fn_light fill_fn,
                             const ND_int iqpt_iBZ, int bz_mode_code);

void compute_and_write_dmats(const char* file_name, const struct WFC* wfcs,
                             const struct Lattice* lattice,
                             const ND_int nph_sym,
                             const struct symmetry* sym_data,
                             const struct ELPH_MPI_Comms* Comm);

void init_kernel(struct kernel_info* kernel);

void set_kernel(const char* kernel_str, struct kernel_info* kernel);
