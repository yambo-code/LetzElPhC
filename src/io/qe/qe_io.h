#pragma once
#include "../ezxml/ezxml.h"
#include "../../elphC.h"
#include "../../common/numerical_func.h"
#include "../../common/data_structs.h"
#include "../../common/parallel.h"
#include "../../symmetries/symmetries.h"
#include "../../common/string_func.h"


void parse_qexml(const char * xml_file, ELPH_float * lat_vec, ELPH_float * alat, \
        char * dim, bool * is_soc_present, ND_int * nmag, ND_int * fft_dims, \
        ND_int * nph_sym, ELPH_float ** ph_sym_mats, ELPH_float ** ph_sym_tau, \
        bool * ph_tim_rev, char ** pseudo_dir, char *** pseudo_pots);

void read_pattern_qe(const char * pat_file, struct Lattice * lattice, \
        ELPH_cmplx * restrict pat_vecs);


void read_dvscf_qe(const char * dvscf_file, struct Lattice * lattice, \
            const ELPH_cmplx * eig, const ELPH_cmplx * pats, \
            ELPH_cmplx * restrict dvscf_out, MPI_Comm commK);

void read_qpts_qe(const char * dyn0_file, ND_int * nqpt_iBZ, \
                ND_int * nqpt_fullBZ, ELPH_float ** qpts);

ND_int read_dyn_qe(const char * dyn_file, struct Lattice * lattice, \
    ELPH_float * restrict qpts, ELPH_float * restrict omega, \
    ELPH_cmplx * restrict pol_vecs);


void get_data_from_qe(struct Lattice * lattice, \
    struct Phonon * phonon, const char * ph_save_dir, \
    char *** pseudo_pots, const struct ELPH_MPI_Comms * Comm);

void get_dvscf_dyn_qe(const char * ph_save_dir, struct Lattice * lattice, \
        ND_int iq_BZ, ELPH_cmplx * eig, ELPH_cmplx * dvscf, ELPH_float * omega_ph, \
        const struct ELPH_MPI_Comms * Comm);






