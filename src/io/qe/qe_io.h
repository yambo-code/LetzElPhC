#pragma once
#include "../ezxml/ezxml.h"
#include "../io.h"


void parse_qexml(const char * xml_file, ELPH_float * alat, char * dim, \
                bool * is_soc_present, ND_int * nmag, ND_int * fft_dims, ND_int * nph_sym, \
                ELPH_float ** ph_sym_mats, ELPH_float ** ph_sym_tau, bool * ph_tim_rev, \
                char ** pseudo_dir, char *** pseudo_pots);

void read_pattern_qe(const char * pat_file, struct Lattice * lattice, \
        ELPH_cmplx * restrict pat_vecs);

void read_dvscf_qe(const char * dvscf_file, struct Lattice * lattice, \
            const ND_int nmag, const ELPH_cmplx * eig, const ELPH_cmplx * pats, \
            ELPH_cmplx * restrict dvscf_out, MPI_Comm commK);

void read_qpts_qe(const char * dyn0_file, ND_int * nqpt_iBZ, \
                ND_int * nqpt_fullBZ, ELPH_float ** qpts);

ND_int read_dyn_qe(const char * dyn_file, struct Lattice * lattice, \
    ELPH_float * restrict qpts, ELPH_float * restrict omega, \
    ELPH_cmplx * restrict pol_vecs);







