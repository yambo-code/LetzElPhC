#pragma once
#include "../ezxml/ezxml.h"
#include "../io.h"

void parse_qexml(const char * xml_file, ELPH_float * alat, \
        ND_int * fft_dims, char ** pseudo_dir, char *** pseudo_pots);


void read_pattern_qe(const char * pat_file, struct Lattice * lattice, \
        ELPH_cmplx * restrict pat_vecs);

void read_dvscf_qe(const char * dvscf_file, struct Lattice * lattice, \
            const ND_int nmag, const ELPH_cmplx * eig, const ELPH_cmplx * pats, \
            ELPH_cmplx * restrict dvscf_out, MPI_Comm commK);


ND_int read_dyn_qe(const char * dyn_file, struct Lattice * lattice, \
    ELPH_float * restrict qpts, ELPH_float * restrict omega, \
    ELPH_cmplx * restrict pol_vecs);







