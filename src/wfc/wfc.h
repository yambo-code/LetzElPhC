/*
This is the header for functions related to wfcs
*/

#pragma once
#include "helpers_wfc.h"


/* su2mat.c */
void SU2mat(const ELPH_float * sym_in, const ND_int nspinor, \
            const bool invert_sym, const bool time_rev, ELPH_cmplx * su2mat);

/* box2sphere.c */
void box2sphere(const ELPH_cmplx * restrict wfcGbox, const ND_int nsets, \
                const ELPH_float * restrict Gvec_crys, const ND_int npw, \
                const ND_int * FFT_dims, ELPH_cmplx * restrict wfcGsphere);

/* sphere2box.c */
void sphere2box(const ELPH_cmplx * restrict wfcGsphere , const ND_int nsets, \
                const ELPH_float * restrict Gvec_crys, const ND_int npw, \
                const ND_int * FFT_dims, ELPH_cmplx * restrict wfcGbox);

void get_wfc_from_pool(struct WFC * wfcs, ND_int ik, ND_int niBZ_tot, \
            MPI_Comm commK, MPI_Comm commQ, struct WFC * wfc_out);

