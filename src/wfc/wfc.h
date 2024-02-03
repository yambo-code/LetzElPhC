/*
This is the header for functions related to wfcs
*/

#pragma once
#include "../elphC.h"
#include "../common/numerical_func.h"
#include "../common/data_structs.h"
#include "../common/parallel.h"
#include "../fft/fft.h"
#include "gsort.h"

void rotateGvecs(const ELPH_float * Gvecs, const ELPH_float * sym, const ND_int ngvecs, 
                const  ELPH_float * lat_vec, const bool inverse, const bool crystal,
                ELPH_float * G0, ELPH_float *  Gvec_out);

void su2rotate(const int nspinor, const ND_int npw, const ND_int nsets, \
            const ELPH_cmplx * restrict su2mat, ELPH_cmplx * restrict wfc);
            

void apply_trans_wfc(const ELPH_float * trans_vec, const ELPH_float * kvec, \
                     const ND_int nsets, const ND_int npw, const ELPH_float * gvecs, \
                     ELPH_cmplx * restrict wfc_G, const bool conjugate);

/* su2mat.c */
void SU2mat(const ELPH_float * sym_in, const ND_int nspinor, \
            const bool invert_sym, const bool time_rev, ELPH_cmplx * su2mat);


void Sort_pw(const ND_int npw_tot, const ND_int npw_loc, const ND_int * fft_dims,\
                const ELPH_float * gvec_in, const ELPH_cmplx * wfc_in, \
                const ND_int nsets, ND_int * npw_loc_out, ND_int * nG_xy,
                int ** gvec_out, ELPH_cmplx ** wfc_out, MPI_Comm mpi_comm);






