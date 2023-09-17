/*
This is the header for functions related to wfcs
*/

#pragma once
#include "../elphC.h"
#include "../common/numerical_func.h"
#include "../common/data_structs.h"
#include "../common/parallel.h"
#include "../io/io.h"

int get_fft_idx(ELPH_float idx_in, int FFT_dimension);

void rotateGvecs(const ELPH_float * Gvecs, const ELPH_float * sym, const ND_int ngvecs, 
                const  ELPH_float * lat_vec, const bool inverse, const bool crystal,
                const ELPH_float * G0, ELPH_float *  Gvec_out);

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


void su2rotate(const int nspinor, const ND_int npw, const ND_int nsets, \
            const ELPH_cmplx * restrict su2mat, ELPH_cmplx * restrict wfc);

void wfcFFT(struct wfcBox * wfcRspace, const ELPH_float * sym, \
        const ELPH_float * tau, const int * ulmvec, const bool tim_rev, \
        const ELPH_float * lat_vec, const ND_int npw_total, \
        const ND_array(Nd_floatS)* Gvecs_in, \
        ND_array(Nd_cmplxS) * wfcG, MPI_Comm mpi_comm);

void wfcinVFFT(ND_array(Nd_cmplxS) * wfcG,  const ELPH_float * sym, 
        const ELPH_float * tau, const int * ulmvec, \
        const bool tim_rev, const ELPH_float * lat_vec, 
        const ND_int npw_total, const ND_array(Nd_floatS)* Gvecs_in, \
        struct wfcBox * wfcRspace, MPI_Comm mpi_comm);
        
// void get_wfc_from_pool(struct WFC * wfcs, ND_int ik, ND_int niBZ_tot, \
//             MPI_Comm commK, MPI_Comm commQ, struct WFC * wfc_out);

alloc_wfcBox(struct wfcBox * buffer, const ND_int * dimensions, \
                const ND_int npw_max_total, unsigned flag, MPI_Comm mpi_comm);

void free_wfcBox(struct wfcBox * buffer);
