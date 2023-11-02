#pragma once
#include "../elphC.h"
#include "../common/numerical_func.h"

/* struct store FFT plans for the forward and back.*/
struct ELPH_fft_plan
{      
    ELPH_cmplx * fft_data ; //
    /* This is a buffer used in AlltoAllv calls.
    Created by the planner so must be destroyed
    */
    ELPH_cmplx * nz_buf;  
    /* This is the main buffer where the z transforms are perfomed //
    Created by the planner so must be destroyed
    */
    ND_int fft_dims[3]; // (Nx,Ny,Nz) //
    ND_int nzloc;       // number of local z elements //
    ND_int align_len;   // simd alignment len

    const int * gvecs  ; // local gVecs //
    ND_int ngvecs_loc;   // local number of Gvecs //

    ND_int nGxyloc;     // number of Gvecs local //
    ND_int nGxy;        // total number of Gvecs local //
    ND_int Gxmin;       // min value of Gx
    ND_int Gxmax;       // max value of Gx
    int * Gxy_total;     // total gxy (nGxy,2) must be in [0,N)
    int * ngxy_z;        // number of z components for each x,y

    int * comm_bufs; // this (4,comm_size), used in alltoallv calls
    /*
    This buffer contains 4 x comm_size arrays 
    0) gxy_counts : (gxys present in each cpu) * nzloc
    1) xy_disp    : (starting index of gxy for each cpu) * nzloc
    2) z_counts   : (Nz present in each cpu) * nGxy_loc
    3) z_disp     : (starting index of Nz for each cpu) *nGxy_loc
    */

    /* forward plans r->G */
    ND_function(FFT_plan, Nd_cmplxS) * fplan_x;  // (naligment plans) for x 
    ND_function(FFT_plan, Nd_cmplxS) * fplan_y;  // (naligment plans) for y
    ND_function(FFT_plan, Nd_cmplxS) * fplan_y1; // (naligment plans) for y
    ND_function(FFT_plan, Nd_cmplxS)   fplan_z; 

    /* backward plans G->r */
    ND_function(FFT_plan, Nd_cmplxS) * bplan_x;  // (naligment plans) for x 
    ND_function(FFT_plan, Nd_cmplxS) * bplan_y;  // (naligment plans) for y
    ND_function(FFT_plan, Nd_cmplxS) * bplan_y1; // (naligment plans) for y
    ND_function(FFT_plan, Nd_cmplxS)   bplan_z; 

    MPI_Comm comm;   // comm
};

/* align.c */
// aligment length finder
ND_int alignment_len(void);

/*planner */
void wfc_plan(struct ELPH_fft_plan * plan, const ND_int ngvecs_loc, const ND_int nzloc, \
            const ND_int nGxyloc, const int * gvecs, const ND_int * fft_dims, \
            unsigned fft_flags, MPI_Comm comm);

void wfc_destroy_plan(struct ELPH_fft_plan * plan);

/* fft */
void fft3D(struct ELPH_fft_plan * plan, const ND_int nsets, \
            ELPH_cmplx * wfcr, ELPH_cmplx * wfcG, const bool conjugate);

/*invfft*/
void invfft3D(struct ELPH_fft_plan * plan, const ND_int nsets, \
            ELPH_cmplx * wfcG, ELPH_cmplx * wfcr, const bool conjugate);


void fwd_transpose(struct ELPH_fft_plan * plan);
void bwd_transpose(struct ELPH_fft_plan * plan);


