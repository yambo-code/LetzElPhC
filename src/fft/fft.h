#pragma once
#include "../elphC.h"


/* struct store FFT plans for the forward and back.*/
struct ELPH_fft_plan
{      
    ELPH_cmplx * fft_data ; //
    /* This is input/output buffer that is used to perform inplane xy transforms.
    This buffer is also used as output in z transforms.
    must be free when the plan is destroyed.
    size = max(nGxy*Nz_loc, nGxyloc*Nz )
    */
    ELPH_cmplx * nz_buf;  
    /* second scratch space for transforms along z.
    this has a size = max(nGxy*Nz_loc, nGxyloc*Nz )
    data must be filled in "nz_buf" buffer and transposed output will be in fft_data
    this is created in the planer and will be freed when plan is destroyed //
    */
    ND_int fft_dims[3]; // (Nx,Ny,Nz) //
    ND_int nzloc;       // number of local z elements //

    const int * gvecs  ; // local gVecs //
    ND_int ngvecs_loc;   // local number of Gvecs //

    ND_int nGxyloc;     // number of Gvecs local //
    ND_int nGxy;        // total number of Gvecs local //
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
    ND_function(FFT_plan, Nd_cmplxS) fplan_xy;
    ND_function(FFT_plan, Nd_cmplxS) fplan_z; 

    /* backward plans G->r */
    ND_function(FFT_plan, Nd_cmplxS) bplan_xy;
    ND_function(FFT_plan, Nd_cmplxS) bplan_z; 

    MPI_Comm comm;   // comm
};



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
