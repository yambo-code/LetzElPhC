/*
This file contains box2sphere function which is called in wfc routines for 
map plane waves in FFT box to sphere. this is used when performing FFT on wfcs 
potentials.
*/
#include "wfc.h"

void box2sphere(const ELPH_cmplx * restrict wfcGbox, const ND_int nsets, \
                const ELPH_float * restrict Gvec_crys, const ND_int npw, \
                const ND_int * FFT_dims, ELPH_cmplx * restrict wfcGsphere)
{   
    /*
    Input : 
    nsets -> number of wave functions/potentials
    wfcGbox -> wave functions in the FFT G-box (nsets,Nx*Ny*Nz)
    npw -> number of planes in G sphere
    Gvec_crys -> miller indices in the G sphere (npw,3)
    FFT_dims -> FFT G-box dims (3) {Nx,Ny,Nz}

    Output :
    wfcGsphere --> wfc/potentials in G sphere (nsets,npw)
    */

    ND_int fft_strides[3] = {FFT_dims[1]*FFT_dims[2],FFT_dims[2],1};
    ND_int nFFTs = FFT_dims[0]*fft_strides[0];

    for (ND_int iset = 0; iset<nsets; ++iset )
    {   
        ELPH_OMP_PAR_FOR_SIMD
        for (ND_int ipw = 0; ipw<npw; ++ipw )
        {   
            int Nx, Ny, Nz;
            /* get the plane wave*/
            const ELPH_float * Gvec_temp = Gvec_crys + 3*ipw ;
            
            /* convert to FFT index i.e [-N/2,N/2) to [0,N)*/
            Nx = get_fft_idx(Gvec_temp[0], FFT_dims[0]); 
            Ny = get_fft_idx(Gvec_temp[1], FFT_dims[1]); 
            Nz = get_fft_idx(Gvec_temp[2], FFT_dims[2]); 

            /*Gvecs must be in the full FFT grid*/
            if (Nx<0 || Ny<0 || Nz <0) error_msg("G vectors and FFT grid are incompatible");
            if (Nx>=(int)FFT_dims[0] || Ny>=(int)FFT_dims[1] || Nz>=(int)FFT_dims[2])
                                                error_msg("Incompatible FFT grid");
            /* set the wfc value */
            wfcGsphere[iset*npw+ipw] = wfcGbox[iset*nFFTs + Nx*fft_strides[0] + \
                                        Ny*fft_strides[1] + Nz*fft_strides[2]];
        } // ipw
    } // iset
}

