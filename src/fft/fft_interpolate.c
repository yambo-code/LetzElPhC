// This file contains functions to perform fourier interpolate

#include "../common/numerical_func.h"
#include "../elphC.h"
#include "fft.h"
// complex.h must be before fftw3.h
#include <complex.h>
#include <fftw3.h>

void fft_interpolate(struct fft_interpolate_plan* iplan)
{
    /*
    data_co/fi : Coarse/fine grid data (will be destroyed after output)
    iplan : interpolation plan
    */

    /*
     * 1) Perform FFT from q->R
     * 2) Pad with zeros on the fine grid
     * 3) Perform inverse FFT on fine grid (R->q)
     */

    const ND_int* fft_dims_co = iplan->fft_dims_co;
    const ND_int* fft_dims_fi = iplan->fft_dims_fi;

    ND_int nfft_co = fft_dims_co[0] * fft_dims_co[1] * fft_dims_co[2];
    ND_int nfft_fi = fft_dims_fi[0] * fft_dims_fi[1] * fft_dims_fi[2];

    ELPH_float norm_fft_co = 1.0 / nfft_co;

    // 1) q -> R FFT
    fftw_fun(execute)(iplan->fft_plan_co);

    // 2) pad the buffer
    //  Zero out the fine interpolation buffer"
    for (ND_int i = 0; i < nfft_fi; ++i)
    {
        iplan->data_fi[i] = 0.0;
    }

    // strides [Ny*Nz,Nz,1]
    ND_int NyNz_co = fft_dims_co[1] * fft_dims_co[2];
    ND_int NyNz_fi = fft_dims_fi[1] * fft_dims_fi[2];

    for (ND_int i = 0; i < nfft_co; ++i)
    {
        // iplan->data_fi[i] = 0.0;
        ND_int ix = i / NyNz_co;
        ND_int iy = (i % NyNz_co) / fft_dims_co[2];
        ND_int iz = (i % NyNz_co) % fft_dims_co[2];

        // convert [0,N)-> [-N/2,N/2)
        ix = get_miller_idx(ix, fft_dims_co[0]);
        iy = get_miller_idx(iy, fft_dims_co[1]);
        iz = get_miller_idx(iz, fft_dims_co[2]);

        // convert back to fft idx in fine grid
        ix = get_fft_idx(ix, fft_dims_fi[0]);
        iy = get_fft_idx(iy, fft_dims_fi[1]);
        iz = get_fft_idx(iz, fft_dims_fi[2]);

        // multiply with the normalizition factor.
        iplan->data_fi[ix * NyNz_fi + iy * fft_dims_fi[2] + iz] =
            norm_fft_co * iplan->data_co[i];
    }

    // 3) invFFT on fine grid
    fftw_fun(execute)(iplan->ifft_plan_fi);

    return;
}
