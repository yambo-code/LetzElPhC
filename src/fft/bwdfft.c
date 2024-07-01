/*
This is a routine to perform 3D FFT in parallel
*/
#include "fft.h"

void invfft3D(struct ELPH_fft_plan* plan, const ND_int nsets, ELPH_cmplx* wfcG,
              ELPH_cmplx* wfcr, const bool conjugate)
{
    /*
    Note: if conjugate is true, then the input is conjugated and then FFT is performed!
    please note that the input wfcG is not alterned even when conjugate is set to true.

    a) sphere to box

    b) inv FFT along Z

    c) Gather Gx, Gy and scatter Nz

    d) inv fft along Gx and Gy
    */

    // first some basic stuff
    ND_int Nx = plan->fft_dims[0];
    ND_int Ny = plan->fft_dims[1];
    ND_int Nz = plan->fft_dims[2];

    ND_int fft_buf_size = Nx * Ny * plan->nzloc;

    ND_int Ny_stride = plan->nzloc * Ny;

    for (ND_int iset = 0; iset < nsets; ++iset)
    {
        ELPH_cmplx* restrict wfcG_tmp = wfcG + iset * plan->ngvecs_loc;
        ELPH_cmplx* wfcr_tmp = wfcr + iset * fft_buf_size;

        ND_int igvec = 0;
        for (ND_int ixy = 0; ixy < plan->nGxyloc; ++ixy)
        {
            ELPH_cmplx* restrict zfft_ptr = plan->nz_buf + ixy * Nz;

            // a) sphere to box
            for (ND_int iz = 0; iz < Nz; ++iz)
            {
                zfft_ptr[iz] = 0;
            }

            for (ND_int ig = 0; ig < plan->ngxy_z[ixy]; ++ig)
            {
                ND_int Gz = plan->gvecs[3 * igvec + 2];
                if (Gz < 0)
                {
                    Gz += Nz;
                }
                zfft_ptr[Gz] = wfcG_tmp[igvec];
                if (conjugate)
                {
                    zfft_ptr[Gz] = conj(zfft_ptr[Gz]);
                }
                ++igvec;
            }
        }

        if (igvec != plan->ngvecs_loc)
        {
            error_msg("Gvec mismatch in bwd execute.");
        }

        // b) perform invfft along z
        fftw_fun(execute)(plan->bplan_z);

        // c) (i) transpose the data
        bwd_transpose(plan); // nz_buf has (NGxy,Nz_loc)

        // zero out the xy buffer and fill the fft_data buffer
        for (ND_int ixyz = 0; ixyz < fft_buf_size; ++ixyz)
        {
            wfcr_tmp[ixyz] = 0;
        }

        for (ND_int ixy = 0; ixy < plan->nGxy; ++ixy)
        {
            ND_int Gx = plan->Gxy_total[2 * ixy];
            ND_int Gy = plan->Gxy_total[2 * ixy + 1];
            // [Gx,Gy,0] =  Gx*Ny*Nzloc + Gy*Nzloc
            ELPH_cmplx* xy_buf = wfcr_tmp + plan->nzloc * (Gy + Gx * Ny);
            // (Nxy,Zloc)
            memcpy(xy_buf, plan->nz_buf + ixy * plan->nzloc,
                   sizeof(ELPH_cmplx) * plan->nzloc);
        }

        // c) (ii) FFT along Y
        for (ND_int ix = 0; ix < plan->fft_dims[0]; ++ix)
        {
            if (!(plan->gx_inGvecs[ix]))
            {
                continue;
            }

            ELPH_cmplx* wfcr_tmp_y = wfcr_tmp + ix * Ny_stride; //(Gx,Ny,Nz_log)
            ND_int iax = fftw_fun(alignment_of)((void*)wfcr_tmp_y);
            iax /= sizeof(ELPH_cmplx);
            fftw_fun(execute_dft)(plan->bplan_y[iax], wfcr_tmp_y, wfcr_tmp_y);
        }

        ND_int ia = fftw_fun(alignment_of)((void*)wfcr_tmp);
        ia /= sizeof(ELPH_cmplx);

        // d) FFT along X
        fftw_fun(execute_dft)(plan->bplan_x[ia], wfcr_tmp, wfcr_tmp);
    }
}
