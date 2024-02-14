/* This is a routine to perform 3D FFT convolution of potential and wavefunction
    i.e FFT(V(r)*psi(r)). */

#include "fft.h"

void fft_convolution3D(struct ELPH_fft_plan* plan, const ND_int nspinor,
                       ND_int nmag, const ELPH_cmplx* Vpotr,
                       const ELPH_cmplx* psir, ELPH_cmplx* wfcG,
                       const bool conjugate)
{
    /*
    a) FFt along X for V(r)*psi(r) followed by fft along z

    if nmag = 1 or 2, "i,si=>si"
    if nmag = 4, "rsxyz,sxyz->rxyz"

    b) scatter Gx,Gy to different cpus to perform fft along x and y.
        *Note that we use scatter instead of alltoallv so that some FFTs
        can be computed while being communicating

    c) Perform FFT along Z

    d) box to sphere
    */

    // basic checks
    if (nspinor == 1)
    {
        nmag = 1;
    }
    if (nmag != 1 && nmag != 4)
    {
        error_msg("Error wrong nmag value.");
    }
    if (nmag == 4 && nspinor != 2)
    {
        error_msg("Error wrong nmag and nspinor combination.");
    }
    if (nspinor > 2)
    {
        error_msg("nspinor must be <= 2.");
    }
    /*
    In case of nmag  = 2 (i.e lsda), we have two V(r) for each spin component.
    so this can be treated as nmag = 1 case
    */

    // first some basic stuff
    ND_int Nx = plan->fft_dims[0];
    ND_int Ny = plan->fft_dims[1];
    ND_int Nz = plan->fft_dims[2];

    ELPH_float norm = Nx * Ny * Nz;
    norm = 1.0 / norm;

    ND_int Ny_stride = plan->nzloc * Ny;

    ND_int fft_buf_size = Nx * Ny * plan->nzloc;

    for (ND_int ispinor = 0; ispinor < nspinor; ++ispinor)
    {
        ELPH_cmplx* restrict wfcG_tmp = wfcG + ispinor * plan->ngvecs_loc;
        const ELPH_cmplx* dV_r = Vpotr;

        // we store V(r)*psi(r) in plan->fft_data and perform FFT
        /* now we do blocking to reuse cache. we divide into Ny blocks and
            perform convolution and FFT */
        if (nmag == 1)
        {
            const ELPH_cmplx* psi_spinor = psir + ispinor * fft_buf_size;
            for (ND_int iy = 0; iy < Ny; ++iy)
            {
                ND_int iyshift = iy * plan->nzloc;
                for (ND_int ix = 0; ix < Nx; ++ix)
                {
                    ND_int ixshift = ix * Ny * plan->nzloc + iyshift;
                    ELPH_cmplx* restrict dvpsi_xyz = plan->fft_data + ixshift;
                    const ELPH_cmplx* psi_xyz = psi_spinor + ixshift;
                    const ELPH_cmplx* dV_xyz = dV_r + ixshift;
                    for (ND_int iz = 0; iz < plan->nzloc; ++iz)
                    {
                        dvpsi_xyz[iz] = psi_xyz[iz] * dV_xyz[iz];
                    }
                }
                // perform the fft along X
                ELPH_cmplx* dvpsi_ptr = plan->fft_data + iyshift;
                ND_int iax = fftw_fun(alignment_of)((void*)(dvpsi_ptr));
                iax /= sizeof(ELPH_cmplx);
                fftw_fun(execute_dft)(plan->cplan_x[iax], dvpsi_ptr, dvpsi_ptr);
            }
        }
        else
        {
            if (ispinor == 1)
            {
                dV_r = Vpotr + 2 * fft_buf_size;
            } // [1,2],[3,4]

            for (ND_int iy = 0; iy < Ny; ++iy)
            {
                ND_int iyshift = iy * plan->nzloc;
                for (ND_int ix = 0; ix < Nx; ++ix)
                {
                    ND_int ixshift = ix * Ny * plan->nzloc + iyshift;
                    ELPH_cmplx* restrict dvpsi_out = plan->fft_data + ixshift;
                    const ELPH_cmplx* psi0 = psir + ixshift;
                    const ELPH_cmplx* psi1 = psir + ixshift + fft_buf_size;

                    const ELPH_cmplx* dV0 = dV_r + ixshift;
                    const ELPH_cmplx* dV1 = dV_r + ixshift + fft_buf_size;
                    for (ND_int iz = 0; iz < plan->nzloc; ++iz)
                    {
                        dvpsi_out[iz] = psi0[iz] * dV0[iz] + psi1[iz] * dV1[iz];
                    }
                }
                // perform the fft along X
                ELPH_cmplx* dvpsi_ptr = plan->fft_data + iyshift;
                ND_int iax = fftw_fun(alignment_of)((void*)(dvpsi_ptr));
                iax /= sizeof(ELPH_cmplx);
                fftw_fun(execute_dft)(plan->cplan_x[iax], dvpsi_ptr, dvpsi_ptr);
            }
        }

        // a) (ii) FFT along Y
        for (ND_int ix = 0; ix < plan->fft_dims[0]; ++ix)
        {
            if (!(plan->gx_inGvecs[ix]))
            {
                continue;
            }
            ELPH_cmplx* wfcr_tmp_y = plan->fft_data + ix * Ny_stride; //(Gx,Ny,Nz_log)
            ND_int iax = fftw_fun(alignment_of)((void*)wfcr_tmp_y);
            iax /= sizeof(ELPH_cmplx);
            fftw_fun(execute_dft)(plan->fplan_y[iax], wfcr_tmp_y, wfcr_tmp_y);
        }

        // b) (i) pack the data for transpose
        for (ND_int ixy = 0; ixy < plan->nGxy; ++ixy)
        {
            ND_int Gx = plan->Gxy_total[2 * ixy];
            ND_int Gy = plan->Gxy_total[2 * ixy + 1];
            // [Gx,Gy,0] =  Gx*Ny*Nzloc + Gy*Nzloc
            ELPH_cmplx* xy_buf = plan->fft_data + plan->nzloc * (Gy + Gx * Ny);
            // (Nxy,Zloc)
            memcpy(plan->nz_buf + ixy * plan->nzloc, xy_buf,
                   sizeof(ELPH_cmplx) * plan->nzloc);
        }
        // b) (ii)  transpose the data
        fwd_transpose(plan); // nz_buf has (NGxy_loc,Nz)

        // c) perform fft along z
        fftw_fun(execute)(plan->fplan_z);

        // d) box to sphere
        ND_int igvec = 0;
        for (ND_int ixy = 0; ixy < plan->nGxyloc; ++ixy)
        {
            ELPH_cmplx* zfft_ptr = plan->nz_buf + ixy * Nz;
            for (ND_int ig = 0; ig < plan->ngxy_z[ixy]; ++ig)
            {
                ND_int Gz = plan->gvecs[3 * igvec + 2];
                if (Gz < 0)
                {
                    Gz += Nz;
                }
                wfcG_tmp[igvec] = zfft_ptr[Gz] * norm;
                if (conjugate)
                {
                    wfcG_tmp[igvec] = conj(wfcG_tmp[igvec]);
                }
                ++igvec;
            }
        }
        if (igvec != plan->ngvecs_loc)
        {
            error_msg("Gvec mismatch in fwd execute. ");
        }
    }
}
