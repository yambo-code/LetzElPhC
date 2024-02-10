/*
This file contains function which adds local part of pseudo potential to
d/dR(Ha + Vxc) i.e induced potential
*/
#include "dvloc.h"

void add_dvscf_qe(ELPH_cmplx* restrict dVscf, const ELPH_cmplx* dVloc,
                  const struct Lattice* lattice)
{
    /*
    nffts = Nx*Ny*Nz_loc

    This function adds dVloc is added to dVscf.
    dVscf_{2x2} = d{V}*I + d{Bx}*sigma_x + d{By}*sigma_y + dBz}*sigma_z
                    =   | d{V} + d{Bz}         d{Bx} - I*d{By} |
                        | d{Bx}+ I*d{By}       d{V} - d{Bz}    |
    where Bx,By,Bz = d{Exc}/dm and V = V_Ha + d{Exc}/dn i.e V = V_Ha + V_Xc

    In QE. dvscf is stored as follows :
    nmag = 1 : dV is stored (non magnetic calc with both collinear and non
    collinear) nmag = 2 : d{V} + d{Bz}  and d{V} - d{Bz} are stored (nspin = 2)
    nmag = 4 : d{V}, d{Bx},d{By},d{Bz} are stored (magnetic, non collinear)

    dVscf -> (nmodes, nmag, nffts)
    dVloc -> (nmodes, nffts)
    */
    const ND_int nmag = lattice->nmag;
    const ND_int nmodes = 3 * lattice->natom;
    const ND_int nffts =
        lattice->fft_dims[0] * lattice->fft_dims[1] * lattice->nfftz_loc;

    ND_int mag_iter = 1;

    if (nmag == 2)
    {
        mag_iter = 2;
    }
    for (ND_int iv = 0; iv < nmodes; ++iv)
    {
        ELPH_cmplx* restrict dVscf_mode = dVscf + iv * nffts * nmag;
        ELPH_cmplx* restrict dVbare = dVloc + iv * nffts;

        /* add local bare part of pseudo potential (i,e dVloc) to V_Ha + V_xc */
        for (ND_int im = 0; im < mag_iter; ++im)
        {
            ELPH_cmplx* restrict temp_ptr = dVscf_mode + im * nffts;
            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int i = 0; i < nffts; ++i)
            {
                temp_ptr[i] += dVbare[i];
            }
        }

        if (nmag == 4)
        {
            /* dvscf_{2x2} = vloc*I + Bx*sigma_x + By*sigma_y + Bz*sigma_z
            | V + Bz         Bx - I*By |
            | Bx+ I*By       V - Bz    |
            where Bx,By,Bz = d{Exc}/dm and Vxc = d{Exc}/dn
            */

            /*
            This is an inplace operation to avoid creating a large buffer
            */

            ELPH_cmplx* restrict Vxc_loc = dVscf_mode + 0 * nffts;
            ELPH_cmplx* restrict Bx_loc = dVscf_mode + 1 * nffts;
            ELPH_cmplx* restrict By_loc = dVscf_mode + 2 * nffts;
            ELPH_cmplx* restrict Bz_loc = dVscf_mode + 3 * nffts;

            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int i = 0; i < nffts; ++i)
            {
                ELPH_cmplx Vxc_t, Bx_t, By_t, Bz_t;
                Vxc_t = Vxc_loc[i];
                Bx_t = Bx_loc[i];
                By_t = By_loc[i];
                Bz_t = Bz_loc[i];

                Vxc_loc[i] = Vxc_t + Bz_t;
                Bx_loc[i] = Bx_t - I * By_t;
                By_loc[i] = Bx_t + I * By_t;
                Bz_loc[i] = Vxc_t - Bz_t;
            }
        }
    }
}
