/*
This file contains function which adds local part of pseudo potential to 
d/dR(Ha + Vxc) i.e induced potential 
*/
#include "dvloc.h"


void add_dvscf(ND_array(Nd_cmplxS) * dVscf, ND_array(Nd_cmplxS) * dVloc)
{   
    /*
    This function adds dVloc is added to dVscf.
    /* dVscf_{2x2} = d{V}*I + d{Bx}*sigma_x + d{By}*sigma_y + dBz}*sigma_z
                    =   | d{V} + d{Bz}         d{Bx} - I*d{By} |
                        | d{Bx}+ I*d{By}       d{V} - d{Bz}    |
    where Bx,By,Bz = d{Exc}/dm and V = V_Ha + d{Exc}/dn i.e V = V_Ha + V_Xc

    In QE. dvscf is stored as follows :
    nmag = 1 : dV is stored (non magnetic calc with both collinear and non collinear)
    nmag = 2 : d{V} + d{Bz}  and d{V} - d{Bz} are stored (nspin = 2)
    nmag = 4 : d{V}, d{Bx},d{By},d{Bz} are stored (magnetic, non collinear)
    
    dVscf -> (nmodes, nmag, nffts_in_this_cpu)
    dVloc -> (nmodes, nffts_in_this_cpu)
    */  
    if (dVloc->dims[1] != dVscf->dims[2]) error_msg("Wrong FFT dimensions in dVscf");

    ND_int nmodes = dVscf->dims[0];
    ND_int nmag   = dVscf->dims[1];
    ND_int mag_iter = 1;
    if (nmag == 2) mag_iter = 2;
    
    for (ND_int iv = 0 ; iv < nmodes; ++iv)
    {
        ELPH_cmplx * restrict dVscf_mode = dVscf->data + iv*dVscf->strides[0];
        ELPH_cmplx * restrict dVbare     = dVloc->data + iv*dVloc->strides[0];

        /* add local bare part of pseudo potential (i,e dVloc) to V_Ha + V_xc */
        for (ND_int im = 0; im<mag_iter; ++im)
        {   
            ELPH_cmplx * restrict temp_ptr = dVscf_mode + (im*dVscf->strides[1]) ; 
            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int i = 0 ; i< dVscf->strides[1]; ++i) temp_ptr[i] += dVbare[i] ;
        }

        if(nmag == 4)
        {
            /* dvscf_{2x2} = vloc*I + Bx*sigma_x + By*sigma_y + Bz*sigma_z
            | V + Bz         Bx - I*By |
            | Bx+ I*By       V - Bz    |
            where Bx,By,Bz = d{Exc}/dm and Vxc = d{Exc}/dn
            */

            /*
            This is an inplace operation to avoid creating a large buffer
            */

            ELPH_cmplx * restrict Vxc_loc = dVscf_mode + (0*dVscf->strides[1]) ;
            ELPH_cmplx * restrict  Bx_loc = dVscf_mode + (1*dVscf->strides[1]) ;
            ELPH_cmplx * restrict  By_loc = dVscf_mode + (2*dVscf->strides[1]) ;
            ELPH_cmplx * restrict  Bz_loc = dVscf_mode + (3*dVscf->strides[1]) ;
            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int i = 0 ; i< dVscf->strides[1]; ++i)
            {
                ELPH_cmplx Vxc_t, Bx_t, By_t, Bz_t;
                Vxc_t = Vxc_loc[i];
                Bx_t  =  Bx_loc[i];
                By_t  =  By_loc[i];
                Bz_t  =  Bz_loc[i];

                Vxc_loc[i] = Vxc_t +   Bz_t;
                Bx_loc[i]  = Bx_t  - I*By_t;
                By_loc[i]  = Bx_t  + I*By_t;
                Bz_loc[i]  = Vxc_t -   Bz_t;
            }
        }
    }
}



