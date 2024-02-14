#include "Vnonloc.h"
/*
This file contains an function which is called in elphNonLocal function that is
used to compute the sandwich with nonlocal potential
*/

void sum_VNL_KKp(ELPH_cmplx* K_ptr, ELPH_cmplx* Kp_ptr, ELPH_cmplx* fcoeff,
                 ND_int nspin, ND_int nbnd, ND_int nspinor,
                 ELPH_cmplx* restrict out)
{
    /*This is a helper function to sum over K and K' in non local part */
    // sum_K_K' V_NL = \sum_{sigma,sigma'} npwK^\sigma[0]*npwKp^\sigma'[1:4]*f +
    // f*npwKp^\sigma'[0]*npwK^\sigma[1:4] ;
    //(4,isn,nb,spnor)  fcoeff is (nspinor, nspinor)

    ELPH_cmplx* K3_ptr = K_ptr + nspin * nbnd * nspinor; //(3,isn,nb,spnor)
    ELPH_cmplx* Kp3_ptr = Kp_ptr + nspin * nbnd * nspinor; //(3,isn,nb,spnor)

    // need to perform Kp3_ptr@f.T@K0_ptr + Kp0_ptr@f.T@K3_ptr -> (3,isn,nk,nkq)

    for (ND_int ix = 0; ix < 3; ++ix)
    {
        for (ND_int ispin = 0; ispin < nspin; ++ispin)
        {
            for (ND_int ib1 = 0; ib1 < nbnd; ++ib1) // Kp (k)
            {
                ND_int strideb1 = nspinor * ib1 + nbnd * nspinor * ispin;
                ELPH_cmplx* ib1_temp0 = Kp_ptr + strideb1;
                ELPH_cmplx* ib1_temp3 = Kp3_ptr + strideb1 + nspin * nbnd * nspinor * ix;

                for (ND_int ib2 = 0; ib2 < nbnd; ++ib2) // K (kq)
                {
                    ND_int strideb2 = nspinor * ib2 + nbnd * nspinor * ispin;
                    ELPH_cmplx* ib2_temp0 = K_ptr + strideb2;
                    ELPH_cmplx* ib2_temp3 = K3_ptr + strideb2 + nspin * nbnd * nspinor * ix;

                    //(3 * nspin, nbnd * nbnd )
                    //(nbnd*nbnd*nspin,nbnd*nbnd,nbnd,1)
                    ELPH_cmplx* restrict temp_out = out + ib2 + ib1 * nbnd + ispin * nbnd * nbnd + ix * nbnd * nbnd * nspin;

                    for (ND_int isp1 = 0; isp1 < nspinor; ++isp1) // for K
                    {
                        for (ND_int isp2 = 0; isp2 < nspinor; ++isp2) // for K'
                        {
                            *temp_out += fcoeff[nspinor * isp1 + isp2] * (ib2_temp3[isp1] * ib1_temp0[isp2] + ib2_temp0[isp1] * ib1_temp3[isp2]);
                        }
                    }
                }
            }
        }
    }
}
