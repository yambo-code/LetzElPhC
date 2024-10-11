/*
This file contains function su2rotate
*/

/*
This function rotates the wave-function in spin space
*/

#include "../common/omp_pragma_def.h"
#include "../elphC.h"
#include "wfc.h"

void su2rotate(const int nspinor, const ND_int npw, const ND_int nsets,
               const ELPH_cmplx* restrict su2mat, ELPH_cmplx* restrict wfc)
{
    /*
    \sum_j su_ij *wfc_njG -> wfc_niG
    Note this is inplace operation, so it will modify the input
    */

    if (nspinor != 2)
    {
        return;
    }

    ND_int stride_wfc = nspinor * npw;  // i.e 2*npw

    for (ND_int iset = 0; iset < nsets; ++iset)
    {
        ELPH_cmplx* restrict wfc_up = wfc + iset * stride_wfc;
        ELPH_cmplx* restrict wfc_dw = wfc + iset * stride_wfc + npw;

        ELPH_OMP_PAR_FOR_SIMD
        for (ND_int ipw = 0; ipw < npw; ++ipw)
        {
            ELPH_cmplx temp = wfc_up[ipw] * su2mat[0] + wfc_dw[ipw] * su2mat[1];
            wfc_dw[ipw] = wfc_up[ipw] * su2mat[2] + wfc_dw[ipw] * su2mat[3];
            wfc_up[ipw] = temp;
        }
    }
}
