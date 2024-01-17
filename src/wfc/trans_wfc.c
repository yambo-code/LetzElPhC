/*
* This file contains function that applies translation operator
 to the wavefunction
*/

#include "wfc.h"

#define dot3_macro(a,b) ((a)[0]*(b)[0] +  (a)[1]*(b)[1] + (a)[2]*(b)[2])

void apply_trans_wfc(const ELPH_float * trans_vec, const ELPH_float * kvec, \
                     const ND_int nsets, const ND_int npw, const ELPH_float * gvecs, \
                     ELPH_cmplx * restrict wfc_G, const bool conjugate)
{                                           
    /*
    This is an in-place function. wfc_G (nsets,npw) will be over written.
    If conjugate is true, the output is conjugated

    T(a) C_G * e^{i(G+k).r} = C_G * e^{i(G+k).r} * e^{-i(G+k).a}
    note : trans_vec, kvec, gvecs can be in crystal/cart coordinates but 
    we must use the same for all. kvec and gvecs must be in 2*pi units 
    when using cart coordinates.
    */

    // first compute the phase from k i.e e^{-ik*a}
    ELPH_cmplx kphase = cexp(-I*2*ELPH_PI*dot3_macro(kvec,trans_vec));
    
    for (ND_int iset = 0; iset < nsets; ++iset)
    {   
        ELPH_cmplx * restrict wfc_G_tmp = wfc_G + iset*npw; 
        
        for (ND_int ig = 0; ig < npw; ++ig)
        {   
            const ELPH_float * gvec_tmp = gvecs + 3*ig;

            ELPH_cmplx gphase = kphase*cexp(-I*2*ELPH_PI*dot3_macro(gvec_tmp,trans_vec)); 
            
            wfc_G_tmp[ig] *= gphase;
            
            if (conjugate)  wfc_G_tmp[ig] = conj(wfc_G_tmp[ig]) ;
        }
    }


}








