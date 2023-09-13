// computes Vloc(r)* psi(r)
// /* Static helper functions */
#include "dvloc.h"

void dvpsi(const ND_int nmag, const ND_int nspinor, const ND_int nFFT, 
            const ELPH_cmplx * restrict dV_r, ELPH_cmplx * restrict psi_r)
{
    /*
    Computes V(r)*psi(r) on box grid
    */

    /*
    if nmag = 1 or 2, "xyz,sxyz->sxyz"
    if nmag = 4, "rsxyz,sxyz->rxyz"

    nFFT = local number of ffts vectors present in this cpu
    */
    if (nmag == 1 || nmag == 2 )
    {
        for (ND_int is = 0; is<nspinor ; ++is)
        {   
            ELPH_cmplx * restrict tmp_ptr = psi_r + is*nFFT ;

            //"xyz,sxyz->sxyz"
            
            ELPH_OMP_PAR_SIMD
            for (ND_int i = 0 ; i < nFFT; ++i ) tmp_ptr[i] *= dV_r[i] ;
        }
    }
    else if (nmag == 4)
    {   
        if (nspinor != 2 ) error_msg("IF nmag = 4 then nspinor must be 2 .");
        
        /* Hard coded for nspinor == 2*/
        //"rsxyz,sxyz->rxyz"
        ELPH_cmplx * restrict psir0    = psi_r ;
        ELPH_cmplx * restrict psir1    = psi_r + nFFT  ;

        const ELPH_cmplx * restrict dV_r0    = dV_r   ;
        const ELPH_cmplx * restrict dV_r1    = dV_r + 1*nFFT  ;
        const ELPH_cmplx * restrict dV_r2    = dV_r + 2*nFFT  ;
        const ELPH_cmplx * restrict dV_r3    = dV_r + 3*nFFT  ;

        ELPH_OMP_PAR_SIMD
        for (ND_int i = 0 ; i < nFFT; ++i )
        {   ELPH_cmplx p_temp;
            p_temp   = dV_r0[i]*psir0[i]  +  dV_r1[i]*psir1[i]  ;
            psir1[i] = dV_r2[i]*psir0[i]  +  dV_r3[i]*psir1[i]  ;
            psir0[i] = p_temp;
        }
    }
    else error_msg("nmag != 4. only npsin = 1 or 2 and noncolin =.true. supported ");

}
