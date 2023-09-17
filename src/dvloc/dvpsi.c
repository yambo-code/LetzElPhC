// computes Vloc(r)* psi(r)
// /* Static helper functions */
#include "dvloc.h"

void dvpsi(const ND_int nmag, const ND_int nspinor, const ND_int nFFT, 
            const ELPH_cmplx * restrict dV_r, const ELPH_cmplx * restrict psi_r, \
            ELPH_cmplx * restrict dVpsi_out)
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
            ELPH_cmplx * restrict tmp_psi = psi_r     + is*nFFT ;
            ELPH_cmplx * restrict tmp_out = dVpsi_out + is*nFFT ;
            //"xyz,sxyz->sxyz"
            
            ELPH_OMP_PAR_SIMD
            for (ND_int i = 0 ; i < nFFT; ++i ) tmp_out[i] = dV_r[i]*tmp_psi[i] ;
        }
    }
    else if (nmag == 4)
    {   
        if (nspinor != 2 ) error_msg("IF nmag = 4 then nspinor must be 2 .");
        
        /* Hard coded for nspinor == 2*/
        //"rsxyz,sxyz->rxyz"
        ELPH_cmplx * restrict psir0    = psi_r ;
        ELPH_cmplx * restrict psir1    = psi_r + nFFT  ;

        ELPH_cmplx * restrict dVpsir0  = dVpsi_out ;
        ELPH_cmplx * restrict dVpsir1  = dVpsi_out + nFFT  ;

        const ELPH_cmplx * restrict dV_r0    = dV_r   ;
        const ELPH_cmplx * restrict dV_r1    = dV_r + 1*nFFT  ;
        const ELPH_cmplx * restrict dV_r2    = dV_r + 2*nFFT  ;
        const ELPH_cmplx * restrict dV_r3    = dV_r + 3*nFFT  ;

        ELPH_OMP_PAR_SIMD
        for (ND_int i = 0 ; i < nFFT; ++i )
        {   
            ELPH_cmplx psi0_t = psir0[i];
            ELPH_cmplx psi1_t = psir1[i];
            dVpsir0[i]   = dV_r0[i]*psi0_t  +  dV_r1[i]*psi1_t  ;
            dVpsir1[i]   = dV_r2[i]*psi0_t  +  dV_r3[i]*psi1_t  ;
        }
    }
    else error_msg("nmag != 4. only npsin = 1 or 2 and noncolin =.true. supported ");

}


