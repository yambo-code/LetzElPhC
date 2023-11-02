/* determines the simd aligment length */

#include "fft.h"

#define ELPH_FFTW_SIMD_LEN 64

ND_int alignment_len(void)
{   
    /* get the aligment len in units of sizeof(ELPH_cmplx) bytes*/
    // return -1 incase of error which means we should turn off simd
    ND_int align_len = ELPH_FFTW_SIMD_LEN+10;
    ND_int temp_len = sizeof(ELPH_cmplx)*ELPH_FFTW_SIMD_LEN;

    ELPH_cmplx * temp = fftw_fun(malloc)(temp_len);

    // first get the alignment value of temp;
    int buf_alignment = fftw_fun(alignment_of)((void *)temp);

    for (int i = 1 ; i < ELPH_FFTW_SIMD_LEN; ++i)
    {
        ELPH_cmplx * ptr = temp+i;
        int ptr_align = fftw_fun(alignment_of)((void *)ptr);
        if (ptr_align == buf_alignment)
        {
            align_len = i ;
            break;
        }
    }
    
    fftw_fun(free)(temp);

    if (align_len >= ELPH_FFTW_SIMD_LEN) error_msg("Error determining simd len");
    return align_len ;
}