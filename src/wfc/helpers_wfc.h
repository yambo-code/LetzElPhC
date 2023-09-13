/*
Header for functions defined in helpers_wfc.c
*/
#pragma once
#include "../elphC.h"


int get_fft_idx(ELPH_float idx_in, int FFT_dimension);

void rotateGvecs(const ELPH_float * Gvecs, const ELPH_float * sym, const ND_int ngvecs, 
                const  ELPH_float * lat_vec, const bool inverse, const bool crystal,
                const ELPH_float * G0, ELPH_float *  Gvec_out);