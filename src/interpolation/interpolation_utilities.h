#pragma once

#include "elphC.h"

void Sorted_qpts_idxs(const ND_int nqpts, ELPH_float* qpts, ND_int* indices);

void rearrange_qpt_grid(const ND_int nqpts, const ELPH_cmplx* in_buf,
                        const ND_int* idx, ELPH_cmplx* out_buf);

void find_qpt_grid(const ND_int nqpts, const ELPH_float* qpts, ND_int* q_grid);

ELPH_float* parse_qpt_entries(const char* filename, ND_int* count_out);

void fft_q2R(ELPH_cmplx* data, const ND_int* qgrid, const ND_int nsets);

void fft_R2q_dvscf(const ELPH_cmplx* dataR, const ELPH_float* qpt_crys,
                   const ND_int* qgrid, const ND_int natom, const ND_int nsets,
                   const ND_int* ws_vecs, const ND_int nws,
                   const ND_int* nws_vecs, ELPH_cmplx* dataq);

void fft_R2q_dyn(const ELPH_cmplx* dataR, const ELPH_float* qpt_crys,
                 const ND_int* qgrid, const ND_int natom, const ND_int* ws_vecs,
                 const ND_int nws, const ND_int* nws_vecs, ELPH_cmplx* dataq);
