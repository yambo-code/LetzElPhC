#pragma once

#include "../elphC.h"

void Sorted_qpts_idxs(const ND_int nqpts, ELPH_float* qpts, ND_int* indices);

void rearrange_qpt_grid(const ND_int nqpts, const ELPH_cmplx* in_buf,
                        const ND_int* idx, ELPH_cmplx* restrict out_buf);

void find_qpt_grid(const ND_int nqpts, const ELPH_float* qpts,
                   ND_int* restrict q_grid);
