#pragma once

#include "elphC.h"

ND_int build_wigner_seitz_vectors(const ND_int* grid,
                                  const ELPH_float* lat_vecs, double eps,
                                  const ELPH_float* rvec_m, ND_int nrvec_m,
                                  const ELPH_float* rvec_n, ND_int nrvec_n,
                                  ND_int** ws_vecs, ND_int** ws_degen);
