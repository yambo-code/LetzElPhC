#pragma once
#include "../elphC.h"

#if defined(ELPH_OMP_PARALLEL_BUILD)
#include <omp.h>
#define ELPH_OMP_PAR_FOR_SIMD _Pragma("omp parallel for simd")

#define ELPH_OMP_PAR_FOR _Pragma("omp parallel for")

#define ELPH_OMP_FOR _Pragma("omp for")

#define ELPH_OMP_PAR_SIMD _Pragma("omp simd")

#define ELPH_OMP_PAR_CRITICAL _Pragma("omp critical")

#define ELPH_OMP_PAR _Pragma("omp parallel")

#define ELPH_OMP_ATOMIC _Pragma("omp atomic")

#define ELPH_OMP_PAR_COLLAPSE_3 _Pragma("omp for collapse(3)")

#define ELPH_OMP_SINGLE _Pragma("omp single")
#else
#define ELPH_OMP_PAR_FOR_SIMD

#define ELPH_OMP_PAR_FOR

#define ELPH_OMP_FOR

#define ELPH_OMP_PAR_SIMD

#define ELPH_OMP_PAR_CRITICAL

#define ELPH_OMP_PAR

#define ELPH_OMP_ATOMIC

#define ELPH_OMP_PAR_COLLAPSE_3

#define ELPH_OMP_SINGLE
#endif
