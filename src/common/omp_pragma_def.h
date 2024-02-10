#pragma once

// Pragma for OMP paralleization

#if defined(ELPH_OMP_PARALLEL_BUILD)
#define ELPH_OMP_PAR_FOR_SIMD _Pragma("omp parallel for simd")

#define ELPH_OMP_PAR_FOR _Pragma("omp parallel for")

#define ELPH_OMP_PAR_SIMD _Pragma("omp simd")

#else
#define ELPH_OMP_PAR_FOR_SIMD

#define ELPH_OMP_PAR_FOR

#define ELPH_OMP_PAR_SIMD
#endif
