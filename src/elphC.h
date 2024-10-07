#pragma once

#if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#ifdef __STDC_NO_COMPLEX__
#error Your compiler does not support C99 complex numbers, Please use a supported compiler.
#endif
#else
#error Your compiler does support C99 standard.
#endif

typedef long long int ND_int;
// always set this to alteast 64 bit integer
#define ELPH_MPI_ND_INT MPI_LONG_LONG_INT

// ===================================
#if defined(COMPILE_ELPH_DOUBLE)
typedef double ELPH_float;
typedef double _Complex ELPH_cmplx;

#define ELPH_MPI_float MPI_DOUBLE
#define ELPH_MPI_cmplx MPI_C_DOUBLE_COMPLEX

#define ELPH_NC4_IO_FLOAT NC_DOUBLE
#else
typedef float ELPH_float;
typedef float _Complex ELPH_cmplx;

#define ELPH_MPI_float MPI_FLOAT
#define ELPH_MPI_cmplx MPI_C_FLOAT_COMPLEX

#define ELPH_NC4_IO_FLOAT NC_FLOAT
#endif
