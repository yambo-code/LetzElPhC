#pragma once
#include <stdbool.h>

#include "elphC.h"

ND_int find_kidx_in_list(ND_int nkpts, const ELPH_float* kpts_list,
                         const ELPH_float* kpt);

void get_KplusQ_idxs(const ND_int Nbz, const ELPH_float* kpoints,
                     const ELPH_float* Q_pt, int* KplusQidxs);

// macro functions
#define dot3_macro(a, b) ((a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2])

#if !defined(MIN)
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#if !defined(MAX)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

#if !defined(ARRAY_LEN)
#define ARRAY_LEN(arr) (sizeof(arr) / sizeof((arr)[0]))
#endif

/* numerical_func.c */
ELPH_float legendre(int l_val, int m_val, ELPH_float x_in);
ELPH_float Ylm(int l_val, int m_val, ELPH_float* vec);
ELPH_float simpson(const ELPH_float* func_vals, const ELPH_float* dx,
                   ND_int npts);

ELPH_float cos_angle_bw_Vec(const ELPH_float* vec1, const ELPH_float* vec2);

void MatVec3f(const ELPH_float* Mat, const ELPH_float* vec, const bool trans,
              ELPH_float* restrict out);

ELPH_cmplx Cmplxdot(const ELPH_cmplx* vec1, const ELPH_cmplx* vec2,
                    const ND_int n);
void normalize_Cmplx_vec(ELPH_cmplx* vec, const ND_int n);
ELPH_float det3x3(const ELPH_float* mat);

void reciprocal_vecs(const ELPH_float* lat_vec, ELPH_float* restrict blat);

void aXpY(const ND_int n, const ELPH_cmplx a, const ELPH_cmplx* X,
          ELPH_cmplx* Y);

void transpose3x3f(const ELPH_float* inmat, ELPH_float* restrict outmat);

void transpose3x3f_inplace(ELPH_float* mat);

ND_int find_maxint(ND_int* in_arr, ND_int nelements);
ELPH_float find_maxfloat(ELPH_float* in_arr, ND_int nelements);

void Gemm3x3f(const ELPH_float* A, const char transA, const ELPH_float* B,
              const char transB, ELPH_float* restrict C);

void matmul_Cmpl2x2(ELPH_cmplx* mat1, ELPH_cmplx* mat2,
                    ELPH_cmplx* restrict out);

int get_fft_idx(ELPH_float idx_in, int FFT_dimension);
ND_int get_miller_idx(ND_int idx_in, ND_int FFT_dimension);

// swap functions
void swap_ints(int* a, int* b);
void swap_floats(ELPH_float* a, ELPH_float* b);

/* spline.c*/
// spline interpolation functions
ELPH_float spline_interpolate(const ELPH_float x, ND_int inear,
                              const ELPH_float* xi, const ELPH_float* yi,
                              const ELPH_float* dy);

void prepare_spline(const ND_int nvals, ELPH_float* xin, ELPH_float* yin,
                    ELPH_float* dy);

void matmul_cmplx(const char TransA, const char TransB, const ELPH_cmplx* arr_A,
                  const ELPH_cmplx* arr_B, ELPH_cmplx* arr_C,
                  const ELPH_cmplx alpha, const ELPH_cmplx beta,
                  const ND_int ldA, const ND_int ldB, const ND_int ldC,
                  const ND_int m, const ND_int n, const ND_int k);

void matmul_float(const char TransA, const char TransB, const ELPH_float* arr_A,
                  const ELPH_float* arr_B, ELPH_float* arr_C,
                  const ELPH_float alpha, const ELPH_float beta,
                  const ND_int ldA, const ND_int ldB, const ND_int ldC,
                  const ND_int m, const ND_int n, const ND_int k);

int diagonalize_hermitian(const char jobz, const char uplo, const ND_int N,
                          const ND_int LDA, ELPH_cmplx* A, ELPH_float* w);

ND_int orthogonal_projection(const ND_int M, const ND_int N, const ND_int LDA,
                             ELPH_float* A, ELPH_float* x0,
                             const ELPH_float tol);

int lsmr_solver(ND_int m, ND_int n,
                void (*matvec)(const int, const double* restrict,
                               double* restrict, void*),
                void* userdata, const double* restrict b, double damp,
                double atol, double btol, double conlim, ND_int maxiter,
                double* restrict x, ND_int* itn_out);
