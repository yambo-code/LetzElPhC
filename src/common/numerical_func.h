#pragma once
#include "../elphC.h"
#include "dtypes.h"
#include "error.h"

#define Function(FUN_NAME, TYPE_SMALL) Function_HIDDEN(FUN_NAME, TYPE_SMALL)
#define Function_HIDDEN(FUN_NAME, TYPE_SMALL) ELPH_##TYPE_SMALL##FUN_NAME

ND_int find_kidx_in_list(ND_int nkpts, const ELPH_float* kpts_list,
                         const ELPH_float* kpt);

void get_KplusQ_idxs(const ND_int Nbz, const ELPH_float* kpoints,
                     const ELPH_float* Q_pt, int* KplusQidxs);

/* numerical_func.c */
ELPH_float legendre(int l_val, int m_val, ELPH_float x_in);
ELPH_float Ylm(int l_val, int m_val, ELPH_float* vec);
ELPH_float simpson(const ELPH_float* restrict func_vals, const ELPH_float* restrict dx,
                   ND_int npts);

ELPH_float cos_angle_bw_Vec(const ELPH_float* vec1, const ELPH_float* vec2);
ELPH_float dotVec3(const ELPH_float* vec1, const ELPH_float* vec2);
void MatVec3f(const ELPH_float* restrict Mat, const ELPH_float* restrict vec,
              const bool trans, ELPH_float* restrict out);
ELPH_cmplx Cmplxdot(const ELPH_cmplx* vec1, const ELPH_cmplx* vec2,
                    const ND_int n);
void normalize_Cmplx_vec(ELPH_cmplx* vec, const ND_int n);
ELPH_float det3x3(const ELPH_float* mat);
void reciprocal_vecs(const ELPH_float* restrict lat_vec,
                     ELPH_float* restrict blat);
void aXpY(const ND_int n, const ELPH_cmplx a, const ELPH_cmplx* restrict X,
          ELPH_cmplx* restrict Y);
void transpose3x3f(const ELPH_float* restrict inmat,
                   ELPH_float* restrict outmat);
ND_int find_maxint(ND_int* in_arr, ND_int nelements);
ELPH_float find_maxfloat(ELPH_float* in_arr, ND_int nelements);
void Gemm3x3f(const ELPH_float* restrict A, const char transA,
              const ELPH_float* restrict B, const char transB,
              ELPH_float* restrict C);

void matmul_Cmpl2x2(ELPH_cmplx* restrict mat1, ELPH_cmplx* restrict mat2,
                    ELPH_cmplx* restrict out);

int get_fft_idx(ELPH_float idx_in, int FFT_dimension);
ND_int get_miller_idx(ND_int idx_in, ND_int FFT_dimension);

/* spline.c*/
// spline interpolation functions
ELPH_float spline_interpolate(const ELPH_float x, ND_int inear,
                              const ELPH_float* restrict xi,
                              const ELPH_float* restrict yi,
                              const ELPH_float* restrict dy);

void prepare_spline(const ND_int nvals, ELPH_float* restrict xin,
                    ELPH_float* restrict yin, ELPH_float* restrict dy);

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
