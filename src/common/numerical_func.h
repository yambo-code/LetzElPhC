#pragma once
#include "../elphC.h"
#include "data_structs.h"

#define Function(FUN_NAME, TYPE_SMALL)                  Function_HIDDEN(FUN_NAME, TYPE_SMALL)
#define Function_HIDDEN(FUN_NAME, TYPE_SMALL)           ELPH_ ## TYPE_SMALL ## FUN_NAME

/* BZ_expand.c */
bool Function(isVECpresent, Nd_cmplxS) ( const ND_array(Nd_floatS) * array,  const ELPH_float * vec, ND_int * idx  );;

void bz_expand(ND_array(Nd_floatS) * ibz_kpts, ND_array(Nd_floatS) * sym_mats, \
            ND_array(Nd_floatS) * lat_vec,ND_array(Nd_floatS) * kpoints, nd_arr_i * kmap);

void get_KplusQ_idxs(ND_array(Nd_floatS) * kpoints, int * KplusQidxs , \
                    ELPH_float * Q_pt, ND_array(Nd_floatS) * lat_vec, bool Qincrystal);



/* numerical_func.c */
ELPH_float legendre(int l_val, int m_val, ELPH_float x_in);
ELPH_float Ylm(int l_val, int m_val, ELPH_float * vec);
ELPH_float simpson(ELPH_float * restrict func_vals, ELPH_float * restrict dx, ND_int npts);
ELPH_float cos_angle_bw_Vec(const ELPH_float * vec1, const ELPH_float * vec2);
ELPH_float dotVec3(const ELPH_float * vec1, const ELPH_float * vec2);
void MatVec3f(const ELPH_float * restrict Mat, const ELPH_float * restrict vec, const bool trans, ELPH_float * restrict out);
ELPH_cmplx Cmplxdot(const ELPH_cmplx * vec1, const ELPH_cmplx * vec2, const ND_int n);
void normalize_Cmplx_vec(ELPH_cmplx * vec, const ND_int n);
ELPH_float det3x3(const ELPH_float * mat);
void reciprocal_vecs(const ELPH_float * restrict lat_vec, ELPH_float * restrict blat);
void aXpY(const ND_int n, const ELPH_cmplx a, const ELPH_cmplx * restrict X, ELPH_cmplx * restrict Y);
void transpose3x3f(const ELPH_float * restrict inmat, ELPH_float * restrict outmat);
ND_int find_maxint(ND_int * in_arr, ND_int nelements);
ELPH_float find_maxfloat(ELPH_float * in_arr, ND_int nelements);
void Gemm3x3f(const ELPH_float * restrict A, const char transA, \
                const ELPH_float * restrict B,  const char transB, \
                ELPH_float * restrict C);


