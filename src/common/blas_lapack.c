/*
This file contails blas wrappers
*/
#define CBLAS_INT int
// note the above header must be called before #include "cblas.h"
#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <tgmath.h>

#include "cblas.h"
#include "elphC.h"
#include "error.h"
#include "numerical_func.h"

#define BLAS_INT_MAX INT_MAX

#if defined(COMPILE_ELPH_DOUBLE)
#define CBLAS_cmplx(FUN_NAME) CBLAS_cmplx_hidden(FUN_NAME)
#define CBLAS_cmplx_hidden(FUN_NAME) cblas_z##FUN_NAME

#define CBLAS_float(FUN_NAME) CBLAS_float_hidden(FUN_NAME)
#define CBLAS_float_hidden(FUN_NAME) cblas_d##FUN_NAME

#define LAPACK_cmplx(FUN_NAME) LAPACK_cmplx_hidden(FUN_NAME)
#define LAPACK_cmplx_hidden(FUN_NAME) z##FUN_NAME##_

#else
#define CBLAS_cmplx(FUN_NAME) CBLAS_cmplx_hidden(FUN_NAME)
#define CBLAS_cmplx_hidden(FUN_NAME) cblas_c##FUN_NAME

#define CBLAS_float(FUN_NAME) CBLAS_float_hidden(FUN_NAME)
#define CBLAS_float_hidden(FUN_NAME) cblas_s##FUN_NAME

#define LAPACK_cmplx(FUN_NAME) LAPACK_cmplx_hidden(FUN_NAME)
#define LAPACK_cmplx_hidden(FUN_NAME) c##FUN_NAME##_
#endif

static enum CBLAS_TRANSPOSE get_gemmn_T(char Trans);

// =================================================================
// NM : Fix me : Check for overflows

void matmul_cmplx(const char TransA, const char TransB, const ELPH_cmplx* arr_A,
                  const ELPH_cmplx* arr_B, ELPH_cmplx* arr_C,
                  const ELPH_cmplx alpha, const ELPH_cmplx beta,
                  const ND_int ldA, const ND_int ldB, const ND_int ldC,
                  const ND_int m, const ND_int n, const ND_int k)
{
    // first check for interger overflows before passing to blas library
    //
    CHECK_OVERFLOW_ERROR(ldA, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(ldB, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(ldC, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(m, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(n, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(k, BLAS_INT_MAX);

    CBLAS_cmplx(gemm)(CblasRowMajor, get_gemmn_T(TransA), get_gemmn_T(TransB),
                      m, n, k, &alpha, arr_A, ldA, arr_B, ldB, &beta, arr_C,
                      ldC);
}

void matmul_float(const char TransA, const char TransB, const ELPH_float* arr_A,
                  const ELPH_float* arr_B, ELPH_float* arr_C,
                  const ELPH_float alpha, const ELPH_float beta,
                  const ND_int ldA, const ND_int ldB, const ND_int ldC,
                  const ND_int m, const ND_int n, const ND_int k)
{
    // first check for interger overflows before passing to blas library
    //
    CHECK_OVERFLOW_ERROR(ldA, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(ldB, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(ldC, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(m, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(n, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(k, BLAS_INT_MAX);

    CBLAS_float(gemm)(CblasRowMajor, get_gemmn_T(TransA), get_gemmn_T(TransB),
                      m, n, k, alpha, arr_A, ldA, arr_B, ldB, beta, arr_C, ldC);
}

int diagonalize_hermitian(const char jobz, const char uplo, const ND_int N,
                          const ND_int LDA, ELPH_cmplx* A, ELPH_float* w)
{
    // Expect in row major format.
    // Note that ith eigen vector is A[i*ldA]
    if (N < 1)
    {
        return 0;
    }

    CHECK_OVERFLOW_ERROR(N, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(LDA, BLAS_INT_MAX);

    CBLAS_INT n = N;
    CBLAS_INT lda = LDA;
    CBLAS_INT info = 0;
    char jobz_ = jobz;
    // we need to swap as we are in row major.
    _Bool is_upper = (uplo == 'U' || uplo == 'u');
    char uplo_ = is_upper ? 'L' : 'U';

    ELPH_cmplx wkopt;
    ELPH_cmplx* work = NULL;
    CBLAS_INT lwork = -1;

    ELPH_float* rwork = malloc(3 * N * sizeof(*rwork));
    if (!rwork)
    {
        return -1;
    }

    // lapack is coloum major, so we need to send the transpose i.e we conj
    for (ND_int i = 0; i < N; ++i)
    {
        for (ND_int j = 0; j < N; ++j)
        {
            // The if condition ensure we only touch upper or lower
            if ((is_upper && j >= i) || (!is_upper && j <= i))
            {
                ND_int idx = i * LDA + j;
                A[idx] = conj(A[idx]);
            }
        }
    }

    LAPACK_cmplx(heev)(&jobz_, &uplo_, &n, A, &lda, w, &wkopt, &lwork, rwork,
                       &info);
    if (info)
    {
        free(rwork);
        return info;
    }

    lwork = (CBLAS_INT)rint(creal(wkopt) * 1.005);

    work = malloc(lwork * sizeof(*work));
    if (!work)
    {
        free(rwork);
        return -1;
    }

    LAPACK_cmplx(heev)(&jobz_, &uplo_, &n, A, &lda, w, work, &lwork, rwork,
                       &info);

    free(work);
    free(rwork);
    return info;
}

// =================================================================
static enum CBLAS_TRANSPOSE get_gemmn_T(char Trans)
{
    if (Trans == 'N')
    {
        return CblasNoTrans;
    }
    else if (Trans == 'T')
    {
        return CblasTrans;
    }
    else if (Trans == 'C')
    {
        return CblasConjTrans;
    }
    else
    {
        error_msg("Can only take 'C' or 'T' or 'N' for Trans input");
        return CblasNoTrans;
    }
}
