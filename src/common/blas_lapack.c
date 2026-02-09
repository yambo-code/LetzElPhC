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

#define LAPACK_float(FUN_NAME) LAPACK_float_hidden(FUN_NAME)
#define LAPACK_float_hidden(FUN_NAME) d##FUN_NAME##_

#else
#define CBLAS_cmplx(FUN_NAME) CBLAS_cmplx_hidden(FUN_NAME)
#define CBLAS_cmplx_hidden(FUN_NAME) cblas_c##FUN_NAME

#define CBLAS_float(FUN_NAME) CBLAS_float_hidden(FUN_NAME)
#define CBLAS_float_hidden(FUN_NAME) cblas_s##FUN_NAME

#define LAPACK_cmplx(FUN_NAME) LAPACK_cmplx_hidden(FUN_NAME)
#define LAPACK_cmplx_hidden(FUN_NAME) c##FUN_NAME##_

#define LAPACK_float(FUN_NAME) LAPACK_float_hidden(FUN_NAME)
#define LAPACK_float_hidden(FUN_NAME) s##FUN_NAME##_
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

ND_int orthogonal_projection(const ND_int M, const ND_int N, const ND_int LDA,
                             ELPH_float* A, ELPH_float* x0,
                             const ELPH_float tol)
{
    /*
     * Projects vector x0 onto the null space of A (Ax = 0).
     * Minimizes |x - x0| subject to Ax = 0.
     * A is M x N matrix in Row-Major order.
     * x0 is length N.
     * result is stored in-place in x0.
     * A is destroyed on output
     */
    if (M == 0 || N == 0)
    {
        return 0;
    }

    CHECK_OVERFLOW_ERROR(M, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(N, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(LDA, BLAS_INT_MAX);

    // 1. Prepare LAPACK dimensions
    // We treat the input Row-Major A (MxN) as Col-Major A^T (NxM).
    // This allows us to compute QR of A^T without physical transpose.
    CBLAS_INT m_lapack = (CBLAS_INT)N;  // Rows of AT
    CBLAS_INT n_lapack = (CBLAS_INT)M;  // Cols of AT
    CBLAS_INT lda = (CBLAS_INT)LDA;
    CBLAS_INT info = 0;

    // Pivot array (must be zero-initialized for free pivoting)
    CBLAS_INT* jpvt = malloc(n_lapack * sizeof(*jpvt));
    if (!jpvt)
    {
        return -1;
    }
    for (ND_int i = 0; i < n_lapack; ++i)
    {
        jpvt[i] = 0;
    }

    // Tau array (scalar factors of reflectors)
    CBLAS_INT min_dim = (m_lapack < n_lapack) ? m_lapack : n_lapack;
    ELPH_float* tau = malloc(min_dim * sizeof(*tau));
    if (!tau)
    {
        free(jpvt);
        return -1;
    }

    // 3. Workspace Query
    ELPH_float work_query;
    CBLAS_INT lwork = -1;

    // Query for geqp3
    LAPACK_float(geqp3)(&m_lapack, &n_lapack, A, &lda, jpvt, tau, &work_query,
                        &lwork, &info);
    work_query *= 1.005;
    CBLAS_INT lwork_qp3 = (CBLAS_INT)work_query;

    // Query for ormqr
    char side = 'L';
    //
    char trans = 'T';  // First pass we use Transpose
    CBLAS_INT one = 1;
    LAPACK_float(ormqr)(&side, &trans, &m_lapack, &one, &min_dim, A, &lda, tau,
                        x0, &m_lapack, &work_query, &lwork, &info);
    work_query *= 1.005;
    CBLAS_INT lwork_orm = (CBLAS_INT)work_query;

    trans = 'N';  // second pass we use normal
    LAPACK_float(ormqr)(&side, &trans, &m_lapack, &one, &min_dim, A, &lda, tau,
                        x0, &m_lapack, &work_query, &lwork, &info);
    work_query *= 1.005;
    CBLAS_INT lwork_orm2 = (CBLAS_INT)work_query;

    // Allocate max required workspace
    lwork = (lwork_qp3 > lwork_orm) ? lwork_qp3 : lwork_orm;
    if (lwork_orm2 > lwork)
    {
        lwork = lwork_orm2;
    }

    ELPH_float* work = malloc(lwork * sizeof(ELPH_float));
    if (!work)
    {
        free(jpvt);
        free(tau);
        return -1;
    }

    // Perform QR with Pivoting on A^T
    // AT becomes [ Q reflectors \ R ]
    LAPACK_float(geqp3)(&m_lapack, &n_lapack, A, &lda, jpvt, tau, work, &lwork,
                        &info);
    if (info != 0)
    {
        free(jpvt);
        free(tau);
        free(work);
        return info;
    }

    // Determine Rank of A
    // Rank is number of non-zero diagonal elements in R
    CBLAS_INT rank = 0;
    for (CBLAS_INT i = 0; i < min_dim; ++i)
    {
        // Diagonal of R is at AT[i + i*lda]
        if (fabs(A[i + i * lda]) > tol)
        {
            rank++;
        }
    }

    // Rotate x0 into basis of Q: y = Q^T * x0
    // x0 is overwritten by y
    trans = 'T';
    LAPACK_float(ormqr)(&side, &trans, &m_lapack, &one, &min_dim, A, &lda, tau,
                        x0, &m_lapack, work, &lwork, &info);

    // Filter: Project onto Null Space
    // The top 'rank' elements correspond to the Row Space (constraints).
    // Set them to zero to satisfy Ax=0.
    for (CBLAS_INT i = 0; i < rank; ++i)
    {
        x0[i] = 0.0;
    }

    // Rotate back to standard basis: x = Q * y_filtered
    trans = 'N';  // No transpose
    LAPACK_float(ormqr)(&side, &trans, &m_lapack, &one, &min_dim, A, &lda, tau,
                        x0, &m_lapack, work, &lwork, &info);

    free(jpvt);
    free(tau);
    free(work);

    return 0;
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
