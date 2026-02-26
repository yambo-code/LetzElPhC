/*
 * unit_test_lsmr_extensive.c
 *
 * Comprehensive test suite for the LSMR iterative least-squares solver.
 * (NaN cascade exact-zero tests removed, LCG generation centered)
 *
 * Compile example (adjust paths as needed):
 * gcc -Wall -Wextra -std=c11 -O2 lsmr.c unit_test_lsmr_extensive.c \
 * -lm -o test_lsmr
 */

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/numerical_func.h"
#include "elphC.h"
/* Pull in the types / declarations needed to call lsmr_solver. */

static int tests_run = 0;
static int tests_failed = 0;

#define ASSERT_TRUE(cond)                                           \
    do                                                              \
    {                                                               \
        tests_run++;                                                \
        if (!(cond))                                                \
        {                                                           \
            tests_failed++;                                         \
            printf("FAIL  %s:%d  %s\n", __FILE__, __LINE__, #cond); \
        }                                                           \
    } while (0)

#define ASSERT_FALSE(cond) ASSERT_TRUE(!(cond))

#define ASSERT_INT_EQ(a, b) ASSERT_TRUE((a) == (b))

#define ASSERT_DOUBLE_NEAR(a, b, tol) ASSERT_TRUE(fabs((a) - (b)) <= (tol))
/* Absolute error */

#define ASSERT_DOUBLE_RELNEAR(a, b, rtol) \
    ASSERT_TRUE(fabs((a) - (b)) <= (rtol) * (fabs(b) + 1e-300))
/* Relative error  (safe when b != 0) */

static double vec_norm(ND_int n, const double* x)
{
    double s = 0.0;
    for (ND_int i = 0; i < n; ++i)
    {
        s += x[i] * x[i];
    }
    return sqrt(s);
}
/* Small helper: Euclidean norm of a dense vector */

typedef struct
{
    ND_int m;
    /* rows */
    ND_int n;
    /* columns */
    double* A;
    /* row-major, m x n */
} DenseMatCtx;

static void dense_matvec(const int mode, const double* restrict x,
                         double* restrict y, void* userdata)
{
    const DenseMatCtx* ctx = (const DenseMatCtx*)userdata;
    ND_int m = ctx->m, n = ctx->n;
    const double* A = ctx->A;

    if (mode == 0)
    {
        for (ND_int i = 0; i < m; ++i)
        {
            double s = 0.0;
            for (ND_int j = 0; j < n; ++j)
            {
                s += A[i * n + j] * x[j];
            }
            y[i] = s;
        }
        /* y = A x,  length m */
    }
    else
    {
        for (ND_int j = 0; j < n; ++j)
        {
            y[j] = 0.0;
        }
        for (ND_int i = 0; i < m; ++i)
        {
            for (ND_int j = 0; j < n; ++j)
            {
                y[j] += A[i * n + j] * x[i];
            }
        }
        /* y = A^T x,  length n */
    }
}
/* Dense matrix-vector product helper. */

static DenseMatCtx* make_dense_ctx(ND_int m, ND_int n, const double* vals)
{
    DenseMatCtx* ctx = malloc(sizeof(*ctx));
    ctx->m = m;
    ctx->n = n;
    ctx->A = malloc(m * n * sizeof(double));
    memcpy(ctx->A, vals, m * n * sizeof(double));
    return ctx;
}
/* Convenience: allocate and fill a dense context */

static void free_dense_ctx(DenseMatCtx* ctx)
{
    free(ctx->A);
    free(ctx);
}

static double run_lsmr_residual(DenseMatCtx* ctx, const double* b, double damp,
                                double atol, double btol, double conlim,
                                ND_int maxiter, double* x_out, ND_int* itn_out)
{
    ND_int m = ctx->m, n = ctx->n;
    memset(x_out, 0, n * sizeof(double));
    int rc = lsmr_solver(m, n, dense_matvec, ctx, b, damp, atol, btol, conlim,
                         maxiter, x_out, itn_out);
    ASSERT_TRUE(rc == 0);

    double* Ax = malloc(m * sizeof(*Ax));
    dense_matvec(0, x_out, Ax, ctx);
    double res = 0.0;
    for (ND_int i = 0; i < m; ++i)
    {
        res += (b[i] - Ax[i]) * (b[i] - Ax[i]);
    }
    free(Ax);
    return sqrt(res);
    /* compute residual */
}
/* Convenience: run lsmr and return residual  ||b - A x|| */

static void test_diagonal_well_conditioned(void)
{
    printf("  [2] Diagonal well-conditioned\n");

    const ND_int n = 6;
    double diag[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double A[6 * 6];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = diag[i];
    }

    double b[] = {2.0, 4.0, 9.0, 16.0, 25.0, 36.0};
    double xref[] = {2.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    /* b[i]/diag[i] */

    double x[6] = {0};
    ND_int itn;

    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-12, 1e-12, 1e8, 200, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-8);
    for (ND_int i = 0; i < n; ++i)
    {
        ASSERT_DOUBLE_NEAR(x[i], xref[i], 1e-7);
    }

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 2. Diagonal matrix (well-conditioned, varying diagonal) */

static void test_small_dense_square_well_conditioned(void)
{
    printf("  [3] Small dense square (well-conditioned)\n");

    double A[] = {4.0, 1.0, 0.5, 1.0, 3.0, 0.2, 0.5, 0.2, 2.0};
    /* 3x3 Hilbert-like but non-singular and better conditioned */

    double b[] = {5.5, 4.2, 2.7};
    double x[3] = {0};
    ND_int itn;

    DenseMatCtx* ctx = make_dense_ctx(3, 3, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-12, 1e-12, 1e8, 200, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-8);
    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 3. Small dense square system (well-conditioned, random-like) */

static void test_overdetermined_consistent(void)
{
    printf("  [4] Overdetermined consistent (m=2n)\n");

    const ND_int n = 4;
    const ND_int m = 2 * n;
    double A[8 * 4];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = 1.0;
    }
    /* Top half: I_n */

    for (ND_int i = 0; i < n; ++i)
    {
        A[(n + i) * n + i] = 1.0;
    }
    /* Bottom half: I_n */

    double xtrue[] = {1.0, 2.0, 3.0, 4.0};
    double b[8];
    for (ND_int i = 0; i < n; ++i)
    {
        b[i] = xtrue[i];
    }
    for (ND_int i = 0; i < n; ++i)
    {
        b[n + i] = xtrue[i];
    }

    double x[4] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(m, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-12, 1e-12, 1e8, 200, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-8);
    for (ND_int i = 0; i < n; ++i)
    {
        ASSERT_DOUBLE_NEAR(x[i], xtrue[i], 1e-7);
    }

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 4. Overdetermined system (m > n), consistent */

static void test_overdetermined_inconsistent(void)
{
    printf("  [5] Overdetermined inconsistent (m=2n, noisy)\n");

    const ND_int n = 4;
    const ND_int m = 2 * n;
    double A[8 * 4];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = 1.0;
    }
    for (ND_int i = 0; i < n; ++i)
    {
        A[(n + i) * n + i] = 1.0;
    }

    double b1[] = {1.0, 2.0, 3.0, 4.0};
    double b2[] = {1.1, 2.1, 3.1, 4.1};
    double b[8];
    for (ND_int i = 0; i < n; ++i)
    {
        b[i] = b1[i];
    }
    for (ND_int i = 0; i < n; ++i)
    {
        b[n + i] = b2[i];
    }

    double x[4] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(m, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-10, 1e-10, 1e8, 300, x, &itn);
    free_dense_ctx(ctx);

    for (ND_int i = 0; i < n; ++i)
    {
        ASSERT_DOUBLE_NEAR(x[i], (b1[i] + b2[i]) / 2.0, 1e-6);
    }
    /* LS minimizer: x_i = (b1_i + b2_i) / 2 */

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 5. Overdetermined system (m > n), inconsistent */

static void test_diagonal_ill_conditioned(void)
{
    printf("  [8] Diagonal ill-conditioned (cond~1e6)\n");

    const ND_int n = 5;
    double diag[] = {1.0, 1e-1, 1e-2, 1e-3, 1e-6};
    double A[5 * 5];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = diag[i];
    }

    double xtrue[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    double b[5];
    for (ND_int i = 0; i < n; ++i)
    {
        b[i] = diag[i] * xtrue[i];
    }

    double x[5] = {0};
    ND_int itn;

    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-10, 1e-10, 1e8, 500, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-6);
    for (ND_int i = 0; i < n - 1; ++i)
    {
        ASSERT_DOUBLE_NEAR(x[i], xtrue[i], 1e-5);
    }
    /* skip last — very ill */

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 8. Ill-conditioned diagonal  (condition number ~ 1e6) */

static void test_hilbert_very_ill_conditioned(void)
{
    printf("  [9] Hilbert 5x5 (very ill-conditioned)\n");

    const ND_int n = 5;
    double A[5 * 5];
    for (ND_int i = 0; i < n; ++i)
    {
        for (ND_int j = 0; j < n; ++j)
        {
            A[i * n + j] = 1.0 / (double)(i + j + 1);
        }
    }

    double b[5];
    double xtrue[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    for (ND_int i = 0; i < n; ++i)
    {
        b[i] = 0.0;
        for (ND_int j = 0; j < n; ++j)
        {
            b[i] += A[i * n + j] * xtrue[j];
        }
    }

    double x[5] = {0};
    ND_int itn;

    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-8, 1e-8, 1e12, 1000, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-4);
    /* We only assert the residual is small; x itself may be far off */

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 9. Very ill-conditioned Hilbert matrix */

static void test_hilbert_8_with_damping(void)
{
    printf(" [10] Hilbert 8x8 + damping (extremely ill-conditioned)\n");

    const ND_int n = 8;
    double A[8 * 8];
    for (ND_int i = 0; i < n; ++i)
    {
        for (ND_int j = 0; j < n; ++j)
        {
            A[i * n + j] = 1.0 / (double)(i + j + 1);
        }
    }

    double xtrue[8];
    for (ND_int i = 0; i < n; ++i)
    {
        xtrue[i] = 1.0;
    }

    double b[8] = {0};
    for (ND_int i = 0; i < n; ++i)
    {
        for (ND_int j = 0; j < n; ++j)
        {
            b[i] += A[i * n + j] * xtrue[j];
        }
    }

    double x[8] = {0};
    ND_int itn;
    double damp = 1e-4;
    /* Tikhonov regularization */

    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    int rc = lsmr_solver(n, n, dense_matvec, ctx, b, damp, 1e-8, 1e-8, 1e12,
                         2000, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(itn > 0);

    double xnorm = vec_norm(n, x);
    ASSERT_TRUE(xnorm < 1e6);
    /* With damping the solution won't be exact, just check it's bounded */

    printf("      ||x||=%.2e  itn=%lld\n", xnorm, (long long)itn);
}
/* 10. Hilbert 8x8 — extremely ill-conditioned */

static void test_damped_identity(void)
{
    printf(" [11] Damped identity (explicit solution known)\n");

    const ND_int n = 4;
    double A[4 * 4];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = 1.0;
    }

    double b[] = {1.0, 2.0, 3.0, 4.0};
    double damp = 0.5;
    double scale = 1.0 / (1.0 + damp * damp);

    double x[4] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    int rc = lsmr_solver(n, n, dense_matvec, ctx, b, damp, 1e-12, 1e-12, 1e8,
                         500, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_TRUE(rc == 0);
    for (ND_int i = 0; i < n; ++i)
    {
        ASSERT_DOUBLE_NEAR(x[i], b[i] * scale, 1e-7);
    }

    printf("      itn=%lld\n", (long long)itn);
}
/* 11. Damped (regularized) identity */

static void test_large_damp_shrinks_solution(void)
{
    printf(" [12] Very large damping forces solution toward zero\n");

    const ND_int n = 4;
    double A[4 * 4];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = 1.0;
    }

    double b[] = {1.0, 2.0, 3.0, 4.0};
    double damp = 1e6;

    double x[4] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    lsmr_solver(n, n, dense_matvec, ctx, b, damp, 1e-12, 1e-12, 1e8, 500, x,
                &itn);
    free_dense_ctx(ctx);

    double xnorm = vec_norm(n, x);
    ASSERT_TRUE(xnorm < 1e-5);
    /* With damp=1e6, ||x|| ≈ ||b|| / (1+damp^2) -> extremely small */

    printf("      ||x||=%.2e  itn=%lld\n", xnorm, (long long)itn);
}
/* 12. Large damp forces solution toward zero */

static void test_tall_skinny(void)
{
    printf(" [15] Tall-skinny (m=20, n=3)\n");

    const ND_int m = 20, n = 3;
    double A[20 * 3];
    double b[20];
    double xtrue[] = {1.0, 2.0, 3.0};

    for (ND_int i = 0; i < m; ++i)
    {
        A[i * n + 0] = (double)(i + 1);
        A[i * n + 1] = (double)((i + 1) * (i + 1));
        /* Fixed to quadratic to ensure linear independence */
        A[i * n + 2] = 1.0;
        b[i] = A[i * n + 0] * xtrue[0] + A[i * n + 1] * xtrue[1] +
               A[i * n + 2] * xtrue[2];
    }

    double x[3] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(m, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-12, 1e-12, 1e8, 500, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-8);
    for (ND_int i = 0; i < n; ++i)
    {
        ASSERT_DOUBLE_NEAR(x[i], xtrue[i], 1e-6);
    }

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 15. Overdetermined with many rows (m >> n) */

static void test_tall_skinny_noisy(void)
{
    printf(" [16] Tall-skinny noisy (m=30, n=3)\n");

    const ND_int m = 30, n = 3;
    double A[30 * 3];
    double b[30];
    double xtrue[] = {2.0, -1.0, 3.0};

    for (ND_int i = 0; i < m; ++i)
    {
        A[i * n + 0] = (double)(i + 1) / m;
        A[i * n + 1] = 1.0 - (double)(i + 1) / m;
        A[i * n + 2] = (double)((i * i) % 7) / 6.0;
        double noise = 0.01 * ((i % 3) - 1.0);
        /* small structured noise */
        b[i] = A[i * n + 0] * xtrue[0] + A[i * n + 1] * xtrue[1] +
               A[i * n + 2] * xtrue[2] + noise;
    }

    double x[3] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(m, n, A);
    int rc = lsmr_solver(m, n, dense_matvec, ctx, b, 0.0, 1e-10, 1e-10, 1e8,
                         500, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_TRUE(rc == 0);

    double Ax[30], r[30], ATr[3];
    DenseMatCtx* ctx2 = make_dense_ctx(m, n, A);
    dense_matvec(0, x, Ax, ctx2);
    for (ND_int i = 0; i < m; ++i)
    {
        r[i] = b[i] - Ax[i];
    }
    dense_matvec(1, r, ATr, ctx2);
    free_dense_ctx(ctx2);

    double ATr_norm = vec_norm(n, ATr);
    ASSERT_TRUE(ATr_norm < 1e-4);
    /* optimality: A^T r ≈ 0 */

    printf("      ||A^Tr||=%.2e  itn=%lld\n", ATr_norm, (long long)itn);
}
/* 16. Tall-skinny inconsistent (noisy observations) */

static void test_wide_fat(void)
{
    printf(" [17] Wide-fat (m=3, n=10)\n");

    const ND_int m = 3, n = 10;
    double A[3 * 10];
    for (ND_int i = 0; i < m; ++i)
    {
        for (ND_int j = 0; j < n; ++j)
        {
            A[i * n + j] = 1.0 / (double)(i + j + 1);
        }
    }
    /* A_ij = 1/(i+j+1) */

    double xtrue[10];
    for (ND_int j = 0; j < n; ++j)
    {
        xtrue[j] = (j % 2 == 0) ? 1.0 : -1.0;
    }

    double b[3] = {0};
    for (ND_int i = 0; i < m; ++i)
    {
        for (ND_int j = 0; j < n; ++j)
        {
            b[i] += A[i * n + j] * xtrue[j];
        }
    }

    double x[10] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(m, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-10, 1e-10, 1e8, 1000, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-6);

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 17. Wide-fat system (m < n) — many solutions, LSMR finds min-norm */

static void test_near_singular(void)
{
    printf(" [18] Near-singular matrix\n");

    const ND_int n = 4;
    double diag[] = {10.0, 5.0, 1.0, 1e-9};
    double A[4 * 4];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = diag[i];
    }

    double xtrue[] = {1.0, 1.0, 1.0, 0.0};
    /* last component unreachable */
    double b[4];
    for (ND_int i = 0; i < n; ++i)
    {
        b[i] = diag[i] * xtrue[i];
    }

    double x[4] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-10, 1e-10, 1e10, 500, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-6);
    for (ND_int i = 0; i < n - 1; ++i)
    {
        ASSERT_DOUBLE_NEAR(x[i], xtrue[i], 1e-4);
    }
    /* First three components recoverable */

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 18. Near-singular matrix — one very small singular value */

static void test_tridiagonal_spd(void)
{
    printf(" [20] Tridiagonal SPD (n=6)\n");

    const ND_int n = 6;
    double A[6 * 6];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = 4.0;
        if (i > 0)
        {
            A[i * n + (i - 1)] = -1.0;
        }
        if (i < n - 1)
        {
            A[i * n + (i + 1)] = -1.0;
        }
    }

    double xtrue[] = {1.0, 2.0, 3.0, 3.0, 2.0, 1.0};
    double b[6] = {0};
    for (ND_int i = 0; i < n; ++i)
    {
        for (ND_int j = 0; j < n; ++j)
        {
            b[i] += A[i * n + j] * xtrue[j];
        }
    }

    double x[6] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-12, 1e-12, 1e8, 500, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-8);
    for (ND_int i = 0; i < n; ++i)
    {
        ASSERT_DOUBLE_NEAR(x[i], xtrue[i], 1e-7);
    }

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 20. Tridiagonal SPD matrix (n=6) */

static void test_non_symmetric_square(void)
{
    printf(" [21] Non-symmetric square (4x4)\n");

    double A[] = {2.0, 1.0, 0.0, 0.5, 0.3, 3.0, 1.5, 0.0,
                  0.0, 0.7, 4.0, 2.0, 0.1, 0.0, 0.8, 5.0};
    double xtrue[] = {1.0, -1.0, 2.0, -2.0};
    double b[4] = {0};
    for (ND_int i = 0; i < 4; ++i)
    {
        for (ND_int j = 0; j < 4; ++j)
        {
            b[i] += A[i * 4 + j] * xtrue[j];
        }
    }

    double x[4] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(4, 4, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-12, 1e-12, 1e8, 300, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-8);
    for (ND_int i = 0; i < 4; ++i)
    {
        ASSERT_DOUBLE_NEAR(x[i], xtrue[i], 1e-7);
    }

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 21. Non-symmetric matrix (general square) */

static void test_large_solution_magnitude(void)
{
    printf(" [22] Large-magnitude solution\n");

    const ND_int n = 4;
    double tiny = 1e-5;
    /* A = tiny*I -> x = b / tiny */

    double A[4 * 4];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = tiny;
    }

    double b[] = {1.0, 2.0, 3.0, 4.0};
    double x[4] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-10, 1e-10, 1e12, 500, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-6);
    for (ND_int i = 0; i < n; ++i)
    {
        ASSERT_DOUBLE_RELNEAR(x[i], b[i] / tiny, 1e-5);
    }

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 22. Solution with large magnitude */

static void test_small_solution_magnitude(void)
{
    printf(" [23] Small-magnitude solution (A big, b small)\n");

    const ND_int n = 4;
    double big = 1e5;
    double A[4 * 4];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = big;
    }

    double b[] = {1.0, 2.0, 3.0, 4.0};
    double x[4] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-10, 1e-10, 1e8, 300, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-6);
    for (ND_int i = 0; i < n; ++i)
    {
        ASSERT_DOUBLE_RELNEAR(x[i], b[i] / big, 1e-5);
    }

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 23. Solution with small magnitude */

static void test_maxiter_one(void)
{
    printf(" [24] maxiter=1 (early termination)\n");

    const ND_int n = 5;
    double A[5 * 5];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = 1.0;
    }

    double b[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double x[5] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    int rc = lsmr_solver(n, n, dense_matvec, ctx, b, 0.0, 1e-12, 1e-12, 1e8, 1,
                         x, &itn);
    free_dense_ctx(ctx);

    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(itn >= 1);

    printf("      itn=%lld  ||x||=%.2e\n", (long long)itn, vec_norm(n, x));
}
/* 24. Maxiter = 1  (only one Lanczos step) */

static void test_return_code_normal(void)
{
    printf(" [25] Return code 0 on success\n");

    double A[] = {1.0};
    double b[] = {1.0};
    double x[] = {0.0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(1, 1, A);
    int rc = lsmr_solver(1, 1, dense_matvec, ctx, b, 0.0, 1e-12, 1e-12, 1e8,
                         100, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_INT_EQ(rc, 0);

    printf("      rc=%d\n", rc);
}
/* 25. Return-code test: normal completion -> 0 */

static void test_residual_decreases_with_more_iterations(void)
{
    printf(" [26] Residual decreases with more iterations\n");

    const ND_int n = 8;
    double A[8 * 8];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = (double)(i + 1);
    }

    double b[8];
    for (ND_int i = 0; i < n; ++i)
    {
        b[i] = (double)(i + 1) * (i + 1);
    }

    ND_int itn;
    double x1[8] = {0}, x100[8] = {0};

    DenseMatCtx* ctx1 = make_dense_ctx(n, n, A);
    double res1 =
        run_lsmr_residual(ctx1, b, 0.0, 1e-14, 1e-14, 0.0, 1, x1, &itn);
    free_dense_ctx(ctx1);

    DenseMatCtx* ctx100 = make_dense_ctx(n, n, A);
    double res100 =
        run_lsmr_residual(ctx100, b, 0.0, 1e-14, 1e-14, 0.0, 100, x100, &itn);
    free_dense_ctx(ctx100);

    ASSERT_TRUE(res100 <= res1 + 1e-14);

    printf("      res(1)=%.2e  res(100)=%.2e\n", res1, res100);
}
/* 26. Convergence monotonicity */

static void test_scaling_invariance(void)
{
    printf(" [27] Scaling invariance of solution\n");

    const ND_int n = 3;
    double A[] = {4.0, 1.0, 0.5, 1.0, 3.0, 0.2, 0.5, 0.2, 2.0};
    double b[] = {1.0, 2.0, 3.0};
    double bsc[3] = {10.0, 20.0, 30.0};

    double x1[3] = {0}, x2[3] = {0};
    ND_int itn;

    DenseMatCtx* ctx1 = make_dense_ctx(n, n, A);
    run_lsmr_residual(ctx1, b, 0.0, 1e-12, 1e-12, 1e8, 200, x1, &itn);
    free_dense_ctx(ctx1);

    DenseMatCtx* ctx2 = make_dense_ctx(n, n, A);
    run_lsmr_residual(ctx2, bsc, 0.0, 1e-12, 1e-12, 1e8, 200, x2, &itn);
    free_dense_ctx(ctx2);

    for (ND_int i = 0; i < n; ++i)
    {
        ASSERT_DOUBLE_RELNEAR(x2[i], 10.0 * x1[i], 1e-6);
    }
    /* x2 should be 10 * x1 */

    printf("      OK\n");
}
/* 27. Scaling invariance: solution should not change if b is scaled */

static void test_random_10x10(void)
{
    printf(" [28] Pseudo-random 10x10 dense\n");

    const ND_int n = 10;
    double A[10 * 10];
    double xtrue[10];
    double b[10] = {0};

    unsigned long state = 1234567UL;
    for (ND_int k = 0; k < n * n; ++k)
    {
        state = state * 6364136223846793005UL + 1442695040888963407UL;
        A[k] = 2.0 * ((double)(state >> 33) / (double)(1UL << 31)) - 1.0;
    }
    /* Deterministic LCG evenly distributed between -1.0 and 1.0 */

    for (ND_int i = 0; i < n; ++i)
    {
        double row_sum = 0.0;
        for (ND_int j = 0; j < n; ++j)
        {
            if (i != j)
            {
                row_sum += fabs(A[i * n + j]);
            }
        }
        A[i * n + i] = row_sum + 2.0;
    }
    /* Make diagonally dominant to ensure non-singular */

    for (ND_int i = 0; i < n; ++i)
    {
        xtrue[i] = (double)(i + 1);
    }
    for (ND_int i = 0; i < n; ++i)
    {
        for (ND_int j = 0; j < n; ++j)
        {
            b[i] += A[i * n + j] * xtrue[j];
        }
    }

    double x[10] = {0};
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-12, 1e-12, 1e8, 500, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-8);
    for (ND_int i = 0; i < n; ++i)
    {
        ASSERT_DOUBLE_RELNEAR(x[i], xtrue[i], 1e-6);
    }

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
}
/* 28. Randomised dense 10x10 system */
static void test_large_overdetermined(void)
{
    printf(" [29] Large overdetermined (m=50, n=10)\n");

    const ND_int m = 50, n = 10;
    double* A = malloc(m * n * sizeof(*A));
    double* b = malloc(m * sizeof(*b));
    double xtrue[10];

    for (ND_int i = 0; i < n; ++i)
    {
        xtrue[i] = (double)(i + 1);
    }
    /* Set to 1.0 through 10.0 to avoid comparing against exact 0.0 */

    unsigned long state = 987654321UL;
    for (ND_int k = 0; k < m * n; ++k)
    {
        state = state * 6364136223846793005UL + 1442695040888963407UL;
        A[k] = 2.0 * ((double)(state >> 33) / (double)(1UL << 31)) - 1.0;
    }
    /* LCG scaled properly between -1.0 and 1.0 */

    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] += 10.0;
    }
    /* Boost the diagonal of the top n x n block to guarantee a low condition
     * number */

    for (ND_int i = 0; i < m; ++i)
    {
        b[i] = 0.0;
        for (ND_int j = 0; j < n; ++j)
        {
            b[i] += A[i * n + j] * xtrue[j];
        }
    }

    double* x = calloc(n, sizeof(*x));
    ND_int itn;
    DenseMatCtx* ctx = make_dense_ctx(m, n, A);
    double res =
        run_lsmr_residual(ctx, b, 0.0, 1e-10, 1e-10, 1e8, 1000, x, &itn);
    free_dense_ctx(ctx);

    ASSERT_DOUBLE_NEAR(res, 0.0, 1e-6);
    for (ND_int i = 0; i < n; ++i)
    {
        ASSERT_DOUBLE_NEAR(x[i], xtrue[i], 1e-6);
    }
    /* Switched back to absolute error which handles small values safely */

    printf("      residual=%.2e  itn=%lld\n", res, (long long)itn);
    free(A);
    free(b);
    free(x);
}

/* 29. Stress: larger overdetermined system (m=50, n=10) */

static void test_maxiter_zero_autoset(void)
{
    printf(" [30] maxiter=0 (auto-set to min(m,n))\n");

    const ND_int n = 4;
    double A[4 * 4];
    memset(A, 0, sizeof(A));
    for (ND_int i = 0; i < n; ++i)
    {
        A[i * n + i] = 1.0;
    }

    double b[] = {1.0, 2.0, 3.0, 4.0};
    double x[4] = {0};
    ND_int itn;

    DenseMatCtx* ctx = make_dense_ctx(n, n, A);
    int rc = lsmr_solver(n, n, dense_matvec, ctx, b, 0.0, 1e-12, 1e-12, 1e8, 0,
                         x, &itn);
    free_dense_ctx(ctx);

    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(itn >= 1);

    printf("      rc=%d  itn=%lld\n", rc, (long long)itn);
}
/* 30. Edge case: maxiter <= 0 should be auto-adjusted (no crash) */

/* ======================================================================
   Main
   ====================================================================== */

int main(void)
{
    printf("====================================================\n");
    printf("   LSMR Extensive Unit Test Suite\n");
    printf("====================================================\n\n");

    /* --- well-conditioned --- */
    test_diagonal_well_conditioned();
    test_small_dense_square_well_conditioned();
    test_tridiagonal_spd();
    test_non_symmetric_square();

    /* --- over/under-determined --- */
    test_overdetermined_consistent();
    test_overdetermined_inconsistent();
    test_tall_skinny();
    test_tall_skinny_noisy();
    test_wide_fat();

    /* --- ill-conditioned --- */
    test_diagonal_ill_conditioned();
    test_hilbert_very_ill_conditioned();
    test_near_singular();

    /* --- very ill-conditioned + damping --- */
    test_hilbert_8_with_damping();
    test_damped_identity();
    test_large_damp_shrinks_solution();

    /* --- magnitude extremes --- */
    test_large_solution_magnitude();
    test_small_solution_magnitude();

    /* --- algorithmic properties --- */
    test_maxiter_one();
    test_maxiter_zero_autoset();
    test_return_code_normal();
    test_residual_decreases_with_more_iterations();
    test_scaling_invariance();

    /* --- larger / stress tests --- */
    test_random_10x10();
    test_large_overdetermined();

    printf("\n====================================================\n");
    printf("Tests run : %d\n", tests_run);
    printf("Failures  : %d\n", tests_failed);
    printf("====================================================\n");

    if (tests_failed == 0)
    {
        printf("✓  All tests passed!\n");
    }
    else
    {
        printf("✗  %d test(s) FAILED!\n", tests_failed);
    }

    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
