// adapated from scipy lsmr
//
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "elphC.h"

static inline void LSMR_vec_set(ND_int n, double* restrict x, double val)
{
    // set a value to all elements in a vector
    for (ND_int i = 0; i < n; ++i)
    {
        x[i] = val;
    }
}

static inline void LSMR_vec_scale(ND_int n, double* restrict x, double alpha)
{
    // vector scaling (x = alpha * x)
    for (ND_int i = 0; i < n; ++i)
    {
        x[i] *= alpha;
    }
}

static inline void LSMR_vec_add_scaled(ND_int n, double* restrict y,
                                       const double* x, double alpha)
{
    // axpy (y = alpha * x + y)
    for (ND_int i = 0; i < n; ++i)
    {
        y[i] += alpha * x[i];
    }
}

static inline double LSMR_vec_norm(ND_int n, const double* x)
{
    // euclidian norm
    double sum_sq = 0.0;
    for (ND_int i = 0; i < n; ++i)
    {
        sum_sq += x[i] * x[i];
    }
    return sqrt(sum_sq);
}

static inline void LSMR_update_x(ND_int n, double* restrict x,
                                 const double* restrict hbar, double zeta,
                                 double rho, double rhobar)
{
    // x = x + (zeta * hbar) / (rho * rhobar)
    double scalar = (zeta) / (rho * rhobar);
    for (ND_int i = 0; i < n; ++i)
    {
        x[i] += scalar * hbar[i];
    }
}

static inline void LSMR_sym_ortho(double a, double b, double* restrict c,
                                  double* restrict s, double* restrict r)
{
    /*
    Stable implementation of Givens rotation.

    Notes
    -----
    The routine 'SymOrtho' was added for numerical stability. This is
    recommended by S.-C. Choi in [1]_.  It removes the unpleasant potential of
    ``1/eps`` in some important places (see, for example text following
    "Compute the next plane rotation Qk" in minres.py).

    References
    ----------
    .. [1] S.-C. Choi, "Iterative Methods for Singular Linear Equations
           and Least-Squares Problems", Dissertation,
           http://www.stanford.edu/group/SOL/dissertations/sou-cheng-choi-thesis.pdf
    */
    if (b == 0.0)
    {
        *c = (a > 0) - (a < 0);
        *s = 0.0;
        *r = fabs(a);
    }
    else if (a == 0.0)
    {
        *c = 0.0;
        *s = (b > 0) - (b < 0);
        *r = fabs(b);
    }
    else if (fabs(b) > fabs(a))
    {
        double tau = a / b;
        *s = ((b > 0) - (b < 0)) / sqrt(1.0 + tau * tau);
        *c = (*s) * tau;
        *r = b / (*s);
    }
    else
    {
        double tau = b / a;
        *c = ((a > 0) - (a < 0)) / sqrt(1.0 + tau * tau);
        *s = (*c) * tau;
        *r = a / (*c);
    }
}

/* -------------------------------------------------------------------------
   LSMR Implementation
   ------------------------------------------------------------------------- */
int lsmr_solver(ND_int m, ND_int n,
                void (*matvec)(const int, const double* restrict,
                               double* restrict, void*),
                void* userdata, const double* restrict b, double damp,
                double atol, double btol, double conlim, ND_int maxiter,
                double* restrict x, ND_int* itn_out)
{
    /**
    * @brief LSMR iterative solver for least-squares problems.
    *
    * Solves the linear system  A@X = b. If the system is inconsistent,
    * the algorithm computes the least-squares solution that minimizes
    * ||b - A x||_2.
    *
    * LSMR is an iterative Krylov subspace method suitable for large and/or
    * sparse matrices. The matrix A may be rectangular (m x n) and is accessed
    * through user-supplied matrix-vector multiplication callbacks.
    *
    * Optionally, a damping parameter allows solving the regularized
    * least-squares problem:
    *
    * ```
      minimize || [ b ] - [  A   ] x ||_2
      ```
    * ```
               || [ 0 ]   [ damp I ]   ||
      ```
    *
    * @param[in] m
    * Number of rows of matrix A.
    *
    * @param[in] n
    * Number of columns of matrix A.
    *
    * @param[in] A_mul
    * Function pointer for forward matrix-vector multiplication.
    * Computes y = A * x (mode = 0) or A^T * x (mode != 0)
    *
    * @param[in] userdata
    * Optional user data passed to A_mul and AT_mul function.
    *
    * @param[in] b
    * Right-hand side vector of length m.
    *
    * @param[in] damp
    * Damping (regularization) parameter. Set to 0.0 to disable
    * regularization.
    *
    * @param[in] atol
    * Absolute stopping tolerance.
    *
    * @param[in] btol
    * Relative stopping tolerance.
    *
    * @param[in] conlim
    * Condition number limit. Iteration stops if the estimated
    * condition number exceeds this value.
    *
    * @param[in] maxiter
    * Maximum number of iterations (typically min(m, n)).
    *
    * @param[in,out] x
    * On input: initial guess for the solution (may be NULL or zeroed).
    * On output: computed solution vector of length n.
    *
    * @param[out] itn_out
    * Number of iterations performed.
    *
    *
    * @return
    * 0 on successful completion, nonzero on error.
    *
    * @note
    * Convergence depends on the conditioning of A and the choice of
    * tolerances. Smaller tolerances may increase iteration count.
    *
    * @see
    * D. C.-L. Fong and M. A. Saunders,
    * "LSMR: An iterative algorithm for sparse least-squares problems",
    * SIAM J. Sci. Comput., 33(5):2950–2971, 2011.
  */

    // Constants
    const double eps = DBL_EPSILON;
    const double ctol = (conlim > 0.0) ? (1.0 / conlim) : 0.0;

    // Adjust maxiter if not set
    if (maxiter <= 0)
    {
        maxiter = (m < n) ? m : n;
    }

    // Allocate working vectors
    double* u = malloc(m * sizeof(*u));
    double* v = malloc(n * sizeof(*v));
    double* h = malloc(n * sizeof(*h));
    double* hbar = malloc(n * sizeof(*hbar));
    double* Ax = malloc(m * sizeof(*Ax));
    double* Atu = malloc(n * sizeof(*Atu));

    if (!u || !v || !h || !hbar || !Ax)
    {
        free(u);
        free(v);
        free(h);
        free(hbar);
        free(Ax);
        free(Atu);
        return -1;
    }

    // Initialize u = b - A*x
    memcpy(u, b, m * sizeof(*u));

    double normb = LSMR_vec_norm(m, b);

    matvec(0, x, Ax, userdata);
    LSMR_vec_add_scaled(m, u, Ax, -1.0);

    double beta = LSMR_vec_norm(m, u);
    double alpha = 0.0;

    if (beta > 0.0)
    {
        LSMR_vec_scale(m, u, 1.0 / beta);
        matvec(1, u, v, userdata);
        alpha = LSMR_vec_norm(n, v);
    }
    else
    {
        LSMR_vec_set(n, v, 0.0);
    }

    if (alpha > 0.0)
    {
        LSMR_vec_scale(n, v, 1.0 / alpha);
    }

    // Initialize variables for 1st iteration
    double zetabar = alpha * beta;
    double alphabar = alpha;
    double rho = 1.0;
    double rhobar = 1.0;
    double cbar = 1.0;
    double sbar = 0.0;

    memcpy(h, v, n * sizeof(*h));
    LSMR_vec_set(n, hbar, 0.0);

    // Initialize variables for estimation of ||r||
    double betadd = beta;
    double betad = 0.0;
    double rhodold = 1.0;
    double tautildeold = 0.0;
    double thetatilde = 0.0;
    double zeta = 0.0;
    double d = 0.0;

    // Initialize variables for estimation of ||A|| and cond(A)
    double normA2 = alpha * alpha;
    double maxrbar = 0.0;
    double minrbar = 0.99 * DBL_MAX;
    double normA = sqrt(normA2);

    ND_int itn = 0;

    // Main iteration loop
    for (itn = 1; itn <= maxiter; ++itn)
    {
        /*  Perform the next step of the bidiagonalization to obtain the
         next  beta, u, alpha, v.  These satisfy the relations
                 beta*u  =  A@v   -  alpha*u,
                alpha*v  =  A'@u  -  beta*v. */

        LSMR_vec_scale(m, u, -alpha);
        matvec(0, v, Ax, userdata);
        LSMR_vec_add_scaled(m, u, Ax, 1.0);

        beta = LSMR_vec_norm(m, u);

        if (beta > 0.0)
        {
            LSMR_vec_scale(m, u, 1.0 / beta);
            LSMR_vec_scale(n, v, -beta);

            matvec(1, u, Atu, userdata);
            LSMR_vec_add_scaled(n, v, Atu, 1.0);

            alpha = LSMR_vec_norm(n, v);

            if (alpha > 0)
            {
                LSMR_vec_scale(n, v, 1.0 / alpha);
            }
        }

        // At this point, beta = beta_{k+1}, alpha = alpha_{k+1}.
        //
        // Construct rotation Qhat_{k,2k+1}
        double chat, shat, alphahat;
        LSMR_sym_ortho(alphabar, damp, &chat, &shat, &alphahat);

        // Use a plane rotation (Q_i) to turn B_i to R_i

        double rhoold = rho;

        // Construct rotation Q_{k,k+1}
        double c, s;
        LSMR_sym_ortho(alphahat, beta, &c, &s, &rho);

        double thetanew = s * alpha;
        alphabar = c * alpha;

        // Use a plane rotation (Qbar_i) to turn R_i^T to R_i^bar

        double rhobarold = rhobar;
        double zetaold = zeta;
        double thetabar = sbar * rho;
        double rhotemp = cbar * rho;

        LSMR_sym_ortho(cbar * rho, thetanew, &cbar, &sbar, &rhobar);

        zeta = cbar * zetabar;
        zetabar = -sbar * zetabar;

        // Update h, h_hat, x.
        double hbar_scale = (-thetabar * rho) / (rhoold * rhobarold);
        LSMR_vec_scale(n, hbar, hbar_scale);
        LSMR_vec_add_scaled(n, hbar, h, 1.0);
        LSMR_update_x(n, x, hbar, zeta, rho, rhobar);

        double h_scale = (-thetanew) / rho;
        LSMR_vec_scale(n, h, h_scale);
        LSMR_vec_add_scaled(n, h, v, 1.0);

        // Estimate of ||r||
        // Apply rotation Qhat_{k,2k+1}.
        double betaacute = chat * betadd;
        double betacheck = -shat * betadd;
        double betahat = c * betaacute;
        betadd = -s * betaacute;

        // Apply rotation Qtilde_{k-1}.
        // betad = betad_{k-1} here.

        double thetatildeold = thetatilde;

        double ctildeold, stildeold, rhotildeold;
        LSMR_sym_ortho(rhodold, thetabar, &ctildeold, &stildeold, &rhotildeold);

        thetatilde = stildeold * rhobar;
        rhodold = ctildeold * rhobar;
        betad = (betad * -stildeold) + (ctildeold * betahat);

        // betad   = betad_k here.
        // rhodold = rhod_k  here.

        tautildeold = ((tautildeold * -thetatildeold) + zetaold) / rhotildeold;
        double taud = (zeta - thetatilde * tautildeold) / rhodold;
        d = d + (betacheck * betacheck);
        double normr = sqrt(d + pow(betad - taud, 2) + pow(betadd, 2));

        // Estimate ||A||
        normA2 = normA2 + (beta * beta);
        normA = sqrt(normA2);
        normA2 = normA2 + (alpha * alpha);

        // Estimate cond(A)
        if (rhobarold > maxrbar)
        {
            maxrbar = rhobarold;
        }
        if (itn > 1)
        {
            if (rhobarold < minrbar)
            {
                minrbar = rhobarold;
            }
        }

        // Test for convergence
        if (itn % 10 == 0)
        {
            double normar = fabs(zetabar);
            double normx = LSMR_vec_norm(n, x);

            double max_val = (maxrbar > rhotemp) ? maxrbar : rhotemp;
            double min_val = (minrbar < rhotemp) ? minrbar : rhotemp;
            double condA = max_val / min_val;

            double test1 = normr / normb;
            double test2 = normar / (normA * normr + eps);
            double test3 = 1.0 / (condA + eps);
            double t1 = test1 / (1.0 + normA * normx / normb);
            double rtol = btol + atol * normA * normx / normb;

            bool stop = ((1.0 + test3 <= 1.0) || (1.0 + test2 <= 1.0) ||
                         (1.0 + t1 <= 1.0) || (test3 <= ctol) ||
                         (test2 <= atol) || (test1 <= rtol));

            if (stop)
            {
                break;
            }
        }
    }

    *itn_out = itn;

    // Cleanup
    free(u);
    free(v);
    free(h);
    free(hbar);
    free(Ax);
    free(Atu);

    return 0;
}
