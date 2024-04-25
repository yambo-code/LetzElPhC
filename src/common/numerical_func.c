/*
This file contains numerical functions used in the code.
These functions are called with high frequency, so need a good optimization
*/
#include "numerical_func.h"

/*
Most of the modern compilers will auto vectorize most parts of the code if
compiled with -O3 -march=native
Additionally we can also use omp_simd to provide hint to the compilers
*/

/* Static function declarations */
static ELPH_float factorial(ND_int n);

/* Function bodies */

/* Function to compute legendre transforms */
ELPH_float legendre(int l_val, int m_val, ELPH_float x_in)
{
    /* This does not contain Condon–Shortley phase
    Eq: 9-12 in
    "Miguel A. Blanco, M. Flórez, M. Bermejo
    Evaluation of the rotation matrices in the basis of real spherical harmonics
    J. Mol. Struct. THEOCHEM 419, 19-27 (1997)"
    */
    /* Order of below if-else conditions is important */
    if (l_val == 0 && m_val == 0)
    {
        return 1.0;
    }

    else if (l_val < 0)
    {
        return legendre(-l_val - 1, m_val, x_in);
    }

    else if (abs(m_val) > l_val)
    {
        return 0.0;
    }

    else if (m_val < 0)
    {
        return legendre(l_val, -m_val, x_in) * pow(-1, -m_val) * factorial(l_val + m_val) / factorial(l_val - m_val);
    }
    else
    {
        ELPH_float legendre_val = 1.0;
        ELPH_float sqrt_xin = sqrt(1 - (x_in * x_in));

        /* P^l_l = (2l-1)*sqrt(1-x^2)*P^{l-1}_{l-1}(x)  */
        for (int im = 1; im <= m_val; ++im)
        {
            legendre_val *= (2 * im - 1) * sqrt_xin;
        }

        if (l_val == m_val)
        {
            return legendre_val;
        }

        ELPH_float legendre_val_mm = legendre_val;
        /* P^l_l+1 = (2l+1)*x*P^l_l*/
        legendre_val *= (2 * m_val + 1) * x_in;

        if (l_val == (m_val + 1))
        {
            return legendre_val;
        }

        /* For all other l,m combination */
        for (int im = m_val + 2; im <= l_val; ++im)
        {
            ELPH_float legendre_temp = legendre_val;
            legendre_val = (2 * im - 1) * x_in * legendre_val - (im + m_val - 1) * legendre_val_mm;
            legendre_val = legendre_val / (im - m_val);
            legendre_val_mm = legendre_temp;
        }
        return legendre_val;
    }
}

/* Compute real spherical Harmonics for vector direction*/
ELPH_float Ylm(int l_val, int m_val, ELPH_float* vec)
{
    /*
    Miguel A. Blanco, M. Flórez, M. Bermejo
    Evaluation of the rotation matrices in the basis of real spherical harmonics
    J. Mol. Struct. THEOCHEM 419, 19-27 (1997)
    */
    if (l_val < 0)
    {
        return 0.0; // error
    }

    ELPH_float cost, phi;
    ELPH_float norm = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    if (norm < ELPH_EPS)
    {   
        if (l_val != 0)
        {
            return 0; 
        }
        else 
        {
            return 1.0/sqrt(4*ELPH_PI);
        }
    }

    cost = vec[2] / norm; // z/r
    phi = atan2(vec[1], vec[0]);

    ELPH_float sh = 0.0;

    /* phi part*/
    if (m_val > 0)
    {
        sh = sqrt(2) * cos(m_val * phi);
    }
    else if (m_val < 0)
    {
        m_val = -m_val;
        sh = sqrt(2) * sin(m_val * phi);
    }
    else
    {
        sh = 1; // 1 for m =0
    }

    /* Multiply with legendre polynomial */
    sh *= legendre(l_val, m_val, cost);

    /* Multiply with normalization factor*/
    sh *= sqrt(factorial(l_val - m_val) * (2 * l_val + 1) / (4 * ELPH_PI * factorial(l_val + m_val)));

    return sh;
}

/* Function for simpson integration*/
ELPH_float simpson(const ELPH_float* restrict func_vals, const ELPH_float* restrict dx,
                   ND_int npts)
{
    /*
    Compute the integral using simpson 1/3 rules i.e /int f(x) dx
    */
    if (npts < 3)
    {
        return 0; // Need atleast 3 points . Return an error instead
    }

    ELPH_float sum = func_vals[0] * dx[0];

/* Odd terms*/
#pragma omp simd reduction(+ : sum)
    for (ND_int i = 1; i < npts - 1; i += 2)
    {
        sum += 4.0 * func_vals[i] * dx[i];
    }

/* even terms*/
#pragma omp simd reduction(+ : sum)
    for (ND_int i = 2; i < npts - 1; i += 2)
    {
        sum += 2.0 * func_vals[i] * dx[i];
    }

    /* add/remove the last point if npts is odd/even */
    if (npts % 2 != 0)
    {
        sum += func_vals[npts - 1] * dx[npts - 1];
    }
    else
    {
        sum -= func_vals[npts - 2] * dx[npts - 2];
    }

    return sum / 3.0;
}

/* Function to compute cosine of angle between two vectors*/
ELPH_float cos_angle_bw_Vec(const ELPH_float* vec1, const ELPH_float* vec2)
{
    /*
    This function returns Cos(w) where w is angle between vec1 and vec2
    */
    ELPH_float norm1 = vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2];
    ELPH_float norm2 = vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2];
    ELPH_float norm = sqrt(norm1 * norm2);
    if (norm < ELPH_EPS)
    {
        return 0;
    }
    ELPH_float dot_12 = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
    return dot_12 / norm;
}

/** Static functions **/

/* Compute factorial */
static ELPH_float factorial(ND_int n)
{
    if (n < 0)
    {
        return 0.0; // error
    }
    if (n == 0 || n == 1)
    {
        return 1.0;
    }
    ELPH_float factorial = 1.0;
    for (ND_int i = 2; i <= n; ++i)
    {
        factorial *= i;
    }
    return factorial;
}

/* Better to inline them ? */
/**This are hardcoded to avoid over heads */
ELPH_float dotVec3(const ELPH_float* vec1, const ELPH_float* vec2)
{
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

void MatVec3f(const ELPH_float* restrict Mat, const ELPH_float* restrict vec,
              const bool trans, ELPH_float* restrict out)
{
    if (trans)
    {
        out[0] = Mat[0] * vec[0] + Mat[3] * vec[1] + Mat[6] * vec[2];
        out[1] = Mat[1] * vec[0] + Mat[4] * vec[1] + Mat[7] * vec[2];
        out[2] = Mat[2] * vec[0] + Mat[5] * vec[1] + Mat[8] * vec[2];
    }
    else
    {
        out[0] = Mat[0] * vec[0] + Mat[1] * vec[1] + Mat[2] * vec[2];
        out[1] = Mat[3] * vec[0] + Mat[4] * vec[1] + Mat[5] * vec[2];
        out[2] = Mat[6] * vec[0] + Mat[7] * vec[1] + Mat[8] * vec[2];
    }
}

ELPH_cmplx Cmplxdot(const ELPH_cmplx* vec1, const ELPH_cmplx* vec2,
                    const ND_int n)
{
    ELPH_cmplx sum = 0.0;

    for (ND_int i = 0; i < n; ++i)
    {
        sum += conj(vec1[i]) * vec2[i];
    }

    return sum;
}

void normalize_Cmplx_vec(ELPH_cmplx* vec, const ND_int n)
{
    ELPH_float norm = sqrt(cabs(Cmplxdot(vec, vec, n)));

    if (norm < ELPH_EPS)
    {
        return;
    }

    for (ND_int i = 0; i < n; ++i)
    {
        vec[i] = vec[i] / norm;
    }
}

ELPH_float det3x3(const ELPH_float* mat)
{
    ELPH_float det = 0;
    det += mat[0] * (mat[4] * mat[8] - mat[5] * mat[7]);
    det -= mat[1] * (mat[3] * mat[8] - mat[5] * mat[6]);
    det += mat[2] * (mat[7] * mat[3] - mat[4] * mat[6]);
    return det;
}

void reciprocal_vecs(const ELPH_float* restrict lat_vec,
                     ELPH_float* restrict blat)
{
    /*
    a[:,i]  are latvecs
    b[:,i]  are blat
    result is multiplied with 2*pi
    */
    ELPH_float det = det3x3(lat_vec);
    if (det < ELPH_EPS)
    {
        error_msg("Inverting singular matrix");
    }
    det = 2 * ELPH_PI / det;
    blat[0] = (lat_vec[4] * lat_vec[8] - lat_vec[7] * lat_vec[5]) * det;
    blat[1] = (lat_vec[5] * lat_vec[6] - lat_vec[3] * lat_vec[8]) * det;
    blat[2] = (lat_vec[3] * lat_vec[7] - lat_vec[6] * lat_vec[4]) * det;
    blat[3] = (lat_vec[2] * lat_vec[7] - lat_vec[1] * lat_vec[8]) * det;
    blat[4] = (lat_vec[0] * lat_vec[8] - lat_vec[2] * lat_vec[6]) * det;
    blat[5] = (lat_vec[1] * lat_vec[6] - lat_vec[0] * lat_vec[7]) * det;
    blat[6] = (lat_vec[1] * lat_vec[5] - lat_vec[2] * lat_vec[4]) * det;
    blat[7] = (lat_vec[2] * lat_vec[3] - lat_vec[0] * lat_vec[5]) * det;
    blat[8] = (lat_vec[0] * lat_vec[4] - lat_vec[1] * lat_vec[3]) * det;
}

void aXpY(const ND_int n, const ELPH_cmplx a, const ELPH_cmplx* restrict X,
          ELPH_cmplx* restrict Y)
{
    /* computes y= aX + y */
    // ELPH_OMP_PAR_FOR_SIMD
    for (ND_int i = 0; i < n; ++i)
    {
        Y[i] += a * X[i];
    }
}

void transpose3x3f(const ELPH_float* restrict inmat,
                   ELPH_float* restrict outmat)
{
    /* Transpose 3x3 matrix */
    outmat[0] = inmat[0];
    outmat[1] = inmat[3];
    outmat[2] = inmat[6];

    outmat[3] = inmat[1];
    outmat[4] = inmat[4];
    outmat[5] = inmat[7];

    outmat[6] = inmat[2];
    outmat[7] = inmat[5];
    outmat[8] = inmat[8];
}

ND_int find_maxint(ND_int* in_arr, ND_int nelements)
{
    ND_int max = in_arr[0];
    for (ND_int imax = 0; imax < nelements; ++imax)
    {
        if (in_arr[imax] > max)
        {
            max = in_arr[imax];
        }
    }
    // printf("Max : %llu \n ", max);
    return max;
}

ELPH_float find_maxfloat(ELPH_float* in_arr, ND_int nelements)
{
    ELPH_float max = in_arr[0];
    for (ND_int imax = 0; imax < nelements; ++imax)
    {
        if (in_arr[imax] > max)
        {
            max = in_arr[imax];
        }
    }
    // printf("Max : %llu \n ", max);
    return max;
}

void Gemm3x3f(const ELPH_float* restrict A, const char transA,
              const ELPH_float* restrict B, const char transB,
              ELPH_float* restrict C)
{
    if (transA == 'N' && transB == 'N')
    {
        C[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
        C[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
        C[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
        C[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
        C[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
        C[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
        C[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
        C[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
        C[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
    }
    else if (transA == 'T' && transB == 'N')
    {
        C[0] = A[0] * B[0] + A[3] * B[3] + A[6] * B[6];
        C[1] = A[0] * B[1] + A[3] * B[4] + A[6] * B[7];
        C[2] = A[0] * B[2] + A[3] * B[5] + A[6] * B[8];
        C[3] = A[1] * B[0] + A[4] * B[3] + A[7] * B[6];
        C[4] = A[1] * B[1] + A[4] * B[4] + A[7] * B[7];
        C[5] = A[1] * B[2] + A[4] * B[5] + A[7] * B[8];
        C[6] = A[2] * B[0] + A[5] * B[3] + A[8] * B[6];
        C[7] = A[2] * B[1] + A[5] * B[4] + A[8] * B[7];
        C[8] = A[2] * B[2] + A[5] * B[5] + A[8] * B[8];
    }
    else if (transA == 'N' && transB == 'T')
    {
        C[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
        C[1] = A[0] * B[3] + A[1] * B[4] + A[2] * B[5];
        C[2] = A[0] * B[6] + A[1] * B[7] + A[2] * B[8];
        C[3] = A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
        C[4] = A[3] * B[3] + A[4] * B[4] + A[5] * B[5];
        C[5] = A[3] * B[6] + A[4] * B[7] + A[5] * B[8];
        C[6] = A[6] * B[0] + A[7] * B[1] + A[8] * B[2];
        C[7] = A[6] * B[3] + A[7] * B[4] + A[8] * B[5];
        C[8] = A[6] * B[6] + A[7] * B[7] + A[8] * B[8];
    }
    else if (transA == 'T' && transB == 'T')
    {
        C[0] = A[0] * B[0] + A[3] * B[1] + A[6] * B[2];
        C[1] = A[0] * B[3] + A[3] * B[4] + A[6] * B[5];
        C[2] = A[0] * B[6] + A[3] * B[7] + A[6] * B[8];
        C[3] = A[1] * B[0] + A[4] * B[1] + A[7] * B[2];
        C[4] = A[1] * B[3] + A[4] * B[4] + A[7] * B[5];
        C[5] = A[1] * B[6] + A[4] * B[7] + A[7] * B[8];
        C[6] = A[2] * B[0] + A[5] * B[1] + A[8] * B[2];
        C[7] = A[2] * B[3] + A[5] * B[4] + A[8] * B[5];
        C[8] = A[2] * B[6] + A[5] * B[7] + A[8] * B[8];
    }
    else
    {
        error_msg("Wrong Trans type");
    }
    return;
}

void matmul_Cmpl2x2(ELPH_cmplx* restrict mat1, ELPH_cmplx* restrict mat2,
                    ELPH_cmplx* restrict out)
{
    out[0] = mat1[0] * mat2[0] + mat1[1] * mat2[2];
    out[1] = mat1[0] * mat2[1] + mat1[1] * mat2[3];
    out[2] = mat1[2] * mat2[0] + mat1[3] * mat2[2];
    out[3] = mat1[2] * mat2[1] + mat1[3] * mat2[3];
}

/* functions related to fft */

ND_int get_miller_idx(ND_int idx_in, ND_int FFT_dimension)
{
    // returns FFT Indices to [-N/2, N/2) if FFT_dimension is even
    // [-(n-1)/2,(n-1)/2 )] of FFT_dimension is odd
    ND_int mid_pnt = (FFT_dimension - 1) / 2 + 1;
    if (idx_in < mid_pnt)
    {
        return idx_in;
    }
    else
    {
        return (idx_in - mid_pnt) - FFT_dimension / 2;
    }
}

/* This function converts miller indices to FFT indices*/
int get_fft_idx(ELPH_float idx_in, int FFT_dimension)
{
    /*
    We need this functions because FFT libraries assume that the
    indices run from [0,N) but where as miller indices are from
    [-N/2,N/1)
    */
    // returns FFT Indices to [0, N-1]
    int temp_idx = rint(idx_in);
    if (temp_idx >= 0)
    {
        return temp_idx;
    }
    else
    {
        return FFT_dimension + temp_idx;
    }
}

ND_int find_kidx_in_list(ND_int nkpts, const ELPH_float* kpts_list,
                         const ELPH_float* kpt)
{
    /*
    kpt and kpts_list must be in crystal coordinates
    return -1 if not found
    */
    ND_int kidx = -1;
    for (ND_int ik = 0; ik < nkpts; ++ik)
    {
        const ELPH_float* ik_vec_tmp = kpts_list + 3 * ik;
        ELPH_float sum = 0;
        for (int i = 0; i < 3; ++i)
        {
            ELPH_float diff_tmp = ik_vec_tmp[i] - kpt[i];
            diff_tmp = diff_tmp - rint(diff_tmp);
            sum += diff_tmp * diff_tmp;
        }
        sum = sqrt(sum);
        if (sum < ELPH_EPS)
        {
            kidx = ik;
            break;
        }
    }
    return kidx;
}

/* This Function finds the K+Q indices in kpoint grid */
void get_KplusQ_idxs(const ND_int Nbz, const ELPH_float* kpoints,
                     const ELPH_float* Q_pt, int* KplusQidxs)
{
    /* returns index of the k+q point in the kpoints*/
    // kpoitsn and  kpoints must be in crystal coordinates

    for (ND_int i = 0; i < Nbz; ++i)
    {
        const ELPH_float* ktemp = kpoints + 3 * i;

        ELPH_float KplusQ[3];
        KplusQ[0] = ktemp[0] + Q_pt[0];
        KplusQ[1] = ktemp[1] + Q_pt[1];
        KplusQ[2] = ktemp[2] + Q_pt[2];

        KplusQidxs[i] = find_kidx_in_list(Nbz, kpoints, KplusQ);

        if (KplusQidxs[i] < 0)
        {
            error_msg(
                "K+Q point cannot be found in kgrid. "
                "The k-grid and q-grid are not commensurate! \n");
        }
    }
}

// some helper functions
void swap_ints(int* restrict a, int* restrict b)
{
    const int c = *b;
    *b = *a;
    *a = c;
}
