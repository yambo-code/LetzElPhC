#include "nd_array_src.h"


/* Linear algebra ND_function (blas, tblis, lapack ...)*/
/****************************************************************************************************/

#if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX) || defined(COMPILE_ND_FLOAT) || defined(COMPILE_ND_DOUBLE)

static enum CBLAS_TRANSPOSE get_gemmn_T(char Trans);

/****************************************************************************************************/
                                /* BLAS and LAPACK based*/


void ND_function(matmul, TYPE_S) (const char TransA, const char TransB, \
        const ND_array(TYPE_S) * nd_arr_A, const ND_array(TYPE_S) * nd_arr_B, \
        ND_array(TYPE_S) * nd_arr_C, const TYPE_L alpha, const TYPE_L beta, \
        const ND_int * A_idx,  const ND_int * B_idx,  const ND_int * C_idx)

{
    /* THis will perform the following matrix multiplication.
    ----- A({A_idx},m,k)@ B({B_idx},k,n) = C({C_idx},m,n)  ------ */

    /* Some checks on data*/
    
    if ( (((nd_arr_A->rank) == NULL )  || ((nd_arr_B->rank) == NULL )) || ((nd_arr_C->rank) == NULL ) ) \
                                                    error_msg("Cannot accept uninitilized array") ; 
    
    if ( (((nd_arr_A->data) == NULL )  || ((nd_arr_B->data) == NULL )) || ((nd_arr_C->data) == NULL ) ) \
                                                    error_msg("Cannot accept NULL array. allocate them beforing passing to matmul");

    /*Find the max BLAS int to check for overflow of indices */
    intmax_t max_len_x;
    ND_int max_len;
    
    for (max_len_x=INTMAX_MAX; (BLAS_INT)max_len_x!=max_len_x; max_len_x/=2);
    
    max_len = max_len_x; 

    if ( ((*(nd_arr_A->rank) <  2)  || (*(nd_arr_B->rank) <  2)) || (*(nd_arr_C->rank) <  2) ) \
                                                    error_msg("Matmul only accepts arrays with atleast rank 2") ;

    if ( (A_idx == NULL) && (*(nd_arr_A->rank)  != 2)) error_msg("Null is passed for A with dim > 2") ;

    if ( (B_idx == NULL) && (*(nd_arr_B->rank)  != 2)) error_msg("Null is passed for B with dim > 2") ;

    if ( (C_idx == NULL) && (*(nd_arr_C->rank)  != 2)) error_msg("Null is passed for C with dim > 2.") ;

    ND_int ldA = (nd_arr_A->dims)[*(nd_arr_A->rank)-1] ;
    ND_int ldB = (nd_arr_B->dims)[*(nd_arr_B->rank)-1];
    ND_int ldC = (nd_arr_C->dims)[*(nd_arr_C->rank)-1];

    ND_int m1, k1, k2, n2, m3,n3 ;

    if (TransA == 'C' || TransA == 'T')
    {
        m1 = (nd_arr_A->dims)[*(nd_arr_A->rank)-1]; /*  m,k--> k,m */
        k1 = (nd_arr_A->dims)[*(nd_arr_A->rank)-2];
    }

    else
    {
        k1 = (nd_arr_A->dims)[*(nd_arr_A->rank)-1]; /*  m,k--> k,m */
        m1 = (nd_arr_A->dims)[*(nd_arr_A->rank)-2];
    }

    if (TransB == 'C' || TransB == 'T')
    {
        k2 = (nd_arr_B->dims)[*(nd_arr_B->rank)-1]; /*  m,k--> k,m */
        n2 = (nd_arr_B->dims)[*(nd_arr_B->rank)-2];
    }

    else
    {
        n2 = (nd_arr_B->dims)[*(nd_arr_B->rank)-1]; /*  m,k--> k,m */
        k2 = (nd_arr_B->dims)[*(nd_arr_B->rank)-2];
    }

    n3 = (nd_arr_C->dims)[*(nd_arr_C->rank)-1]; /*  m,k--> k,m */
    m3 = (nd_arr_C->dims)[*(nd_arr_C->rank)-2];
    
    if ( ((k1 != k2) || (m1 !=m3)) || (n2 !=n3)  )
    {
        fprintf(stdout, "# [ Error !!!] : Incompatible matrix dimensions. Matrix product not possible ("PRI_NdInt", "PRI_NdInt") \
        , ("PRI_NdInt", "PRI_NdInt") -> ("PRI_NdInt", "PRI_NdInt")! \n", m1, k1, k2, n2, m3,n3);
        error_msg("Incompatible matrix dimensions. Matrix product not possible");
    }

    ND_int idx_A = 0, idx_B = 0, idx_C = 0; 

    for (ND_int i = 0; i< *(nd_arr_A->rank) -  2; ++i)
    {
        if ( A_idx[i] >= (nd_arr_A->dims)[i] ) error_msg("Provided indices of A matrix for nd_matmul ND_function are out of bounds");

        idx_A = idx_A + ((nd_arr_A->strides)[i] * A_idx[i]) ; 
    }

    // error_msg("");

    for (ND_int i = 0; i< *(nd_arr_B->rank) -  2; ++i)
    {
        if ( B_idx[i] >= (nd_arr_B->dims)[i] ) error_msg("Provided indices of B matrix for nd_matmul ND_function are out of bounds");

        idx_B = idx_B + ((nd_arr_B->strides)[i] * B_idx[i]) ; 
    }

    for (ND_int i = 0; i< *(nd_arr_C->rank) -  2; ++i){
        if ( C_idx[i] >= (nd_arr_C->dims)[i] ) error_msg("Provided indices of C matrix for nd_matmul ND_function are out of bounds");

        idx_C = idx_C + ((nd_arr_C->strides)[i] * C_idx[i]) ; 
    }
    /* call the gemm*/
    // 
    if ( (((((m1 >= max_len ||  n2 >= max_len) ||  k1 >= max_len) ||  ldA >= max_len) ||  ldB >= max_len) ||  ldC >= max_len)  ) \
                                            error_msg("BLAS indices Overflow. Compile the blas library with higher bit indices!");
    
    // Call blas gemm

    BLAS_CALL(gemm , TYPE_S) (CblasRowMajor, get_gemmn_T(TransA), get_gemmn_T(TransB), (BLAS_INT) m1 , (BLAS_INT) n2, (BLAS_INT) k1 ,\
                BLAS_POINTER(alpha), (nd_arr_A->data) + idx_A, (BLAS_INT) ldA,  (nd_arr_B->data) + idx_B,  (BLAS_INT) ldB,  BLAS_POINTER(beta), (nd_arr_C->data) + idx_C,  (BLAS_INT) ldC);
}


/** Matmul Expert version.*/
void ND_function(matmulX, TYPE_S) (const char TransA, const char TransB, const TYPE_L * arr_A, const TYPE_L * arr_B, TYPE_L * arr_C, \
                const TYPE_L alpha, const TYPE_L beta, const ND_int ldA, const ND_int ldB, const ND_int ldC, \
                const ND_int m, const ND_int n, const ND_int k)

{
    /* THis will perform the following matrix multiplication. Expert version
    ----- Simply Calls Gemm ND_function from Cblas  ------ */

    /* Some checks on data*/
    

    /*Find the max BLAS int to check for overflow of indices */
    intmax_t max_len_x;
    ND_int max_len;
    
    for (max_len_x=INTMAX_MAX; (BLAS_INT)max_len_x!=max_len_x; max_len_x/=2);
    
    max_len = max_len_x; 

    // 
    if ( m >= max_len ||  n >= max_len ||  k >= max_len  ) \
        error_msg("BLAS indices Overflow. Compile the blas library with higher bit indices!");
    
    // Call blas gemm

    BLAS_CALL(gemm , TYPE_S) (CblasRowMajor, get_gemmn_T(TransA), get_gemmn_T(TransB), (BLAS_INT) m , (BLAS_INT) n, (BLAS_INT) k ,\
                BLAS_POINTER(alpha), arr_A, (BLAS_INT) ldA,  arr_B,  (BLAS_INT) ldB,  BLAS_POINTER(beta), arr_C,  (BLAS_INT) ldC);    
}


/* ND_function OF MATMUL to identify transpose/conjuate*/
static enum CBLAS_TRANSPOSE get_gemmn_T(char Trans)
{
    if (Trans == 'N')      return CblasNoTrans ;
    else if (Trans == 'T') return CblasTrans;
    else if (Trans == 'C') return CblasConjTrans;
    else 
    {
        error_msg("Can only take 'C' or 'T' or 'N' for Trans input");
        return CblasNoTrans ;
    }
}
#endif
