#include "nd_array_src.h"


/* Linear algebra ND_function (blas, tblis, lapack ...)*/
/****************************************************************************************************/

#if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX) || defined(COMPILE_ND_FLOAT) || defined(COMPILE_ND_DOUBLE)

static enum CBLAS_TRANSPOSE get_gemmn_T(char Trans);

#if defined(COMPILE_ND_TBLIS)
static void nd_free_tblis(tblis_tensor * A_arr);
static void ND_function(to_tblis, TYPE_S)(const ND_array(TYPE_S) * nd_arr_A, tblis_tensor * A_arr);
static void ND_function(to_tblis_scaled, TYPE_S)(const ND_array(TYPE_S) * nd_arr_A, tblis_tensor * A_arr, const TYPE_L alpha);
static void check_string_dims(char * str_A, char * str_B, ND_int * dims_A,  ND_int * dims_B);
static void check_all_indices(char * str_A, char * str_B);
static void get_einsum_excluded_idxs(const char * str_A, const char* str_B, const char* str_C, char * out_str);
static void check_str_char_unq(char * inp_str, char* err_msg);
#endif
/****************************************************************************************************/
                                /* BLAS and LAPACK based*/


void ND_function(matmul, TYPE_S) (const char TransA, const char TransB, const ND_array(TYPE_S) * nd_arr_A, const ND_array(TYPE_S) * nd_arr_B, ND_array(TYPE_S) * nd_arr_C,
                const TYPE_L alpha, const TYPE_L beta, const ND_int * A_idx,  const ND_int * B_idx,  const ND_int * C_idx)

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

/****************************************************************************************************/
                            /* TBLIS based*/
#if defined(COMPILE_ND_TBLIS)
/* Sum ND_function*/
void ND_function(sum, TYPE_S) (char * str_A, char * str_C, ND_array(TYPE_S) * nd_arrA, ND_array(TYPE_S) * nd_arrC, const TYPE_L alpha, const TYPE_L beta)
{
    //
    tblis_tensor A_arr, C_arr ;

    if (str_A == NULL || str_C == NULL) error_msg("sum ND_function indices cannot take NULL pointer") ;

    if (((nd_arrA->rank) == NULL )  || ((nd_arrC->rank) == NULL )) \
                                        error_msg("Cannot accept uninitilized array in sum ND_function") ; 
    
    if (((nd_arrA->data) == NULL ) || ((nd_arrC->data) == NULL )) \
                        error_msg("Cannot accept NULL data array. allocate them beforing passing to sum ND_function");

    if ((strlen(str_A) != (size_t) (nd_arrA->rank)[0]) || (strlen(str_C) != (size_t) (nd_arrC->rank)[0])) \
                            error_msg("Rank of input arrays doesn't match with given indices in sum ND_function");

    if (strlen(str_A) == 0) error_msg("Rank of input array for sum ND_function cannot be 0.");

    check_all_indices(str_A,str_C);//check if the final indices are present in input arrays
    check_str_char_unq(str_C, "Cannot have identical indices in output subscripts in sum ND_function");

    check_string_dims(str_A, str_A, nd_arrA->dims, nd_arrA->dims); // check is repeated indices are consistant in same array
    
    check_string_dims(str_C, str_C, nd_arrC->dims, nd_arrC->dims);

    check_string_dims(str_A, str_C, nd_arrA->dims, nd_arrC->dims); // check is repeated indices are consistant in other arrays
    
    //to_tblis_scaled
    ND_function(to_tblis_scaled, TYPE_S)(nd_arrA, &A_arr, alpha);

    ND_function(to_tblis_scaled, TYPE_S)(nd_arrC, &C_arr, beta);

    tblis_tensor_add(NULL, NULL, &A_arr, str_A, &C_arr, str_C);
    
    //
    nd_free_tblis(&A_arr);

    nd_free_tblis(&C_arr);

}


/* Einsum ND_function*/
void ND_function(einsum, TYPE_S) (char * einsum_indices, ND_array(TYPE_S) * nd_arrA, ND_array(TYPE_S) * nd_arrB, ND_array(TYPE_S) * nd_arrC, const TYPE_L alpha, const TYPE_L beta)
{
    //

    if (einsum_indices == NULL) error_msg("einsum indices cannot take NULL pointer") ;

    if ( (((nd_arrA->rank) == NULL )  || ((nd_arrB->rank) == NULL )) || ((nd_arrC->rank) == NULL ) ) \
                                                    error_msg("Cannot accept uninitilized array in einsum") ; 
    
    if ( (((nd_arrA->data) == NULL )  || ((nd_arrB->data) == NULL )) || ((nd_arrC->data) == NULL ) ) \
                        error_msg("Cannot accept NULL data array. allocate them beforing passing to einsum");
    
    tblis_tensor A_arr, B_arr, C_arr ;

    ND_int einsum_indices_len = strlen(einsum_indices) ; 
    
    if (einsum_indices_len>400) error_msg("Max length of Einsum indices string is 400");


    char * str_A = malloc(sizeof(char) * 5*400);
    char * str_B = str_A+400 ;
    char * str_C = str_A+800 ; 
    char * str_AB_cat = str_A+1200 ; // stores str_A + str_B
    //
    //
    ND_int ldaA = 0, ldaB = 0, ldaC = 0 ; 
    ND_int token_counter = 0;

    for (ND_int i =0 ; i<einsum_indices_len; ++i)
    {
        if (isspace(einsum_indices[i]) == 0)
        {
            if (token_counter == 0)
            {
                if (einsum_indices[i] == '-' || einsum_indices[i+1] == '>' ) error_msg("Incompatible einsum indices provided") ;

                else if (einsum_indices[i] == ',') token_counter++ ; 

                else
                {
                    str_A[ldaA] = einsum_indices[i];
                    ldaA++ ;
                }
            }
            else if (token_counter == 1)
            {
                if (einsum_indices[i] == ',') error_msg("Incompatible einsum indices provided");
                //
                else if (einsum_indices[i] == '-')
                {
                    if (i <einsum_indices_len-1)
                    {
                        if (einsum_indices[i+1] != '>') error_msg("Incompatible einsum indices provided");
                    }
                    else error_msg("Incompatible einsum indices provided");
                }
                //
                else if (einsum_indices[i] == '>')
                {
                    if (i>0)
                    {
                        if (einsum_indices[i-1] != '-') error_msg("Incompatible einsum indices provided");

                        else token_counter++ ;
                    }
                    else error_msg("Incompatible einsum indices provided");
                }
                //
                else{
                    str_B[ldaB] = einsum_indices[i];
                    ldaB++ ; 
                }
            }
            else if (token_counter == 2)
            {
                if ((einsum_indices[i] == '>' || einsum_indices[i] == '-') || einsum_indices[i] == ',') error_msg("Incompatible einsum indices provided");

                else
                {
                    str_C[ldaC] = einsum_indices[i];
                    ldaC++;
                }
                
            }
            else error_msg("Incompatible einsum indices provided");

        }
    }

    str_C[ldaC] = '\0'; str_B[ldaB] = '\0'; str_A[ldaA] = '\0';

    if (token_counter != 2) error_msg("Incompatible einsum indices provided");

    if (  ((strlen(str_A) != (size_t) (nd_arrA->rank)[0]) || (strlen(str_B) != (size_t) (nd_arrB->rank)[0]) ) || (strlen(str_C) != (size_t) (nd_arrC->rank)[0])) \
                            error_msg("Rank of input arrays doesn't match with given indices");

    if (  strlen(str_A) == 0 || strlen(str_B) == 0 ) error_msg("Rank of input arrays cannot be 0. Use scale/sum ND_function instead");

    strcpy(str_AB_cat,str_A);
    strcat(str_AB_cat,str_B);

    check_all_indices(str_AB_cat,str_C);//check if the final indices are present in input arrays
    check_str_char_unq(str_C, "Cannot have identical indices in output subscripts in einsum ND_function");

    check_string_dims(str_A, str_A, nd_arrA->dims, nd_arrA->dims); // check is repeated indices are consistant in same array
    check_string_dims(str_B, str_B, nd_arrB->dims, nd_arrB->dims);
    check_string_dims(str_C, str_C, nd_arrC->dims, nd_arrC->dims);

    check_string_dims(str_A, str_B, nd_arrA->dims, nd_arrB->dims); // check is repeated indices are consistant in other arrays
    check_string_dims(str_A, str_C, nd_arrA->dims, nd_arrC->dims);
    check_string_dims(str_B, str_C, nd_arrB->dims, nd_arrC->dims);

    //to_tblis_scaled
    if (  ND_function(size, TYPE_S) (nd_arrA)  < ND_function(size, TYPE_S) (nd_arrB) )
    {
        ND_function(to_tblis_scaled, TYPE_S)(nd_arrA, &A_arr, alpha);
        ND_function(to_tblis, TYPE_S)(nd_arrB, &B_arr);
    }
    else
    {
        ND_function(to_tblis, TYPE_S)(nd_arrA, &A_arr);
        ND_function(to_tblis_scaled, TYPE_S)(nd_arrB, &B_arr, alpha);
    }

    ND_function(to_tblis_scaled, TYPE_S)(nd_arrC, &C_arr, beta);

    /* 
        Find the indices of final array which contracted but only present in one array 
        for ex ijk, ijl ->ij k and l are contracted but only present in one array. 
        This is because, einsum only contracts correctly when contracted indices are present on both inputs.
    */
    
    get_einsum_excluded_idxs(str_A, str_B, str_C, str_AB_cat); // now str_AB_cat is used to ger missing indices

    if (strlen(str_AB_cat) == 0)
    {   /* Don only einsum if the str_AB_cat is empty */
        tblis_tensor_mult(NULL, NULL, &A_arr, str_A, &B_arr, str_B, &C_arr, str_C);
    }
    else
    {
        /* do einsum + reduction*/
        /*
            Get the dims of str_AB_cat indices -> create memory -> einsum -> add -> free
        */

        char * einsum_dummy_str = malloc(sizeof(char) * 2*400);

        strcpy(einsum_dummy_str,str_C);

        strcat(einsum_dummy_str,str_AB_cat);

        /* get the dimensions of the array of the indices in einsum_dummy_str */
        ND_int C_temp_rank =  strlen(einsum_dummy_str);

        ND_int * C_temp_dims = malloc(C_temp_rank * sizeof(ND_int));

        // get the dimensions;

        ND_int str_len_A , str_len_B ; 

        str_len_A = strlen(str_A);
        str_len_B = strlen(str_B);

        for (ND_int i_count =0 ; i_count < C_temp_rank ; ++i_count )
        {
            for (ND_int j_count =0; j_count <str_len_A ; ++j_count )
            {
                if ( str_A[j_count] == einsum_dummy_str[i_count]  ) C_temp_dims[i_count] = (nd_arrA->dims)[j_count] ; 
            }

            for (ND_int j_count =0; j_count <str_len_B ; ++j_count )
            {
                if ( str_B[j_count] == einsum_dummy_str[i_count]  ) C_temp_dims[i_count] = (nd_arrB->dims)[j_count] ; 
            }

        }

        ND_array(TYPE_S) nd_C_temp ; 

        tblis_tensor C_temp_blis;

        ND_function(init, TYPE_S) (&nd_C_temp, C_temp_rank, C_temp_dims); // initiate temp array

        ND_function(malloc, TYPE_S) (&nd_C_temp); // allocate memory

        ND_function(set_all, TYPE_S) (&nd_C_temp,  (TYPE_L) 0.0f ); // set all elements to zero

        ND_function(to_tblis, TYPE_S)(&nd_C_temp, &C_temp_blis);   // map to tblis tensor object
        
        tblis_tensor_mult(NULL, NULL, &A_arr, str_A, &B_arr, str_B, &C_temp_blis, einsum_dummy_str); // einsum first

        tblis_tensor_add(NULL, NULL, &C_temp_blis, einsum_dummy_str, &C_arr, str_C); // contraction

        nd_free_tblis(&C_temp_blis);    // free tblis tensor object

        ND_function(free, TYPE_S) (&nd_C_temp); // free the memory

        ND_function(uninit, TYPE_S) (&nd_C_temp); //initiate temp array

        free(C_temp_dims);

        free(einsum_dummy_str);

    }

    //
    free(str_A);
    nd_free_tblis(&A_arr);
    nd_free_tblis(&B_arr);

    nd_free_tblis(&C_arr);

}

#endif

/************************************* INTERNAL STATIC HELPER ND_functionS ******************************************************/

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


#if defined(COMPILE_ND_TBLIS)
static void check_string_dims(char * str_A, char * str_B, ND_int * dims_A,  ND_int * dims_B)
{   
    /* This ND_function checks if the string indices are consistant with the dims*/
    /* Before calling this ND_function, you must check the ranks of A and B match with string lenghts*/
    if (str_A == NULL || str_B == NULL) error_msg("Cannot pass NULL strings to check_string_dims");

    char * str_found_ptr;
    ND_int index_str_A;
    //
    ND_int str_len_B ;
    //
    str_len_B = strlen(str_B);
    //
    for (ND_int i = 0; i < str_len_B; ++i )
    {
        str_found_ptr = strchr(str_A, str_B[i]);
        if (str_found_ptr != NULL)
        {   
            index_str_A = (str_found_ptr - str_A);
            if (dims_A[index_str_A] != dims_B[i]) error_msg("Dimensions are inconsistant with indices provided");
        }
    }

}


static void check_all_indices(char * str_A, char * str_B)
{   
    /* This ND_function checks if the all string indices of B are present in A are consistant with the dims*/
    /* To pass, mutiple strings, concatnate and pass to str_A*/
    if (str_A == NULL || str_B == NULL) error_msg("Cannot pass NULL strings to check_all_indices");

    char * str_found_ptr;
    ND_int str_len_B ;

    str_len_B = strlen(str_B);

    for (ND_int i = 0; i < str_len_B; ++i )
    {
        str_found_ptr = strchr(str_A, str_B[i]);
        if (str_found_ptr == NULL) error_msg("Final indices are not present in the initial indices");
    }

}

/* ND_function to check if there is any repeatition in output indices*/
static void check_str_char_unq(char * inp_str, char* err_msg)
{
    int freq[256];
    
    for (ND_int i =0; i<256; ++i) freq[i] = 0 ;

    int ascii_to_int ;

    ND_int str_len = strlen(inp_str);

    for (ND_int i =0; i<str_len; ++i)
    {
        ascii_to_int = (int) inp_str[i];

        if (freq[ascii_to_int]>0) error_msg(err_msg) ; //"Cannot have identical indices in output subscripts");

        else ++freq[ascii_to_int];
    }
}

static void get_einsum_excluded_idxs(const char * str_A, const char* str_B, const char* str_C, char * out_str)
{   /* The ND_function computes elements that are present in one string but not in other.
        length of out_str is set to  str_A + str_B+ str_C + 1 to be safe. These indices must
        be appended to einsum output indices and later contract them */
    
    //
    int freqA[256];
    int freqB[256];
    int freqC[256];

    
    for (ND_int i =0; i<256; ++i) freqA[i] = 0 ;
    for (ND_int i =0; i<256; ++i) freqB[i] = 0 ;
    for (ND_int i =0; i<256; ++i) freqC[i] = 0 ;

    int ascii_to_int ;

    ND_int str_lenA = strlen(str_A);
    ND_int str_lenB = strlen(str_B);
    ND_int str_lenC = strlen(str_C);

    for (ND_int i =0; i<str_lenA; ++i)
    {
        ascii_to_int = (int) str_A[i];

        freqA[ascii_to_int]++;
    }

    for (ND_int i =0; i<str_lenB; ++i)
    {
        ascii_to_int = (int) str_B[i];

        freqB[ascii_to_int]++;
    }

    for (ND_int i =0; i<str_lenC; ++i)
    {
        ascii_to_int = (int) str_C[i];

        freqC[ascii_to_int]++;
    }

    ND_int str_counter = 0 ;

    for (int i =0; i<256; ++i)
    {
        if ( ((freqA[i] >0 && freqB[i] ==0 ) || ( freqA[i] == 0 && freqB[i] >0  )) && (freqC[i] == 0) )
        {   
            out_str[str_counter] = (char) i ; 
            str_counter++ ; 
        }
    }
    out_str[str_counter] = '\0' ; 

}


/* Convert nd_array to tblis array*/
static void ND_function(to_tblis, TYPE_S)(const ND_array(TYPE_S) * nd_arr_A, tblis_tensor * A_arr)
{   
    /*Must be freed after using*/

    if (nd_arr_A->rank == NULL) error_msg("Cannot convert a nd_array to blis array");
    if (nd_arr_A->data == NULL) error_msg("Cannot convert a nd_array to blis array");

    /* Get max values to check for overflow*/
    intmax_t max_len, max_stride;

    for (max_len=INTMAX_MAX; (len_type)max_len!=max_len; max_len/=2);
    for (max_stride=INTMAX_MAX; (stride_type)max_stride!=max_stride; max_stride/=2);
    
    
    ND_int rank = (nd_arr_A->rank)[0];
    len_type * blis_dims         =  malloc(sizeof(len_type)*rank) ;
    stride_type * blis_strides =  malloc(sizeof(stride_type)*rank) ;
    
    for (ND_int i =0; i<rank; ++i )
    {
        if ( ((nd_arr_A->dims)[i] >= (ND_int)max_len  ) || ((nd_arr_A->strides)[i] >= (ND_int)max_stride ) ) \
                error_msg("Indices overflow. Complile tblis with larger bit indices for len_type and stride_type");

        else
        {
            blis_dims[i] =    (len_type) (nd_arr_A->dims)[i] ;
            blis_strides[i] = (stride_type) (nd_arr_A->strides)[i] ;
        }
    }

    TBLIS_CALL(init_tensor,TYPE_S) (A_arr, rank , blis_dims, nd_arr_A->data, blis_strides);

    // free(blis_dims);
    // free(blis_strides);
}



/* Convert nd_array to tblis array with scaling*/
static void ND_function(to_tblis_scaled, TYPE_S)(const ND_array(TYPE_S) * nd_arr_A, tblis_tensor * A_arr, const TYPE_L alpha)
{
    /*Must be freed after using*/

    if (nd_arr_A->rank == NULL) error_msg("Cannot convert a nd_array to blis array");
    if (nd_arr_A->data == NULL) error_msg("Cannot convert a nd_array to blis array");

    /* Get max values to check for overflow*/
    intmax_t max_len, max_stride;

    for (max_len=INTMAX_MAX; (len_type)max_len!=max_len; max_len/=2);
    for (max_stride=INTMAX_MAX; (stride_type)max_stride!=max_stride; max_stride/=2);
    
    
    ND_int rank = (nd_arr_A->rank)[0];
    len_type * blis_dims         =  malloc(sizeof(len_type)*rank) ;
    stride_type * blis_strides =  malloc(sizeof(stride_type)*rank) ;
    
    for (ND_int i =0; i<rank; ++i )
    {
        if ( ((nd_arr_A->dims)[i] >= (ND_int)max_len  ) || ((nd_arr_A->strides)[i] >= (ND_int)max_stride ) ) \
                error_msg("Indices overflow. Complile tblis with larger bit indices for len_type and stride_type");
        
        else
        {
            blis_dims[i] =     (nd_arr_A->dims)[i] ;
            blis_strides[i] =  (nd_arr_A->strides)[i] ;
        }
    }

    TBLIS_CALL(init_tensor_scaled,TYPE_S) (A_arr, alpha, rank , blis_dims, nd_arr_A->data, blis_strides);

    // free(blis_dims);
    // free(blis_strides);
}



/* Free tblis arrays*/
static void nd_free_tblis(tblis_tensor * A_arr)
{   /*
    Note: This ND_function only free, dimensions and strinds, It DOES NOT free the main data.
    */
    if ( (A_arr->len == NULL )|| (A_arr->stride == NULL)) error_msg("tblis_tensor already freed");
    //
    else
    {
        free(A_arr->len);
        free(A_arr->stride);
        A_arr->len = NULL ;
        A_arr->stride = NULL ;
    }
}
#endif

#endif


/* Level -1 Blas*/

/* Level -2 Blas*/

/* Level -3 Blas*/
/* Gemm. i-> j-->k  C(.....,i,k) = A(....,i,j)*B(.....,j,k) 
    C(...,m,n) = A(....,m,k)@B(....,k,n),
    "(2,3,:,:) (3,5,:,:) -> (4,5,:,:)", m, n, k
*/
