#include "nd_array_src.h"
#include <complex.h> 
#include <string.h>


#if defined(COMPILE_ND_DOUBLE_COMPLEX)
    #define FFTW_CALL(FUN_NAME) FFTW_CALL_HIDDEN(FUN_NAME)
    #define FFTW_CALL_HIDDEN(FUN_NAME) fftw_##FUN_NAME
#elif defined(COMPILE_ND_SINGLE_COMPLEX)
    #define FFTW_CALL(FUN_NAME) FFTW_CALL_HIDDEN(FUN_NAME)
    #define FFTW_CALL_HIDDEN(FUN_NAME) fftwf_##FUN_NAME
#endif

/*For now only complex FFTs*/
#if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX) //|| defined(COMPILE_ND_FLOAT) || defined(COMPILE_ND_DOUBLE)

static bool is_idx_present(ND_int idx, const ND_int * idx_arr, ND_int len_arr)
{   
    for (ND_int i = 0; i < len_arr; ++i)
    {
        if (idx == idx_arr[i]) return true;
    }
    return false;
}

/******* These are FFTW versions of malloc and free so that we get aligned array*/
void ND_function(FFT_malloc, TYPE_S) (ND_array(TYPE_S) * nd_arr_in)
{
    /* Allocated memory for the data. Must call free ND_function.
    */
    if (nd_arr_in->owner == true) error_msg("Overwriting already allocated array (i.e, Owner array). Possible memory leak !") ;

    ND_int size_arr = 1;
    //
    if (*(nd_arr_in->rank) ==  0)
    {
        size_arr =   1 ;
    }

    else
    {
        size_arr = *(nd_arr_in->dims) * *(nd_arr_in->strides);
    }
    //

    nd_arr_in->owner = true ; /* Set this to true as it owns the data created*/
    /* Create data pointer*/
    nd_arr_in->data = FFTW_CALL(malloc)(size_arr * sizeof(TYPE_L));

    if (nd_arr_in->data == NULL) error_msg("Failed to allocate data for nd_malloc array") ;
}


void ND_function(FFT_calloc, TYPE_S) (ND_array(TYPE_S) * nd_arr_in)
{
    /* Allocated memory for the data. Must call free ND_function.
    */
    ND_function(FFT_malloc, TYPE_S) (nd_arr_in);
    size_t size_arr = sizeof(TYPE_L) * ND_function(size, TYPE_S) (nd_arr_in);
    memset(nd_arr_in->data,0,size_arr);
}



/* Free ND_function for FFTW_malloc array */
void ND_function(FFT_free, TYPE_S) (ND_array(TYPE_S) * nd_arr_in)
{
    /* 
    Must be called when nd_malloc_? or nd_calloc_? ND_functions are called, else memory leak !
    */
    if (nd_arr_in->owner)
    {   
        if (nd_arr_in->data == NULL) error_msg("Trying to free array with NULL data");
        FFTW_CALL(free)(nd_arr_in->data);
        nd_arr_in->data = NULL ; 
        nd_arr_in->owner = false;
    }
}

/* Destroy array created with for FFT_malloc */
void ND_function(FFT_destroy, TYPE_S) (ND_array(TYPE_S) * nd_arr_in)
{
    /* simply free + uninit */
    ND_function(FFT_free, TYPE_S) (nd_arr_in);
    ND_function(uninit, TYPE_S) (nd_arr_in);
}


/** this is quick fft routine */
void ND_function(fft, TYPE_S) (const ND_array(TYPE_S) * nd_arr_A, ND_array(TYPE_S) * nd_arr_B, const ND_int rank, \
                               const ND_int * in_idx, int exp_sign, bool normalize)
{
    FFTW_CALL(plan) fftplan;

    TYPE_L Norm_factor = 0;
    ND_function(fft_planner, TYPE_S) (nd_arr_A, nd_arr_B, rank, in_idx, exp_sign, \
                                        &Norm_factor, FFTW_ESTIMATE, &fftplan);
    
    if (fftplan == NULL) error_msg("Failed to create the FFT plan"); 

    FFTW_CALL(execute)(fftplan);

    FFTW_CALL(destroy_plan)(fftplan);

    if (normalize)
    {
        Norm_factor = 1.0/Norm_factor ; 
        ND_int size_arr = ND_function(size,TYPE_S)(nd_arr_A);
        for (ND_int i = 0; i<size_arr; ++i) nd_arr_B->data[i] = nd_arr_B->data[i]*Norm_factor ;
    }

    return ;
}



/* Create a plan for the arrays */
void ND_function(fft_planner, TYPE_S) (const ND_array(TYPE_S) * nd_arr_A, ND_array(TYPE_S) * nd_arr_B, const ND_int rank, \
                        const ND_int * in_idx, int exp_sign, TYPE_L * norm_out, unsigned flag, FFTW_CALL(plan) * fftplan)
{
    /* This will perform the ndimentional (ndim) FFT/invFFT */
    
    /* Some checks on data */
    // if (normalize) printf("Using normalize");
    if (( nd_arr_A->rank == NULL )  || ( nd_arr_B->rank == NULL )) \
                                                    error_msg("Cannot accept uninitilized array") ; 
    
    if ( (nd_arr_A->data == NULL )  || (nd_arr_B->data == NULL ) ) \
                    error_msg("Cannot accept NULL array. allocate them beforing passing to matmul");

    if (nd_arr_A->rank[0] != nd_arr_B->rank[0]) error_msg("Rank mismatch between input and output arrays");
    if (rank > nd_arr_A->rank[0]) error_msg("FFT Rank cannot be greater than Dimension of array"); 

    for (ND_int i = 0 ; i<nd_arr_A->rank[0]; ++i)
    {   
        if (nd_arr_A->dims[i] != nd_arr_B->dims[i]) \
            error_msg("Dimension mismatch between input and output arrays");  
    }

    
    int howmany_rank =  nd_arr_A->rank[0]-rank ;
    int FFT_rank = rank;
    
    FFTW_CALL(iodim64) * dims = malloc(rank*sizeof(FFTW_CALL(iodim64)));

    FFTW_CALL(iodim64) * howmany_dims = malloc(howmany_rank*sizeof(FFTW_CALL(iodim64))) ;

    if (dims == NULL || howmany_dims == NULL)
    {   
        error_msg("Error allocating FFT dimension arrays");  
        return ;
    }

    ND_int  dimCnt =0;
    ND_int hdimCnt =0;

    TYPE_L Norm_factor = 1.0;

    for (ND_int i = 0 ; i< nd_arr_A->rank[0]; ++i)
    {   
        if (is_idx_present(i,in_idx,rank))
        {
            dims[dimCnt].n  =  nd_arr_A->dims[i] ;
            dims[dimCnt].is =  nd_arr_A->strides[i];
            dims[dimCnt].os =  nd_arr_B->strides[i];
            ++dimCnt ;
            Norm_factor *= nd_arr_A->dims[i];
        }
        else
        {
            howmany_dims[hdimCnt].n  =  nd_arr_A->dims[i] ;
            howmany_dims[hdimCnt].is =  nd_arr_A->strides[i];
            howmany_dims[hdimCnt].os =  nd_arr_B->strides[i];
            ++hdimCnt;
        }
    }

    /*Debug*/
    if (dimCnt != rank || hdimCnt != (nd_arr_A->rank[0]-rank) )
    {   
        error_msg("Error allocating FFT dimension arrays");  
        return ;
    }

    *fftplan = FFTW_CALL(plan_guru64_dft)(FFT_rank, dims,howmany_rank, howmany_dims, \
                nd_arr_A->data, nd_arr_B->data, exp_sign, flag);
    
    free(dims); free(howmany_dims);
    
    if (norm_out != NULL) *norm_out = Norm_factor ; 

    return ;
}

void ND_function(fft_execute_plan, TYPE_S) (FFTW_CALL(plan) fftplan)
{
    FFTW_CALL(execute)(fftplan);
}


void ND_function(fft_destroy_plan, TYPE_S) (FFTW_CALL(plan) fftplan)
{
    FFTW_CALL(destroy_plan)(fftplan);
}



int ND_function(fft_export_wisdom_to_file, TYPE_S) (const char *filename)
{
    return FFTW_CALL(export_wisdom_to_filename)(filename);
}


int ND_function(fft_import_wisdom_to_file, TYPE_S) (const char *filename)
{
    return FFTW_CALL(import_wisdom_from_filename)(filename);
}


void ND_function(fft_forget_wisdom, TYPE_S) (void)
{
    FFTW_CALL(forget_wisdom)();
}

void ND_function(fft_cleanup, TYPE_S) (void)
{
    FFTW_CALL(cleanup)();
}

/*Threads*/
int ND_function(fft_init_threads, TYPE_S) (void)
{
    return FFTW_CALL(init_threads)();
}

void ND_function(fft_plan_with_nthreads, TYPE_S) (int nthreads)
{
    FFTW_CALL(plan_with_nthreads)(nthreads);
}

void ND_function(fft_cleanup_threads, TYPE_S) (void)
{
    FFTW_CALL(cleanup_threads)();
}

#endif