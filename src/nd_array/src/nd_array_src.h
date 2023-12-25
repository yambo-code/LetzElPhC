//Compile with C99 or above
// Written by Muralidhar Nalabothula July 2022
/**
 * @brief A simple library for ND arrays and perform blas operations on them. Read and write I/O is 
 * performed using netCDF library
 * Must link some blas and netCDF
 * Contains ND_functions for float(s), double(d), float complex (c), double complex (z)
 * array object nd_arr(##TYPE) Ex: for double complex nd_arrz
 * all indices, counters must be in ND_int
 */

#ifdef __STDC_NO_COMPLEX__
#error Your compiler does not support C99 complex numbers, Please use a supported compiler.
#endif

#ifndef GENERATE_HEADERS // this is used to generate headers
#pragma once
#include <stdio.h>      
#include <complex.h>
#include <stddef.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdbool.h> 
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <netcdf.h>
#if defined(COMPILE_ND_TBLIS)
    #include <tblis/tblis.h>
#endif
#include <stdint.h>
#include "nd_error.h"
#include <fftw3.h>
//
#include "common_def.h" 

#ifdef ND_MKL_BLAS
    typedef MKL_INT BLAS_INT;
#else
    #ifdef USE64BITINT
        typedef int64_t BLAS_INT;
    #else
        typedef int BLAS_INT; /* Internal type used for blas indices. set this according to the blas library (32 or 64 bit) */   
    #endif
#endif

#define CBLAS_INT BLAS_INT

#ifdef ND_MKL_BLAS
    #define lapack_complex_float    float _Complex
    #define lapack_complex_double   double _Complex
    #include <mkl.h>
#elif defined ND_DEFAULT_BLAS_HEADER
    #include "cblas.h"
#else
    #define LAPACK_COMPLEX_C99
    #include <cblas.h>
//    #include <lapacke.h> // for now no lapack is needed
#endif

//


#define ERR(e) {nd_error_msg(nc_strerror(e), "netcdf function");}

#define BLAS_CALL(FUN_NAME, TYPE_SMALL)         BLAS_CALL_HIDDEN(FUN_NAME, TYPE_SMALL)
#define BLAS_CALL_HIDDEN(FUN_NAME, TYPE_SMALL)  cblas_ ## TYPE_SMALL ## FUN_NAME

#define TBLIS_CALL(FUN_NAME, TYPE_SMALL)         TBLIS_CALL_HIDDEN(FUN_NAME, TYPE_SMALL)
#define TBLIS_CALL_HIDDEN(FUN_NAME, TYPE_SMALL)  tblis_ ## FUN_NAME ## _ ## TYPE_SMALL 

#else // GENERATE_HEADERS
#define ND_function(FUN_NAME, TYPE_SMALL)         ND_function_HIDDEN(FUN_NAME, TYPE_SMALL)
#define ND_function_HIDDEN(FUN_NAME, TYPE_SMALL)  nd_ ## FUN_NAME ## _ ## TYPE_SMALL

#define ND_array(TYPE_SMALL)         ND_array_HIDDEN(TYPE_SMALL)
#define ND_array_HIDDEN(TYPE_SMALL)  nd_arr ## _ ## TYPE_SMALL
#endif // GENERATE_HEADERS


#if defined(COMPILE_ND_FLOAT)
    #define TYPE_S s    /* */
    #define TYPE_L float /* type */
    #define NetCDF_IO_TYPE NC_FLOAT
    #define NetCDF_FUN_TYPE float

#elif defined(COMPILE_ND_DOUBLE)
    #define TYPE_S d    /* */
    #define TYPE_L double /* type */
    #define NetCDF_IO_TYPE NC_DOUBLE
    #define NetCDF_FUN_TYPE double

#elif defined(COMPILE_ND_SINGLE_COMPLEX)
    #define TYPE_S c    /* */
    #define TYPE_L float complex /* type */
    #define NetCDF_IO_TYPE NC_FLOAT
    #define NetCDF_FUN_TYPE float

#elif defined(COMPILE_ND_DOUBLE_COMPLEX)
    #define TYPE_S z    /* */
    #define TYPE_L double complex /* type */
    #define NetCDF_IO_TYPE NC_DOUBLE
    #define NetCDF_FUN_TYPE double

#elif defined(COMPILE_ND_INT)
    #define TYPE_S i    /* */
    #define TYPE_L int /* type */
    #define NetCDF_IO_TYPE NC_INT
    #define NetCDF_FUN_TYPE int

#endif


#if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)
    #define BLAS_POINTER(VARIABLE_NAME) &VARIABLE_NAME    /* for C and Z, pointer is passed for constant in blas routines */
#else
    #define BLAS_POINTER(VARIABLE_NAME) VARIABLE_NAME    /*  for S and D, value is passed for constant in blas routines */
#endif

/*********************************************************************************************************************************/


typedef struct ND_array(TYPE_S) {
    TYPE_L * data;        // store data pointer
    ND_int * rank;    // rank of the tensor
    ND_int * dims;    // pointer to dims array
    ND_int * strides; // pointer to strides of an array (in elements and not in bytes)
    bool owner; /* set this to true if data belongs to this array. if it is referened then false. 
                when array is free data is not freed if owner is false */
} ND_array(TYPE_S);



/**************************************************** alloc.c ND_functions **********************************************************/
/*********************************************************************************************************************************/

void ND_function(init, TYPE_S) (ND_array(TYPE_S) * nd_arr_in, const ND_int rank, const ND_int * dimensions);

void ND_function(uninit, TYPE_S) (ND_array(TYPE_S) * nd_arr_in);

void ND_function(malloc, TYPE_S) (ND_array(TYPE_S) * nd_arr_in);

void ND_function(calloc, TYPE_S) (ND_array(TYPE_S) * nd_arr_in);

void ND_function(free, TYPE_S) (ND_array(TYPE_S) * nd_arr_in);

void ND_function(destroy, TYPE_S) (ND_array(TYPE_S) * nd_arr_in);

void ND_function(init_tranpose, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, const ND_int * order, ND_array(TYPE_S) * nd_arr_out);

void ND_function(init_slice, TYPE_S) (const ND_int * start_idx, const ND_int * end_idx, const ND_int * step_idx, \
                            const ND_array(TYPE_S) * nd_arr_in, ND_array(TYPE_S) * nd_arr_out );

void ND_function(init_strip_dims, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, const ND_int n_dims_strip, ND_array(TYPE_S) * nd_arr_out);

/**************************************************** array.c ND_functions **********************************************************/
/*********************************************************************************************************************************/

TYPE_L * ND_function(ele, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, const ND_int * dimensions);

ND_int ND_function(size, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in);

void ND_function(set_all, TYPE_S) (ND_array(TYPE_S) * nd_arr_in, const TYPE_L set_constant );

void ND_function(reshape, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, ND_array(TYPE_S) * nd_arr_out);

void ND_function(strip_dims, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, const ND_int n_dims_strip, \
                                    const ND_int * stripped_idxs, ND_array(TYPE_S) * nd_arr_out);

void ND_function(slice, TYPE_S) (const ND_int * start_idx, const ND_int * end_idx, const ND_int * step_idx, \
                            const ND_array(TYPE_S) * nd_arr_in, ND_array(TYPE_S) * nd_arr_out );

void ND_function(tranpose, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, const ND_int * order, ND_array(TYPE_S) * nd_arr_out);

void ND_function(copy, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, ND_array(TYPE_S) * nd_arr_out);


/**************************************************** netcdf_io.c ND_functions ******************************************************/
/*********************************************************************************************************************************/


void ND_function(read, TYPE_S) (const char* file_name, const char* var_name, ND_array(TYPE_S) * nd_arr_in);

void ND_function(read_sub, TYPE_S) (const char* file_name, const char* var_name, ND_array(TYPE_S) * nd_arr_in, ...);

void ND_function(write, TYPE_S) (const char* file_name, const char* var_name, const ND_array(TYPE_S) * nd_arr_in, char ** dim_names, size_t * chunksize);

#if defined(COMPILE_ONCE)
void NC_open_file(const char* file_name, char mode, int * ncid);

void NC_close_file(const int ncid); 
#endif

void ND_function(readVar, TYPE_S) (const int ncid, const char* var_name, ND_array(TYPE_S) * nd_arr_in);
void ND_function(readVar_sub, TYPE_S) (const int ncid, const char* var_name, ND_array(TYPE_S) * nd_arr_in, ...);
void ND_function(writeVar, TYPE_S) (const int ncid, const char* var_name, const ND_array(TYPE_S) * nd_arr_in, char ** dim_names, size_t * chunksize);
/*************************************************** linalg.c ND_functions **********************************************************/
/*********************************************************************************************************************************/

#if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX) || defined(COMPILE_ND_FLOAT) || defined(COMPILE_ND_DOUBLE)

void ND_function(matmul, TYPE_S) (const char TransA, const char TransB, const ND_array(TYPE_S) * nd_arr_A, const ND_array(TYPE_S) * nd_arr_B, \
                                ND_array(TYPE_S) * nd_arr_C, const TYPE_L alpha, const TYPE_L beta, const ND_int * A_idx,  \
                                                                        const ND_int * B_idx,  const ND_int * C_idx);

/** Matmul Expert version.*/
void ND_function(matmulX, TYPE_S) (const char TransA, const char TransB, const TYPE_L * arr_A, const TYPE_L * arr_B, TYPE_L * arr_C, \
                const TYPE_L alpha, const TYPE_L beta, const ND_int ldA, const ND_int ldB, const ND_int ldC, \
                const ND_int m, const ND_int n, const ND_int k);

#if defined(COMPILE_ND_TBLIS)
void ND_function(sum, TYPE_S) (char * str_A, char * str_C, ND_array(TYPE_S) * nd_arrA, ND_array(TYPE_S) * nd_arrC, const TYPE_L alpha, const TYPE_L beta);

void ND_function(einsum, TYPE_S) (char * einsum_indices, ND_array(TYPE_S) * nd_arrA, ND_array(TYPE_S) * nd_arrB, ND_array(TYPE_S) * nd_arrC, \
                                                                                                const TYPE_L alpha, const TYPE_L beta);
#endif

#endif



/*************************************************** nd_fft.c functions (FFT) ****************************************************/
/*********************************************************************************************************************************/
/** these are fftw malloc and free. returns aligned memory to make use of simd fftw*/
#if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX) //|| defined(COMPILE_ND_FLOAT) || defined(COMPILE_ND_DOUBLE)


#if defined(COMPILE_ND_DOUBLE_COMPLEX)
typedef fftw_plan ND_function(FFT_plan, TYPE_S);
#elif defined(COMPILE_ND_SINGLE_COMPLEX)
typedef fftwf_plan ND_function(FFT_plan, TYPE_S);
#endif

void ND_function(FFT_malloc, TYPE_S) (ND_array(TYPE_S) * nd_arr_in);
void ND_function(FFT_calloc, TYPE_S) (ND_array(TYPE_S) * nd_arr_in);

void ND_function(FFT_free, TYPE_S) (ND_array(TYPE_S) * nd_arr_in);
void ND_function(FFT_destroy, TYPE_S) (ND_array(TYPE_S) * nd_arr_in);

void ND_function(fft, TYPE_S) (const ND_array(TYPE_S) * nd_arr_A, ND_array(TYPE_S) * nd_arr_B, const ND_int rank, \
                               const ND_int * in_idx, int exp_sign, bool normalize);


void ND_function(fft_planner, TYPE_S) (const ND_array(TYPE_S) * nd_arr_A, ND_array(TYPE_S) * nd_arr_B, const ND_int rank, \
                        const ND_int * in_idx, int exp_sign, TYPE_L * norm_out, unsigned flag, ND_function(FFT_plan, TYPE_S) * fftplan);


void ND_function(fft_execute_plan, TYPE_S) (ND_function(FFT_plan, TYPE_S) fftplan);

void ND_function(fft_destroy_plan, TYPE_S) (ND_function(FFT_plan, TYPE_S) fftplan);

int ND_function(fft_export_wisdom_to_file, TYPE_S) (const char *filename);

int ND_function(fft_import_wisdom_to_file, TYPE_S) (const char *filename);

void ND_function(fft_forget_wisdom, TYPE_S) (void);

void ND_function(fft_cleanup, TYPE_S) (void);

int ND_function(fft_init_threads, TYPE_S) (void);

void ND_function(fft_plan_with_nthreads, TYPE_S) (int nthreads);

void ND_function(fft_cleanup_threads, TYPE_S) (void);
#endif



