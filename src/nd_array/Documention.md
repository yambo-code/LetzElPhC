# Documentation


Please look into the test.c file in test folder to get an idea of this library usage

**?** denote s/d/c/z/i . For example nd_array_s stores data with float type
> s -> *float type*
> d -> *double type*
> c -> *float complex type*
> z -> *double complex type*
> i -> *int type*

## Basic array type in nd_array library
The basic array object in this library is nd_array_? . For example, to nd_array_c denotes a complex array. The array type is defined as follow:
``` c
typedef struct nd_arr_? {
    type * data;
    ND_indices * rank;
    ND_indices * dims;
    ND_indices * strides;
    bool owner;

} nd_arr_?;
```

In the code, the nd_array_? type has to go through following phases
> - define the object
> - initiate the object
> - allocate/assign data to the created object
> - use the object ....
> - free the data (if allocated). If assgined no need to free the oject
> - uninitialize the object


Below are the functions that are defined for these objects:
**Note that all the functions are pass by reference w.r.t to nd_array_? types**

# Functions
----
Function to initiate the array:
``` c
void nd_init_? (nd_arr_? * nd_arr_in, const ND_indices rank, const ND_indices * dimensions);
```
----
Function to destroy the array:
``` c
void nd_uninit_? (nd_arr_? * nd_arr_in);
``` 
----
Function to allocate data the array:
``` c
void nd_malloc_? (nd_arr_? * nd_arr_in);
``` 
----
Function to allocate data to array with all bits set to zero bit initiated:
``` c
void nd_calloc_? (nd_arr_? * nd_arr_in);
``` 
----
free the data if allocated
``` c
void nd_free_? (nd_arr_? * nd_arr_in);
``` 
----
Same as nd_init_ but specially for initalizing arrays for transposing array (to make it easier)
``` c
void nd_init_tranpose_? (const nd_arr_? * nd_arr_in, const ND_indices * order, nd_arr_? * nd_arr_out);
``` 
----
Same as nd_init_ but specially for initalizing arrays for slicing (to make it easier)
``` c
void nd_init_slice_? (const ND_indices * start_idx, const ND_indices * end_idx, const ND_indices * step_idx, const nd_arr_? * nd_arr_in, nd_arr_? * nd_arr_out );
``` 
----
Same as nd_init_ but specially for initalizing arrays for slipping the dims (to make it easier)
``` c
void nd_init_strip_dims_? (const nd_arr_? * nd_arr_in, const ND_indices n_dims_strip, nd_arr_? * nd_arr_out);
``` 
----
get pointer to that element
``` c
float * nd_ele_? (const nd_arr_? * nd_arr_in, const ND_indices * dimensions);
``` 
----
get size of array
``` c
ND_indices nd_size_? (const nd_arr_? * nd_arr_in);
``` 
----
set all elements to some constant
``` c
void nd_set_all_? (nd_arr_? * nd_arr_in, const float set_constant );
``` 
----
reshape the array
``` c
void nd_reshape_? (const nd_arr_? * nd_arr_in, nd_arr_? * nd_arr_out);
```
----
Strinp the first n dimesntions. No data is to be created is created. returns view
``` c
void nd_strip_dims_? (const nd_arr_? * nd_arr_in, const ND_indices n_dims_strip, const ND_indices * stripped_idxs, nd_arr_? * nd_arr_out);
``` 
----
creates the slice: you must allocate the data and then pass
``` c
void nd_slice_? (const ND_indices * start_idx, const ND_indices * end_idx, const ND_indices * step_idx, const nd_arr_? * nd_arr_in, nd_arr_? * nd_arr_out );
``` 
----
creates the transpose: you must allocate the data and then pass
``` c
void nd_tranpose_? (const nd_arr_? * nd_arr_in, const ND_indices * order, nd_arr_? * nd_arr_out);
``` 
----
deep compy from one to ther array
``` c
void nd_copy_? (const nd_arr_? * nd_arr_in, nd_arr_? * nd_arr_out);
``` 
----
Read data from netCDF. **Note that DONOT create data and pass to this function. This function will automatically create data for you**
``` c
void nd_read_? (const char* file_name, const char* var_name, nd_arr_? * nd_arr_in);
``` 
----
Read subarray data from netCDF. **Note that DONOT create data and pass to this function. This function will automatically create data for you** .
Few MACROS:
Unline slice function, we pass ndim arrays with start, end element, step. USe **ND_END** for end and use **ND_ALL** to read entire dimension 
``` c
void nd_read_sub_? (const char* file_name, const char* var_name, nd_arr_? * nd_arr_in, ...);
// nd_read_sub_c("nc.temp", "exc_elph", &read_sub_arr, ND_ALL, nd_idx{1,3,1}, nd_idx{173,356,7}, nd_idx{324,ND_END,9} );  == [: , 1:3:1, 173:356:7, 324::9]

``` 
----
Write data from netCDF.
``` c
void nd_write_? (const char* file_name, const char* var_name, const nd_arr_? * nd_arr_in, char ** dim_names, size_t * chunk_size);
``` 
If chunk_size is set to NULL, then it is stored as single Contiguous array. Note, for when writing complex array, the size of chunk_size array must be (ndim+1) as there is additional rea/imag dimension.

----
Matmul
> 'N' Normal
> 'T' for transpose
> 'C' Conjugate  + Transpose

A_idx, B_idx, C_idx are indices of first (n-2) dimensions
``` c
void nd_matmul_? (const char TransA, const char TransB, const nd_arr_? * nd_arr_A, const nd_arr_? * nd_arr_B, nd_arr_? * nd_arr_C, const double complex alpha, const double complex beta, const ND_indices * A_idx, const ND_indices * B_idx, const ND_indices * C_idx);

Ex: nd_matmul_c('T', 'C', &write_arr, &write_arr, &mat_prod, 1.0f + 0.0f*I, 0.0f + 0.0f*I, nd_idx{0,1},  nd_idx{0,2}, NULL /*nd_idx{},  nd_idx{}*/);
``` 
----
Sum + transpose 
``` c
void nd_sum_? (char * str_A, char * str_C, nd_arr_? * nd_arrA, nd_arr_? * nd_arrC, const double complex alpha, const double complex beta);
// Ex: nd_sum_c("ijl","jl",&einsum_test,&sum_test,1.0f + 0.0f*I, 0.0f + 0.0f*I);
``` 
----
Einsum
``` c
void nd_einsum_? (char * einsum_indices, nd_arr_? * nd_arrA, nd_arr_? * nd_arrB, nd_arr_? * nd_arrC, const double complex alpha, const double complex beta);
// Ex: nd_einsum_c("ijkx,ijky->",&write_arr,&write_arr,&einsum_test2,1.0f + 0.0f*I, 0.0f + 0.0f*I);
``` 
----



