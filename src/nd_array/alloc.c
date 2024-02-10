#include "nd_array_src.h"

/********************************************************************************************************************************************/
/********************************************************* INITIALIZATION  ND_functionS ********************************************************/


void ND_function(init, TYPE_S) (ND_array(TYPE_S) * nd_arr_in, const ND_int rank, const ND_int * dimensions)
{
    /* first initiate the ND_int array. Only stirdes, rank and dim arrays are created
    */
    ND_int size_arr = 1;

    nd_arr_in->owner = false ; /* When initialized, set the data to false.*/ 

    ND_int * rank_array = malloc((1 + 2*rank) * sizeof(ND_int)); /* (2*rank + 1) elements*/

    if (rank_array == NULL) 
        error_msg("Failed to allocate dimensions for nd_malloc array");
    
    nd_arr_in->rank = rank_array;
    //
    if (rank == 0)
    {
        nd_arr_in->dims = NULL;
        nd_arr_in->strides = NULL;
    }
    else
    {
        nd_arr_in->dims = rank_array+1;
        nd_arr_in->strides = rank_array+1+rank;
    }
    rank_array[0] = rank;
    /* Get all the args of type ND_int*/
    
    if (rank != 0) memcpy(rank_array+1,dimensions,sizeof(ND_int)*rank);
    /* set strides */ 
    for (ND_int i = 0; i < rank; ++i)
    {
        rank_array[2*rank - i] = size_arr;
        size_arr = size_arr*rank_array[rank - i];
    }
    //
    nd_arr_in->data = NULL ; 
    //
}


void ND_function(uninit, TYPE_S) (ND_array(TYPE_S) * nd_arr_in)
{
    /* Note that rank, dims, strides are all created with single malloc with rank being first pointer
        So Only free data and rank pointers
    */
    if (nd_arr_in->rank == NULL ) 
        error_msg("Array already uninitialized");

    else if (nd_arr_in->owner) 
        error_msg("Free the array before uninitialization");

    else 
    {
        free(nd_arr_in->rank);
        nd_arr_in->rank = NULL ; 
    }
}



/********************************************************************************************************************************************/
/********************************************************* MALLOC/CALLOC/FREE  ND_functionS ****************************************************/

/* malloc ND_function*/
void ND_function(malloc, TYPE_S) (ND_array(TYPE_S) * nd_arr_in)
{
    /* Allocated memory for the data. Must call free ND_function.
    */
    if (nd_arr_in->owner == true) 
        error_msg("Overwriting already allocated array (i.e, Owner array). Possible memory leak !") ;

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
    nd_arr_in->data = malloc(size_arr * sizeof(TYPE_L));

    if (nd_arr_in->data == NULL) error_msg("Failed to allocate data for nd_malloc array") ;
}



/* calloc ND_function*/
void ND_function(calloc, TYPE_S) (ND_array(TYPE_S) * nd_arr_in)
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
    nd_arr_in->data = calloc(size_arr, sizeof(TYPE_L));

    if (nd_arr_in->data == NULL) error_msg("Failed to allocate data for nd_malloc array") ;
}




/* Free ND_function */
void ND_function(free, TYPE_S) (ND_array(TYPE_S) * nd_arr_in)
{
    /* 
    Must be called when nd_malloc_? or nd_calloc_? ND_functions are called, else memory leak !
    */
    if (nd_arr_in->owner)
    {   
        if (nd_arr_in->data == NULL) error_msg("Trying to free array with NULL data");
        free(nd_arr_in->data);
        nd_arr_in->data = NULL ; 
        nd_arr_in->owner = false;
    }
}

/* Destroy array */
void ND_function(destroy, TYPE_S) (ND_array(TYPE_S) * nd_arr_in)
{
    /* simply free + uninit */
    ND_function(free, TYPE_S) (nd_arr_in);
    ND_function(uninit, TYPE_S) (nd_arr_in);
}




