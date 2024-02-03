#include "nd_array_src.h"

/* Array operation ND_function */
/****************************************************************************************************/

/* ND_function to get element pointer of an array */
TYPE_L * ND_function(ele, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, const ND_int * dimensions)
{
    /*returns the pointer to the particular element . so we can set and get the elements*/
    //
    if (nd_arr_in->rank == NULL) 
        error_msg("uninitialized input array for nd_ele ND_function");
    if (nd_arr_in->data == NULL) 
        error_msg("NULL data array received for nd_ele ND_function");

    ND_int index_arr = 0, arg_val; 

    for(ND_int i = 0; i < *(nd_arr_in->rank); ++i)
    {
        arg_val = dimensions[i];

        if (arg_val < ((nd_arr_in->dims)[i]) ) index_arr = index_arr + arg_val * ((nd_arr_in->strides)[i]);

        else 
            error_msg("Array out of bound");
    }

    return (nd_arr_in->data) + index_arr;
}





/* ND_function to get Size of an array .*/
ND_int ND_function(size, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in)
{   
    if (nd_arr_in->rank == NULL) 
        error_msg("uninitialized input array for nd_size ND_function");

    ND_int size ; 

    if (*(nd_arr_in->rank) ==  0) size =   1 ;
    
    else size = *(nd_arr_in->dims) * *(nd_arr_in->strides);
    

    return size; 
}


/* ND_function to set all  elements of an array to constant*/
void ND_function(set_all, TYPE_S) (ND_array(TYPE_S) * nd_arr_in, const TYPE_L set_constant )
{   
    if (nd_arr_in->rank == NULL) 
        error_msg("uninitialized input array for nd_ele ND_function");
    if (nd_arr_in->data == NULL) 
        error_msg("NULL data array received for nd_ele ND_function");

    const ND_int arr_in_size = ND_function(size, TYPE_S)(nd_arr_in);
    for (ND_int i = 0 ; i < arr_in_size ; ++i ) 
    {
        (nd_arr_in->data)[i] = set_constant ;
    }

}


/* ND_function to deepcopy from one array to other*/
void ND_function(copy, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, ND_array(TYPE_S) * nd_arr_out)
{

    if (nd_arr_in->rank == NULL) 
        error_msg("uninitialized input array for copy ND_function");
    if (nd_arr_in->data == NULL) 
        error_msg("NULL data array received for copy ND_function");

    if (nd_arr_out->rank == NULL) 
        error_msg("uninitialized output array for copy ND_function");
    if (nd_arr_out->data == NULL) 
        error_msg("NULL data array received for copy ND_function");

    if (ND_function(size, TYPE_S) (nd_arr_in) != ND_function(size, TYPE_S) (nd_arr_out)) 
        error_msg("Incompatible size between input and copy array");

    if (*(nd_arr_in->rank) != *(nd_arr_out->rank)) 
        error_msg("Incompatible rank between input and copy array");


    /* check dims*/
    for ( ND_int i = 0 ; i < *(nd_arr_in->rank); ++i)
    {
        if (  (nd_arr_in->dims)[i] != (nd_arr_out->dims)[i] ) 
            error_msg("Incompatible dims for input and copy arrays.");

    }
    memcpy(nd_arr_out->data, nd_arr_in->data, sizeof(TYPE_L)* ND_function(size, TYPE_S) (nd_arr_in));
    memcpy(nd_arr_out->rank, nd_arr_in->rank, sizeof(ND_int) * ( (1 + 2* (nd_arr_in->rank)[0] ) ));
}



