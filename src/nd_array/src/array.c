#include "nd_array_src.h"

static void ND_function(slice_internal, TYPE_S) ( const ND_int * restrict start_idx, const ND_int * restrict end_idx, const ND_int * restrict step_idx, \
        const ND_int * restrict stride_F, const ND_int * restrict stride_S, ND_int idx_F,  ND_int idx_S, const ND_int ndim, \
        const ND_int idim, const TYPE_L * restrict arrF, TYPE_L * restrict arrS);




/* Array operation ND_function */
/****************************************************************************************************/

/* ND_function to get element pointer of an array */
TYPE_L * ND_function(ele, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, const ND_int * dimensions)
{
    /*returns the pointer to the particular element . so we can set and get the elements*/
    //
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for nd_ele ND_function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for nd_ele ND_function");

    ND_int index_arr = 0, arg_val; 

    for(ND_int i = 0; i < *(nd_arr_in->rank); ++i)
    {
        arg_val = dimensions[i];

        if (arg_val < ((nd_arr_in->dims)[i]) ) index_arr = index_arr + arg_val * ((nd_arr_in->strides)[i]);

        else error_msg("Array out of bound");
    }

    return (nd_arr_in->data) + index_arr;
}





/* ND_function to get Size of an array .*/
ND_int ND_function(size, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in)
{   
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for nd_size ND_function");

    ND_int size ; 

    if (*(nd_arr_in->rank) ==  0) size =   1 ;
    
    else size = *(nd_arr_in->dims) * *(nd_arr_in->strides);
    

    return size; 
}


/* ND_function to set all  elements of an array to constant*/
void ND_function(set_all, TYPE_S) (ND_array(TYPE_S) * nd_arr_in, const TYPE_L set_constant )
{   
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for nd_ele ND_function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for nd_ele ND_function");

    const ND_int arr_in_size = ND_function(size, TYPE_S)(nd_arr_in);
    for (ND_int i = 0 ; i < arr_in_size ; ++i ) 
    {
        (nd_arr_in->data)[i] = set_constant ;
    }

}

/* ND_function to get reshape an array .*/

void ND_function(reshape, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, ND_array(TYPE_S) * nd_arr_out)
{

    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for reshape ND_function");
    if (nd_arr_out->rank == NULL) error_msg("uninitialized output array for reshape ND_function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for reshape ND_function");
    if (nd_arr_out->data != NULL && nd_arr_out->owner) \
        error_msg("Cannot pass array with data to reshape ND_function. This leads to memory leak");

    //nd_arr_out->owner = false ;

    /*check size*/
    if ( ND_function(size, TYPE_S) (nd_arr_in) !=  ND_function(size, TYPE_S) (nd_arr_out) ) \
                                    error_msg("Size mismatch of reshaped array and original array");

    nd_arr_out->data = nd_arr_in->data;

}


/* ND_function to strip off the first n dimensions */
void ND_function(strip_dims, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, const ND_int n_dims_strip, const ND_int * stripped_idxs, ND_array(TYPE_S) * nd_arr_out)
{
    /* 
        Strips first n dims: 
        Ex: 
        for 5 dim array A, when we strip first two dims we get A[i,j,:,:,:]. n_dims_strip=2, stripped_idxs={i,j}
     */
    
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for strip_dims ND_function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for strip_dims ND_function");

    if (nd_arr_out->rank == NULL) error_msg("uninitialized output array for strip_dims ND_function");
    if (nd_arr_out->data != NULL && nd_arr_out->owner) \
        error_msg("Cannot pass array with data to strip_dims ND_function. This leads to memory leak");

    //nd_arr_out->owner = false ;

    if (n_dims_strip >= (nd_arr_in->rank)[0]) error_msg("Number of dims to stip must be less than the rank of original array");

    ND_int sub_tensor_rank = (nd_arr_in->rank)[0]-n_dims_strip ;
    
    if ( *(nd_arr_out->rank) != sub_tensor_rank) error_msg("Rank mismatch in nd_strip_dims ND_function") ;

    for (ND_int i = 1; i <= sub_tensor_rank; ++i)
    {
        if ((nd_arr_out->dims)[sub_tensor_rank - i] != (nd_arr_in->dims)[(nd_arr_in->rank)[0] - i]) error_msg("Dimension mismatch in nd_strip_dims ND_function") ; 
    }

    ND_int index_arr = 0, arg_val; 
    
    for(ND_int i = 0; i < n_dims_strip; ++i)
    {
        arg_val = stripped_idxs[i];

        if (arg_val < ((nd_arr_in->dims)[i]) )
        {
            index_arr = index_arr + arg_val * ((nd_arr_in->strides)[i]);
        }
        else error_msg("Array out of bound");
    }

    nd_arr_out->data = nd_arr_in->data + index_arr;
}


/* Availble for user :) */
/* ND_function to slice a subtensor and copy to new tensor*/
void ND_function(slice, TYPE_S) (const ND_int * start_idx, const ND_int * end_idx, const ND_int * step_idx, \
                            const ND_array(TYPE_S) * nd_arr_in, ND_array(TYPE_S) * nd_arr_out )
{
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for slice ND_function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for slice ND_function");

    if (nd_arr_out->rank == NULL) error_msg("uninitialized output array for slice ND_function");
    if (nd_arr_out->data == NULL) error_msg("NULL data array received for slice ND_function");

    ND_int slice_rank =  0;
    ND_int arr_in_rank = *(nd_arr_in->rank);

    ND_int * sliced_dims_temp = malloc(3*arr_in_rank * sizeof(ND_int)); // create both temp and sliced and strides
    ND_int * sliced_dims = sliced_dims_temp + arr_in_rank ; 
    ND_int * out_strides = sliced_dims_temp + 2*arr_in_rank ;
    
    for (ND_int i =0 ; i < arr_in_rank; i++ )
    {

        if ((end_idx[i] <= start_idx[i] || end_idx[i] > (nd_arr_in->dims)[i]) || start_idx[i] >= (nd_arr_in->dims)[i] ) \
                error_msg("Slicing index out of bound or ending index is larger than the starting index when slicing array.");

        else
        {
            sliced_dims_temp[i] = 1 + (end_idx[i]-start_idx[i] - 1)/step_idx[i];

            if (sliced_dims_temp[i] > 1)
            {
                sliced_dims[slice_rank] = sliced_dims_temp[i];
                slice_rank++;
            }
        }
    }

    ND_int slice_size = 1;

    for (ND_int i = 1; i < arr_in_rank+1; ++i)
    {
        out_strides[arr_in_rank - i] = slice_size;
        slice_size = slice_size*sliced_dims_temp[arr_in_rank-i];
    }

    if (slice_rank == 0 || *(nd_arr_in->rank) == 0 ) error_msg("Ranks of sliced array or input array cannot be zero. Use nd_ele ND_function");

    if (slice_size !=  ND_function(size, TYPE_S) (nd_arr_out))
    {
        fprintf(stdout, "# [ Error !!!] : Sliced array size is "PRI_NdInt" but the given sliced array size is "PRI_NdInt"  \n",\
            slice_size, ND_function(size, TYPE_S) (nd_arr_out) );
        error_msg("Incompatible sized in nd_slice ND_function");
    }

    if (slice_rank !=  *(nd_arr_out->rank) ) error_msg("Ranks of inputed sliced array are wrong");

    for (ND_int i = 0; i < slice_rank; ++i)
    {
        if (sliced_dims[i] != (nd_arr_out->dims)[i] ) error_msg("Dimensions of inputed sliced array are wrong");
    }
    
    ND_function(slice_internal, TYPE_S) (start_idx, end_idx,step_idx,nd_arr_in->strides,out_strides,  0,\
       0, *(nd_arr_in->rank),  0, nd_arr_in->data, nd_arr_out->data);

    free(sliced_dims_temp);
}



// void ND_function(swap_owner, TYPE_S) ( ND_array(TYPE_S) * nd_arr_new, ND_array(TYPE_S) * nd_arr_old )
//{
//     nd_arr_new->owner = true;
//     nd_arr_old->owner = false;
// }

/* ND_function to transpose copy from one array to other. Same as like numpy.transpose()*/
void ND_function(tranpose, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, const ND_int * order, ND_array(TYPE_S) * nd_arr_out)
{
    /*Note that you must input and out put array. No array is created. it just sets elements.
    This is not very cache friendly way !
    */
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for transpose ND_function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for transpose ND_function");

    if (nd_arr_out->rank == NULL) error_msg("uninitialized output array for transpose ND_function");
    if (nd_arr_out->data == NULL) error_msg("NULL data array received for transpose ND_function");

    if (ND_function(size, TYPE_S) (nd_arr_in) != ND_function(size, TYPE_S) (nd_arr_out)) \
                                error_msg("Incompatible size between input and transposed array"); //check dim

    if (*(nd_arr_in->rank) != *(nd_arr_out->rank)) error_msg("Incompatible size between input and transposed array"); // check rank

    // check if all transpose indices are in range
    for ( ND_int i = 0 ; i < *(nd_arr_in->rank); ++i)
    {
        if ( order[i] >= *(nd_arr_in->rank) ) error_msg("Transpose indices must be less than rank.");
    }

    /* check dims*/
    for ( ND_int i = 0 ; i < *(nd_arr_in->rank); ++i)
    {
        if (  (nd_arr_in->dims)[order[i]] != (nd_arr_out->dims)[i] ) error_msg("Incompatible transpose order.");
    }

    ND_int rank = *(nd_arr_in->rank );
    ND_int * idx_array = malloc(( 1 + rank ) * sizeof(ND_int)); /*last value of index is remainder*/

    /* set the data*/
    ND_int out_idx ;

    //
    for (ND_int i = 0; i < ND_function(size, TYPE_S) (nd_arr_in); ++i)
    {
        idx_array[rank] = i;
        out_idx = 0;

        for (ND_int j = 0; j < rank; ++j)
        {
            idx_array[j] = idx_array[rank]/(nd_arr_in->strides)[j] ;
            idx_array[rank] = idx_array[rank] % (nd_arr_in->strides)[j] ;
        }
        
        for (ND_int j = 0; j < rank; ++j)
        {
            out_idx = out_idx + idx_array[order[j]] * (nd_arr_out->strides)[j];
        }

        (nd_arr_out->data)[out_idx] = (nd_arr_in->data)[i] ; 
    }

    free(idx_array);

}




/* ND_function to deepcopy from one array to other*/
void ND_function(copy, TYPE_S) (const ND_array(TYPE_S) * nd_arr_in, ND_array(TYPE_S) * nd_arr_out)
{

    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for copy ND_function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for copy ND_function");

    if (nd_arr_out->rank == NULL) error_msg("uninitialized output array for copy ND_function");
    if (nd_arr_out->data == NULL) error_msg("NULL data array received for copy ND_function");

    if (ND_function(size, TYPE_S) (nd_arr_in) != ND_function(size, TYPE_S) (nd_arr_out)) error_msg("Incompatible size between input and copy array");

    if (*(nd_arr_in->rank) != *(nd_arr_out->rank)) error_msg("Incompatible rank between input and copy array");


    /* check dims*/
    for ( ND_int i = 0 ; i < *(nd_arr_in->rank); ++i)
    {
        if (  (nd_arr_in->dims)[i] != (nd_arr_out->dims)[i] ) error_msg("Incompatible dims for input and copy arrays.");

    }
    memcpy(nd_arr_out->data, nd_arr_in->data, sizeof(TYPE_L)* ND_function(size, TYPE_S) (nd_arr_in));
    memcpy(nd_arr_out->rank, nd_arr_in->rank, sizeof(ND_int) * ( (1 + 2* (nd_arr_in->rank)[0] ) ));
}


/************************************************************************ STATIC ND_functionS **********************************************************************/

/* Not availble for user */

static void ND_function(slice_internal, TYPE_S) ( const ND_int * restrict start_idx, const ND_int * restrict end_idx, const ND_int * restrict step_idx, \
        const ND_int * restrict stride_F, const ND_int * restrict stride_S, ND_int idx_F,  ND_int idx_S, const ND_int ndim, \
        const ND_int idim, const TYPE_L * restrict arrF, TYPE_L * restrict arrS)
{
    
    if (idim == ndim) arrS[idx_S] = arrF[idx_F] ;
    //
    else if ( (idim == ndim-1) && (step_idx[idim] ==  1) )
    {
        memcpy(arrS+idx_S, arrF + idx_F + (start_idx[idim] * stride_F[idim]), (end_idx[idim]-start_idx[idim])*sizeof(TYPE_L));
    }
    //
    else 
    {
        ND_int idx_F1, idx_S1;

        for (ND_int i = start_idx[idim] ; i < end_idx[idim] ; i = i + step_idx[idim] ){
            if (idim ==  0)
            {
                idx_F =  0  ; 
                idx_S =  0  ; 
            }
            idx_F1 = idx_F + (i * stride_F[idim]) ;
            idx_S1 = idx_S + (((i - start_idx[idim] )/step_idx[idim]) * stride_S[idim]) ;
            ND_function(slice_internal, TYPE_S)(start_idx, end_idx, step_idx, stride_F, stride_S, idx_F1, idx_S1, ndim, idim +  1, arrF, arrS) ; 
        }
    }
}