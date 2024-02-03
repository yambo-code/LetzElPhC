#include "nd_array_src.h"


/* netCDF4 I/O ND_function */
/****************************************************************************************************/


/* ND_function to read from netCDF */
void ND_function(read, TYPE_S) (const char* file_name, const char* var_name, ND_array(TYPE_S) * nd_arr_in)
{   
    /* Note this ND_function create the ND_array(TYPE_S), so ND_function(free, TYPE_S) () must be called to free the memory
        DO NOT pass a pointer which already has data. this leads to memory leak
        Ex: ND_function(read, TYPE_S) ("ndb.BS_elph", "exc_elph", &temp_array);
    */
    //
    if (nd_arr_in->data != NULL) error_msg("Input array is uninitialized or has data."); // check if there is data or uninitialized

    if (nd_arr_in->rank != NULL) free(nd_arr_in->rank); // if not null, free the memeory

    //
    int ncid, var_id, retval, nc_rank; // variables 

    if ((retval = nc_open(file_name, NC_NOWRITE, &ncid))) ERR(retval); // open file

    if ((retval = nc_inq_varid(ncid, var_name, &var_id))) ERR(retval); // get the id of the req variable

    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, &nc_rank, NULL, NULL ))) ERR(retval); // get rank

    int * dim_ids             = malloc(nc_rank*sizeof(int));
    ND_int * nd_dims      = malloc(nc_rank*sizeof(ND_int));
    size_t * nd_dims_temp     = malloc(nc_rank*sizeof(size_t));

    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, NULL, dim_ids, NULL ))) ERR(retval); // get dims
    //
    for (ND_int i = 0; i < ( (ND_int)nc_rank); ++i)
    {
        if ((retval = nc_inq_dimlen(ncid, dim_ids[i], nd_dims_temp + i))) ERR(retval);
        nd_dims[i] =  nd_dims_temp[i] ;
    }

    #if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)

        if ((nc_rank < 1) || nd_dims[ (nc_rank-1)] != 2 ) error_msg("Cannot convert a real to complex array.") ; 

        ND_function(init, TYPE_S) (nd_arr_in, ( (nc_rank-1)), nd_dims); 

    #else

        ND_function(init, TYPE_S) (nd_arr_in, ( nc_rank), nd_dims); 

    #endif

    ND_function(malloc, TYPE_S) (nd_arr_in);  // this must be free outside else memory leak

    if ((retval = nc_get_var(ncid, var_id, nd_arr_in->data))) ERR(retval); //get data in floats

    if ((retval = nc_close(ncid))) ERR(retval); // close the file

    // free all temp arrays
    free(dim_ids);
    free(nd_dims);
    free(nd_dims_temp);
}

/* Read only subsection of array*/
void ND_function(read_sub, TYPE_S) (const char* file_name, const char* var_name, ND_array(TYPE_S) * nd_arr_in, ...)
{   
    /* Note this ND_function create the ND_array(TYPE_S), so ND_function(free, TYPE_S) () must be called to free the memory
        DO NOT pass a pointer which already has data. this leads to memory leak
        Ex: ND_function(read, TYPE_S) ("ndb.BS_elph", "exc_elph", &temp_array, ...);
    */
    //

    if (nd_arr_in->data != NULL) error_msg("Input array is uninitialized or has data."); // check if there is data or uninitialized

    if (nd_arr_in->rank != NULL) free(nd_arr_in->rank); // if not null, free the memeory


    //
    int ncid, var_id, retval, nc_rank; // variables 

    if ((retval = nc_open(file_name, NC_NOWRITE, &ncid))) ERR(retval); // open file

    if ((retval = nc_inq_varid(ncid, var_name, &var_id))) ERR(retval); // get the id of the req variable

    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, &nc_rank, NULL, NULL ))) ERR(retval); // get rank



    int * dim_ids             = malloc(nc_rank*sizeof(int));
    ND_int * nd_dims      = malloc(nc_rank*sizeof(ND_int));
    size_t * nd_dims_temp     = malloc(nc_rank*sizeof(size_t));

    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, NULL, dim_ids, NULL ))) ERR(retval); // get dims
    //
    for (ND_int i = 0; i < ( (ND_int)nc_rank); ++i)
    {
        if ((retval = nc_inq_dimlen(ncid, dim_ids[i], nd_dims_temp + i))) ERR(retval);

        nd_dims[i] =  nd_dims_temp[i] ;
    }


    size_t * startp = malloc( 2 * nc_rank * sizeof(size_t) );

    size_t * countp = startp + nc_rank ;
    
    ptrdiff_t * stridep = malloc(nc_rank * sizeof(ptrdiff_t));

    ND_int * temp_ptr ; 

    /*** va_args section **/

    va_list valist;

    va_start(valist, nd_arr_in);

    #if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)

    if ((nc_rank < 1) || nd_dims_temp[nc_rank-1] != 2 ) error_msg("Cannot convert a real to complex array. The dimension of last must be 2 ") ; 

    for (size_t i = 0; i < (size_t) nc_rank - 1; i++)

    #else
    
    for (size_t i = 0; i < (size_t) nc_rank; i++)
    
    #endif
    {   /* [0] -> start, [1] -> end, [2] - step*/

        temp_ptr = va_arg(valist,  ND_int * );

        if (temp_ptr == ND_ALL)
        {
            startp[i] = 0 ; 

            countp[i] = nd_dims_temp[i] ; 

            stridep[i] = 1 ;
        }
        else
        {
            startp[i] = temp_ptr[0];

            if (temp_ptr[2] == 0) error_msg("NetCDF read error - Strides must be postive integers.") ; 

            stridep[i] = temp_ptr[2] ;

            if ((size_t) temp_ptr[1] > nd_dims_temp[i] ) error_msg("NetCDF read error - End indices are out of bound.") ; 

            if (temp_ptr[1] == 0 ) countp[i]  = 1 + (nd_dims_temp[i] - temp_ptr[0] - 1)/(stridep[i]) ;
            
            else
            {   
                if (temp_ptr[0] >= temp_ptr[1] ) error_msg("NetCDF read error - start indices must be \
                                    less than end indices (expect when end indices are not zeros).") ;

                countp[i]  =  1 + (temp_ptr[1] - temp_ptr[0] - 1)/(stridep[i])  ; 
            }
            
        }

        nd_dims[i] =  countp[i] ; 


    }
	
    va_end(valist);

    /*** end va_args **/

    #if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)

        startp[nc_rank-1] =  0 ; 

        countp[nc_rank-1] =  2 ; 

        stridep[nc_rank-1]=  1 ;

        ND_function(init, TYPE_S) (nd_arr_in, ( (nc_rank-1)), nd_dims); 
        
    #else

        ND_function(init, TYPE_S) (nd_arr_in, ( nc_rank), nd_dims); 

    #endif

    ND_function(malloc, TYPE_S) (nd_arr_in);  // this must be free outside else memory leak

    if ((retval = nc_get_vars(ncid, var_id,  startp, countp, stridep, nd_arr_in->data))) ERR(retval); //get data in floats

    if ((retval = nc_close(ncid))) ERR(retval); // close the file

    free(startp);
    free(stridep);
    free(dim_ids);
    free(nd_dims);
    free(nd_dims_temp);

}




/* ND_function to write to netCDF */
void ND_function(write, TYPE_S) (const char* file_name, const char* var_name, const ND_array(TYPE_S) * nd_arr_in, char ** dim_names, size_t * chunksize)
{
    /* Dumps the ND_array(TYPE_S) to file with filemane (file_name) and dimenstion name (dim_names)
        varible names
        Ex: ND_function(write, TYPE_S) ("nc.temp", "elph_ex", &temp_array, (char * [5]) {"nq", "modes", "Sf", "Si", "re_im"});
    */
    if (nd_arr_in->rank == NULL) error_msg("Cannot write a uninitialized array");

    if (nd_arr_in->data == NULL) error_msg("Cannot write a empty array");

    int retval; // error iD 

    int ncid, varid; // var iDs

    if ((retval = nc_create(file_name, NC_CLOBBER | NC_NETCDF4 , &ncid))) ERR(retval); // create file

    #if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)

    int * dimids = malloc((*(nd_arr_in->rank) + 1) * sizeof(int));
    //
    for (ND_int i =0 ; i< (*(nd_arr_in->rank)); ++i )
    {
        if ((retval = nc_def_dim(ncid, dim_names[i],  (nd_arr_in->dims)[i], dimids+i))) ERR(retval);   // get ids for the omega dimensions
    }

    if ((retval = nc_def_dim(ncid, "re_im",  2, dimids+ (*(nd_arr_in->rank))  ))) ERR(retval);   // get ids for the omega dimensions
    //

    if ((retval = nc_def_var(ncid, var_name, NetCDF_IO_TYPE, (int) (1 + *(nd_arr_in->rank)), dimids, &varid))) ERR(retval); // Writing the data with name Raman_tensor
    
    #else
    
    int * dimids = malloc((*(nd_arr_in->rank)) * sizeof(int));
    //
    for (ND_int i =0 ; i< (*(nd_arr_in->rank)); ++i )
    {
        if ((retval = nc_def_dim(ncid, dim_names[i],  (nd_arr_in->dims)[i], dimids+i))) ERR(retval);   // get ids for the omega dimensions
    }
    //
    if ((retval = nc_def_var(ncid, var_name, NetCDF_IO_TYPE, (int) (*(nd_arr_in->rank)), dimids, &varid))) ERR(retval); // Writing the data with name Raman_tensor
    
    #endif

    // Set chunking if there
    if (chunksize == NULL)
    {
        if ((retval = nc_def_var_chunking(ncid, varid, NC_CONTIGUOUS, NULL))) ERR(retval);
    }
    else
    {
        if ((retval = nc_def_var_chunking(ncid, varid, NC_CHUNKED, chunksize))) ERR(retval);
    }
    

    if ((retval = nc_enddef(ncid))) ERR(retval);

    if ((retval =  nc_put_var(ncid, varid, nd_arr_in->data))) ERR(retval);

    if ((retval = nc_close(ncid))) ERR(retval);

    free(dimids);

    
} 



#if defined(COMPILE_ONCE)
void NC_open_file(const char* file_name, char mode, int * ncid)
{   
    int retval ;

    if (mode == 'r')       retval = nc_open  (file_name, NC_NOWRITE,                ncid) ;
    else if (mode == 'a')  retval = nc_open  (file_name, NC_WRITE,                  ncid) ;
    else if (mode == 'w')  retval = nc_create(file_name, NC_NOCLOBBER | NC_NETCDF4, ncid) ;
    else if (mode == 'W')  retval = nc_create(file_name, NC_CLOBBER   | NC_NETCDF4, ncid) ;

    else {
        printf("Only mode = 'w' (write); 'W' (write with overwriting the file); \
                'a' (append); 'r' (read) are supported \n");
        exit(EXIT_FAILURE);
    }

    if (retval != NC_NOERR)
    {   
        printf("Error opening/creating the netCDF file %s \n", file_name);
        ERR(retval); 
    }

}


void NC_close_file(const int ncid)
{   
    int retval ;
    retval = nc_close(ncid);
    if (retval != NC_NOERR)
    {   
        printf("Error closing the netCDF file. ncid do not correspond to any open file \n");
        ERR(retval); 
    }
}
#endif

/* ND_function to read from netCDF */
void ND_function(readVar, TYPE_S) (const int ncid, const char* var_name, ND_array(TYPE_S) * nd_arr_in)
{   
    /* Note this ND_function create the ND_array(TYPE_S), so ND_function(free, TYPE_S) () must be called to free the memory
        DO NOT pass a pointer which already has data. this leads to memory leak
        Ex: ND_function(read, TYPE_S) ("ndb.BS_elph", "exc_elph", &temp_array);
    */
    //
    if (nd_arr_in->data != NULL) error_msg("Input array is uninitialized or has data."); // check if there is data or uninitialized

    if (nd_arr_in->rank != NULL) free(nd_arr_in->rank); // if not null, free the memeory

    //
    int var_id, retval, nc_rank; // variables 

    if ((retval = nc_inq(ncid, NULL, NULL, NULL, NULL))) ERR(retval); // check the ncid

    if ((retval = nc_inq_varid(ncid, var_name, &var_id))) ERR(retval); // get the id of the req variable

    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, &nc_rank, NULL, NULL ))) ERR(retval); // get rank

    int * dim_ids             = malloc(nc_rank*sizeof(int));
    ND_int * nd_dims      = malloc(nc_rank*sizeof(ND_int));
    size_t * nd_dims_temp     = malloc(nc_rank*sizeof(size_t));

    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, NULL, dim_ids, NULL ))) ERR(retval); // get dims
    //
    for (ND_int i = 0; i < ( (ND_int)nc_rank); ++i)
    {
        if ((retval = nc_inq_dimlen(ncid, dim_ids[i], nd_dims_temp + i))) ERR(retval);
        nd_dims[i] =  nd_dims_temp[i] ;
    }

    #if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)

        if ((nc_rank < 1) || nd_dims[ (nc_rank-1)] != 2 ) error_msg("Cannot convert a real to complex array.") ; 

        ND_function(init, TYPE_S) (nd_arr_in, ( (nc_rank-1)), nd_dims); 

    #else

        ND_function(init, TYPE_S) (nd_arr_in, ( nc_rank), nd_dims); 

    #endif

    ND_function(malloc, TYPE_S) (nd_arr_in);  // this must be free outside else memory leak

    if ((retval = nc_get_var(ncid, var_id, nd_arr_in->data))) ERR(retval); //get data in floats

    // free all temp arrays
    free(dim_ids);
    free(nd_dims);
    free(nd_dims_temp);
}

/* Read only subsection of array*/
void ND_function(readVar_sub, TYPE_S) (const int ncid, const char* var_name, ND_array(TYPE_S) * nd_arr_in, ...)
{   
    /* Note this ND_function create the ND_array(TYPE_S), so ND_function(free, TYPE_S) () must be called to free the memory
        DO NOT pass a pointer which already has data. this leads to memory leak
        Ex: ND_function(read, TYPE_S) ("ndb.BS_elph", "exc_elph", &temp_array, ...);
    */
    //

    if (nd_arr_in->data != NULL) error_msg("Input array is uninitialized or has data."); // check if there is data or uninitialized

    if (nd_arr_in->rank != NULL) free(nd_arr_in->rank); // if not null, free the memeory


    //
    int var_id, retval, nc_rank; // variables 

    if ((retval = nc_inq(ncid, NULL, NULL, NULL, NULL))) ERR(retval); // check the ncid

    if ((retval = nc_inq_varid(ncid, var_name, &var_id))) ERR(retval); // get the id of the req variable

    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, &nc_rank, NULL, NULL ))) ERR(retval); // get rank



    int * dim_ids             = malloc(nc_rank*sizeof(int));
    ND_int * nd_dims          = malloc(nc_rank*sizeof(ND_int));
    size_t * nd_dims_temp     = malloc(nc_rank*sizeof(size_t));

    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, NULL, dim_ids, NULL ))) ERR(retval); // get dims
    //
    for (ND_int i = 0; i < ( (ND_int)nc_rank); ++i)
    {
        if ((retval = nc_inq_dimlen(ncid, dim_ids[i], nd_dims_temp + i))) ERR(retval);

        nd_dims[i] =  nd_dims_temp[i] ;
    }


    size_t * startp = malloc( 2 * nc_rank * sizeof(size_t) );

    size_t * countp = startp + nc_rank ;
    
    ptrdiff_t * stridep = malloc(nc_rank * sizeof(ptrdiff_t));

    ND_int * temp_ptr ; 

    /*** va_args section **/

    va_list valist;

    va_start(valist, nd_arr_in);

    #if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)

    if ((nc_rank < 1) || nd_dims_temp[nc_rank-1] != 2 ) error_msg("Cannot convert a real to complex array. The dimension of last must be 2 ") ; 

    for (size_t i = 0; i < (size_t) nc_rank - 1; i++)

    #else
    
    for (size_t i = 0; i < (size_t)nc_rank; i++)
    
    #endif
    {   /* [0] -> start, [1] -> end, [2] - step*/

        temp_ptr = va_arg(valist,  ND_int * );

        if (temp_ptr == ND_ALL)
        {
            startp[i] = 0 ; 

            countp[i] = nd_dims_temp[i] ; 

            stridep[i] = 1 ;
        }
        else
        {
            startp[i] = temp_ptr[0];

            if (temp_ptr[2] == 0) error_msg("NetCDF read error - Strides must be postive integers.") ; 

            stridep[i] = temp_ptr[2] ;

            if ((size_t) temp_ptr[1] > nd_dims_temp[i] ) error_msg("NetCDF read error - End indices are out of bound.") ; 

            if (temp_ptr[1] == 0 ) countp[i]  = 1 + (nd_dims_temp[i] - temp_ptr[0] - 1)/(stridep[i]) ;
            
            else
            {   
                if (temp_ptr[0] >= temp_ptr[1] ) error_msg("NetCDF read error - start indices must be \
                                    less than end indices (expect when end indices are not zeros).") ;

                countp[i]  =  1 + (temp_ptr[1] - temp_ptr[0] - 1)/(stridep[i])  ; 
            }
            
        }

        nd_dims[i] =  countp[i] ; 


    }
	
    va_end(valist);

    /*** end va_args **/

    #if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)

        startp[nc_rank-1] =  0 ; 

        countp[nc_rank-1] =  2 ; 

        stridep[nc_rank-1]=  1 ;

        ND_function(init, TYPE_S) (nd_arr_in, ( (nc_rank-1)), nd_dims); 
        
    #else

        ND_function(init, TYPE_S) (nd_arr_in, ( nc_rank), nd_dims); 

    #endif

    ND_function(malloc, TYPE_S) (nd_arr_in);  // this must be free outside else memory leak

    bool use_vara = true ;

    for (size_t i = 0; i < (size_t) nc_rank ; i++)
    {
        if (stridep[i] != 1) 
        {
            use_vara = false ;
            break ;
        }
    }
    if (use_vara)
    {   
        if ((retval = nc_get_vara(ncid, var_id,  startp, countp, nd_arr_in->data))) ERR(retval); //get data in floats
    }
    else 
    {   
        if ((retval = nc_get_vars(ncid, var_id,  startp, countp, stridep, nd_arr_in->data))) ERR(retval); //get data in floats
    }
    free(startp);
    free(stridep);
    free(dim_ids);
    free(nd_dims);
    free(nd_dims_temp);

}




/* ND_function to write to netCDF */
void ND_function(writeVar, TYPE_S) (const int ncid, const char* var_name, const ND_array(TYPE_S) * nd_arr_in, char ** dim_names, size_t * chunksize)
{
    /* Dumps the ND_array(TYPE_S) to file with filemane (file_name) and dimenstion name (dim_names)
        varible names
        Ex: ND_function(write, TYPE_S) ("nc.temp", "elph_ex", &temp_array, (char * [5]) {"nq", "modes", "Sf", "Si", "re_im"});
    */
    if (nd_arr_in->rank == NULL) error_msg("Cannot write a uninitialized array");

    if (nd_arr_in->data == NULL) error_msg("Cannot write a empty array");

    int retval; // error iD 

    int varid; // var iDs

    if ((retval = nc_inq(ncid, NULL, NULL, NULL, NULL))) ERR(retval); // check the ncid

    #if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)

    int * dimids = malloc((*(nd_arr_in->rank) + 1) * sizeof(int));
    //
    for (ND_int i =0 ; i< (*(nd_arr_in->rank)); ++i )
    {   
        if ((retval = nc_inq_dimid(ncid, dim_names[i], dimids+i)))
        {
            if ((retval = nc_def_dim(ncid, dim_names[i],  (nd_arr_in->dims)[i], dimids+i))) ERR(retval);   // get ids for the omega dimensions
        }
        /* Check is dim len are consistant */
        size_t temp_dim ;
        if ((retval = nc_inq_dimlen(ncid,dimids[i],&temp_dim))) ERR(retval);
        if ((ND_int)temp_dim != (nd_arr_in->dims)[i]) error_msg("Dimension name already exist with different dimension length");
    }
    if ((retval = nc_inq_dimid(ncid, "re_im", dimids+ (*(nd_arr_in->rank)) )))
    {
        if ((retval = nc_def_dim(ncid, "re_im",  2, dimids+ (*(nd_arr_in->rank))  ))) ERR(retval);   // get ids for the omega dimensions
    }//

    if ((retval = nc_def_var(ncid, var_name, NetCDF_IO_TYPE, (int) (1 + *(nd_arr_in->rank)), dimids, &varid))) ERR(retval); // Writing the data with name Raman_tensor
    
    #else
    
    int * dimids = malloc((*(nd_arr_in->rank)) * sizeof(int));
    //
    for (ND_int i =0 ; i< (*(nd_arr_in->rank)); ++i )
    {   
        if ((retval = nc_inq_dimid(ncid, dim_names[i], dimids+i)))
        {
            if ((retval = nc_def_dim(ncid, dim_names[i],  (nd_arr_in->dims)[i], dimids+i))) ERR(retval);   // get ids for the omega dimensions
        }
        /* Check is dim len are consistant */
        size_t temp_dim ;
        if ((retval = nc_inq_dimlen(ncid,dimids[i],&temp_dim))) ERR(retval);
        if ((ND_int)temp_dim != (nd_arr_in->dims)[i]) error_msg("Dimension name already exist with different dimension length");
    }
    //
    if ((retval = nc_def_var(ncid, var_name, NetCDF_IO_TYPE, (int) (*(nd_arr_in->rank)), dimids, &varid))) ERR(retval); // Writing the data with name Raman_tensor
    
    #endif

    // Set chunking if there
    if (chunksize == NULL)
    {
        if ((retval = nc_def_var_chunking(ncid, varid, NC_CONTIGUOUS, NULL))) ERR(retval);
    }
    else
    {
        if ((retval = nc_def_var_chunking(ncid, varid, NC_CHUNKED, chunksize))) ERR(retval);
    }
    

    if ((retval = nc_enddef(ncid))) ERR(retval);

    if ((retval =  nc_put_var(ncid, varid, nd_arr_in->data))) ERR(retval);

    free(dimids);

    
} 



