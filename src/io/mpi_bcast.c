#include "io.h"
#include <limits.h>


// float
void Bcast_ND_arrayFloat(ND_array(Nd_floatS) * array_in, bool alloc_mem, int root, MPI_Comm comm)
{
    if (alloc_mem)
    {
        ND_int dims[24];
        ND_int rank;

        ND_int * dim_ptr;
        ND_int * rank_ptr;

        int my_rank, comm_size, mpi_error;

        MPI_Comm_rank(comm, &my_rank);
        MPI_Comm_size(comm, &comm_size);

        if (my_rank == root)
        {
            dim_ptr =  array_in->dims;
            rank_ptr = array_in->rank;
            if (*rank_ptr > 24 ) error_msg("Max dimension for Bcast is limited to 24");
        }
        else
        {
            dim_ptr = dims ;
            rank_ptr = &rank;
        }
        // Bcast rank and dims
        mpi_error = MPI_Bcast(rank_ptr, 1, ELPH_MPI_ND_INT, root, comm );
        mpi_error = MPI_Bcast(dim_ptr, *rank_ptr, ELPH_MPI_ND_INT, root, comm );
        
        if (my_rank != root)
        {
            ND_function(init, Nd_floatS) (array_in, *rank_ptr, dim_ptr);
            ND_function(malloc, Nd_floatS) (array_in);
        }
    }

    // Bcast data
    ND_int size_arr = ND_function(size, Nd_floatS) (array_in);

    /* Check if size is less than int max */

    if (size_arr < ( ((ND_int)INT_MAX) -10 ) ) 
    {
        mpi_error = MPI_Bcast(array_in->data, size_arr, ELPH_MPI_float, root, comm);
    }
    else 
    {   
        ND_int max_val_int = ((ND_int)INT_MAX) -10  ;
        ND_int nloops = size_arr/max_val_int ;
        ND_int rem = size_arr%max_val_int ;

        for (ND_int iloop = 0; iloop < nloops; ++iloop)
        {   
            mpi_error = MPI_Bcast(array_in->data + iloop*max_val_int,  max_val_int, ELPH_MPI_float, root, comm);
        }
        if (rem != 0 ) 
        {
            mpi_error = MPI_Bcast(array_in->data + nloops*max_val_int, rem, ELPH_MPI_float, root, comm);
        }
    }

}


// complex
void Bcast_ND_arrayCmplx(ND_array(Nd_cmplxS) * array_in, bool alloc_mem, int root, MPI_Comm comm)
{
    if (alloc_mem)
    {
        ND_int dims[24];
        ND_int rank;

        ND_int * dim_ptr;
        ND_int * rank_ptr;

        int my_rank, comm_size, mpi_error;

        MPI_Comm_rank(comm, &my_rank);
        MPI_Comm_size(comm, &comm_size);

        if (my_rank == root)
        {
            dim_ptr =  array_in->dims;
            rank_ptr = array_in->rank;
            if (*rank_ptr > 24 ) error_msg("Max dimension for Bcast is limited to 24");
        }
        else
        {
            dim_ptr = dims ;
            rank_ptr = &rank;
        }
        // Bcast rank and dims
        mpi_error = MPI_Bcast(rank_ptr, 1, ELPH_MPI_ND_INT, root, comm );
        mpi_error = MPI_Bcast(dim_ptr, *rank_ptr, ELPH_MPI_ND_INT, root, comm );
        
        if (my_rank != root)
        {
            ND_function(init, Nd_cmplxS) (array_in, *rank_ptr, dim_ptr);
            ND_function(malloc, Nd_cmplxS) (array_in);
        }
    }

    // Bcast data
    ND_int size_arr = ND_function(size, Nd_cmplxS) (array_in);

    /* Check if size is less than int max */
    if (size_arr < ( ((ND_int)INT_MAX) -10 ) ) 
    {
        mpi_error = MPI_Bcast(array_in->data, size_arr, ELPH_MPI_cmplx, root, comm);
    }
    else 
    {   
        ND_int max_val_int = ((ND_int)INT_MAX) -10  ;
        ND_int nloops = size_arr/max_val_int ;
        ND_int rem = size_arr%max_val_int ;

        for (ND_int iloop = 0; iloop < nloops; ++iloop)
        {   
            mpi_error = MPI_Bcast(array_in->data + iloop*max_val_int,  max_val_int, ELPH_MPI_cmplx, root, comm);
        }
        if (rem != 0 ) 
        {
            mpi_error = MPI_Bcast(array_in->data + nloops*max_val_int, rem, ELPH_MPI_cmplx, root, comm);
        }
    }
}



void Bcast_wfc(struct WFC * wfc, bool alloc_mem, int root, MPI_Comm comm)
{
    /*
    Bcast wfc
    */
    Bcast_ND_arrayCmplx(wfc->wfc, alloc_mem, root, comm);
    // gvecs
    Bcast_ND_arrayFloat(wfc->gvec, alloc_mem, root, comm);
    //Fk
    Bcast_ND_arrayFloat(wfc->Fk, alloc_mem, root, comm);

    // npw_total
    mpi_error = MPI_Bcast(&(wfc->npw_total), 1, ELPH_MPI_ND_INT, root, comm );
    //npw_loc
    mpi_error = MPI_Bcast(&(wfc->npw_loc),   1, ELPH_MPI_ND_INT, root, comm );
}

