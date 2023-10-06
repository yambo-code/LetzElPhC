/*
This file contains functions for FFT planers
*/
#include "wfc.h"


void alloc_wfcBox(struct wfcBox * buffer, const ND_int * dimensions, \
                const ND_int npw_max_total, unsigned flag, MPI_Comm mpi_comm)
{
    /*
    This function creates a buffer(s) array(s) and generates a FFT plans
    */
    /*
    dimensions = array[6]
    Note that dimensions are same for all cpus i.e (nspin,nbnds,nspinor,Nx,Ny,Nz)
    This code will internally distribute plane and FFT vectors over the nodes
    */
    int mpi_error, my_rank, Comm_size;
    mpi_error = MPI_Comm_size(mpi_comm, &Comm_size);
    mpi_error = MPI_Comm_rank(mpi_comm, &my_rank);

    const ND_int * FFT_dims = dimensions+3;
    memcpy(buffer->FFT_dimensions,FFT_dims,sizeof(ND_int)*3);

    ND_int nFFT = FFT_dims[0]*FFT_dims[1]*FFT_dims[2];

    ND_int nffts_per_core = nFFT/Comm_size;
    ND_int nffts_rem      = nFFT%Comm_size;

    /*
    cross check this with the inputs FIX ME 
    */
    ND_int nffts_inthis_cpu = nffts_per_core;
    if (my_rank < nffts_rem) ++nffts_inthis_cpu;

    ND_int nsets = dimensions[0]*dimensions[1]*dimensions[2]; 
    // total number of sets of p.ws i.e nspin*nbnd*nspinor

    int nset_per_cpu = nsets/Comm_size; //
    int nset_rem     = nsets%Comm_size; // any remaining
    int nset_inthis_cpu = nset_per_cpu;
    if (my_rank < nset_rem) ++nset_inthis_cpu;

    ND_int dim_Buffer[4] = {dimensions[0],dimensions[1],dimensions[2],nffts_inthis_cpu};
    ND_int dim_FFTBuf[4] = {nset_inthis_cpu,FFT_dims[0],FFT_dims[1],FFT_dims[2]};

    buffer->nset_loc = nset_inthis_cpu;

    if (nset_inthis_cpu == 0) 
    {   
        // although the std specicifies we can pass 0 to malloc. 
        // should be removed in future.
        dim_FFTBuf[0] =1;
        dim_FFTBuf[1] =1;
        dim_FFTBuf[2] =1;
        dim_FFTBuf[3] =1;
    }

    int pw_max_per_core = npw_max_total/Comm_size;
    int pw_max_rem      = npw_max_total%Comm_size;

    int pw_max_cpu = pw_max_per_core + 10;
    

    ND_function(init, Nd_cmplxS) (&(buffer->Buffer), 4, dim_Buffer);
    ND_function(init, Nd_cmplxS) (&(buffer->Buffer_temp), 4, dim_Buffer);
    ND_function(init, Nd_cmplxS) (&(buffer->FFTBuf), 4, dim_FFTBuf);

    ND_function(malloc, Nd_cmplxS) (&(buffer->Buffer));
    ND_function(malloc, Nd_cmplxS) (&(buffer->Buffer_temp));
    ND_function(FFT_calloc, Nd_cmplxS) (&(buffer->FFTBuf));

    buffer->ft_plan = malloc(sizeof(struct fft_plans)*(nset_inthis_cpu+1));

    if (nset_inthis_cpu != 0)
    {   
        // allocate fft and invFFT plans
        const ND_int in_idx[3] = {0,1,2}; // 3D fft is performed 
        ND_array(Nd_cmplxS) fft_buft[1];
        ND_function(init_strip_dims, Nd_cmplxS) (&(buffer->FFTBuf), 1, fft_buft);
        
        for (ND_int iset = 0 ; iset <nset_inthis_cpu ; ++iset )
        {   
            /* create a slice [iset,...]*/
            ND_function(strip_dims, Nd_cmplxS) (&(buffer->FFTBuf), 1, \
                            nd_idx{iset}, fft_buft);
            /* fft planner */
            ND_function(fft_planner, Nd_cmplxS) (fft_buft, fft_buft, 3, \
                    in_idx, -1, &(buffer->norm), flag, &(buffer->ft_plan[iset].fft_plan));

            if (buffer->ft_plan[iset].fft_plan == NULL) error_msg("FFT Plan creation failed \n");
            /* inv fft planner */
            ND_function(fft_planner, Nd_cmplxS) (fft_buft, fft_buft, 3, \
                            in_idx, 1, NULL, flag, &(buffer->ft_plan[iset].inVfft_plan));

            if (buffer->ft_plan[iset].inVfft_plan == NULL) error_msg("inVFFT Plan creation failed \n");
        }
        ND_function(uninit, Nd_cmplxS) (fft_buft);
    }
    
    buffer->BufGsphere = malloc(sizeof(ELPH_cmplx)*(nset_inthis_cpu*npw_max_total + 1));
    // 1 was added above just that we do not pass 0 to malloc (passing 0 malloc is fine though).
    buffer->Gvecs      = malloc(sizeof(ELPH_float)*3*npw_max_total);
    buffer->Gvecs_loc  = malloc(sizeof(ELPH_float)*3*pw_max_cpu);
    int* counts_send   = malloc(sizeof(int)*4*Comm_size);
    
    if (counts_send == NULL) error_msg("Failed to allocated mpi FFT buffer array \n");
    
    buffer->comm_buffer = counts_send;

    return ;
}


void free_wfcBox(struct wfcBox * buffer)
{
    /*
    This function free the wfcBox datastructures
    */
    free(buffer->comm_buffer);
    free(buffer->Gvecs);
    free(buffer->Gvecs_loc);
    free(buffer->BufGsphere);
    ND_function(destroy, Nd_cmplxS) (&(buffer->Buffer));
    ND_function(destroy, Nd_cmplxS) (&(buffer->Buffer_temp));
    ND_function(FFT_destroy, Nd_cmplxS) (&(buffer->FFTBuf));
    for  (ND_int i = 0; i<buffer->nset_loc; ++i)
    {
        ND_function(fft_destroy_plan, Nd_cmplxS) (buffer->ft_plan[i].fft_plan);
        ND_function(fft_destroy_plan, Nd_cmplxS) (buffer->ft_plan[i].inVfft_plan);
    }
    free(buffer->ft_plan);

}

