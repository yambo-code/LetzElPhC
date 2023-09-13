/*
This file contains functions for FFT planers
*/
#include "wfc.h"


void alloc_wfcBox(struct wfcBox * buffer, const ND_int rank, const ND_int * dimensions, \
                const ND_int npw_max_total, 
                const ND_int nFFT_rank, const ND_int * in_idx, unsigned flag, \
                MPI_Comm mpi_comm)
{
    /*
    This function creates a buffer(s) array(s) and generates a FFT plans
    */
    /*
    Note that dimensions are same for all cpus i.e (nspin,nbnds,nspinor,Nx,Ny,Nz)
    This code will internally distribute plane and FFT vectors over the nodes
    */
    int mpi_error, my_rank, Comm_size;
    mpi_error = MPI_Comm_size(mpi_comm, &Comm_size);
    mpi_error = MPI_Comm_rank(mpi_comm, &my_rank);

    const ND_int * FFT_dims = dimensions+3;

    ND_int nFFT = FFT_dims[0]*FFT_dims[1]*FFT_dims[2];

    int mpi_error;
    // get the total cpus in comm and rank of each cpu
    int my_rank, Comm_size;
    mpi_error = MPI_Comm_size(mpi_comm, &Comm_size);
    mpi_error = MPI_Comm_rank(mpi_comm, &my_rank);

    ND_int nffts_per_core = nFFT/Comm_size;
    ND_int nffts_rem      = nFFT%Comm_size;

    /*
    cross check this with the inputs FIX ME 
    */
    ND_int nffts_inthis_cpu = nffts_per_core;
    if (my_rank < nffts_rem) ++nffts_inthis_cpu;

    int nset_per_cpu = nsets/Comm_size; //
    int nset_rem     = nsets%Comm_size; // any remaining
    int nset_inthis_cpu = nset_per_cpu;
    if (my_rank < nset_rem) ++nset_inthis_cpu;

    ND_int dim_Buffer[4] = {dimensions[0],dimensions[1],dimensions[2],nffts_inthis_cpu};
    ND_int dim_FFTBuf[4] = {nset_inthis_cpu,FFT_dims[0],FFT_dims[1],FFT_dims[2]};
    if (nset_inthis_cpu == 0) 
    {
        dim_FFTBuf[0] =1;
        dim_FFTBuf[1] =1;
        dim_FFTBuf[2] =1;
        dim_FFTBuf[3] =1;
    }

    int pw_per_core = npw_max_total/Comm_size;
    int pw_rem      = npw_max_total%Comm_size;

    int pw_max_cpu = pw_per_core + 10;
    

    ND_function(init, Nd_cmplxS) (&(buffer->Buffer), 4, dim_Buffer);
    ND_function(init, Nd_cmplxS) (&(buffer->FFTBuf), 4, dim_FFTBuf);

    ND_function(malloc, Nd_cmplxS) (&(buffer->Buffer));
    ND_function(FFT_calloc, Nd_cmplxS) (&(buffer->FFTBuf));

    ND_function(fft_planner, Nd_cmplxS) (&(buffer->FFTBuf), &(buffer->FFTBuf), nFFT_rank, \
                    in_idx, -1, &(buffer->norm), flag, &(buffer->fft_plan));

    if (buffer->fft_plan == NULL) error_msg("FFT Plan creation failed \n");

    ND_function(fft_planner, Nd_cmplxS) (&(buffer->FFTBuf), &(buffer->FFTBuf), nFFT_rank, \
                    in_idx, 1, NULL, flag, &(buffer->inVfft_plan));

    if (buffer->inVfft_plan == NULL) error_msg("inVFFT Plan creation failed \n");


    buffer->BufGsphere = malloc(sizeof(ELPH_cmplx)*nset_inthis_cpu*npw_max_total);
    buffer->Gvecs      = malloc(sizeof(ELPH_float)*3*npw_max_total);
    buffer->Gvecs_loc  = malloc(sizeof(ELPH_float)*3*pw_max_cpu);

    int* counts_send = malloc(sizeof(int)*4*Comm_size);
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
    free(buffer->Gvecs_loc)
    free(buffer->BufGsphere)
    ND_function(destroy, Nd_cmplxS) (&(buffer->Buffer));
    ND_function(FFT_destroy, Nd_cmplxS) (&(buffer->FFTBuf));
    ND_function(fft_destroy_plan, Nd_cmplxS) (buffer->fft_plan);
    ND_function(fft_destroy_plan, Nd_cmplxS) (buffer->inVfft_plan);

}