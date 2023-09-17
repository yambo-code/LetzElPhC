/*
This file contains functions that performs invFFT on local potentials
*/
#include "dvloc.h"
/**/
void VlocinVFFT(ND_array(Nd_cmplxS) * VlocPot, struct Lattice * lattice, MPI_Comm mpi_comm)
{   
    /*

    1) Transpose the data (where A,B,C are bands*spin*spinor)
    | A_fft1 B_fft1 C_fft1 ... |               | A_fft1 A_fft2 A_fft3 ... |
    | A_fft2 B_fft2 C_fft2 ... |   ----->      | B_fft1 B_fft2 B_fft3 ... |
    | A_fft3 B_fft3 C_fft3 ... |  AlltoAllv    | C_fft1 C_fft2 C_fft3 ... |
    | .....   .....  ..... ... |               | .....   .....  ..... ... |

    2) Perform the nset FFTs on different availble cpus 
    (! Note this might leave some process empty. )

    3) Transpose back the data

    | A_fft1 A_fft2 A_fft3 ... |                 | A_fft1 B_fft1 C_fft1 ... | 
    | B_fft1 B_fft2 B_fft3 ... |    ----->       | A_fft2 B_fft2 C_fft2 ... |
    | C_fft1 C_fft2 C_fft3 ... |   AlltoAllv     | A_fft3 B_fft3 C_fft3 ... |
    | .....   .....  ..... ... |                 | .....   .....  ..... ... |
    */
    int mpi_error;
    // get the total cpus in comm and rank of each cpu
    int my_rank, Comm_size;
    mpi_error = MPI_Comm_size(mpi_comm, &Comm_size);
    mpi_error = MPI_Comm_rank(mpi_comm, &my_rank);

    ND_int nsets            = VlocPot->dims[0]; 
    const ND_int * FFT_dims = lattice->fft_dims;
    ND_int nFFT             = FFT_dims[0]*FFT_dims[1]*FFT_dims[2]; 

    int nset_per_cpu = nsets/Comm_size; //
    int nset_rem     = nsets%Comm_size; // any remaining

    int nset_inthis_cpu = nset_per_cpu;
    if (my_rank < nset_rem) ++nset_inthis_cpu;

    int nffts_per_core = nFFT/Comm_size;
    int nffts_rem      = nFFT%Comm_size;

    int nffts_inthis_cpu = nffts_per_core;
    if (my_rank < nffts_rem) ++nffts_inthis_cpu;

    ELPH_cmplx * VlocPot_in = VlocPot->data;

    ND_array(Nd_cmplxS) VlocPot_loc[1];

    ND_int dim_FFTBuf[4] = {nset_inthis_cpu,FFT_dims[0],FFT_dims[1],FFT_dims[2]};
    if (nset_inthis_cpu == 0) 
    {   
        // although the std specicifies we can pass 0 to malloc. 
        // should be removed in future.
        dim_FFTBuf[0] =1;
        dim_FFTBuf[1] =1;
        dim_FFTBuf[2] =1;
        dim_FFTBuf[3] =1;
    }
    ND_function(init, Nd_cmplxS) (VlocPot_loc, 4, dim_FFTBuf);
    ND_function(FFT_calloc, Nd_cmplxS) (VlocPot_loc);

    ELPH_cmplx * VlocPot_fft_loc = VlocPot_loc->data;
    // 
    // create the FFT plan. Note planner must be created before filling the data
    ND_function(FFT_plan, Nd_cmplxS) inVfft_plan;

    if (nset_inthis_cpu != 0 )
    {   
        /* perform the FFT */
        const ND_int in_idx[3] = {1,2,3}; // fft is performed over last three dimensions 
        
        ND_function(fft_planner, Nd_cmplxS) (VlocPot_loc, VlocPot_loc, 3, in_idx, 1, NULL, \
                                                FFTW_MEASURE, &(inVfft_plan));

        if (inVfft_plan == NULL) error_msg("inVFFT Plan creation failed \n");

    }
    else inVfft_plan = NULL ;

    int* counts_send = malloc(sizeof(int)*4*Comm_size);
    if(counts_send == NULL) error_msg("allocation of mpi_buffer arrays failed");
    int* displacements_send = counts_send +   Comm_size ;
    int* counts_recv        = counts_send + 2*Comm_size ;
    int* displacements_recv = counts_send + 3*Comm_size ;

    /* zero out buffers */
    // memset ?
    for (ND_int i = 0 ; i<Comm_size; ++i) counts_send[i]        = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) displacements_send[i] = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) counts_recv[i]        = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) displacements_recv[i] = 0; 
    //
    /*First sent nset_per_cpu batches
    */
    int disp_rectemp = 0;
    for (ND_int i = 0 ; i<Comm_size; ++i)
    {
        counts_send[i] = nffts_inthis_cpu;
        displacements_send[i] = i*nffts_inthis_cpu ;
        // this is number of pw per code
        counts_recv[i] = nffts_per_core;
        if (i < nffts_rem) ++counts_recv[i]; 

        displacements_recv[i] = disp_rectemp;
        disp_rectemp += counts_recv[i];
    }

    for (int iset =0 ; iset <nset_per_cpu; ++ iset)
    {   
        mpi_error = MPI_Alltoallv(VlocPot_in + iset*Comm_size*nffts_inthis_cpu, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, VlocPot_fft_loc+iset*nFFT, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm);
    }
    /*
    Scatter the remainder sets
    */
    if (nset_rem != 0)
    {   
        int input_shift = nset_per_cpu*Comm_size*nffts_inthis_cpu;
        int out_shift = nset_per_cpu*nFFT;

        if (my_rank >= nset_rem) out_shift = 0;

        for (ND_int i = 0 ; i<Comm_size; ++i) counts_send[i]        = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) displacements_send[i] = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) counts_recv[i]        = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) displacements_recv[i] = 0; 

        for (ND_int i = 0 ; i<nset_rem; ++i)
        {               
            counts_send[i] = nffts_inthis_cpu;
            displacements_send[i] = i*nffts_inthis_cpu ;
        }

        if (my_rank < nset_rem)
        {   
            disp_rectemp = 0;
            for (ND_int i = 0 ; i<Comm_size; ++i)
            {
                counts_recv[i] = nffts_per_core;
                if (i < nffts_rem) ++counts_recv[i]; 

                displacements_recv[i] = disp_rectemp;
                disp_rectemp += counts_recv[i];
            }
        }
        
        mpi_error = MPI_Alltoallv(VlocPot_in + input_shift, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, VlocPot_fft_loc + out_shift, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm);
    }
    
    /* Now every cpu has nset_per_cpu + 0/1 VlocPots. we do a invfft */
    if (nset_inthis_cpu != 0 )
    {   
        /* perform the FFT */
        ND_function(fft_execute_plan, Nd_cmplxS) (inVfft_plan);

        /* destroy the plan */
        ND_function(fft_destroy_plan, Nd_cmplxS) (inVfft_plan);

    }

    /* Finally scatter back i.e Transpose back the data */
    for (ND_int i = 0 ; i<Comm_size; ++i) counts_send[i]        = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) displacements_send[i] = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) counts_recv[i]        = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) displacements_recv[i] = 0; 

    int disp_sendtemp = 0;
    for (ND_int i = 0 ; i<Comm_size; ++i)
    {
        counts_send[i] = nffts_per_core;
        if(i<nffts_rem) ++counts_send[i];
        displacements_send[i] = disp_sendtemp ;
        disp_sendtemp += counts_send[i];
        counts_recv[i] = nffts_inthis_cpu ;
        displacements_recv[i] = i*nffts_inthis_cpu;
    }
    /*
    send the sets
    */
    for (int iset =0 ; iset <nset_per_cpu; ++ iset)
    {   
        mpi_error = MPI_Alltoallv(VlocPot_fft_loc + iset*nFFT, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, VlocPot_in + iset*nffts_inthis_cpu*Comm_size, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm);
    }

    /* scatter the remainder sets*/
    if (nset_rem != 0)
    {   

        int input_shift = nset_per_cpu*nFFT;
        int out_shift = nset_per_cpu*nffts_inthis_cpu*Comm_size;
        if (my_rank >= nset_rem) input_shift = 0;
        
        for (ND_int i = 0 ; i<Comm_size; ++i) counts_send[i]        = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) displacements_send[i] = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) counts_recv[i]        = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) displacements_recv[i] = 0; 

        for (ND_int i = 0 ; i<nset_rem; ++i)
        {               
            counts_recv[i] = nffts_inthis_cpu;
            displacements_recv[i] = i*nffts_inthis_cpu ;
        }

        if (my_rank < nset_rem)
        {   
            disp_sendtemp = 0;
            for (ND_int i = 0 ; i<Comm_size; ++i)
            {
                counts_send[i] = nffts_per_core;
                if(i<nffts_rem) ++counts_send[i];
                displacements_send[i] = disp_sendtemp ;
                disp_sendtemp += counts_send[i];
            }
        }
        
        mpi_error = MPI_Alltoallv(VlocPot_fft_loc + input_shift, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, VlocPot_in + out_shift, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm);
    }
    
    /* free the memory */
    ND_function(FFT_destroy, Nd_cmplxS) (VlocPot_loc);
    free(counts_send);

}

