/*
This file contains functions that performs FFT on wfc/potentials
*/
#include "wfc.h"
/**/
void wfcinVFFT(ND_array(Nd_cmplxS) * wfcG,  const ELPH_float * sym, 
        const ELPH_float * tau, const int * ulmvec, \
        const bool tim_rev, const ELPH_float * lat_vec, 
        const ND_int npw_total, const ND_array(Nd_floatS)* Gvecs_in, \
        struct wfcBox * wfcRspace, MPI_Comm mpi_comm)
{   
    /*
    0) Get the gvectors 

    1) Transpose the data (where A,B,C are bands*spin*spinor)
    | A_pw1 B_pw1 C_pw1 ... |               | A_pw1 A_pw2 A_pw3 ... |
    | A_pw2 B_pw2 C_pw2 ... |   ----->      | B_pw1 B_pw2 B_pw3 ... |
    | A_pw3 B_pw3 C_pw3 ... |  AlltoAllv    | C_pw1 C_pw2 C_pw3 ... |
    | ..... ..... ..... ... |               | ..... ..... ..... ... |

    2) Perform the spin*band*spinor FFTs on different availble cpus 
    (! Note this might leave some process empty. )

    3) sphere to box 
    
    4) Perform ffts on each cpu 

    5) Transpose back the data

    | A_fft1 A_fft2 A_fft3 ... |                 | A_fft1 B_fft1 C_fft1 ... | 
    | B_fft1 B_fft2 B_fft3 ... |    ----->       | A_fft2 B_fft2 C_fft2 ... |
    | C_fft1 C_fft2 C_fft3 ... |   AlltoAllv     | A_fft3 B_fft3 C_fft3 ... |
    | .....   .....  ..... ... |                 | .....   .....  ..... ... |
    */
    
    /* First get the gvectors */

    ELPH_float * Gvecs_loc = wfcRspace->Gvecs_loc; // this is buffer to store local pws
    ELPH_float * Gtemp     = wfcRspace->Gvecs;      // this is buffer to store total pws
    
    ND_int loc_pw    = wfcG->dims[3];                      
    // number of plane waves in each local processor // this should be less than the INT_MAX 
    
    ND_int nsets     = wfcG->dims[0]*wfcG->dims[1]*wfcG->dims[2]; 
    // total number of sets of p.ws i.e nspin*nbnd*nspinor
    
    const ND_int * FFT_dims = wfcRspace->FFT_dimensions;

    ND_int nFFT      = FFT_dims[0]*FFT_dims[1]*FFT_dims[2]; 
    // product of fft dimensions

    ELPH_float G0[3] = {-ulmvec[0],-ulmvec[1],-ulmvec[2]};

    rotateGvecs(Gvecs_in->data, sym, loc_pw, lat_vec, false, true, G0, Gvecs_loc);

    ELPH_cmplx * wfc_pw_in   = wfcG->data;              // (nspin,nbnd,nspinor,loc_pw)
    ELPH_cmplx * wfc_pw_loc  = wfcRspace->BufGsphere ;  // (nsets_per_cpu,npw)
    ELPH_cmplx * wfc_fft_loc = wfcRspace->FFTBuf.data;  // (nsets_per_cpu,,Nx,Ny,Ny)
    ELPH_cmplx * wfc_fft_in  = wfcRspace->Buffer.data;  // (nspin,nbnd,nspinor,nffts_percpu)

    int mpi_error;
    // get the total cpus in comm and rank of each cpu
    int my_rank, Comm_size;
    mpi_error = MPI_Comm_size(mpi_comm, &Comm_size);
    mpi_error = MPI_Comm_rank(mpi_comm, &my_rank);

    // 
    int* counts_send = wfcRspace->comm_buffer;
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
    int nset_per_cpu = nsets/Comm_size; //
    int nset_rem     = nsets%Comm_size; // any remaining

    int pw_per_core = npw_total/Comm_size;
    int pw_rem      = npw_total%Comm_size;

    /** cross check if number of sets in this cpu and plane waves are consistant with inputs */
    int nset_inthis_cpu = nset_per_cpu;
    if (my_rank < nset_rem) ++nset_inthis_cpu;

    int temp_pw = pw_per_core;
    if (my_rank < pw_rem) ++temp_pw ;


    int disp_rectemp = 0;
    for (ND_int i = 0 ; i<Comm_size; ++i)
    {
        // this is number of pw per code
        counts_recv[i] = pw_per_core;
        if (i < pw_rem) ++counts_recv[i]; 

        counts_recv[i] *= 3;

        displacements_recv[i] = disp_rectemp;
        disp_rectemp += counts_recv[i];
    }
    /* get all the gvecs */
    MPI_Allgatherv(Gvecs_loc, 3*loc_pw, ELPH_MPI_float, Gtemp, counts_recv, \
                    displacements_recv, ELPH_MPI_float, mpi_comm);
    /*
    First sent nset_per_cpu batches
    */
    for (ND_int i = 0 ; i<Comm_size; ++i) counts_send[i]        = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) displacements_send[i] = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) counts_recv[i]        = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) displacements_recv[i] = 0; 

    disp_rectemp = 0;
    for (ND_int i = 0 ; i<Comm_size; ++i)
    {
        counts_send[i] = loc_pw;
        displacements_send[i] = i*loc_pw ;
        // this is number of pw per code
        counts_recv[i] = pw_per_core;
        if (i < pw_rem) ++counts_recv[i]; 

        displacements_recv[i] = disp_rectemp;
        disp_rectemp += counts_recv[i];
    }
    
    for (int iset =0 ; iset <nset_per_cpu; ++ iset)
    {   
        mpi_error = MPI_Alltoallv(wfc_pw_in + iset*Comm_size*loc_pw, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, wfc_pw_loc+iset*npw_total, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm);
    }
    /*
    Scatter the remainder sets
    */
    if (nset_rem != 0)
    {   
        int input_shift = nset_per_cpu*Comm_size*loc_pw;
        int out_shift = nset_per_cpu*npw_total;

        if (my_rank >= nset_rem) out_shift = 0;

        for (ND_int i = 0 ; i<Comm_size; ++i) counts_send[i]        = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) displacements_send[i] = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) counts_recv[i]        = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) displacements_recv[i] = 0; 

        for (ND_int i = 0 ; i<nset_rem; ++i)
        {               
            counts_send[i] = loc_pw;
            displacements_send[i] = i*loc_pw ;
        }

        if (my_rank < nset_rem)
        {   
            disp_rectemp = 0;
            for (ND_int i = 0 ; i<Comm_size; ++i)
            {
                counts_recv[i] = pw_per_core;
                if (i < pw_rem) ++counts_recv[i]; 

                displacements_recv[i] = disp_rectemp;
                disp_rectemp += counts_recv[i];
            }
        }
        
        mpi_error = MPI_Alltoallv(wfc_pw_in + input_shift, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, wfc_pw_loc + out_shift, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm);
    }
    
    
    /* Now we have the data */ 
    /* we FFT and scatter back i.e Transpose back the data */
    int nffts_per_core = nFFT/Comm_size;
    int nffts_rem      = nFFT%Comm_size;

    int nffts_inthis_cpu = nffts_per_core;
    if (my_rank < nffts_rem) ++nffts_inthis_cpu;


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

    for (int iset =0 ; iset <nset_per_cpu; ++ iset)
    {   
        ELPH_cmplx * restrict wpwloc  = wfc_pw_loc + iset*npw_total;
        ELPH_cmplx * restrict wfftloc = wfc_fft_loc + iset*nFFT;
        
        sphere2box(wpwloc , 1, Gtemp, npw_total, FFT_dims, wfftloc);
        /* perform the FFT */
        ND_function(fft_execute_plan, Nd_cmplxS) (wfcRspace->ft_plan[iset].inVfft_plan);
        /* send them back */
        mpi_error = MPI_Ialltoallv(wfftloc, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, wfc_fft_in + iset*nffts_inthis_cpu*Comm_size, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm, \
                    &(wfcRspace->ft_plan[iset].request) );
    }

    MPI_Request req_rem;

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
        

        /* The last set is anyways blocking */
        if (my_rank < nset_rem)
        {   
            ND_int iset = nset_per_cpu;
            ELPH_cmplx * restrict wpwloc  = wfc_pw_loc + iset*npw_total;
            ELPH_cmplx * restrict wfftloc = wfc_fft_loc + iset*nFFT;
        
            sphere2box(wpwloc , 1, Gtemp, npw_total, FFT_dims, wfftloc);
            /* perform the FFT */
            ND_function(fft_execute_plan, Nd_cmplxS) (wfcRspace->ft_plan[iset].inVfft_plan);
        }
        mpi_error = MPI_Ialltoallv(wfc_fft_loc + input_shift, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, wfc_fft_in + out_shift, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm,&req_rem);
    }

    /* Wait for all alltoallv */

    for (int iset =0 ; iset <nset_per_cpu; ++ iset)
    {   
        MPI_Wait(&(wfcRspace->ft_plan[iset].request), MPI_STATUS_IGNORE);
    }
    if (nset_rem != 0) MPI_Wait(&req_rem, MPI_STATUS_IGNORE);

    /*
    Now rotate the wfcs in spin space
    */
    ND_int nspinor = wfcG->dims[2];
    ELPH_cmplx su2mat[4];

    SU2mat(sym, nspinor, false, tim_rev, su2mat);

    // rotate the wavefunction
    su2rotate(nspinor, nffts_inthis_cpu, wfcG->dims[0]*wfcG->dims[1], su2mat, wfc_fft_in);
    if (tim_rev)
    {   
        ND_int size_wfc_in = ND_function(size, Nd_cmplxS) (&(wfcRspace->Buffer));
        for (ND_int i = 0; i<size_wfc_in; ++i) wfc_fft_in[i] = conj(wfc_fft_in[i]);
    }


}