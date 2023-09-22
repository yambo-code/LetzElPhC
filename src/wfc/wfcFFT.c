#include "wfc.h"

void wfcFFT(struct wfcBox * wfcRspace, const ELPH_float * sym, \
        const ELPH_float * tau, const int * ulmvec, const bool tim_rev, \
        const ELPH_float * lat_vec, const ND_int npw_total, \
        const ND_array(Nd_floatS)* Gvecs_in, \
        ND_array(Nd_cmplxS) * wfcG, MPI_Comm mpi_comm)
{
    /*
    
    1) Transpose the data (where A,B,C are bands*spin*spinor)
    | A_fft1 B_fft1 C_fft1 ... |               | A_fft1 A_fft2 A_fft3 ... |
    | A_fft2 B_fft2 C_fft2 ... |   ----->      | B_fft1 B_fft2 B_fft3 ... |
    | A_fft3 B_fft3 C_fft3 ... |  AlltoAllv    | C_fft1 C_fft2 C_fft3 ... |
    | .....   .....  ..... ... |               | .....   .....  ..... ... |
    2) Perform the spin*band*spinor FFTs on different availble cpus 
    (! Note this might leave some process empty. )

    3) Perform ffts on each cpu
    
    4) box to sphere

    5) Transpose back the data

    
    | A_pw1 A_pw2 A_pw3 ... |             | A_pw1 B_pw1 C_pw1 ... | 
    | B_pw1 B_pw2 B_pw3 ... |   ----->    | A_pw2 B_pw2 C_pw2 ... | 
    | C_pw1 C_pw2 C_pw3 ... |  AlltoAllv  | A_pw3 B_pw3 C_pw3 ... | 
    | ..... ..... ..... ... |             | ..... ..... ..... ... | 
    */

    ELPH_float * Gvecs_loc = wfcRspace->Gvecs_loc; // this is buffer to store local pws
    ELPH_float * Gtemp     = wfcRspace->Gvecs;      // this is buffer to store total pws
    
    ND_int loc_pw    = wfcG->dims[3];                      
    // number of plane waves in each local processor // this should be less than the INT_MAX 
    
    ND_int nsets     = wfcG->dims[0]*wfcG->dims[1]*wfcG->dims[2]; 
    // total number of sets of p.ws i.e nspin*nbnd*nspinor
    
    const ND_int * FFT_dims = wfcRspace->FFT_dimensions;

    ND_int nFFT      = FFT_dims[0]*FFT_dims[1]*FFT_dims[2]; 
    // product of fft dimensions

    ELPH_float G0[3] = {ulmvec[0],ulmvec[1],ulmvec[2]};

    rotateGvecs(Gvecs_in->data, sym, loc_pw, lat_vec, true, true, G0, Gvecs_loc);

    ELPH_cmplx * wfc_pw_in   = wfcG->data;              // (nspin,nbnd,nspinor,loc_pw)
    ELPH_cmplx * wfc_pw_loc  = wfcRspace->BufGsphere ;  // (nsets_per_cpu,npw)
    ELPH_cmplx * wfc_fft_loc = wfcRspace->FFTBuf.data;  // (nsets_per_cpu,,Nx,Ny,Ny)
    ELPH_cmplx * wfc_fft_in  = wfcRspace->Buffer.data;  // (nspin,nbnd,nspinor,nffts_percpu)


    int mpi_error;
    // get the total cpus in comm and rank of each cpu
    int my_rank, Comm_size;
    mpi_error = MPI_Comm_size(mpi_comm, &Comm_size);
    mpi_error = MPI_Comm_rank(mpi_comm, &my_rank);


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

    int nffts_per_core = nFFT/Comm_size;
    int nffts_rem      = nFFT%Comm_size;

    /*
    cross check this with the inputs FIX ME 
    */
    int nffts_inthis_cpu = nffts_per_core;
    if (my_rank < nffts_rem) ++nffts_inthis_cpu;

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
        counts_send[i] = nffts_inthis_cpu ;
        displacements_send[i] = i*nffts_inthis_cpu;  //displacements_send

        counts_recv[i] = nffts_per_core;
        if(i<nffts_rem)  ++counts_recv[i];
        displacements_recv[i] = disp_rectemp ;

        disp_rectemp += counts_recv[i];
        
    }
    /*
    send the sets
    */
    for (int iset =0 ; iset <nset_per_cpu; ++ iset)
    {   
        mpi_error = MPI_Alltoallv(wfc_fft_in + iset*nffts_inthis_cpu*Comm_size, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, wfc_fft_loc + iset*nFFT, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm);
    }

    /* scatter the remainder sets*/

    if (nset_rem != 0)
    {   
        int input_shift = nset_per_cpu*nffts_inthis_cpu*Comm_size;
        int out_shift = nset_per_cpu*nFFT;

        if (my_rank >= nset_rem) out_shift = 0;

        for (ND_int i = 0 ; i<Comm_size; ++i) counts_send[i]        = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) displacements_send[i] = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) counts_recv[i]        = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) displacements_recv[i] = 0; 

        for (ND_int i = 0 ; i<nset_rem; ++i)
        {               
            counts_send[i] = nffts_inthis_cpu ;
            displacements_send[i] = i*nffts_inthis_cpu;  //displacements_send
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
        
        mpi_error = MPI_Alltoallv(wfc_fft_in + input_shift, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, wfc_fft_loc + out_shift, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm);
    }
    
    /* Now every cpu has nset_per_cpu + 0/1 wfcs. we do a fft */

    if (nset_inthis_cpu != 0 )
    {   
        ND_function(fft_execute_plan, Nd_cmplxS) (wfcRspace->fft_plan);
        ELPH_cmplx fft_norm = 1.0/nFFT;
        ND_int size_in = ND_function(size, Nd_cmplxS) (&(wfcRspace->FFTBuf));
        /* Normalize */
        ELPH_OMP_PAR_FOR_SIMD
        for (ND_int i = 0 ; i<size_in; ++i) wfc_fft_loc[i] *= fft_norm;

        /* box2sphere */
        box2sphere(wfc_fft_loc, nset_inthis_cpu, \
                Gtemp, npw_total, FFT_dims, wfc_pw_loc);
    }


    for (ND_int i = 0 ; i<Comm_size; ++i) counts_send[i]        = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) displacements_send[i] = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) counts_recv[i]        = 0; 
    for (ND_int i = 0 ; i<Comm_size; ++i) displacements_recv[i] = 0; 

    /*
    First sent nset_per_cpu batches
    */
    int disp_sendtemp = 0;
    for (ND_int i = 0 ; i<Comm_size; ++i)
    {
        counts_send[i] = pw_per_core;
        if(i<pw_rem) ++counts_send[i];
        displacements_send[i] = disp_sendtemp ;
        disp_sendtemp += counts_send[i];
        counts_recv[i] = loc_pw ;
        displacements_recv[i] = i*loc_pw;
    }
    /*
    send the sets
    */
    for (int iset =0 ; iset <nset_per_cpu; ++ iset)
    {   
        mpi_error = MPI_Alltoallv(wfc_pw_loc + iset*npw_total, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, wfc_pw_in+iset*Comm_size*loc_pw, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm);
    }
    /*
    Scatter the remainder sets
    */
    if (nset_rem != 0)
    {   
        int input_shift = nset_per_cpu*npw_total;
        int out_shift = nset_per_cpu*loc_pw*Comm_size;
        if (my_rank >= nset_rem) input_shift = 0;

        for (ND_int i = 0 ; i<Comm_size; ++i) counts_send[i]        = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) displacements_send[i] = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) counts_recv[i]        = 0; 
        for (ND_int i = 0 ; i<Comm_size; ++i) displacements_recv[i] = 0; 

        for (ND_int i = 0 ; i<nset_rem; ++i)
        {               
            counts_recv[i] = loc_pw ;
            displacements_recv[i] = i*loc_pw;
        }

        if (my_rank < nset_rem)
        {   
            disp_sendtemp = 0;
            for (ND_int i = 0 ; i<Comm_size; ++i)
            {
                counts_send[i] = pw_per_core;
                if(i<pw_rem) ++counts_send[i];
                displacements_send[i] = disp_sendtemp ;
                disp_sendtemp += counts_send[i];
            }
        }
        
        mpi_error = MPI_Alltoallv(wfc_pw_loc + input_shift, counts_send, \
                    displacements_send, ELPH_MPI_cmplx, wfc_pw_in + out_shift, \
                    counts_recv, displacements_recv, ELPH_MPI_cmplx, mpi_comm);
    }


    /*
    Now rotate the wfcs in spin space
    */
    ND_int nspinor = wfcG->dims[2];
    ELPH_cmplx su2mat[4];

    SU2mat(sym, nspinor, false, tim_rev, su2mat);

    su2rotate(nspinor, loc_pw, wfcG->dims[0]*wfcG->dims[1], su2mat, wfc_pw_in);
    /*
    conjugate in case of time reversal
    */
    if (tim_rev)
    {   
        ND_int size_wfc_in = ND_function(size, Nd_cmplxS) (wfcG);
        for (ND_int i = 0; i<size_wfc_in; ++i) wfc_pw_in[i] = conj(wfc_pw_in[i]);
    }

}