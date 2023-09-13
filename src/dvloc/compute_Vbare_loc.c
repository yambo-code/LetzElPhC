#include "dvloc.h"

/* Compute the electron phonon matrix elements i.e sandwich for Local part of KS potential */
static void elphLocal(const ELPH_float * qpt, struct WFC * wfcs, struct Lattice * lattice, \
                int ikq, int ik, int kqsym, int ksym, ND_array(Nd_cmplxS) * dVlocr, \
                MPI_Comm commK, MPI_Comm commQ,  struct wfcBox * wfcRspace, ELPH_cmplx * elph_kq)
{
    /* Computes <S2*k2 | dV_{q}local | S1*k1>
    Note that the inputs kvectors must full the following condition S2*k2 = S1*k1 + q + ulmveckq
    Input : . k+q and k wave functions, (nspin, bands, nspinor, ng) 
            . kqvec, kvec are kvectors in (cart units ) for wfc_kq and wfc_k
            . Gkq, Gk  are Gvectors for wfc_kq and wfc_k,
            . ulmveckq is shift that is applyed to k+q vector i.e Skq*kq - S*k - q
            . symkq and symk are symmetric operations to be applied on k+q and k wave functions
            . taukq and tauk are fractional translation
            . Gkq and Gk are Gvectors for k+q and k wf respectively
            . timerev(k)/(kq) -- true if sym(k/kq) is time rev
            . dVlocr -> change in local KS potential i.e dVlocal + dVscf (nmag, Nx, Ny, Nz) in for one mode
            . qpt -> qpoint in crystal units
    output : (nu, nspin, mk, nk+q) // electron-phonon mat in mode basis
    */
    
    int mpi_error;
    int krank;
    mpi_error = MPI_Comm_rank(commK, &krank);
    /*
    First we get the wfcs.
    */
    struct WFC wfcs_k_q, wfcs_k; // Note these are allocated by 
    // get_wfc_from_pool and must be freed in this function

    /* create buffers which must be destroyed in the end */
    ND_array(Nd_cmplxS) wfc_kq[1],wfc_k[1];
    ND_array(Nd_floatS) Gkq[1], Gk[1];
    ND_array(Nd_floatS) Fkq[1], Fk[1];

    wfcs_k_q.wfc = wfc_kq;
    wfcs_k.wfc   = wfc_k;
    wfcs_k_q.gvec = Gkq;
    wfcs_k.gvec   = Gk;
    wfcs_k_q.Fk = Fkq;
    wfcs_k.Fk   = Fk;

    /* Now get the wfcs for k+q and k */
    get_wfc_from_pool(wfcs, ikq, lattice->kpt_iredBZ->dims[0], commK, commQ, &wfcs_k_q);
    get_wfc_from_pool(wfcs, ik, lattice->kpt_iredBZ->dims[0], commK, commQ, &wfcs_k);
    
    /*initialization and setup */
    const ELPH_float * kqvec            = ND_function(ele,Nd_floatS)(lattice->kpt_iredBZ, nd_idx{ikq,0});
    const ELPH_float * kvec             = ND_function(ele,Nd_floatS)(lattice->kpt_iredBZ, nd_idx{ik,0}); 
    ELPH_float * lat_vec                = lattice->alat_vec->data;
    const ELPH_float * symkq            = ND_function(ele,Nd_floatS)(lattice->sym_mat, nd_idx{kqsym,0,0}) ;
    const ELPH_float * symk             = ND_function(ele,Nd_floatS)(lattice->sym_mat, nd_idx{ksym,0,0}) ;
    const ELPH_float * taukq            = ND_function(ele,Nd_floatS)(lattice->frac_trans, nd_idx{kqsym,0});
    const ELPH_float * tauk             = ND_function(ele,Nd_floatS)(lattice->frac_trans, nd_idx{ksym,0});
    const bool timerevkq                = lattice->time_rev_array[kqsym];
    const bool timerevk                 = lattice->time_rev_array[ksym];

    //(nspin,nbnds,nspinor,npw)

    ND_int nspin, nbndsk, nbndskq, nspinor, npwkq, npwk ;
    nspin   = wfc_k->dims[0] ;
    nbndsk  = wfc_k->dims[1] ;
    nbndskq = wfc_kq->dims[1];
    nspinor = wfc_k->dims[2] ;
    npwkq   = wfc_kq->dims[3] ;
    npwk    = wfc_k->dims[3] ;
    ND_int * FFT_dims = lattice->fft_dims;
    
    int ulmveckq[3]; // ulmveckq is shift that is applyed to k+q vector i.e -(Skq*kq - S*k - q)
    ELPH_float tempSkq[3] = {0,0,0}; //S2*k2
    ELPH_float tempSk[3] = {0,0,0};  // S1*k1

    MatVec3f(symkq, kqvec, false, tempSkq);
    MatVec3f(symk,  kvec,  false,  tempSk);
    for (int xi = 0 ; xi<3 ; ++xi ) tempSkq[xi] -= tempSk[xi] ;
    MatVec3f(lat_vec,  tempSkq,  true,  tempSk);
    for (int xi = 0 ; xi<3 ; ++xi )  ulmveckq[xi] = -(int)rint(tempSk[xi]-qpt[xi]);

    if (!trivial_phase)
    {
        kqvec = NULL ;
        kvec = NULL ;
    }
    ND_array(Nd_cmplxS) dVpsiG, wfc_SkqG;
    ND_array(Nd_floatS) GSkq_rot ; // Gvectors for rotated wavefunction i.e S2*k2


    ND_int elph_buffer_len = nbndsk*nbndskq*nspin ; 

    ELPH_OMP_PAR_FOR_SIMD
    for (ND_int i =0 ; i<elph_buffer_len; ++i) elph_kq[i] = 0.0 ;


    ND_function(init,Nd_cmplxS)(&dVpsiG, 4, nd_idx{nspin,nbndsk,nspinor,npwkq});
    ND_function(calloc,Nd_cmplxS)(&dVpsiG);
    
    /* Rotated wavefunction k+q*/
    ND_function(init,Nd_cmplxS)(&wfc_SkqG, wfc_kq->rank[0], wfc_kq->dims);
    ND_function(calloc,Nd_cmplxS)(&wfc_SkqG);
    /* G vectors for the above rotated wavefunction */
    ND_function(init,Nd_floatS)(&GSkq_rot, Gkq->rank[0], Gkq->dims); 
    ND_function(calloc,Nd_floatS)(&GSkq_rot);

    const ELPH_float Imat3[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}; // 3x3 identity matrix 

    const ELPH_float zero_temp[3] = {0.0,0.0,0.0} ;
    /*Get the S2K2 wavefunction */
    wfcG_rotate(wfc_kq, Gkq, kqvec, symkq, taukq, ulmveckq, false, false, timerevkq, &GSkq_rot, &wfc_SkqG);

    ND_int nmag; /* nmag = 1 for non magnetic and = 2/4 for spin polarized/magnetic systems */
    /* nmag is the spinor dimension of the change in potential */
    nmag = dVlocr->dims[0];

    const ND_int nFFT_tot = FFT_dims[0]*FFT_dims[1]*FFT_dims[2]; 

    // Now compute dV(r)*psi(r)

    ELPH_OMP_PAR
    {
    struct wfcBox wfcr_buffer ;    
    //(nspin, bands, nspinor, Nx, Ny, Nz )
    ELPH_OMP_PAR_CRITICAL
    {   
        /* must be in critical, not thread safe*/
        alloc_wfcBox(&wfcr_buffer,  6, nd_idx{1, 1, lattice->nspinor, FFT_dims[0], \
                            FFT_dims[1], FFT_dims[2]}, 3, nd_idx{3,4,5}, true, FFTW_MEASURE);
    }
    ND_array(Nd_cmplxS) * dVpsi = &(wfcr_buffer.inBuf);

    ND_array(Nd_cmplxS) wfc_k_temp, dVpsiG_temp;
    
    ELPH_OMP_PAR_CRITICAL
    {
        ND_function(init,Nd_cmplxS)(&wfc_k_temp,  4, nd_idx{1,1,nspinor,npwk});
        ND_function(init,Nd_cmplxS)(&dVpsiG_temp, 4, nd_idx{1,1,nspinor,npwkq});
    }

    for (ND_int is = 0 ; is <nspin; ++is )
    {   
        ELPH_cmplx *  dV_r ;
        if (nmag == 2) dV_r = dVlocr->data + is*nFFT_tot ;
        else           dV_r = dVlocr->data ;

        ELPH_OMP_FOR
        for (ND_int ibnd = 0 ; ibnd <nbndsk ; ++ibnd)
        {   
            wfc_k_temp.data = wfc_k->data  + is*wfc_k->strides[0] + ibnd*wfc_k->strides[1];
            dVpsiG_temp.data = dVpsiG.data + is*dVpsiG.strides[0] + ibnd*dVpsiG.strides[1];
            /* Get Sk in real space*/
            wfcinVFFT(&wfc_k_temp, kvec, symk, tauk, NULL, false, false, timerevk, lat_vec, Gk, &wfcr_buffer);
            
            ELPH_OMP_PAR_CRITICAL
            if (!wfcr_buffer.inplace) ND_function(copy,Nd_cmplxS)(&wfcr_buffer.outBuf,dVpsi);

            /* this is inplace */
            dVlocPsi(nmag, nspinor, nFFT_tot, dV_r, dVpsi->data);
            
            /***/
            wfcFFT(&wfcr_buffer, zero_temp, Imat3, zero_temp, NULL, false, false, false, lat_vec, &GSkq_rot, &dVpsiG_temp);
        }
    }
    ELPH_OMP_PAR_CRITICAL
    {   
        /* must be in critical, not thread safe*/
        ND_function(uninit,Nd_cmplxS)(&wfc_k_temp);
        ND_function(uninit,Nd_cmplxS)(&dVpsiG_temp);
        free_wfcBox(&wfcr_buffer);
    }
    }

    // Compute the sandwich
    for (ND_int is = 0 ; is <nspin; ++is )
    {   /*msg, nsg -> mn */
        ND_function(matmulX, Nd_cmplxS) ('N', 'C', dVpsiG.data + is*(dVpsiG.strides[0]), \
        wfc_SkqG.data + is*(wfc_SkqG.strides[0]), elph_kq+ is*nbndsk*nbndskq, 1.0, 0.0, dVpsiG.strides[1], \
        wfc_SkqG.strides[1], nbndskq, nbndsk, nbndskq, dVpsiG.strides[1]);
    }
    
    // Free stuff
    
    /*free wfc buffers */
    ND_function(destroy,Nd_cmplxS)(&wfc_tempkq);
    ND_function(destroy,Nd_cmplxS)(&wfc_tempk);
    ND_function(destroy,Nd_floatS)(&gvec_tempkq);
    ND_function(destroy,Nd_floatS)(&gvec_tempk);
    ND_function(destroy,Nd_floatS)(&Fk_tempkq);
    ND_function(destroy,Nd_floatS)(&Fk_tempk);
    
}






