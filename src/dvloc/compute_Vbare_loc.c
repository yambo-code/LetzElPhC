#include "dvloc.h"

/* Compute the electron phonon matrix elements i.e sandwich for Local part of KS potential */
void elphLocal(const ELPH_float * qpt, struct WFC * wfcs, struct Lattice * lattice, \
                int ikq, int ik, int kqsym, int ksym, ND_array(Nd_cmplxS) * dVlocr, \
                MPI_Comm commK, struct wfcBox * wfcRspace, ELPH_cmplx * elph_kq)
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
            . dVlocr -> change in local KS potential i.e dVlocal + dVscf (nmodes, nmag, nffts_in_this_cpu) in for one mode
            . qpt -> qpoint in crystal units
    output : (nu, nspin, mk, nk+q) // electron-phonon mat in mode basis (only the master process writes)
    */
    
    int mpi_error;
    int krank;
    mpi_error = MPI_Comm_rank(commK, &krank);
    /*
    First we get the wfcs. */

    ND_array(Nd_cmplxS) wfc_kq[1], dVpsiG[1];
    ND_array(Nd_floatS) Gkq[1];

    ND_int npwkq_total = (wfcs+ikq)->npw_total;
    ND_int npwk_total =  (wfcs+ik)->npw_total;

    ND_array(Nd_cmplxS) * wfc_k = (wfcs+ik)->wfc;
    ND_array(Nd_floatS) * Gk    = (wfcs+ik)->gvec;
    /*initialization and setup */
    ND_function(init,Nd_cmplxS)(wfc_kq, (wfcs+ikq)->wfc->rank[0],  (wfcs+ikq)->wfc->dims);
    ND_function(malloc,Nd_cmplxS)(wfc_kq);

    ND_function(init,Nd_cmplxS)(dVpsiG, (wfcs+ikq)->wfc->rank[0],  (wfcs+ikq)->wfc->dims);
    ND_function(malloc,Nd_cmplxS)(dVpsiG);

    ND_function(init,Nd_floatS)(Gkq,  (wfcs+ikq)->gvec->rank[0], (wfcs+ikq)->gvec->dims);
    ND_function(malloc,Nd_floatS)(Gkq);
    
    // copy k+q wfcs and gvecs
    ND_function(copy, Nd_cmplxS) ((wfcs+ikq)->wfc,wfc_kq);  
    ND_function(copy, Nd_floatS) ((wfcs+ikq)->gvec,Gkq);

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

    ND_int nspin, nbndsk, nbndskq, nspinor, nmodes, npwkq, npwk ;
    nspin   = wfc_k->dims[0] ;
    nbndsk  = wfc_k->dims[1] ;
    nbndskq = wfc_kq->dims[1];
    nspinor = wfc_k->dims[2] ;
    npwkq   = Gkq->dims[1] ;
    npwk    = Gk->dims[1] ;
    nmodes  = dVlocr->dims[0];

    ELPH_float ulmveckq[3]; // ulmveckq is shift that is applyed to k+q vector i.e -(Skq*kq - S*k - q)
    ELPH_float tempSkq[3] = {0,0,0}; //S2*k2
    ELPH_float tempSk[3] = {0,0,0};  // S1*k1

    MatVec3f(symkq, kqvec, false, tempSkq);
    MatVec3f(symk,  kvec,  false,  tempSk);
    for (int xi = 0 ; xi<3 ; ++xi ) tempSkq[xi] -= tempSk[xi] ;

    // convert qpt to cartisian coord
    // tempSkq is Skq*kq - S*k
    {
        ELPH_float blat[9];
        reciprocal_vecs(lat_vec, blat); // note this has 2*pi//
        MatVec3f(blat,  qpt,  false,  tempSk);
        for (int xi = 0 ; xi<3 ; ++xi ) tempSk[xi] = tempSk[xi]/(2*ELPH_PI) ;
    }

    for (int xi = 0 ; xi<3 ; ++xi )  ulmveckq[xi] = -(tempSkq[xi]-tempSk[xi]);

    /* rotate gvectors and shift */
    // C'_{S*G - G0} = C_{G}
    for (ND_int ipw= 0 ; ipw <npwkq ; ++ipw)
    {   
        ELPH_float * restrict GPtr = Gkq->data + 3*ipw;
        ELPH_float tempG[3]={0,0,0};
        MatVec3f(symkq, GPtr, false, tempG);

        GPtr[0] = tempG[0]-ulmveckq[0] ;
        GPtr[1] = tempG[1]-ulmveckq[1] ;
        GPtr[2] = tempG[2]-ulmveckq[2] ;
    }
    // rotate the wave function in spin space 
    ELPH_cmplx su2kq[4]  = {1,0,0,1};
    /* Get SU(2) matrices for spinors*/
    SU2mat(symkq,  nspinor,  false,  timerevkq,  su2kq);
    /* Apply spinors to wfcs */
    su2rotate(nspinor, npwkq,  nspin*nbndskq, su2kq,  wfc_kq->data);
    // in case of time rev, we have to complex conj the wavefunction, 
    // which is done at sandwiching.
    
    ND_int elph_buffer_len = nbndsk*nbndskq*nspin ; 

    
    if (krank ==0)
    {   
        ELPH_OMP_PAR_FOR_SIMD // FIX ME, only master node?
        for (ND_int i =0 ; i<elph_buffer_len; ++i) elph_kq[i] = 0.0 ;
    }
    /* temporart buffer arrays */
    const ELPH_float Imat3[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}; // 3x3 identity matrix 
    const ELPH_float zero_temp[3] = {0.0,0.0,0.0} ;
    int zero3[3] = {0,0,0};
    
    ND_int nmag = dVlocr->dims[1]; /* nmag = 1 for non magnetic and = 2/4 for spin polarized/magnetic systems */
    /* nmag is the spinor dimension of the change in potential */
    if (nmag == 2 && nspin != 2 ) error_msg("Incompatible dvscf ");
    if (nmag == 2 && nspinor != 1 ) error_msg("Incompatible dvscf ");
    ///
    //-------
    
    // create a temporary buffer to store local el-ph mat elements 
    ELPH_cmplx * elph_kq_mn  = calloc(nbndsk*nbndskq, sizeof(ELPH_cmplx));
    /* Now Get Sk in real space*/
    
    wfcinVFFT(wfc_k,symk,tauk,zero3,timerevk,lat_vec, npwk_total, Gk, wfcRspace, commK); 
    // store the wfc in the other buffer
    ND_function(copy, Nd_cmplxS) (&(wfcRspace->Buffer), &(wfcRspace->Buffer_temp));
    
    ND_array(Nd_cmplxS) wfc_Sk_real = wfcRspace->Buffer_temp;
    ND_array(Nd_cmplxS) dv_psi_r    = wfcRspace->Buffer ;

    /* Compute dVpsi in G space and compute the sandwich */        
    for (ND_int iv =0 ; iv <nmodes; ++iv)
    {   
        /*
        Compute dv*psi
        */
        ELPH_cmplx * dv_nu = dVlocr->data + dVlocr->strides[0]*iv ; // (nmodes, nmag, nffts_in_this_cpu)
        for (ND_int is = 0 ; is <nspin; ++is)
        {   
            ELPH_cmplx * psi_r_spin = wfc_Sk_real.data  + is*wfc_Sk_real.strides[0];
            ELPH_cmplx * dV_r = dv_nu +  dVlocr->strides[1]*is ; // only in nspin = 2 case, nmag represent nspin dimension

            ELPH_cmplx * dv_psi = dv_psi_r.data + is*dv_psi_r.strides[0];
            
            dvpsi(nbndsk, nmag, nspinor, dVlocr->strides[1], dV_r, psi_r_spin, dv_psi);
        }
        /***/
        // get back to G space 
        wfcFFT(wfcRspace, Imat3, zero_temp, zero3, false, lat_vec, npwkq_total, Gkq, dVpsiG, commK);
        // Compute the sandwich
        char blas_char = 'C';
        if (timerevkq) blas_char = 'T'; // we perform the time reversal conjugation (if any) here now

        for (ND_int is = 0 ; is <nspin; ++is )
        {   /* msg, nsg -> mn */
            // FIX ME donot forget time reversal here. // wfc_kq->data
            // elph_kq_mn+ (iv*nspin+is)*nbndsk*nbndskq
            ND_function(matmulX, Nd_cmplxS) ('N', blas_char, dVpsiG->data + is*(dVpsiG->strides[0]), \
            wfc_kq->data + is*(wfc_kq->strides[0]), elph_kq_mn , 1.0, 0.0, dVpsiG->strides[1], \
            wfc_kq->strides[1], nbndskq, nbndsk, nbndskq, dVpsiG->strides[1]);
            // reduce the electron phonon matrix elements
            ELPH_cmplx * elph_sum_buf ;
            ELPH_cmplx temp_sum = 0 ; // dummy
            if (krank ==0) elph_sum_buf = elph_kq + (iv*nspin+is)*nbndsk*nbndskq;
            else elph_sum_buf = &temp_sum;
            MPI_Reduce(elph_kq_mn, elph_sum_buf , nbndsk*nbndskq, ELPH_MPI_cmplx, MPI_SUM, 0, commK);
        }
    }
    // Free stuff
    free(elph_kq_mn);
    // Free stuff

    /*free wfc buffers */
    ND_function(destroy,Nd_cmplxS)(wfc_kq);
    ND_function(destroy,Nd_cmplxS)(dVpsiG);
    ND_function(destroy,Nd_floatS)(Gkq);
}






