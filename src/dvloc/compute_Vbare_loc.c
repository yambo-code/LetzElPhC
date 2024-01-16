#include "dvloc.h"

/* Compute the electron phonon matrix elements i.e sandwich for Local part of KS potential */
void elphLocal(const ELPH_float * qpt, struct WFC * wfcs, struct Lattice * lattice, \
                int ikq, int ik, int kqsym, int ksym, ND_array(Nd_cmplxS) * dVlocr, \
                MPI_Comm commK, ELPH_cmplx * elph_kq)
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
            . dVlocr -> change in local KS potential i.e dVlocal + dVscf (nmodes, nmag, Nx,Ny,Nz-loc) in for one mode
            . qpt -> qpoint in crystal units
    output : (nu, nspin, mk, nk+q) // electron-phonon mat in mode basis (only the master process writes)
    */
    
    int mpi_error;
    int krank;
    mpi_error = MPI_Comm_rank(commK, &krank);
    /*
    First we get the wfcs. */

    ND_array(Nd_cmplxS) * wfc_kq = (wfcs+ikq)->wfc;;
    ND_array(Nd_floatS) * Gkq    = (wfcs+ikq)->gvec;;

    ND_array(Nd_cmplxS) * wfc_k = (wfcs+ik)->wfc;
    ND_array(Nd_floatS) * Gk    = (wfcs+ik)->gvec;

    ND_int npwkq_total = (wfcs+ikq)->npw_total;
    ND_int npwk_total =  (wfcs+ik)->npw_total;
    /*initialization and setup */

    ELPH_float * gSkq_buf = calloc(3*Gkq->dims[0],sizeof(ELPH_float));
    ELPH_float * gSk_buf =  calloc(3*Gk->dims[0] ,sizeof(ELPH_float));

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
    npwkq   = Gkq->dims[0] ;
    npwk    = Gk->dims[0] ;
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
    // Skq+G = Sk + q => G = Sk+q-Skq
    for (int xi = 0 ; xi<3 ; ++xi )  ulmveckq[xi] = (tempSkq[xi]-tempSk[xi]);

    rotateGvecs(Gkq->data, symkq, npwkq, lat_vec, false, true, ulmveckq, gSkq_buf);
    rotateGvecs(Gk->data, symk, npwk, lat_vec, false, true, NULL, gSk_buf);

    // rotate the wave function in spin space 
    ELPH_cmplx su2kq[4]  = {1,0,0,1};
    ELPH_cmplx su2k[4]  = {1,0,0,1};
    /* Get SU(2) matrices for spinors*/
    SU2mat(symkq,  nspinor,  false,  timerevkq,  su2kq);
    SU2mat(symk,   nspinor,  false,  timerevk,   su2k);
    
    // in case of time rev, we have to complex conj the wavefunction, 
    // which is done at sandwiching.
    
    /* scatter the wfc and gvecs */
    int * gvecSGkq ; // gvecs after rearragement // k+q
    int * gvecSGk ;  // k
    ELPH_cmplx * wfcSkq; // wfc after rearragement // k+q
    ELPH_cmplx * wfcSk; // k

    ND_int nGxySkq, nGxySk;

    // npwkq and npwk are overwritten by number of gvecs in gvecSGkq and gvecSGk respectively.
    Sort_pw(npwkq_total, npwkq, lattice->fft_dims, gSkq_buf , wfc_kq->data, \
                nspin*nspinor*nbndskq, &npwkq, &nGxySkq, &gvecSGkq, &wfcSkq, commK);

    Sort_pw(npwk_total, npwk, lattice->fft_dims, gSk_buf , wfc_k->data, \
                nspin*nspinor*nbndsk, &npwk, &nGxySk, &gvecSGk, &wfcSk, commK);
                
    /* Apply spinors to wfcs */
    su2rotate(nspinor, npwkq,  nspin*nbndskq, su2kq,  wfcSkq);

    free(gSkq_buf);
    free(gSk_buf);

    ND_array(Nd_cmplxS) wfcSk_r[1];
    ND_array(Nd_cmplxS) dVpsiG[1];
    // s,b,sp,nx,ny,nz
    ND_function(init, Nd_cmplxS) (wfcSk_r, 6, nd_idx{nspin,nbndsk,nspinor, \
                                lattice->fft_dims[0], lattice->fft_dims[1],lattice->nfftz_loc});

    ND_function(malloc, Nd_cmplxS) (wfcSk_r);

    struct ELPH_fft_plan fft_plan; 

    // create plan for Sk
    wfc_plan(&fft_plan, npwk, lattice->nfftz_loc, nGxySk, gvecSGk, lattice->fft_dims, FFTW_MEASURE, commK);

    /* FFT Sk wave function in real space */
    for (ND_int iset = 0 ; iset < (nspin*nbndsk) ; ++iset )
    {
        // rotate spinor wfc
        ELPH_cmplx * wfcSk_tmp = wfcSk + iset*nspinor*npwk;
        su2rotate(nspinor, npwk,  1, su2k, wfcSk_tmp);

        ELPH_cmplx * wfcSkr_tmp = wfcSk_r->data + iset*wfcSk_r->strides[1] ;

        invfft3D(&fft_plan, nspinor, wfcSk_tmp, wfcSkr_tmp, timerevk);
        
    }

    // free some buffers
    wfc_destroy_plan(&fft_plan);
    free(gvecSGk) ;
    free(wfcSk);


    // create plan for dvSpi
    wfc_plan(&fft_plan, npwkq, lattice->nfftz_loc, nGxySkq, gvecSGkq, lattice->fft_dims, FFTW_MEASURE, commK);

    ND_function(init, Nd_cmplxS) (dVpsiG, 4, nd_idx{nspin,nbndsk,nspinor, npwkq});
    ND_function(malloc, Nd_cmplxS) (dVpsiG);

    ND_int elph_buffer_len = nbndsk*nbndskq*nspin ; 
    if (krank ==0)
    {   
        ELPH_OMP_PAR_FOR_SIMD // FIX ME, only master node?
        for (ND_int i =0 ; i<elph_buffer_len; ++i) elph_kq[i] = 0.0 ;
    }
    
    ND_int nmag = dVlocr->dims[1]; /* nmag = 1 for non magnetic and = 2/4 for spin polarized/magnetic systems */
    /* nmag is the spinor dimension of the change in potential */
    if (nmag == 2 && nspin != 2 ) error_msg("Incompatible dvscf ");
    if (nmag == 2 && nspinor != 1 ) error_msg("Incompatible dvscf ");
    ///
    //-------
    
    // create a temporary buffer to store local el-ph mat elements 
    ELPH_cmplx * elph_kq_mn  = calloc(nbndsk*nbndskq, sizeof(ELPH_cmplx));
    /* Now Get Sk in real space*/


    /* Compute dVpsi in G space and compute the sandwich */        
    for (ND_int iv =0 ; iv <nmodes; ++iv)
    {   
        /*
        Compute dv*psi
        */
        ELPH_cmplx * dv_nu = dVlocr->data + dVlocr->strides[0]*iv ; // (nmodes, nmag, nffts_in_this_cpu)
        for (ND_int is = 0 ; is <nspin; ++is)
        {   
            ELPH_cmplx * psi_r_spin = wfcSk_r->data  + is*(wfcSk_r->strides[0]);
            ELPH_cmplx * dV_r = dv_nu +  dVlocr->strides[1]*is ; // only in nspin = 2 case, nmag represent nspin dimension
            
            // ND_function(init, Nd_cmplxS) (dVpsiG, 4, nd_idx{nspin,nbndsk,nspinor, npwkq});
            for (ND_int ibnd = 0 ; ibnd < nbndsk; ++ibnd)
            {
                /* compute the convolution FFT(dV(r)*psi(r))*/
                ELPH_cmplx * dV_psiG_ptr = dVpsiG->data + is*dVpsiG->strides[0] + ibnd*dVpsiG->strides[1];
                fft_convolution3D(&fft_plan, nspinor, nmag, dV_r, psi_r_spin+ibnd*wfcSk_r->strides[1], dV_psiG_ptr, false);
            }
            
            // Compute the sandwich
            char blas_char = 'C';
            if (timerevkq) blas_char = 'T'; // we perform the time reversal conjugation (if any) here now
        
            /* msg, nsg -> mn */
            // FIX ME donot forget time reversal here. // wfc_kq->data
            // elph_kq_mn+ (iv*nspin+is)*nbndsk*nbndskq
            ND_function(matmulX, Nd_cmplxS) ('N', blas_char, dVpsiG->data + is*dVpsiG->strides[0], \
            wfcSkq + is*(dVpsiG->strides[0]), elph_kq_mn , 1.0, 0.0, dVpsiG->strides[1], \
            dVpsiG->strides[1], nbndskq, nbndsk, nbndskq, dVpsiG->strides[1]);
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
    wfc_destroy_plan(&fft_plan);
    ND_function(destroy, Nd_cmplxS) (dVpsiG);
    ND_function(destroy, Nd_cmplxS) (wfcSk_r);

    free(gvecSGkq) ;
    free(wfcSkq);
    // Free stuff

}






