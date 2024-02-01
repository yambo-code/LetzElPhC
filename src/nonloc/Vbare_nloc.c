/*
This routine computes the non local part to 
bare electron-phonon mat elements.
*/
#include "Vnonloc.h"

/*
**** Yambo stores <K|X^a_{lm}>sqrt(E_l) with out spherical harmonics i.e P^a_l(K) = Sqrt(|E_l|)*F^a_l(K)

**** dV(K,K')/dR_a = \sum_{a,l,m}  -1j*[sigh(E_l)* P^a_l(K) * P^a_l(K')]* [exp(-1j*(K-K')*R_a) * (K-K') ]* Y^l_m(K) * (Y^l_m(K'))^\dagger

*/


/*********** Function bodies ************/

void add_elphNonLocal(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo, \
                    int ikq, int ik, int kqsym, int ksym,  ND_array(Nd_cmplxS) * eigVec, \
                    ELPH_cmplx * elph_kq_mn, const struct ELPH_MPI_Comms * Comm)
{
    /*
    Compute <psi_K| dV_nl/dtau |psi_K'> 
    where K' = k and K = k+q

    In this function variables with Kp/K represent K'/K i.e k/k+q respectively

    (Wavefunctions) wfc_K, wfc_Kp = (ispin, nbnd, nspinor,npw)

    fCoeffs ()  # (natom_types,l*j) of Nd arrays with dim (2l+1, 2l+1, nspinor, nspinor)

    Fkq --> Kb projectors (nl,natom_types,ngx)
    
    SymK, Symkp, symmetry operators on K and K'
    
    atom_pos --> atomic positions in cart coordinates
    
    Gvecs --> ngx Gvectors in cart coordinates

    eigVec : eigen vectors , (nu,atom,3)
    
    elph_kq_mn --> Output <k+q| dV_nl/dt|k>. !! Warning . this must be initialized else Undefined behaviour !
    *** The non local part is added to existing value of elph_kq_mn

    Note : Pass only un-rotated wave functions, this functions uses symmetries internally.
    
    General comments: 
        This functions will never explicitly construct full V(K,K') but loops and computes the 
        matrix elements on fly.
    */
    int mpi_error;

    /*
    First we get the wfcs.
    */

    ND_array(Nd_cmplxS) wfc_K[1], wfc_Kp[1];
    ND_array(Nd_floatS) GvecK[1], GvecKp[1];


    /*initialization and setup */
    ND_function(init,Nd_cmplxS)(wfc_K, (wfcs+ikq)->wfc->rank[0],  (wfcs+ikq)->wfc->dims);
    ND_function(init,Nd_cmplxS)(wfc_Kp,(wfcs+ik)->wfc->rank[0],(wfcs+ik)->wfc->dims);
    ND_function(malloc,Nd_cmplxS)(wfc_K);
    ND_function(malloc,Nd_cmplxS)(wfc_Kp);

    ND_function(init,Nd_floatS)(GvecK,  (wfcs+ikq)->gvec->rank[0], (wfcs+ikq)->gvec->dims);
    ND_function(init,Nd_floatS)(GvecKp, (wfcs+ik)->gvec->rank[0] ,   (wfcs+ik)->gvec->dims);
    ND_function(malloc,Nd_floatS)(GvecK);
    ND_function(malloc,Nd_floatS)(GvecKp);


    // copy wavefunctions
    ND_function(copy, Nd_cmplxS) ((wfcs+ikq)->wfc, wfc_K);
    ND_function(copy, Nd_cmplxS) ((wfcs+ik)->wfc, wfc_Kp);
    // copy gvecs
    ND_function(copy, Nd_floatS) ((wfcs+ikq)->gvec, GvecK);
    ND_function(copy, Nd_floatS) ((wfcs+ik)->gvec, GvecKp);


    ND_array(Nd_floatS)* FK  = (wfcs+ikq)->Fk;
    ND_array(Nd_floatS)* FKp = (wfcs+ik)->Fk;
    
    //printf("Debug-%d \n",1);
    ELPH_float * Kvec  = ND_function(ele,Nd_floatS)(lattice->kpt_iredBZ, nd_idx{ikq,0});
    ELPH_float * Kpvec = ND_function(ele,Nd_floatS)(lattice->kpt_iredBZ, nd_idx{ik,0}); 
    
    ND_array(Nd_floatS) * PP_table = pseudo->PP_table;
    
    const int lmax = pseudo->lmax;
    
    ND_array(Nd_cmplxS) * fCoeff = pseudo->fCoeff;
    ND_array(Nd_floatS)* Fsign = pseudo->Fsign;
    
    const ELPH_float * tauK  = ND_function(ele,Nd_floatS)(lattice->frac_trans, nd_idx{kqsym,0});
    const ELPH_float * tauKp = ND_function(ele,Nd_floatS)(lattice->frac_trans, nd_idx{ksym,0});
    
    ELPH_float * SymK  = ND_function(ele,Nd_floatS)(lattice->sym_mat, nd_idx{kqsym,0,0}) ;
    ELPH_float * SymKp = ND_function(ele,Nd_floatS)(lattice->sym_mat, nd_idx{ksym,0,0}) ;

    const bool timerevK  = lattice->time_rev_array[kqsym];
    const bool timerevKp = lattice->time_rev_array[ksym];

    ND_array(Nd_floatS)* atom_pos     = lattice->atomic_pos;
    const int * atom_type             = lattice->atom_type;

    ND_int nspinor = lattice->nspinor ;
    ND_int natom_types = pseudo->ntype; // number of atomic types
    ND_int natom = atom_pos->dims[0];
    ND_int npwK, npwKp;
    npwK  =  GvecK->dims[0];
    npwKp = GvecKp->dims[0];

    ND_int nl_max = (lmax + 1)*(lmax + 1); // we compute for all (lmax+1)^2 projectors
    

    // From here, real stuff starts 
    /* first rotate G vectors. The data is over written on existing gvecs */
    ELPH_OMP_PAR_FOR_SIMD
    for (ND_int ipw= 0 ; ipw <npwK ; ++ipw)
    {   
        ELPH_float * GPtr    = GvecK->data + 3*ipw;
        ELPH_float tempG[3];
        tempG[0] = Kvec[0]+ GPtr[0] ;
        tempG[1] = Kvec[1]+ GPtr[1] ;
        tempG[2] = Kvec[2]+ GPtr[2] ;
        MatVec3f(SymK, tempG, false, GPtr);
    }
    // for K'  // Pragma omp for 
    ELPH_OMP_PAR_FOR_SIMD
    for (ND_int ipw= 0 ; ipw <npwKp ; ++ipw)
    {   
        ELPH_float * GPtr    = GvecKp->data + 3*ipw;  
        ELPH_float tempG[3];
        tempG[0] = Kpvec[0]+ GPtr[0] ;
        tempG[1] = Kpvec[1]+ GPtr[1] ;
        tempG[2] = Kpvec[2]+ GPtr[2] ;
        MatVec3f(SymKp, tempG, false, GPtr);
    }

    /* ----------- */
    /* Now pre compute Ylm(K) for 0-lmax */
    /*
    The idea behind computing Ylm before is that, the calls for Ylm are reduced
    One could store these for every wavefunction and apply wigner D matrices for rotation. 
    But we always compute on fly as for large k points, these would add to more memory footprint 
    per core with less performance gain
    */
    ELPH_float * YlmK  = malloc( npwK*nl_max*sizeof(ELPH_float)); // (nl_max,npwK)
    ELPH_float * YlmKp = malloc(npwKp*nl_max*sizeof(ELPH_float)); // These are real spherical harmonics

    for (int il =0 ; il<=lmax; ++il)
    {   
        for (int im=0; im<= 2*il; ++im)
        {   //
            int m = im-il ;

            ND_int ilim_idx = (il*il) + im ;
            
            ELPH_float * YlmKtemp  = YlmK  + ilim_idx*npwK;
            ELPH_float * YlmKptemp = YlmKp + ilim_idx*npwKp;

            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int ipw = 0 ; ipw < npwK; ++ipw )
            {   
                ELPH_float * GrotPtr = GvecK->data + 3*ipw ; 
                YlmKtemp[ipw] = Ylm(il, m, GrotPtr);
            }
            
            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int ipw = 0 ; ipw < npwKp; ++ipw )
            {   
                ELPH_float * GrotPtr = GvecKp->data + 3*ipw ; 
                YlmKptemp[ipw] = Ylm(il, m, GrotPtr);
            }
        }
    }
    
    /* ----- */
    /* (ispin, nbnd, nspinor,npw) */
    ND_int nspin, nbndK ;
    nspin  = wfc_K->dims[0] ;
    nbndK  = wfc_K->dims[1] ;
    
    ELPH_cmplx su2K[4]  = {1,0,0,1};
    ELPH_cmplx su2Kp[4] = {1,0,0,1};

    /* Get SU(2) matrices for spinors*/
    SU2mat(SymK,  nspinor,  false,  timerevK,  su2K);
    SU2mat(SymKp, nspinor,  false,  timerevKp, su2Kp);
    
    ND_int nsets = nspin*nbndK;
    ELPH_float kzero[3] = {0,0,0};

    /* Apply spinors and fractional translation to wfcs */
    for (ND_int iset = 0; iset < nsets ; ++iset)
    {   
        ELPH_cmplx * wfc_tmp = wfc_K->data + iset*nspinor*npwK;
        // su2 rotate
        su2rotate(nspinor, npwK,  1, su2K,  wfc_tmp);
        // apply fractional translation
        // note GvecK is k+G so we set kvec to 0
        apply_trans_wfc(tauK, kzero, nspinor, npwK, GvecK->data, wfc_tmp, false);               
    }

    for (ND_int iset = 0; iset < nsets ; ++iset)
    {   
        ELPH_cmplx * wfc_tmp = wfc_Kp->data + iset*nspinor*npwKp;
        // su2 rotate
        su2rotate(nspinor, npwKp, 1, su2Kp, wfc_tmp);
        // apply fractional translation
        // note GvecK is k+G so we set kvec to 0
        apply_trans_wfc(tauKp, kzero, nspinor, npwKp, GvecKp->data, wfc_tmp, false);               
    }
    ///
    
    /* Buffer arrays */
    ND_int nltimesj = PP_table->dims[0];

    /*temporary beta_ia buffers */
    // ((2*lmax+1)*nproj_max,4,nspin*nspinor*nbndK)
    const ND_int bandbuffer_stride = nspin*nbndK*nspinor*4;
    ELPH_cmplx * bandbufferK  = malloc( nltimesj*(2*lmax+1)*bandbuffer_stride*sizeof(ELPH_cmplx));
    ELPH_cmplx * bandbufferKp = malloc( nltimesj*(2*lmax+1)*bandbuffer_stride*sizeof(ELPH_cmplx));
    

    ND_int temp_len = nltimesj*(2*lmax+1);
    ELPH_cmplx * betaK  = malloc(4*temp_len*npwK*sizeof(ELPH_cmplx));// buffer for beta and K*beta (pw,4)
    ELPH_cmplx * betaKp = malloc(4*temp_len*npwKp*sizeof(ELPH_cmplx));// buffer for beta' and K'*beta'
    
    // (natom, 3, nspin, mk, nk+q)
    /* buffer to store elph mat elements in cart coordinates */
    ND_int elph_buffer_len = natom * 3 * nbndK * nbndK * nspin ;
    ELPH_cmplx * elph_buffer;
    if (Comm->commK_rank == 0) elph_buffer  = malloc(elph_buffer_len * sizeof(ELPH_cmplx)); 
    ND_int elph_buffer_stride = 3*nbndK * nbndK * nspin;
    // FIXED TILL HERE
    // zero the buffers
    if (Comm->commK_rank == 0)
    {
        for (ND_int i = 0 ; i<elph_buffer_len ; ++i) elph_buffer[i] = 0.0 ;
    }
    for (ND_int i = 0 ; i<nltimesj*(2*lmax+1)*bandbuffer_stride ; ++i) bandbufferK[i] = 0.0 ; 
    for (ND_int i = 0 ; i<nltimesj*(2*lmax+1)*bandbuffer_stride ; ++i) bandbufferKp[i] = 0.0 ; 
    for (ND_int i = 0 ; i<4*temp_len*npwK ; ++i)  betaK[i] = 0.0 ; 
    for (ND_int i = 0 ; i<4*temp_len*npwKp ; ++i) betaKp[i] = 0.0 ; 

    /* Now compute betas */
    for (ND_int ia = 0; ia < natom ; ++ia)
    {   
        ND_int itype  = atom_type[ia];
        
        ELPH_float * tau = ND_function(ele,Nd_floatS)(atom_pos,nd_idx{ia,0}); // atomic position in cart
        /* First compute betas for each atom */
        
        ND_int idxK =0 ;  // counter for nltimesj*(2*lmax+1) i.e l+m for K
        ND_int idxKp =0 ; // counter for nltimesj*(2*lmax+1) i.e l+m for K'
        //  idxK and idxKp should be same
        for (ND_int lidx = 0 ; lidx < nltimesj ; ++lidx)
        {   
            int l = rint( *ND_function(ele,Nd_floatS)(PP_table, nd_idx{lidx,itype,0})  - 1);
            if (l < 0) continue; // skip fake entries
                
            ELPH_float   Kbsign  = *ND_function(ele,Nd_floatS)(Fsign,nd_idx{lidx,itype});
            ELPH_float * FKtemp  = ND_function(ele,Nd_floatS)(FK,  nd_idx{lidx,itype,0});
            ELPH_float * FKptemp = ND_function(ele,Nd_floatS)(FKp, nd_idx{lidx,itype,0});
            // nltimesj*(2*lmax+1)*4*npw_split; 

            for (ND_int im1 =0 ; im1 <= 2*l ; ++im1)
            {   
                ELPH_float * YlmKtemp  = YlmK  + (l*l+im1)*npwK;
                ELPH_cmplx * betaK_temp = betaK + idxK*4*npwK;

                ELPH_cmplx * restrict betaK0 = betaK_temp;
                ELPH_cmplx * restrict betaK1 = betaK_temp + npwK;
                ELPH_cmplx * restrict betaK2 = betaK_temp + 2*npwK;
                ELPH_cmplx * restrict betaK3 = betaK_temp + 3*npwK;
                /*** WARNING !! DO NOT PARALLELIZE LOOPS except this !! */
                //ELPH_OMP_PAR_FOR_SIMD
                for (ND_int ipw = 0 ; ipw < npwK; ++ipw )
                {   
                    ELPH_float * GrotPtr = GvecK->data + 3*ipw ; 
                    // tau.G
                    ELPH_float tau_dotG = GrotPtr[0]*tau[0] + GrotPtr[1]*tau[1] + GrotPtr[2]*tau[2] ;
                    betaK0[ipw] = FKtemp[ipw]*cexp(-I*2*ELPH_PI*tau_dotG)*YlmKtemp[ipw];
                    betaK1[ipw] = betaK0[ipw]*GrotPtr[0]*Kbsign;
                    betaK2[ipw] = betaK0[ipw]*GrotPtr[1]*Kbsign;
                    betaK3[ipw] = betaK0[ipw]*GrotPtr[2]*Kbsign; // Kbsign is the sign coming from F.T of projectors 
                }
                ++idxK;
            }
            
            for (ND_int im2 =0 ; im2 <= 2*l ; ++im2)
            {   
                ELPH_float * YlmKptemp = YlmKp + (l*l+im2)*npwKp ;
                ELPH_cmplx * betaKp_temp = betaKp + idxKp*4*npwKp;

                ELPH_cmplx * restrict betaKp0 = betaKp_temp;
                ELPH_cmplx * restrict betaKp1 = betaKp_temp + npwKp;
                ELPH_cmplx * restrict betaKp2 = betaKp_temp + 2*npwKp;
                ELPH_cmplx * restrict betaKp3 = betaKp_temp + 3*npwKp;
                //
                /*** WARNING !! DO NOT PARALLELIZE LOOPS except this !! */
                //ELPH_OMP_PAR_FOR_SIMD
                for (ND_int ipw = 0 ; ipw < npwKp; ++ipw )
                {   // K' has -ve sign
                    ELPH_float * GrotPtr = GvecKp->data + 3*ipw; 
                    // tau.G
                    ELPH_float tau_dotG = GrotPtr[0]*tau[0] + GrotPtr[1]*tau[1] + GrotPtr[2]*tau[2] ;
                    betaKp0[ipw] =  FKptemp[ipw]*cexp(I*2*ELPH_PI*tau_dotG)*conj(YlmKptemp[ipw]) ;
                    betaKp1[ipw] = -betaKp0[ipw]*GrotPtr[0]*Kbsign;
                    betaKp2[ipw] = -betaKp0[ipw]*GrotPtr[1]*Kbsign;
                    betaKp3[ipw] = -betaKp0[ipw]*GrotPtr[2]*Kbsign;
                }
                ++idxKp;
            }
        }
        
        char blasK = 'C'; char blasKp = 'T';
        // conjugate if symmetry operation is time rev 
        if (timerevK)  blasK  = 'T';
        if (timerevKp) blasKp = 'C';

        /* matmul with wfcs to get betas*/ 
        // (4, nspin*nbndK*nspinor);  (lj,4,pw)@(nspin,nbnd,spinor,npw)->(lj,4,nspin,nbnd,spinor)
        ND_function(matmulX, Nd_cmplxS) ('N', blasK, betaK, wfc_K->data , bandbufferK, \
                    1.0, 0.0, npwK, npwK, nspin*nspinor*nbndK, 4*idxK, nspin*nspinor*nbndK, npwK);

        ND_function(matmulX, Nd_cmplxS) ('N', blasKp, betaKp, wfc_Kp->data , bandbufferKp, \
                    1.0, 0.0, npwKp, npwKp, nspin*nspinor*nbndK, 4*idxKp, nspin*nspinor*nbndK, npwKp);
        
        if(idxKp != idxK) error_msg("something wrong with number of projectors for K and K' ");

        ND_int reduce_count  = 4*idxK*nspin*nspinor*nbndK; 
        
        
        ND_int max_int_val = ((ND_int)INT_MAX) -10 ;
        // this generally doesn't overflow. but better to check
        if (reduce_count > max_int_val) error_msg("int overflow in MPI_reduce function");
        
        if (Comm->commK_rank == 0) mpi_error = MPI_Reduce(MPI_IN_PLACE, bandbufferKp, reduce_count, ELPH_MPI_cmplx, MPI_SUM, 0, Comm->commK);
        else mpi_error = MPI_Reduce(bandbufferKp, bandbufferKp, reduce_count, ELPH_MPI_cmplx, MPI_SUM, 0, Comm->commK);

        if (Comm->commK_rank == 0) mpi_error = MPI_Reduce(MPI_IN_PLACE, bandbufferK, reduce_count, ELPH_MPI_cmplx, MPI_SUM, 0, Comm->commK);
        else mpi_error = MPI_Reduce(bandbufferK, bandbufferK, reduce_count, ELPH_MPI_cmplx, MPI_SUM, 0, Comm->commK);

        
        // now reduce bandbufferK and bandbufferKp i.e perform the sum
        //Comm->commK

        /* The below section will be run only by master core of each Kpool */
        if (Comm->commK_rank == 0)
        {
            /* Now compute non local contribution to elph matrix  elements*/
            ELPH_cmplx * elph_buffer_temp = elph_buffer + ia*elph_buffer_stride;

            ND_int il_counter = 0 ;
            for (ND_int lidx = 0 ; lidx < nltimesj ; ++lidx)
            {
                int l = rint( *ND_function(ele,Nd_floatS)(PP_table, nd_idx{lidx,itype,0})  - 1);
                //int j = rint( *ND_function(ele,Nd_floatS)(PP_table, nd_idx{lidx,itype,1})); // Careful, this is 2j
        
                if (l < 0) continue; // skip fake entries

                ND_array(Nd_cmplxS) * ifCoeff  = fCoeff + natom_types*lidx + itype ;

                for (ND_int im1 =0 ; im1 <= 2*l ; ++im1)
                {   
                    ELPH_cmplx * betaPsi_K = bandbufferK + (il_counter+im1)*bandbuffer_stride;

                    for (ND_int im2 =0 ; im2 <= 2*l ; ++im2)
                    {   
                        ELPH_cmplx * betaPsi_Kp = bandbufferKp + (il_counter+im2)*bandbuffer_stride;

                        if (nspinor == 1 && im1 != im2) continue;

                        ELPH_cmplx temp_flmm = 1.0 ;
                        ELPH_cmplx * flmm = &temp_flmm ;

                        if (nspinor != 1)
                        {   
                            ND_int * t_strds = ifCoeff->strides;
                            flmm =  ifCoeff->data + t_strds[0]*im1 + t_strds[1]*im2; 
                        }
                        /*** WARNING !! DO NOT PARALLELIZE LOOPS !! */
                        // perform the summation over K,K'
                        // sum_K_K' V_NL = \sum_{sigma,sigma'} npwK^\sigma[0]*npwKp^\sigma'[1:4]*f + f*npwKp^\sigma'[0]*npwK^\sigma[1:4] ;
                        sum_VNL_KKp(betaPsi_K, betaPsi_Kp, flmm, nspin, nbndK , nspinor, elph_buffer_temp); // this is not thread safe
                    }
                }
                // update counter
                il_counter = il_counter + 2*l+1;
            }
        }
        mpi_error = MPI_Barrier(Comm->commK);
    }
    if (Comm->commK_rank == 0) 
    {
        /* Convert to mode basis */
        ND_int nmodes = eigVec->dims[0];
        ND_int elph_stride = nspin*nbndK*nbndK;

        ELPH_cmplx pre_facNL = -2*ELPH_PI*I*ELPH_e2; // this is a prefactor from VKL, but we multiply it here
        /* 2*pi from Gvecs, -I from derivate with ion pos, and e^2 for Ha->Ry. 
        Note that the factor (4*pi)^2/V is included in output of yambo*/

        // !! WARNING : elph_kq_mn and elph_buffer are defined only on master process.
        //(nu,atom,3) , (natom, 3, nspin, mk, nk+q)
        ND_function(matmulX, Nd_cmplxS) ('N', 'N', eigVec->data, elph_buffer, elph_kq_mn, \
                pre_facNL, 1.0, natom*3, elph_stride, elph_stride, nmodes, elph_stride, natom*3);
        free(elph_buffer);
    }

    mpi_error = MPI_Barrier(Comm->commK);

    free(bandbufferK);
    free(bandbufferKp);
    free(betaK);
    free(betaKp);
    free(YlmK);
    free(YlmKp);

    /*free wfc buffers */
    ND_function(destroy,Nd_cmplxS)(wfc_K);
    ND_function(destroy,Nd_cmplxS)(wfc_Kp);
    ND_function(destroy,Nd_floatS)(GvecK);
    ND_function(destroy,Nd_floatS)(GvecKp);


} // end of function


