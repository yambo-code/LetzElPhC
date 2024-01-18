#include "symmetries.h"
#include "../wfc/wfc.h"


void electronic_reps(const struct WFC * wfcs, const struct Lattice * lattice, \
    const ELPH_float * Rsym_mat,  const ELPH_float * tauR, \
    const bool tim_revR, const ND_int ikBZ, ELPH_cmplx * Dkmn_rep, MPI_Comm commK)
{
    /*
    This is a function to compute representation matrices for the given symmetry and k point
    Dkmn_rep (<b Rk |U(R) | k,a>): (nspin, b_{bnd} , a_{bnd}) 
    */

    int my_rank, Comm_size, mpi_error;
    mpi_error = MPI_Comm_size(commK, &Comm_size);
    mpi_error = MPI_Comm_rank(commK, &my_rank);

    // compute the Rk vector and find it in the list of k-points
    ELPH_float Rk_vec[3] = {0,0,0};
    ELPH_float Rk_tmp[3] = {0,0,0};
    MatVec3f(Rsym_mat, lattice->kpt_fullBZ->data + 3*ikBZ, false, Rk_tmp);
    // convert to crystal coordinates
    MatVec3f(lattice->alat_vec->data,Rk_tmp, true, Rk_vec);
    
    // now find the index of the rotated k point
    ND_int iRkBZ = -1;
    // find the index
    for (ND_int ik = 0; ik < lattice->kpt_fullBZ_crys->dims[0]; ++ik)
    {
        ELPH_float * ik_vec_tmp = lattice->kpt_fullBZ_crys->data + 3*ik ;
        ELPH_float sum = 0;
        for (int i = 0; i < 3 ; ++i)
        {   
            ELPH_float diff_tmp = ik_vec_tmp[i]-Rk_vec[i];
            diff_tmp = diff_tmp-rint(diff_tmp);
            sum += diff_tmp*diff_tmp;
        }
        sum = sqrt(sum);
        if (sum < ELPH_EPS)
        {
            iRkBZ = ik;
            break;
        }
    }

    if (iRkBZ < 0) error_msg("Rotated k point not found. Either Wrong Phonon symmetry or using non-uniform kgrid");

    /*
    $\langle Rk,b |U(R) | k,a \rangle = \langle Sym_2*k_2,b | \big(U(R) U(Sym_1) | k_1,a\rangle\big) $
    where Rk = Sym_2*k_2  and k = Sym_1*k_1  (upto to a G-vector)
    */

    /*
    Suppose R belongs to the same set of symmetries used to expand kpoints then
    ik1 = ik2. This is true if the point group is symmorphic. 
    in case of non-symmorphic, yambo uses a lesser group, so in some cases (if phonon symmetries are non-symmorphic) 
    it is possible ik1 != ik2. 
    We code for the most general case i.e assuming the possibility of ik1 != ik2
    */
    // First get the corresponding iBZ point and symmetry for k and Rk
    const int ik1      = *(lattice->kmap->data + ikBZ*2)      ;
    const int iSym1     = *(lattice->kmap->data + ikBZ*2 + 1)  ;

    const int ik2      = *(lattice->kmap->data + iRkBZ*2)      ;
    const int iSym2     = *(lattice->kmap->data + iRkBZ*2 + 1)  ;

    const ELPH_float * Sym1  = lattice->sym_mat->data + 9*iSym1;
    const ELPH_float * Sym2  = lattice->sym_mat->data + 9*iSym2;
    
    const ELPH_float * k1_vec = lattice->kpt_iredBZ->data + 3*ik1;
    const ELPH_float * k2_vec = lattice->kpt_iredBZ->data + 3*ik2;

    const bool tr1 = lattice->time_rev_array[iSym1];
    const bool tr2 = lattice->time_rev_array[iSym2];

    const ELPH_float * tau1  = ND_function(ele,Nd_floatS)(lattice->frac_trans, nd_idx{iSym1,0});
    const ELPH_float * tau2  = ND_function(ele,Nd_floatS)(lattice->frac_trans, nd_idx{iSym2,0});

    ELPH_float tau1_crys[3], tauR_crys[3], tau2_crys[3]; // tau1, tauR, tau2 in crystal coordinates 
    ELPH_float kcrys[3], Rkcrys[3]; // k and R*k vectors in crystal coordinates


    // compute S1*k1 in crystal coordinates
    MatVec3f(Sym1, k1_vec, false, Rkcrys); // we used Rkcrys as tmp storage
    MatVec3f(lattice->alat_vec->data, Rkcrys, true, kcrys);  // convert to crystal coordinates
    // start of small scope 1
    {
        ELPH_float blat[9];
        reciprocal_vecs(lattice->alat_vec->data, blat); // note this has 2*pi//
        for (int xi = 0 ; xi<9 ; ++xi ) blat[xi] = blat[xi]/(2*ELPH_PI);
        
        // convert tau to crystal coordinates
        MatVec3f(blat,  tau1,  true,  tau1_crys);
        MatVec3f(blat,  tauR,  true,  tauR_crys);
        MatVec3f(blat,  tau2,  true,  tau2_crys);
    } // end of small scope 1

    // Get the wfcs in iBZ
    //
    const ELPH_cmplx * wfc_k1 = (wfcs+ik1)->wfc->data ;
    const ELPH_cmplx * wfc_k2 = (wfcs+ik2)->wfc->data ;

    const ELPH_float * gvecs_k1 = (wfcs+ik1)->gvec->data ;
    const ELPH_float * gvecs_k2 = (wfcs+ik2)->gvec->data ;

    const ND_int npw_k1_loc = (wfcs+ik1)->npw_loc ;
    const ND_int npw_k2_loc = (wfcs+ik2)->npw_loc ;

    const ND_int npw_k1_total = (wfcs+ik1)->npw_total ;
    const ND_int npw_k2_total = (wfcs+ik2)->npw_total ;


    // compute the SU(2) mats for spinor rotation
    ELPH_cmplx SU2_S1[4]={1,0,0,1};
    ELPH_cmplx SU2_R[4]={1,0,0,1}; 
    ELPH_cmplx SU2_S2[4]={1,0,0,1}; 
    
    SU2mat(Sym1, lattice->nspinor, false, tr1, SU2_S1); 
    SU2mat(Rsym_mat, lattice->nspinor, false, tim_revR, SU2_R); 
    SU2mat(Sym2, lattice->nspinor, false, tr2, SU2_S2); 

    // compute the rotated gvecs in crystal coordinates
    ELPH_float * G_S1k1    = calloc(3*npw_k1_loc, sizeof(ELPH_float)); // S1*k1 gvecs
    ELPH_float * G_RS1k1   = calloc(3*npw_k1_loc, sizeof(ELPH_float)); // R*S1*k1 gvecs
    ELPH_float * G_S2k2    = calloc(3*npw_k2_loc, sizeof(ELPH_float)); // S2*k2 gvecs

    // compute the ulm vec i.e S2K2 + G = R*S1*k1 = > G = R*S1*k1-S2K2 
    // C'_G-G0 = C_G. we need to add -G0 = S2*k2-R*S1*k1;

    ELPH_float SymRS1[9] ={0,0,0, 0,0,0, 0,0,0}; // R@S1 matrix
    Gemm3x3f(Rsym_mat, 'N', Sym1,  'N', SymRS1);

    // Compute the ulmvec i.e -G0 = S2*k2-R*S1*k1
    ELPH_float ulm_vec[3] = {0,0,0}; // (in cart)
    MatVec3f(Sym2, k2_vec, false, ulm_vec); // first compute and store S2K2 in ulm_vec
    // we store Rk = R*S1*k1 in Rk_vec (in cart)
    MatVec3f(SymRS1, k1_vec, false, Rk_vec);
    MatVec3f(lattice->alat_vec->data, Rk_vec, true, Rkcrys); // R*S1*k1 in crystal units
    
    for (int i = 0; i<3; ++i) ulm_vec[i] -= Rk_vec[i] ;

    // rotate G vectors and out put them in crystal coordinates
    // S1*G
    rotateGvecs(gvecs_k1,   Sym1, npw_k1_loc, lattice->alat_vec->data, false, true, NULL,  G_S1k1);
    // R*S1*G
    rotateGvecs(gvecs_k1, SymRS1, npw_k1_loc, lattice->alat_vec->data, false, true, NULL, G_RS1k1);
    // S2*G
    rotateGvecs(gvecs_k2,   Sym2, npw_k2_loc, lattice->alat_vec->data, false, true, ulm_vec, G_S2k2);
    // Now ulm_vec is in crys coordinates as rotateGvecs internally converted it from cart to crys

    // Now we need to rearrange the wavefunctions. We need to find the indices of G_S2k2 in G_RS1k1

    // First get G_RS1k1 and G_S2k2 on the root node
    ND_int * idx_arr = NULL ; // This is the array that maps G_S2k2 to G_RS1k1; (only allocated on root)
    ELPH_float * G_S2k2_root_all = NULL; // all gvectors collected on root process for sorting
    ELPH_float * G_RS1k1_root_all = NULL; // all gvectors collected on root process for sorting

    // mpi buffers
    int * counts = NULL; 
    int * disp = NULL;

    int * counts2 = NULL; 
    int * disp2 = NULL;

    if (my_rank == 0)
    {   
        G_RS1k1_root_all = malloc(sizeof(ELPH_float)*npw_k1_total*3);
        G_S2k2_root_all  = malloc(sizeof(ELPH_float)*npw_k2_total*3);
        idx_arr          = malloc(sizeof(ND_int)*npw_k2_total);
        counts           = malloc(4*sizeof(int)*Comm_size);
        disp             = counts +   Comm_size ;
        counts2          = counts + 2*Comm_size ;
        disp2            = counts + 3*Comm_size ;
        

        if (counts   == NULL) error_msg("Failed to allocate comm array");
        if (G_S2k2_root_all == NULL)   error_msg("Failed to allocate gvec123 array");
        if (G_RS1k1_root_all == NULL)   error_msg("Failed to allocate gvec1 array");
        if (idx_arr == NULL)     error_msg("Failed to allocate indices array");
    }
    
    // collect R*S1*G on root
    int pw_loc_int = 3*npw_k1_loc;
    mpi_error = MPI_Gather(&pw_loc_int, 1, MPI_INT, counts, 1, MPI_INT, 0, commK);
    if (my_rank == 0)
    {
        int disp_tmp = 0;
        for (int i = 0 ; i<Comm_size; ++i)
        {
            disp[i] = disp_tmp;
            disp_tmp += counts[i];
        }
    }
    // gather 
    MPI_Gatherv(G_RS1k1, pw_loc_int, ELPH_MPI_float, \
        G_RS1k1_root_all, counts, disp, ELPH_MPI_float, 0, commK);

    // collect S2*G on root
    pw_loc_int = 3*npw_k2_loc;
    mpi_error = MPI_Gather(&pw_loc_int, 1, MPI_INT, counts2, 1, MPI_INT, 0, commK);
    if (my_rank == 0)
    {
        int disp_tmp = 0;
        for (int i = 0 ; i<Comm_size; ++i)
        {
            disp2[i] = disp_tmp;
            disp_tmp += counts2[i];
        }
    }
    // gather 
    MPI_Gatherv(G_S2k2, pw_loc_int, ELPH_MPI_float, \
        G_S2k2_root_all, counts2, disp2, ELPH_MPI_float, 0, commK);
    
    // get the indices
    if (my_rank == 0)
    {   
        find_gvecs_idxs(npw_k2_total, G_S2k2_root_all, npw_k1_total, G_RS1k1_root_all, idx_arr);
        // free some space
        free(G_S2k2_root_all);
        free(G_RS1k1_root_all);

        // divide mpi buffers by 3 so that we can use them for mpi communication routines
        for (int i = 0 ; i<Comm_size; ++i)
        {
            counts[i]   /= 3 ;
            disp[i]     /= 3 ;     
            counts2[i]  /= 3 ;  
            disp2[i]    /= 3 ;  
        }
    }
    
    // now we map the wavefunctions
    // allocate memory for rearranged space
    ELPH_cmplx * wfc_k2_root = NULL ; // gather ik2 wavefunctions on root for sorting
    ELPH_cmplx * wfc_k2_sort_root = NULL ; // store sorted ik2 on root

    if (my_rank == 0)
    {   
        wfc_k2_root      = malloc(sizeof(ELPH_cmplx)*npw_k2_total);
        wfc_k2_sort_root = malloc(sizeof(ELPH_cmplx)*npw_k1_total);
    }
    
    ND_int nsets         = lattice->nspin*lattice->nbnds; // nbands * nspin
    ND_int npw_spinor_k1 = lattice->nspinor*npw_k1_loc; // nspinor * npw

    ELPH_cmplx *  wfc_RS1k = malloc(nsets*npw_spinor_k1*sizeof(ELPH_cmplx)); // R*Sym1*k1 wfc
    if (wfc_RS1k == NULL) error_msg("Allocation of sorted local buffer RS1k failed");
    ELPH_cmplx *  wfc_S2k2 = malloc(nsets*npw_spinor_k1*sizeof(ELPH_cmplx)); // Sym2*k2 wfc
    if (wfc_S2k2 == NULL) error_msg("Allocation of sorted local buffer S2k2 failed");
    
    // create a tmp buffer
    ELPH_cmplx * Dkmn_rep_tmp = calloc(lattice->nbnds*lattice->nbnds,sizeof(ELPH_cmplx));
    if (Dkmn_rep_tmp == NULL) error_msg("Allocation of local buffer Dkmn failed");

    // Now rearrage the wavefunctin, and compute the sandwitch
    for (ND_int iset =0; iset<nsets; ++iset)
    {   
        const ELPH_cmplx * wfc_k2_tmp = wfc_k2 + iset*lattice->nspinor*npw_k2_loc;
        ELPH_cmplx * wfc_S2k2_tmp  = wfc_S2k2 + iset*npw_spinor_k1;
        ELPH_cmplx * wfc_RS1k1_tmp = wfc_RS1k + iset*npw_spinor_k1;

        // copy k1 buffer to wfc_RS1k
        memcpy(wfc_RS1k1_tmp, wfc_k1+iset*npw_spinor_k1, sizeof(ELPH_cmplx)*npw_spinor_k1);
        // apply SU(S1)
        su2rotate(lattice->nspinor, npw_k1_loc, 1, SU2_S1, wfc_RS1k1_tmp);
        // apply fractional translation and conjugate if symmetry is time reversal.
        apply_trans_wfc(tau1_crys, kcrys, lattice->nspinor, npw_k1_loc, G_S1k1, wfc_RS1k1_tmp, tr1);
        // note apply_trans_wfc can output conjugate if last parameter is set to true

        // now apply SU(R)
        su2rotate(lattice->nspinor, npw_k1_loc, 1, SU2_R, wfc_RS1k1_tmp);
        // apply fractional translation and conjugate if symmetry is time reversal.
        apply_trans_wfc(tauR_crys, Rkcrys, lattice->nspinor, npw_k1_loc, G_RS1k1, wfc_RS1k1_tmp, tim_revR);
        // note apply_trans_wfc can output conjugate if last parameter is set to true
        
        // now sort S2K2 wavefunction
        for (ND_int ispinor=0; ispinor<lattice->nspinor; ++ispinor)
        {
            // gather the wfc on root process
            mpi_error = MPI_Gatherv(wfc_k2_tmp + ispinor*npw_k2_loc, npw_k2_loc, ELPH_MPI_cmplx, \
                            wfc_k2_root, counts2, disp2, ELPH_MPI_cmplx, 0, commK);
            if(my_rank == 0)
            {
                //rearrange
                // 
                for (ND_int ig = 0; ig<npw_k1_total; ++ig) wfc_k2_sort_root[ig] = 0;
                for (ND_int ii = 0; ii < npw_k2_total; ++ii)
                {
                    ND_int idx_tmp = idx_arr[ii];
                    if (idx_tmp < 0) continue; // set the missing ones to 0
                    wfc_k2_sort_root[idx_tmp] = wfc_k2_root[ii];
                }
            }
            // scatter back the wfc to each process
            mpi_error = MPI_Scatterv(wfc_k2_sort_root, counts, disp, ELPH_MPI_cmplx, \
                        wfc_S2k2_tmp + ispinor*npw_k1_loc, npw_k1_loc, ELPH_MPI_cmplx, 0, commK);
        }

        // apply SU(S2) on k2
        su2rotate(lattice->nspinor, npw_k1_loc, 1, SU2_S2, wfc_S2k2_tmp);
        // conjugate if Sym2 is time reversal
        // also we need to conjugate for sandwich
        // SO no conjugation if Sym2 time revesal and conj if not
        // translate and conjugate(if not time reversal)
        apply_trans_wfc(tau2_crys, Rkcrys, lattice->nspinor, npw_k1_loc, G_RS1k1, wfc_S2k2_tmp, !tr2);
        // note apply_trans_wfc can output conjugate if last parameter is set to true
    }

    // compute the sandwitch
    // <S2*k2 | RS1k1>
    for (ND_int ispin = 0; ispin < lattice->nspin; ++ispin)
    {
        // compute the sandwitch. COnjugation of the left wfc is already done above
        ND_function(matmulX, Nd_cmplxS) ('N', 'T', wfc_S2k2+ispin*lattice->nbnds*npw_spinor_k1, \
                wfc_RS1k + ispin*lattice->nbnds*npw_spinor_k1, Dkmn_rep_tmp, 1.0, 0.0, npw_spinor_k1, npw_spinor_k1, \
                lattice->nbnds, lattice->nbnds, lattice->nbnds, npw_spinor_k1);
        // (nba,pw)@ (nbn,pw)^C
        // reduce to root node
        ELPH_cmplx * Dkmn_ptr = NULL;
        if (my_rank == 0) Dkmn_ptr = Dkmn_rep + ispin*lattice->nbnds*lattice->nbnds;
        MPI_Reduce(Dkmn_rep_tmp, Dkmn_ptr, lattice->nbnds*lattice->nbnds, ELPH_MPI_cmplx, MPI_SUM, 0, commK);
    }

    free(G_S1k1);
    free(G_S2k2);
    free(G_RS1k1);
    free(wfc_RS1k);
    free(wfc_S2k2);
    free(Dkmn_rep_tmp);

    if (my_rank == 0)
    {
        free(wfc_k2_root);
        free(wfc_k2_sort_root);
        free(counts);
        free(idx_arr);  
    }

}














