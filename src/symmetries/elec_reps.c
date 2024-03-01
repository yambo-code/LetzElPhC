#include "../wfc/wfc.h"
#include "symmetries.h"

void electronic_reps(const struct WFC* wfcs, const struct Lattice* lattice,
                     const ELPH_float* Rsym_mat, const ELPH_float* tauR,
                     const bool tim_revR, const ND_int ikBZ,
                     ELPH_cmplx* Dkmn_rep, const struct ELPH_MPI_Comms* Comm)
{
    /*
    This is a function to compute representation matrices for the given symmetry
    and k point Dkmn_rep (<b Rk |U(R) | k,a>): (nspin, b_{bnd} , a_{bnd})
    */

    int mpi_error;

    const ELPH_float* alat = lattice->alat_vec;
    const ELPH_float* blat = lattice->blat_vec;
    // note this has 2*pi in blat//

    // compute the Rk vector and find it in the list of k-points
    ELPH_float Rk_vec[3] = { 0, 0, 0 };
    ELPH_float Rk_tmp[3] = { 0, 0, 0 };
    MatVec3f(Rsym_mat, lattice->kpt_fullBZ + 3 * ikBZ, false, Rk_tmp);
    // convert to crystal coordinates
    MatVec3f(alat, Rk_tmp, true, Rk_vec);

    ND_int nkpts_BZ = lattice->nkpts_BZ;
    // now find the index of the rotated k point
    ND_int iRkBZ = find_kidx_in_list(nkpts_BZ, lattice->kpt_fullBZ_crys, Rk_vec);

    if (iRkBZ < 0)
    {
        error_msg(
            "Rotated k point not found. \
    Either Wrong Phonon symmetry or using non-uniform kgrid");
    }

    /*
    $\langle Rk,b |U(R) | k,a \rangle = \langle Sym_2*k_2,b | \big(U(R) U(Sym_1)
    | k_1,a\rangle\big) $ where Rk = Sym_2*k_2  and k = Sym_1*k_1  (upto to a
    G-vector)
    */

    /*
    Suppose R belongs to the same set of symmetries used to expand kpoints then
    ik1 = ik2. This is true if the point group is symmorphic.
    in case of non-symmorphic, yambo uses a lesser group, so in some cases
    (if phonon symmetries are non-symmorphic)
    it is possible ik1 != ik2.
    We code for the most general case i.e assuming the possibility of ik1 != ik2
    */
    // First get the corresponding iBZ point and symmetry for k and Rk
    const int ik1 = lattice->kmap[ikBZ * 2];
    const int iSym1 = lattice->kmap[ikBZ * 2 + 1];

    const int ik2 = lattice->kmap[iRkBZ * 2];
    const int iSym2 = lattice->kmap[iRkBZ * 2 + 1];

    const ELPH_float* Sym1 = lattice->syms[iSym1].Rmat;
    const ELPH_float* tau1 = lattice->syms[iSym1].tau;
    const bool tr1 = lattice->syms[iSym1].time_rev;

    const ELPH_float* Sym2 = lattice->syms[iSym2].Rmat;
    const ELPH_float* tau2 = lattice->syms[iSym2].tau;
    const bool tr2 = lattice->syms[iSym2].time_rev;

    const ELPH_float* k1_vec = lattice->kpt_iredBZ + 3 * ik1;
    const ELPH_float* k2_vec = lattice->kpt_iredBZ + 3 * ik2;

    ELPH_float tau1_crys[3], tauR_crys[3], tau2_crys[3];
    // tau1, tauR, tau2 in crystal coordinates
    ELPH_float kcrys[3], Rkcrys[3];
    // k and R*k vectors in crystal coordinates

    // compute S1*k1 in crystal coordinates
    MatVec3f(Sym1, k1_vec, false, Rkcrys);
    // we used Rkcrys as tmp storage
    MatVec3f(alat, Rkcrys, true, kcrys);
    // convert to crystal coordinates

    // convert tau to crystal coordinates
    MatVec3f(blat, tau1, true, tau1_crys);
    MatVec3f(blat, tauR, true, tauR_crys);
    MatVec3f(blat, tau2, true, tau2_crys);
    // remove 2*pi from taus
    for (int xi = 0; xi < 3; ++xi)
    {
        tau1_crys[xi] /= (2 * ELPH_PI);
        tauR_crys[xi] /= (2 * ELPH_PI);
        tau2_crys[xi] /= (2 * ELPH_PI);
    }
    // Get the wfcs in iBZ
    //
    const ELPH_cmplx* wfc_k1 = (wfcs + ik1)->wfc;
    const ELPH_cmplx* wfc_k2 = (wfcs + ik2)->wfc;

    const ELPH_float* gvecs_k1 = (wfcs + ik1)->gvec;
    const ELPH_float* gvecs_k2 = (wfcs + ik2)->gvec;

    const ND_int npw_k1_loc = (wfcs + ik1)->npw_loc;
    const ND_int npw_k2_loc = (wfcs + ik2)->npw_loc;

    const ND_int npw_k1_total = (wfcs + ik1)->npw_total;
    const ND_int npw_k2_total = (wfcs + ik2)->npw_total;

    // compute the SU(2) mats for spinor rotation
    ELPH_cmplx SU2_S1[4] = { 1, 0, 0, 1 };
    ELPH_cmplx SU2_R[4] = { 1, 0, 0, 1 };
    ELPH_cmplx SU2_S2[4] = { 1, 0, 0, 1 };

    SU2mat(Sym1, lattice->nspinor, false, tr1, SU2_S1);
    SU2mat(Rsym_mat, lattice->nspinor, false, tim_revR, SU2_R);
    SU2mat(Sym2, lattice->nspinor, false, tr2, SU2_S2);

    // compute the rotated gvecs in crystal coordinates
    ELPH_float* G_S1k1 = calloc(3 * npw_k1_loc, sizeof(ELPH_float)); // S1*k1 gvecs
    CHECK_ALLOC(G_S1k1);

    ELPH_float* G_RS1k1 = calloc(3 * npw_k1_loc, sizeof(ELPH_float)); // R*S1*k1 gvecs
    CHECK_ALLOC(G_RS1k1);

    ELPH_float* G_S2k2 = calloc(3 * npw_k2_loc, sizeof(ELPH_float)); // S2*k2 gvecs
    CHECK_ALLOC(G_S2k2);

    // compute the ulm vec i.e S2K2 + G = R*S1*k1 = > G = R*S1*k1-S2K2
    // C'_G-G0 = C_G. we need to add -G0 = S2*k2-R*S1*k1;

    ELPH_float SymRS1[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 }; // R@S1 matrix
    Gemm3x3f(Rsym_mat, 'N', Sym1, 'N', SymRS1);

    // Compute the ulmvec i.e -G0 = S2*k2-R*S1*k1
    ELPH_float ulm_vec[3] = { 0, 0, 0 }; // (in cart)
    MatVec3f(Sym2, k2_vec, false, ulm_vec);
    // first compute and store S2K2 in ulm_vec
    // we store Rk = R*S1*k1 in Rk_vec (in cart)
    MatVec3f(SymRS1, k1_vec, false, Rk_vec);

    MatVec3f(alat, Rk_vec, true, Rkcrys);
    // R*S1*k1 in crystal units

    for (int i = 0; i < 3; ++i)
    {
        ulm_vec[i] -= Rk_vec[i];
    }

    // rotate G vectors and out put them in crystal coordinates
    // S1*G
    rotateGvecs(gvecs_k1, Sym1, npw_k1_loc, alat, false, true, NULL, G_S1k1);
    // R*S1*G
    rotateGvecs(gvecs_k1, SymRS1, npw_k1_loc, alat, false, true, NULL, G_RS1k1);
    // S2*G
    rotateGvecs(gvecs_k2, Sym2, npw_k2_loc, alat, false, true, ulm_vec, G_S2k2);
    // Now ulm_vec is in crys coordinates as rotateGvecs internally converted it
    // from cart to crys

    // Now we need to rearrange the wavefunctions. We need to find the indices
    // of G_S2k2 in G_RS1k1

    // First get G_RS1k1 and G_S2k2 on the root node
    ND_int* idx_arr = NULL;
    // This is the array that maps G_S2k2 to G_RS1k1; // (only allocated on
    // root)
    ELPH_float* G_S2k2_root_all = NULL;
    // all gvectors collected on root process for sorting
    ELPH_float* G_RS1k1_root_all = NULL;
    // all gvectors collected on root process for sorting

    // mpi buffers
    int* counts = NULL;
    int* disp = NULL;

    int* counts2 = NULL;
    int* disp2 = NULL;

    if (Comm->commK_rank == 0)
    {
        G_RS1k1_root_all = malloc(sizeof(ELPH_float) * npw_k1_total * 3);
        CHECK_ALLOC(G_RS1k1_root_all);

        G_S2k2_root_all = malloc(sizeof(ELPH_float) * npw_k2_total * 3);
        CHECK_ALLOC(G_S2k2_root_all);

        idx_arr = malloc(sizeof(ND_int) * npw_k2_total);
        CHECK_ALLOC(idx_arr);

        counts = malloc(4 * sizeof(int) * Comm->commK_size);
        CHECK_ALLOC(counts);
        
        disp = counts + Comm->commK_size;
        counts2 = counts + 2 * Comm->commK_size;
        disp2 = counts + 3 * Comm->commK_size;

    }

    // collect R*S1*G on root
    int pw_loc_int = 3 * npw_k1_loc;
    mpi_error = MPI_Gather(&pw_loc_int, 1, MPI_INT, counts, 1, MPI_INT, 0, Comm->commK);
    MPI_error_msg(mpi_error);

    if (Comm->commK_rank == 0)
    {
        int disp_tmp = 0;
        for (int i = 0; i < Comm->commK_size; ++i)
        {
            disp[i] = disp_tmp;
            disp_tmp += counts[i];
        }
    }
    // gather
    MPI_Gatherv(G_RS1k1, pw_loc_int, ELPH_MPI_float, G_RS1k1_root_all, counts,
                disp, ELPH_MPI_float, 0, Comm->commK);

    // collect S2*G on root
    pw_loc_int = 3 * npw_k2_loc;
    mpi_error = MPI_Gather(&pw_loc_int, 1, MPI_INT, counts2, 1, MPI_INT, 0,
                           Comm->commK);
    MPI_error_msg(mpi_error);

    if (Comm->commK_rank == 0)
    {
        int disp_tmp = 0;
        for (int i = 0; i < Comm->commK_size; ++i)
        {
            disp2[i] = disp_tmp;
            disp_tmp += counts2[i];
        }
    }
    // gather
    MPI_Gatherv(G_S2k2, pw_loc_int, ELPH_MPI_float, G_S2k2_root_all, counts2,
                disp2, ELPH_MPI_float, 0, Comm->commK);

    // get the indices
    if (Comm->commK_rank == 0)
    {
        find_gvecs_idxs(npw_k2_total, G_S2k2_root_all, npw_k1_total,
                        G_RS1k1_root_all, idx_arr);
        // free some space
        free(G_S2k2_root_all);
        free(G_RS1k1_root_all);

        // divide mpi buffers by 3 so that we can use them for mpi communication
        // routines
        for (int i = 0; i < Comm->commK_size; ++i)
        {
            counts[i] /= 3;
            disp[i] /= 3;
            counts2[i] /= 3;
            disp2[i] /= 3;
        }
    }

    // now we map the wavefunctions
    // allocate memory for rearranged space
    ELPH_cmplx* wfc_k2_root = NULL; // gather ik2 wavefunctions on root for sorting
    ELPH_cmplx* wfc_k2_sort_root = NULL; // store sorted ik2 on root

    if (Comm->commK_rank == 0)
    {
        wfc_k2_root = malloc(sizeof(ELPH_cmplx) * npw_k2_total);
        CHECK_ALLOC(wfc_k2_root);

        wfc_k2_sort_root = malloc(sizeof(ELPH_cmplx) * npw_k1_total);
        CHECK_ALLOC(wfc_k2_sort_root);
    }

    ND_int nsets = lattice->nspin * lattice->nbnds; // nbands * nspin
    ND_int npw_spinor_k1 = lattice->nspinor * npw_k1_loc; // nspinor * npw

    ELPH_cmplx* wfc_RS1k = malloc(nsets * npw_spinor_k1 * sizeof(ELPH_cmplx)); // R*Sym1*k1 wfc
    CHECK_ALLOC(wfc_RS1k);

    ELPH_cmplx* wfc_S2k2 = malloc(nsets * npw_spinor_k1 * sizeof(ELPH_cmplx)); // Sym2*k2 wfc
    CHECK_ALLOC(wfc_S2k2);

    // create a tmp buffer
    ELPH_cmplx* Dkmn_rep_tmp = calloc(lattice->nbnds * lattice->nbnds, sizeof(ELPH_cmplx));
    CHECK_ALLOC(Dkmn_rep_tmp);

    // Now rearrage the wavefunctin, and compute the sandwitch
    for (ND_int iset = 0; iset < nsets; ++iset)
    {
        const ELPH_cmplx* wfc_k2_tmp = wfc_k2 + iset * lattice->nspinor * npw_k2_loc;
        ELPH_cmplx* wfc_S2k2_tmp = wfc_S2k2 + iset * npw_spinor_k1;
        ELPH_cmplx* wfc_RS1k1_tmp = wfc_RS1k + iset * npw_spinor_k1;

        // copy k1 buffer to wfc_RS1k
        memcpy(wfc_RS1k1_tmp, wfc_k1 + iset * npw_spinor_k1,
               sizeof(ELPH_cmplx) * npw_spinor_k1);
        // apply SU(S1)
        su2rotate(lattice->nspinor, npw_k1_loc, 1, SU2_S1, wfc_RS1k1_tmp);
        // apply fractional translation and conjugate if symmetry is time
        // reversal.
        apply_trans_wfc(tau1_crys, kcrys, lattice->nspinor, npw_k1_loc, G_S1k1,
                        wfc_RS1k1_tmp, tr1);
        // note apply_trans_wfc can output conjugate if last parameter is set to
        // true

        // now apply SU(R)
        su2rotate(lattice->nspinor, npw_k1_loc, 1, SU2_R, wfc_RS1k1_tmp);
        // apply fractional translation and conjugate if symmetry is time
        // reversal.
        apply_trans_wfc(tauR_crys, Rkcrys, lattice->nspinor, npw_k1_loc,
                        G_RS1k1, wfc_RS1k1_tmp, tim_revR);
        // note apply_trans_wfc can output conjugate if last parameter is set to
        // true

        // now sort S2K2 wavefunction
        for (ND_int ispinor = 0; ispinor < lattice->nspinor; ++ispinor)
        {
            // gather the wfc on root process
            mpi_error = MPI_Gatherv(
                wfc_k2_tmp + ispinor * npw_k2_loc, npw_k2_loc, ELPH_MPI_cmplx,
                wfc_k2_root, counts2, disp2, ELPH_MPI_cmplx, 0, Comm->commK);
            MPI_error_msg(mpi_error);

            if (Comm->commK_rank == 0)
            {
                // rearrange
                //
                for (ND_int ig = 0; ig < npw_k1_total; ++ig)
                {
                    wfc_k2_sort_root[ig] = 0;
                }
                for (ND_int ii = 0; ii < npw_k2_total; ++ii)
                {
                    ND_int idx_tmp = idx_arr[ii];
                    if (idx_tmp < 0)
                    {
                        continue; // set the missing ones to 0
                    }
                    wfc_k2_sort_root[idx_tmp] = wfc_k2_root[ii];
                }
            }
            // scatter back the wfc to each process
            mpi_error = MPI_Scatterv(wfc_k2_sort_root, counts, disp, ELPH_MPI_cmplx,
                                     wfc_S2k2_tmp + ispinor * npw_k1_loc, npw_k1_loc,
                                     ELPH_MPI_cmplx, 0, Comm->commK);
            MPI_error_msg(mpi_error);
        }

        // apply SU(S2) on k2
        su2rotate(lattice->nspinor, npw_k1_loc, 1, SU2_S2, wfc_S2k2_tmp);
        // conjugate if Sym2 is time reversal
        // also we need to conjugate for sandwich
        // SO no conjugation if Sym2 time revesal and conj if not
        // translate and conjugate(if not time reversal)
        apply_trans_wfc(tau2_crys, Rkcrys, lattice->nspinor, npw_k1_loc,
                        G_RS1k1, wfc_S2k2_tmp, !tr2);
        // note apply_trans_wfc can output conjugate if last parameter is set to
        // true
    }

    // compute the sandwitch
    // <S2*k2 | RS1k1>
    for (ND_int ispin = 0; ispin < lattice->nspin; ++ispin)
    {
        // compute the sandwitch. COnjugation of the left wfc is already done
        // above
        matmul_cmplx(
            'N', 'T', wfc_S2k2 + ispin * lattice->nbnds * npw_spinor_k1,
            wfc_RS1k + ispin * lattice->nbnds * npw_spinor_k1, Dkmn_rep_tmp,
            1.0, 0.0, npw_spinor_k1, npw_spinor_k1, lattice->nbnds,
            lattice->nbnds, lattice->nbnds, npw_spinor_k1);
        // (nba,pw)@ (nbn,pw)^C
        // reduce to root node
        ELPH_cmplx* Dkmn_ptr = NULL;
        if (Comm->commK_rank == 0)
        {
            Dkmn_ptr = Dkmn_rep + ispin * lattice->nbnds * lattice->nbnds;
        }
        MPI_Reduce(Dkmn_rep_tmp, Dkmn_ptr, lattice->nbnds * lattice->nbnds,
                   ELPH_MPI_cmplx, MPI_SUM, 0, Comm->commK);
    }

    free(G_S1k1);
    free(G_S2k2);
    free(G_RS1k1);
    free(wfc_RS1k);
    free(wfc_S2k2);
    free(Dkmn_rep_tmp);

    if (Comm->commK_rank == 0)
    {
        free(wfc_k2_root);
        free(wfc_k2_sort_root);
        free(counts);
        free(idx_arr);
    }
}
