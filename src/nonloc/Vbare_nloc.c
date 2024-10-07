/*
This routine computes the non local part to
bare electron-phonon mat elements.
*/
#include "../common/constants.h"
#include "../common/dtypes.h"
#include "../common/error.h"
#include "../common/numerical_func.h"
#include "../common/omp_pragma_def.h"
#include "../elphC.h"
#include "../wfc/wfc.h"
#include "Vnonloc.h"
#include "fcoeff.h"
#include <complex.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/*
**** Yambo stores <K|X^a_{lm}>sqrt(E_l) with out spherical harmonics i.e
P^a_l(K) = Sqrt(|E_l|)*F^a_l(K)

**** dV(K,K')/dR_a = \sum_{a,l,m}  -1j*[sigh(E_l)* P^a_l(K) * P^a_l(K')]*
[exp(-1j*(K-K')*R_a) * (K-K') ]* Y^l_m(K) * (Y^l_m(K'))^\dagger

*/

/*********** Function bodies ************/

void add_elphNonLocal(struct WFC* wfcs, struct Lattice* lattice,
                      struct Pseudo* pseudo, int ikq, int ik, int kqsym,
                      int ksym, ELPH_cmplx* eigVec, ELPH_cmplx* elph_kq_mn,
                      const struct ELPH_MPI_Comms* Comm)
{
    /*
    Compute <psi_K| dV_nl/dtau |psi_K'>
    where K' = k and K = k+q

    In this function variables with Kp/K represent K'/K i.e k/k+q respectively

    (Wavefunctions) wfc_K, wfc_Kp = (ispin, nbnd, nspinor,npw)

    fCoeffs ()  # (ntype,l*j) of Nd arrays with dim (2l+1, 2l+1, nspinor,
    nspinor)

    Fkq --> Kb projectors (nl,ntype,ngx)

    SymK, Symkp, symmetry operators on K and K'

    atom_pos --> atomic positions in cart coordinates

    Gvecs --> ngx Gvectors in cart coordinates

    eigVec : eigen vectors , (nu,atom,3)

    elph_kq_mn --> Output <k+q| dV_nl/dt|k>. !! Warning . this must be
    initialized else Undefined behaviour !
    *** The non local part is added to existing value of elph_kq_mn

    Note : Pass only un-rotated wave functions, this functions uses symmetries
    internally.

    General comments:
        This functions will never explicitly construct full V(K,K') but loops
    and computes the matrix elements on fly. We do not precompute the beta
    projectors as they can take lot of memory. Neverthless, it is not very slow
    */
    int mpi_error;

    /*
    First we get the wfcs.
    */

    ND_int nspin = lattice->nspin;
    ND_int nbnds = lattice->nbnds;
    ND_int nspinor = lattice->nspinor;
    ND_int ntype = pseudo->ntype; // number of atomic types
    ND_int natom = lattice->natom;
    ND_int npwK = (wfcs + ikq)->npw_loc;
    ND_int npwKp = (wfcs + ik)->npw_loc;

    ELPH_cmplx* wfc_K; // K = S*k+q wfc
    ELPH_cmplx* wfc_Kp; // K^\prime = S*k

    ELPH_float* GvecK; // gvecs for S*k+q wfc
    ELPH_float* GvecKp; // gvecs for S*k wfc

    // wfc structure (nspin,nbnd,nspinor,npw_loc)

    /*initialization and setup */
    ND_int wfcK_buf_size = nspin * nbnds * nspinor * npwK;
    ND_int wfcKp_buf_size = nspin * nbnds * nspinor * npwKp;

    wfc_K = malloc(sizeof(ELPH_cmplx) * wfcK_buf_size);
    CHECK_ALLOC(wfc_K);
    wfc_Kp = malloc(sizeof(ELPH_cmplx) * wfcKp_buf_size);
    CHECK_ALLOC(wfc_Kp);

    GvecK = malloc(sizeof(ELPH_float) * 3 * npwK);
    CHECK_ALLOC(GvecK);
    GvecKp = malloc(sizeof(ELPH_float) * 3 * npwKp);
    CHECK_ALLOC(GvecKp);

    // copy wavefunctions
    memcpy(wfc_K, (wfcs + ikq)->wfc, sizeof(ELPH_cmplx) * wfcK_buf_size);
    memcpy(wfc_Kp, (wfcs + ik)->wfc, sizeof(ELPH_cmplx) * wfcKp_buf_size);
    // copy gvecs
    memcpy(GvecK, (wfcs + ikq)->gvec, sizeof(ELPH_float) * 3 * npwK);
    memcpy(GvecKp, (wfcs + ik)->gvec, sizeof(ELPH_float) * 3 * npwKp);
    // we do not need to copy Fk's as we donot alter them
    const ELPH_float* FK = (wfcs + ikq)->Fk;
    const ELPH_float* FKp = (wfcs + ik)->Fk;

    // printf("Debug-%d \n",1);
    const ELPH_float* Kvec = lattice->kpt_iredBZ + 3 * ikq;
    const ELPH_float* Kpvec = lattice->kpt_iredBZ + 3 * ik;

    const ELPH_float* PP_table = pseudo->PP_table;
    ELPH_cmplx** fCoeff = pseudo->fCoeff;
    const ELPH_float* Fsign = pseudo->Fsign;
    const int lmax = pseudo->lmax;
    const ND_int nl_max = (lmax + 1) * (lmax + 1);
    // we compute for all (lmax+1)^2 projectors

    const struct symmetry* symm_K = lattice->syms + kqsym;
    const struct symmetry* symm_Kp = lattice->syms + ksym;

    const ELPH_float* tauK = symm_K->tau;
    const ELPH_float* tauKp = symm_Kp->tau;

    const ELPH_float* SymK = symm_K->Rmat;
    const ELPH_float* SymKp = symm_Kp->Rmat;

    const bool timerevK = symm_K->time_rev;
    const bool timerevKp = symm_Kp->time_rev;

    const ELPH_float* atom_pos = lattice->atomic_pos;
    const int* atom_type = lattice->atom_type;

    // From here, real stuff starts
    /* first rotate G vectors. The data is over written on existing gvecs */
    ELPH_OMP_PAR_FOR_SIMD
    for (ND_int ipw = 0; ipw < npwK; ++ipw)
    {
        ELPH_float* GPtr = GvecK + 3 * ipw;
        ELPH_float tempG[3];
        tempG[0] = Kvec[0] + GPtr[0];
        tempG[1] = Kvec[1] + GPtr[1];
        tempG[2] = Kvec[2] + GPtr[2];
        MatVec3f(SymK, tempG, false, GPtr);
    }
    // for K'  // Pragma omp for
    ELPH_OMP_PAR_FOR_SIMD
    for (ND_int ipw = 0; ipw < npwKp; ++ipw)
    {
        ELPH_float* GPtr = GvecKp + 3 * ipw;
        ELPH_float tempG[3];
        tempG[0] = Kpvec[0] + GPtr[0];
        tempG[1] = Kpvec[1] + GPtr[1];
        tempG[2] = Kpvec[2] + GPtr[2];
        MatVec3f(SymKp, tempG, false, GPtr);
    }

    /* ----------- */
    /* Now pre compute Ylm(K) for 0-lmax */
    /*
    The idea behind computing Ylm before is that, the calls for Ylm are reduced
    One could store these for every wavefunction and apply wigner D matrices for
    rotation. But we always compute on fly as for large k points, these would
    add to more memory footprint per core with less performance gain
    */
    ELPH_float* YlmK = malloc(npwK * nl_max * sizeof(ELPH_float)); // (nl_max,npwK)
    CHECK_ALLOC(YlmK);
    ELPH_float* YlmKp = malloc(npwKp * nl_max * sizeof(ELPH_float)); // These are real spherical harmonics
    CHECK_ALLOC(YlmKp);

    for (int il = 0; il <= lmax; ++il)
    {
        for (int im = 0; im <= 2 * il; ++im)
        { //
            int m = im - il;

            ND_int ilim_idx = (il * il) + im;

            ELPH_float* YlmKtemp = YlmK + ilim_idx * npwK;
            ELPH_float* YlmKptemp = YlmKp + ilim_idx * npwKp;

            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int ipw = 0; ipw < npwK; ++ipw)
            {
                ELPH_float* GrotPtr = GvecK + 3 * ipw;
                YlmKtemp[ipw] = Ylm(il, m, GrotPtr);
            }

            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int ipw = 0; ipw < npwKp; ++ipw)
            {
                ELPH_float* GrotPtr = GvecKp + 3 * ipw;
                YlmKptemp[ipw] = Ylm(il, m, GrotPtr);
            }
        }
    }

    /* ----- */
    /* (ispin, nbnd, nspinor,npw) */

    ELPH_cmplx su2K[4] = { 1, 0, 0, 1 };
    ELPH_cmplx su2Kp[4] = { 1, 0, 0, 1 };

    /* Get SU(2) matrices for spinors*/
    SU2mat(SymK, nspinor, false, timerevK, su2K);
    SU2mat(SymKp, nspinor, false, timerevKp, su2Kp);

    ND_int nsets = nspin * nbnds;
    ELPH_float kzero[3] = { 0, 0, 0 };

    /* Apply spinors and fractional translation to wfcs */
    for (ND_int iset = 0; iset < nsets; ++iset)
    {
        ELPH_cmplx* wfc_tmp = wfc_K + iset * nspinor * npwK;
        // su2 rotate
        su2rotate(nspinor, npwK, 1, su2K, wfc_tmp);
        // apply fractional translation
        // note GvecK is k+G so we set kvec to 0
        apply_trans_wfc(tauK, kzero, nspinor, npwK, GvecK, wfc_tmp, false);
    }

    for (ND_int iset = 0; iset < nsets; ++iset)
    {
        ELPH_cmplx* wfc_tmp = wfc_Kp + iset * nspinor * npwKp;
        // su2 rotate
        su2rotate(nspinor, npwKp, 1, su2Kp, wfc_tmp);
        // apply fractional translation
        // note GvecK is k+G so we set kvec to 0
        apply_trans_wfc(tauKp, kzero, nspinor, npwKp, GvecKp, wfc_tmp, false);
    }
    ///

    /* Buffer arrays */
    ND_int nltimesj = pseudo->nltimesj;

    /*temporary beta_ia buffers */
    // ((2*lmax+1)*nproj_max,4,nspin*nspinor*nbnds)
    const ND_int bandbuffer_stride = nspin * nbnds * nspinor * 4;
    ELPH_cmplx* bandbufferK = malloc(nltimesj * (2 * lmax + 1) * bandbuffer_stride * sizeof(ELPH_cmplx));
    CHECK_ALLOC(bandbufferK);
    ELPH_cmplx* bandbufferKp = malloc(nltimesj * (2 * lmax + 1) * bandbuffer_stride * sizeof(ELPH_cmplx));
    CHECK_ALLOC(bandbufferKp);

    ND_int temp_len = nltimesj * (2 * lmax + 1);
    ELPH_cmplx* betaK = malloc(4 * temp_len * npwK * sizeof(ELPH_cmplx)); // buffer for beta and K*beta (pw,4)
    CHECK_ALLOC(betaK);
    ELPH_cmplx* betaKp = malloc(4 * temp_len * npwKp * sizeof(ELPH_cmplx)); // buffer for beta' and K'*beta'
    CHECK_ALLOC(betaKp);

    // (natom, 3, nspin, mk, nk+q)
    /* buffer to store elph mat elements in cart coordinates */
    ND_int elph_buffer_len = natom * 3 * nbnds * nbnds * nspin;
    ELPH_cmplx* elph_buffer = NULL;
    if (Comm->commK_rank == 0)
    {
        elph_buffer = malloc(elph_buffer_len * sizeof(ELPH_cmplx));
        CHECK_ALLOC(elph_buffer);
    }
    ND_int elph_buffer_stride = 3 * nbnds * nbnds * nspin;
    // FIXED TILL HERE
    // zero the buffers
    if (Comm->commK_rank == 0)
    {
        for (ND_int i = 0; i < elph_buffer_len; ++i)
        {
            elph_buffer[i] = 0.0;
        }
    }
    for (ND_int i = 0; i < nltimesj * (2 * lmax + 1) * bandbuffer_stride; ++i)
    {
        bandbufferK[i] = 0.0;
    }
    for (ND_int i = 0; i < nltimesj * (2 * lmax + 1) * bandbuffer_stride; ++i)
    {
        bandbufferKp[i] = 0.0;
    }
    for (ND_int i = 0; i < 4 * temp_len * npwK; ++i)
    {
        betaK[i] = 0.0;
    }
    for (ND_int i = 0; i < 4 * temp_len * npwKp; ++i)
    {
        betaKp[i] = 0.0;
    }

    /* Now compute betas */
    for (ND_int ia = 0; ia < natom; ++ia)
    {
        const ND_int itype = atom_type[ia];

        const ELPH_float* tau = atom_pos + 3 * ia; // atom_pos[ia,:]
        /* First compute betas for each atom */

        ND_int idxK = 0; // counter for nltimesj*(2*lmax+1) i.e l+m for K
        ND_int idxKp = 0; // counter for nltimesj*(2*lmax+1) i.e l+m for K'
        // idxK and idxKp should be same
        for (ND_int lidx = 0; lidx < nltimesj; ++lidx)
        {
            const int l = rint(PP_table[ntype * 3 * lidx + itype * 3] - 1); // PP_table[lidx,itype,0]

            if (l < 0)
            {
                continue;
            } // skip fake entries

            const ELPH_float Kbsign = Fsign[lidx * ntype + itype]; // Fsign[lidx,itype]
            const ELPH_float* FKtemp = FK + (lidx * ntype + itype) * npwK; // FK[lidx,itype,:]
            const ELPH_float* FKptemp = FKp + (lidx * ntype + itype) * npwKp; // FKp[lidx,itype,:]

            // nltimesj*(2*lmax+1)*4*npw_split;

            for (ND_int im1 = 0; im1 <= 2 * l; ++im1)
            {
                ELPH_float* YlmKtemp = YlmK + (l * l + im1) * npwK;
                ELPH_cmplx* betaK_temp = betaK + idxK * 4 * npwK;

                ELPH_cmplx* restrict betaK0 = betaK_temp;
                ELPH_cmplx* restrict betaK1 = betaK_temp + npwK;
                ELPH_cmplx* restrict betaK2 = betaK_temp + 2 * npwK;
                ELPH_cmplx* restrict betaK3 = betaK_temp + 3 * npwK;
                /*** WARNING !! DO NOT PARALLELIZE LOOPS except this !! */
                // ELPH_OMP_PAR_FOR_SIMD
                for (ND_int ipw = 0; ipw < npwK; ++ipw)
                {
                    ELPH_float* GrotPtr = GvecK + 3 * ipw;
                    // tau.G
                    ELPH_float tau_dotG = GrotPtr[0] * tau[0] + GrotPtr[1] * tau[1] + GrotPtr[2] * tau[2];
                    betaK0[ipw] = FKtemp[ipw] * cexp(-I * 2 * ELPH_PI * tau_dotG) * YlmKtemp[ipw];
                    betaK1[ipw] = betaK0[ipw] * GrotPtr[0] * Kbsign;
                    betaK2[ipw] = betaK0[ipw] * GrotPtr[1] * Kbsign;
                    betaK3[ipw] = betaK0[ipw] * GrotPtr[2] * Kbsign; // Kbsign is the sign coming from F.T
                                                                     // of projectors
                }
                ++idxK;
            }

            for (ND_int im2 = 0; im2 <= 2 * l; ++im2)
            {
                ELPH_float* YlmKptemp = YlmKp + (l * l + im2) * npwKp;
                ELPH_cmplx* betaKp_temp = betaKp + idxKp * 4 * npwKp;

                ELPH_cmplx* restrict betaKp0 = betaKp_temp;
                ELPH_cmplx* restrict betaKp1 = betaKp_temp + npwKp;
                ELPH_cmplx* restrict betaKp2 = betaKp_temp + 2 * npwKp;
                ELPH_cmplx* restrict betaKp3 = betaKp_temp + 3 * npwKp;
                //
                /*** WARNING !! DO NOT PARALLELIZE LOOPS except this !! */
                // ELPH_OMP_PAR_FOR_SIMD
                for (ND_int ipw = 0; ipw < npwKp; ++ipw)
                { // K' has -ve sign
                    ELPH_float* GrotPtr = GvecKp + 3 * ipw;
                    // tau.G
                    ELPH_float tau_dotG = GrotPtr[0] * tau[0] + GrotPtr[1] * tau[1] + GrotPtr[2] * tau[2];
                    betaKp0[ipw] = FKptemp[ipw] * cexp(I * 2 * ELPH_PI * tau_dotG) * conj(YlmKptemp[ipw]);
                    betaKp1[ipw] = -betaKp0[ipw] * GrotPtr[0] * Kbsign;
                    betaKp2[ipw] = -betaKp0[ipw] * GrotPtr[1] * Kbsign;
                    betaKp3[ipw] = -betaKp0[ipw] * GrotPtr[2] * Kbsign;
                }
                ++idxKp;
            }
        }

        char blasK = 'C';
        char blasKp = 'T';
        // conjugate if symmetry operation is time rev
        if (timerevK)
        {
            blasK = 'T';
        }
        if (timerevKp)
        {
            blasKp = 'C';
        }

        /* matmul with wfcs to get betas*/
        // (4, nspin*nbnds*nspinor);
        // (lj,4,pw)@(nspin,nbnd,spinor,npw)->(lj,4,nspin,nbnd,spinor)
        matmul_cmplx('N', blasK, betaK, wfc_K, bandbufferK, 1.0, 0.0, npwK,
                     npwK, nspin * nspinor * nbnds, 4 * idxK,
                     nspin * nspinor * nbnds, npwK);

        matmul_cmplx('N', blasKp, betaKp, wfc_Kp, bandbufferKp, 1.0, 0.0, npwKp,
                     npwKp, nspin * nspinor * nbnds, 4 * idxKp,
                     nspin * nspinor * nbnds, npwKp);

        if (idxKp != idxK)
        {
            error_msg(
                "something wrong with number of projectors for K and K' ");
        }

        ND_int reduce_count = 4 * idxK * nspin * nspinor * nbnds;

        ND_int max_int_val = ((ND_int)INT_MAX) - 10;
        // this generally doesn't overflow. but better to check
        if (reduce_count > max_int_val)
        {
            error_msg("int overflow in MPI_reduce function");
        }

        if (Comm->commK_rank == 0)
        {
            mpi_error = MPI_Reduce(MPI_IN_PLACE, bandbufferKp, reduce_count,
                                   ELPH_MPI_cmplx, MPI_SUM, 0, Comm->commK);
            MPI_error_msg(mpi_error);
        }
        else
        {
            mpi_error = MPI_Reduce(bandbufferKp, bandbufferKp, reduce_count,
                                   ELPH_MPI_cmplx, MPI_SUM, 0, Comm->commK);
            MPI_error_msg(mpi_error);
        }

        if (Comm->commK_rank == 0)
        {
            mpi_error = MPI_Reduce(MPI_IN_PLACE, bandbufferK, reduce_count,
                                   ELPH_MPI_cmplx, MPI_SUM, 0, Comm->commK);
            MPI_error_msg(mpi_error);
        }
        else
        {
            mpi_error = MPI_Reduce(bandbufferK, bandbufferK, reduce_count,
                                   ELPH_MPI_cmplx, MPI_SUM, 0, Comm->commK);
            MPI_error_msg(mpi_error);
        }

        // now reduce bandbufferK and bandbufferKp i.e perform the sum
        // Comm->commK

        /* The below section will be run only by master core of each Kpool */
        if (Comm->commK_rank == 0)
        {
            /* Now compute non local contribution to elph matrix  elements*/
            ELPH_cmplx* elph_buffer_temp = elph_buffer + ia * elph_buffer_stride;

            ND_int il_counter = 0;
            for (ND_int lidx = 0; lidx < nltimesj; ++lidx)
            {
                int l = rint(PP_table[ntype * 3 * lidx + itype * 3] - 1); // PP_table[lidx,itype,0]
                int j = rint(PP_table[ntype * 3 * lidx + itype * 3 + 1]); // Careful, this is 2j

                if (l < 0)
                {
                    continue;
                } // skip fake entries

                ELPH_float soc_fac = 1.0;
                // incase relativistic pseudo is used in a non-SOC calc, we need to do a j-avg
                if (j != 0 && !lattice->is_soc_present)
                {
                    soc_fac = (j + 1.0) / (4.0 * l + 2.0);
                }

                int two_lp1 = 2 * l + 1;

                ELPH_cmplx* ifCoeff = NULL;
                if (fCoeff)
                {
                    ifCoeff = fCoeff[ntype * lidx + itype];
                }
                for (ND_int im1 = 0; im1 <= 2 * l; ++im1)
                {
                    ELPH_cmplx* betaPsi_K = bandbufferK + (il_counter + im1) * bandbuffer_stride;

                    for (ND_int im2 = 0; im2 <= 2 * l; ++im2)
                    {
                        ELPH_cmplx* betaPsi_Kp = bandbufferKp + (il_counter + im2) * bandbuffer_stride;

                        if (!ifCoeff && im1 != im2)
                        {
                            continue;
                        } // if no soc, diagonal in m's and spins

                        ELPH_cmplx temp_flmm[4] = { soc_fac, 0.0, 0.0, soc_fac };
                        // in case of nspinor = 1, only 1st element is read, for
                        // nspinor =2 it is 2x2 identity.
                        ELPH_cmplx* flmm = temp_flmm;
                        if (ifCoeff)
                        {
                            flmm = ifCoeff + (two_lp1 * im1 + im2) * nspinor * nspinor;
                        }
                        /*** WARNING !! DO NOT PARALLELIZE LOOPS !! */
                        // perform the summation over K,K'
                        // sum_K_K' V_NL = \sum_{sigma,sigma'}
                        // npwK^\sigma[0]*npwKp^\sigma'[1:4]*f +
                        // f*npwKp^\sigma'[0]*npwK^\sigma[1:4] ;
                        sum_VNL_KKp(
                            betaPsi_K, betaPsi_Kp, flmm, nspin, nbnds, nspinor,
                            elph_buffer_temp); // this is not thread safe
                    }
                }
                // update counter
                il_counter = il_counter + 2 * l + 1;
            }
        }
        mpi_error = MPI_Barrier(Comm->commK);
        MPI_error_msg(mpi_error);
    }
    if (Comm->commK_rank == 0)
    {
        /* Convert to mode basis */
        ND_int nmodes = lattice->nmodes;
        ND_int elph_stride = nspin * nbnds * nbnds;

        ELPH_cmplx pre_facNL = -2 * ELPH_PI * I * ELPH_e2; // this is a prefactor from VKL, but we multiply it here
        /* 2*pi from Gvecs, -I from derivate with ion pos, and e^2 for Ha->Ry.
        Note that the factor (4*pi)^2/V is included in output of yambo*/

        // !! WARNING : elph_kq_mn and elph_buffer are defined only on master
        // process.
        //(nu,atom,3) , (natom, 3, nspin, mk, nk+q)
        matmul_cmplx('N', 'N', eigVec, elph_buffer, elph_kq_mn, pre_facNL, 1.0,
                     natom * 3, elph_stride, elph_stride, nmodes, elph_stride,
                     natom * 3);
        free(elph_buffer);
    }

    mpi_error = MPI_Barrier(Comm->commK);
    MPI_error_msg(mpi_error);

    free(bandbufferK);
    free(bandbufferKp);
    free(betaK);
    free(betaKp);
    free(YlmK);
    free(YlmKp);
    free(wfc_K);
    free(wfc_Kp);
    free(GvecK);
    free(GvecKp);

} // end of function
