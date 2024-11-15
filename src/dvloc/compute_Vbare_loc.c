#include <stdbool.h>
#include <stdlib.h>

#include "../common/ELPH_timers.h"
#include "../common/constants.h"
#include "../common/dtypes.h"
#include "../common/error.h"
#include "../common/numerical_func.h"
#include "../common/omp_pragma_def.h"
#include "../elphC.h"
#include "../fft/fft.h"
#include "../wfc/wfc.h"
#include "dvloc.h"

/* Compute the electron phonon matrix elements i.e sandwich for Local part of KS
 * potential */
void elphLocal(const ELPH_float* qpt, struct WFC* wfcs, struct Lattice* lattice,
               int ikq, int ik, int kqsym, int ksym, ELPH_cmplx* dVlocr,
               const struct ELPH_MPI_Comms* Comm, ELPH_cmplx* elph_kq)
{
    /* Computes <S2*k2 | dV_{q}local | S1*k1>
    Note that the inputs kvectors must full the following condition S2*k2 =
    S1*k1 + q + ulmveckq Input : . k+q and k wave functions, (nspin, bands,
    nspinor, ng) . kqvec, kvec are kvectors in (cart units ) for wfc_kq and
    wfc_k . Gkq, Gk  are Gvectors for wfc_kq and wfc_k, .

    ulmveckq is shift that is applyed to k+q vector i.e Skq*kq - S*k - q .

    symkq and symk are symmetric operations to be applied on k+q and k wave
    functions .

    taukq and tauk are fractional translation .

    Gkq and Gk are Gvectors for k+q and k wfs respectively .

    dVlocr -> change in local KS potential i.e dVlocal + dVscf (nmodes, nmag,
    Nx,Ny,Nz-loc) in for one mode . qpt -> qpoint in crystal units output : (nu,
    nspin, mk, nk+q) // electron-phonon mat in mode basis (only the master
    process writes)
    */

    ELPH_start_clock("elph local");

    const ND_int nspin = lattice->nspin;
    const ND_int nbnds = lattice->nbnds;
    const ND_int nspinor = lattice->nspinor;
    const ND_int nmodes = lattice->nmodes;
    const ND_int nmag = lattice->nmag;
    /* nmag = 1 for non magnetic and = 2/4 for spin
    polarized/magnetic systems */

    const ND_int nfft_loc =
        lattice->fft_dims[0] * lattice->fft_dims[1] * lattice->nfftz_loc;

    const ELPH_float* lat_vec = lattice->alat_vec;
    const ELPH_float* blat = lattice->blat_vec;
    // note this has 2*pi//

    // get the wfcs.  (nspin,nbnds,nspinor,npw)
    const ELPH_cmplx* wfc_kq = (wfcs + ikq)->wfc;
    const ELPH_float* Gkq = (wfcs + ikq)->gvec;
    ND_int npwkq = (wfcs + ikq)->npw_loc;
    // no const for the npwkq as sort_pw will change them

    const ELPH_cmplx* wfc_k = (wfcs + ik)->wfc;
    const ELPH_float* Gk = (wfcs + ik)->gvec;
    ND_int npwk = (wfcs + ik)->npw_loc;
    // no const for the npwk as sort_pw will change them

    const ND_int npwkq_total = (wfcs + ikq)->npw_total;
    const ND_int npwk_total = (wfcs + ik)->npw_total;
    /*initialization and setup */

    ELPH_float* gSkq_buf = calloc(3 * npwkq, sizeof(ELPH_float));
    CHECK_ALLOC(gSkq_buf);

    ELPH_float* gSk_buf = calloc(3 * npwk, sizeof(ELPH_float));
    CHECK_ALLOC(gSk_buf);

    /*initialization and setup */
    const ELPH_float* kqvec = lattice->kpt_iredBZ + 3 * ikq;
    const ELPH_float* kvec = lattice->kpt_iredBZ + 3 * ik;

    const ELPH_float* symkq = lattice->syms[kqsym].Rmat;
    const ELPH_float* symk = lattice->syms[ksym].Rmat;

    const ELPH_float* taukq = lattice->syms[kqsym].tau;
    const ELPH_float* tauk = lattice->syms[ksym].tau;

    const bool timerevkq = lattice->syms[kqsym].time_rev;
    const bool timerevk = lattice->syms[ksym].time_rev;

    ELPH_float ulmveckq[3];
    // ulmveckq is shift that is applyed to k+q vector
    // i.e -(Skq*kq - S*k - q)
    ELPH_float tempSkq[3] = {0, 0, 0};  // S2*k2
    ELPH_float tempSk[3] = {0, 0, 0};   // S1*k1
    ELPH_float taukq_crys[3], tauk_crys[3];
    // fractional translation in crystal coordinates
    ELPH_float kvecSkq[3], kvecSk[3];
    // Sk + q and Sk in crystal coordinates

    MatVec3f(symkq, kqvec, false, tempSkq);
    MatVec3f(symk, kvec, false, tempSk);

    MatVec3f(lat_vec, tempSk, true, kvecSk);
    for (int xi = 0; xi < 3; ++xi)
    {
        kvecSkq[xi] = kvecSk[xi] + qpt[xi];
    }
    // note that tempSkq = kvecSkq + G

    for (int xi = 0; xi < 3; ++xi)
    {
        tempSkq[xi] -= tempSk[xi];
    }

    // convert qpt to cartisian coord
    // tempSkq is Skq*kq - S*k
    MatVec3f(blat, qpt, false, tempSk);
    // store qpt( in cart units ) in tempSk
    for (int xi = 0; xi < 3; ++xi)
    {
        tempSk[xi] = tempSk[xi] / (2 * ELPH_PI);
    }

    // convert tau to crystal coordinates
    MatVec3f(blat, taukq, true, taukq_crys);
    for (int xi = 0; xi < 3; ++xi)
    {
        taukq_crys[xi] = taukq_crys[xi] / (2 * ELPH_PI);
    }

    MatVec3f(blat, tauk, true, tauk_crys);
    for (int xi = 0; xi < 3; ++xi)
    {
        tauk_crys[xi] = tauk_crys[xi] / (2 * ELPH_PI);
    }

    // Skq+G = Sk + q => G = Sk+q-Skq
    for (int xi = 0; xi < 3; ++xi)
    {
        ulmveckq[xi] = (tempSkq[xi] - tempSk[xi]);
    }
    // note in the above we store -G0 as we add -G0 to all gvectors

    rotateGvecs(Gkq, symkq, npwkq, lat_vec, false, true, ulmveckq, gSkq_buf);
    rotateGvecs(Gk, symk, npwk, lat_vec, false, true, NULL, gSk_buf);

    // rotate the wave function in spin space
    ELPH_cmplx su2kq[4] = {1, 0, 0, 1};
    ELPH_cmplx su2k[4] = {1, 0, 0, 1};
    /* Get SU(2) matrices for spinors*/
    SU2mat(symkq, nspinor, false, timerevkq, su2kq);
    SU2mat(symk, nspinor, false, timerevk, su2k);

    // in case of time rev, we have to complex conj the wavefunction,
    // which is done at sandwiching.

    /* scatter the wfc and gvecs */
    int* gvecSGkq;       // gvecs after rearragement // k+q
    int* gvecSGk;        // k
    ELPH_cmplx* wfcSkq;  // wfc after rearragement // k+q
    ELPH_cmplx* wfcSk;   // k

    ND_int nGxySkq, nGxySk;

    // Note : npwkq and npwk are overwritten by number of gvecs in gvecSGkq and
    // gvecSGk respectively. Sort_pw internally allocates buffer that must be
    // freed out side of the function
    Sort_pw(npwkq_total, npwkq, lattice->fft_dims, gSkq_buf, wfc_kq,
            nspin * nspinor * nbnds, &npwkq, &nGxySkq, &gvecSGkq, &wfcSkq,
            Comm->commK);

    Sort_pw(npwk_total, npwk, lattice->fft_dims, gSk_buf, wfc_k,
            nspin * nspinor * nbnds, &npwk, &nGxySk, &gvecSGk, &wfcSk,
            Comm->commK);

    // Note that we need to free gvecSGkq, gvecSGk, wfcSkq, wfcSk in this
    // function when no longer need

    free(gSkq_buf);
    free(gSk_buf);

    ND_int max_npw = ((npwkq > npwk) ? npwkq : npwk);
    /* Gvecs are stored in ints but we need it in float type, so we
    allocate a gvec array to store gvecs in floats (needed for translation
    routine) */
    ELPH_float* gvecs_float_tmp = malloc(sizeof(ELPH_float) * 3 * max_npw);
    CHECK_ALLOC(gvecs_float_tmp);

    for (ND_int ipw = 0; ipw < (3 * npwkq); ++ipw)
    {
        gvecs_float_tmp[ipw] = gvecSGkq[ipw];
    }

    /* Apply spinors and fractional translation to k+q wfc */
    for (ND_int iset = 0; iset < (nspin * nbnds); ++iset)
    {
        ELPH_cmplx* wfcSkq_tmp = wfcSkq + iset * nspinor * npwkq;
        // apply su2 rotation
        su2rotate(nspinor, npwkq, 1, su2kq, wfcSkq_tmp);
        // apply fractional translation
        apply_trans_wfc(taukq_crys, kvecSkq, nspinor, npwkq, gvecs_float_tmp,
                        wfcSkq_tmp, false);
    }

    ELPH_cmplx* wfcSk_r =
        malloc(sizeof(ELPH_cmplx) * nspin * nbnds * nspinor * nfft_loc);
    CHECK_ALLOC(wfcSk_r);

    struct ELPH_fft_plan fft_plan;

    // create plan for Sk
    wfc_plan(&fft_plan, npwk, lattice->nfftz_loc, nGxySk, gvecSGk,
             lattice->fft_dims, FFTW_MEASURE, Comm->commK);

    for (ND_int ipw = 0; ipw < (3 * npwk); ++ipw)
    {
        gvecs_float_tmp[ipw] = gvecSGk[ipw];
    }

    /* FFT Sk wave function to get it in real space */
    for (ND_int iset = 0; iset < (nspin * nbnds); ++iset)
    {
        // rotate spinor wfc
        ELPH_cmplx* wfcSk_tmp = wfcSk + iset * nspinor * npwk;
        su2rotate(nspinor, npwk, 1, su2k, wfcSk_tmp);

        // apply fractional translation
        apply_trans_wfc(tauk_crys, kvecSk, nspinor, npwk, gvecs_float_tmp,
                        wfcSk_tmp, false);

        ELPH_cmplx* wfcSkr_tmp = wfcSk_r + iset * nspinor * nfft_loc;

        // conjugate (only incase of timerevk is true) and then perform inverse
        // FFT.
        invfft3D(&fft_plan, nspinor, wfcSk_tmp, wfcSkr_tmp, timerevk);
        // Note that incase of time reversal symmetry,
        // inverse fft of the conjugate of input is performed.
        // Inputs are never altered in invfft3D function !
    }

    // free some buffers
    wfc_destroy_plan(&fft_plan);
    free(gvecSGk);
    free(wfcSk);
    free(gvecs_float_tmp);

    // create plan for dvSpi
    wfc_plan(&fft_plan, npwkq, lattice->nfftz_loc, nGxySkq, gvecSGkq,
             lattice->fft_dims, FFTW_MEASURE, Comm->commK);

    ELPH_cmplx* dVpsiG =
        malloc(sizeof(ELPH_cmplx) * nspin * nbnds * nspinor * npwkq);
    CHECK_ALLOC(dVpsiG);

    if (Comm->commK_rank == 0)
    {
        ND_int elph_buffer_len = nmodes * nbnds * nbnds * nspin;

        ELPH_OMP_PAR_FOR_SIMD
        for (ND_int i = 0; i < elph_buffer_len; ++i)
        {
            elph_kq[i] = 0.0;
        }
    }

    /* nmag is the spinor dimension of the change in potential */
    if (nmag == 2 && nspin != 2)
    {
        error_msg("Incompatible dvscf ");
    }
    if (nmag == 2 && nspinor != 1)
    {
        error_msg("Incompatible dvscf ");
    }
    ///
    //-------

    // create a temporary buffer to store local el-ph mat elements
    ELPH_cmplx* elph_kq_mn = calloc(nbnds * nbnds, sizeof(ELPH_cmplx));
    CHECK_ALLOC(elph_kq_mn);
    /* Now Get Sk in real space*/

    /* Compute dVpsi in G space and compute the sandwich */
    for (ND_int iv = 0; iv < nmodes; ++iv)
    {
        /*
        Compute dv*psi
        */
        ELPH_cmplx* dv_nu = dVlocr + nmag * nfft_loc * iv;
        // (nmodes, nmag, nfft_loc)
        for (ND_int is = 0; is < nspin; ++is)
        {
            ELPH_cmplx* psi_r_spin = wfcSk_r + is * nbnds * nspinor * nfft_loc;
            ELPH_cmplx* dV_r = dv_nu + nfft_loc * is;
            // only in nspin = 2 case, nmag represent nspin dimension

            for (ND_int ibnd = 0; ibnd < nbnds; ++ibnd)
            {
                /* compute the convolution FFT(dV(r)*psi(r))*/
                ELPH_cmplx* dV_psiG_ptr = dVpsiG +
                                          is * nbnds * nspinor * npwkq +
                                          ibnd * nspinor * npwkq;

                fft_convolution3D(&fft_plan, nspinor, nmag, dV_r,
                                  psi_r_spin + ibnd * nspinor * nfft_loc,
                                  dV_psiG_ptr, false);
            }

            // Compute the sandwich
            char blas_char = 'C';
            if (timerevkq)
            {
                blas_char = 'T';
                // we perform the time reversal conjugation for K+q wfc (if any)
                // here now
            }

            /* msg, nsg -> mn */
            // FIX ME donot forget time reversal here. // wfc_kq
            // elph_kq_mn+ (iv*nspin+is)*nbnds*nbnds
            matmul_cmplx('N', blas_char, dVpsiG + is * nbnds * nspinor * npwkq,
                         wfcSkq + is * (nbnds * nspinor * npwkq), elph_kq_mn,
                         1.0, 0.0, nspinor * npwkq, nspinor * npwkq, nbnds,
                         nbnds, nbnds, nspinor * npwkq);
            // reduce the electron phonon matrix elements
            ELPH_cmplx* elph_sum_buf;
            ELPH_cmplx temp_sum = 0;  // dummy
            if (Comm->commK_rank == 0)
            {
                elph_sum_buf = elph_kq + (iv * nspin + is) * nbnds * nbnds;
            }
            else
            {
                elph_sum_buf = &temp_sum;
            }
            int mpi_error = MPI_Reduce(elph_kq_mn, elph_sum_buf, nbnds * nbnds,
                                       ELPH_MPI_cmplx, MPI_SUM, 0, Comm->commK);
            MPI_error_msg(mpi_error);
        }
    }
    // Free stuff
    free(elph_kq_mn);
    wfc_destroy_plan(&fft_plan);
    free(dVpsiG);
    free(wfcSk_r);

    free(gvecSGkq);
    free(wfcSkq);
    ELPH_stop_clock("elph local");
    // Free stuff
}
