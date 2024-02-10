/*
 * This file contains functions that expand electron-phonon matrix
 * elements in q space
 *
 */
#include "symmetries.h"

void elph_q_rotate(const ELPH_cmplx* elph_mat_q, const struct Lattice* lattice,
                   const ELPH_cmplx* Dmats, const ELPH_float* symS,
                   const bool tim_revS, const ELPH_cmplx fac,
                   const ELPH_float* qpt, ELPH_cmplx* restrict elph_mat_Sq)
{
    /*
    WARNING !! : Dmats must be for S^-1 symmetry and not for S
    * qpt is unrotated q-point in crystal coordinates
    We compute gmn(S^-1k)*D_nn'(S^-1,k) * D^*_mm'(S^-1,k+Sq)
    in einsum notation
    'sij, sxy, six->sjy' Dmats(k), conj(Dmats(k + S*q)), elph_mat_q

    "fac" is any factor that is multiplies to output

    elph_mat_Sq and elph_mat_q : ((k, nmodes, nspin, nbands, nbands))
    Dmats = (k, nspin, nbands, nbands)
    */

    // First, we need to find k+Sq and S^-1*k indices

    const ELPH_float* lat_vec = lattice->alat_vec->data;

    ELPH_float Sqvec[3] = {0, 0, 0};

    ELPH_float Sinv_crys[9];  // lat^T@sym^T
    Gemm3x3f(lat_vec, 'T', symS, 'T', Sinv_crys);

    {  // small scope
        ELPH_float blat[9];
        ELPH_float qcart[3] = {0, 0, 0};
        reciprocal_vecs(lat_vec, blat);     // note this has 2*pi//
        MatVec3f(blat, qpt, false, qcart);  // get qpt in cart units
        for (int xi = 0; xi < 3; ++xi)
        {
            qcart[xi] = qcart[xi] / (2 * ELPH_PI);
        }
        // compute Sq
        MatVec3f(symS, qcart, false, Sqvec);
        // convert Sqvec in crystal
        MatVec3f(lat_vec, Sqvec, true, qcart);
        for (int xi = 0; xi < 3; ++xi)
        {
            Sqvec[xi] = qcart[xi];
        }
    }  // end of scope

    ND_int nkpts_BZ = lattice->kpt_fullBZ_crys->dims[0];

    ND_int nmodes = 3 * lattice->atomic_pos->dims[0];
    ND_int nbnds = lattice->nbnds;
    ND_int nspin = lattice->nspin;

    ND_int elph_stride = nmodes * nspin * nbnds * nbnds;
    ND_int dmat_stride = nspin * nbnds * nbnds;  // nspin, nbands, nbands

    ELPH_cmplx* D_tmp = calloc(nbnds * nbnds, sizeof(ELPH_cmplx));

    for (ND_int ik = 0; ik < nkpts_BZ; ++ik)
    {
        const ELPH_float* kcrys = lattice->kpt_fullBZ_crys->data + 3 * ik;
        const ELPH_float* kcart = lattice->kpt_fullBZ->data + 3 * ik;

        ELPH_float kplusSq[3] = {0, 0, 0};  // k + S*q in crystal units
        ELPH_float Sinvk[3] = {0, 0, 0};    // S^-1*k in crystal units

        for (int xi = 0; xi < 3; ++xi)
        {
            kplusSq[xi] = kcrys[xi] + Sqvec[xi];
        }
        MatVec3f(Sinv_crys, kcart, false, Sinvk);

        ND_int ikplusSq = find_kidx_in_list(
            nkpts_BZ, lattice->kpt_fullBZ_crys->data, kplusSq);
        ND_int iSinvk =
            find_kidx_in_list(nkpts_BZ, lattice->kpt_fullBZ_crys->data, Sinvk);

        // now find the indices of S^-1*k and k + S*q in kpoints
        if (ikplusSq < 0 || iSinvk < 0)
        {
            error_msg(
                "Cannot find k + Sq or S^-1*k in kpoint grid. \
        This is be due to incommensurate q-grid or wrong phonon symmetries.");
        }

        ELPH_cmplx* restrict elph_Sq_ktmp = elph_mat_Sq + ik * elph_stride;
        const ELPH_cmplx* elph_q_ktmp = elph_mat_q + iSinvk * elph_stride;

        const ELPH_cmplx* D_k = Dmats + ik * dmat_stride;
        const ELPH_cmplx* D_k_Sq = Dmats + ikplusSq * dmat_stride;

        for (ND_int ispin = 0; ispin < nspin; ++ispin)
        {
            const ELPH_cmplx* D_k_spin = D_k + ispin * nbnds * nbnds;
            const ELPH_cmplx* D_k_Sq_spin = D_k_Sq + ispin * nbnds * nbnds;

            for (ND_int imode = 0; imode < nmodes; ++imode)
            {
                ND_int shift_tmp =
                    nbnds * nbnds *
                    (ispin + nspin * imode);  // nmodes, nspin, nbands, nbands

                ELPH_cmplx* restrict elph_Sq_tmp = elph_Sq_ktmp + shift_tmp;
                const ELPH_cmplx* elph_q_tmp = elph_q_ktmp + shift_tmp;

                for (ND_int ib = 0; ib < (nbnds * nbnds); ++ib)
                {
                    elph_Sq_tmp[ib] = 0;
                }

                // we need to compute 'ij, xy, ix->jy' Dmats(k), conj(Dmats(k +
                // S*q)), elph_mat_q

                // Dmats(k + S*q)^C @ elph_mat_q^T = D_tmp -> yx,xi->yi
                ND_function(matmulX, Nd_cmplxS)(
                    'C', 'T', D_k_Sq_spin, elph_q_tmp, D_tmp, 1.0, 0.0, nbnds,
                    nbnds, nbnds, nbnds, nbnds, nbnds);

                // Dmats^T@D_tmp^T -> ji,iy->jy
                ND_function(matmulX, Nd_cmplxS)(
                    'T', 'T', D_k_spin, D_tmp, elph_Sq_tmp, fac, 0.0, nbnds,
                    nbnds, nbnds, nbnds, nbnds, nbnds);
                /*
                    if symmetry is time reversal we need to conjugate. This is
                   because we change the application of U from right to left
                   thereby picking up a conjugation (due to anti-linearity) i.e
                    <k + Sq| (U dV U^\dagger | k> ) = {(<k + Sq|U) (dV U^\dagger
                   | k>)}^* where U is time reversal.
                */
                if (tim_revS)
                {
                    for (ND_int ib = 0; ib < (nbnds * nbnds); ++ib)
                    {
                        elph_Sq_tmp[ib] = conj(elph_Sq_tmp[ib]);
                    }
                }
            }
        }
    }

    free(D_tmp);
}
