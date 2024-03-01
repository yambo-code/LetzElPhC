/*
 * This file contains functions that expand electron-phonon matrix
 * elements in q space
 *
 */
#include "symmetries.h"

void elph_q_rotate(const ELPH_cmplx* Dmats_l,
                   const ELPH_cmplx* elph_mat_q,
                   const ELPH_cmplx* Dmats_r,
                   const struct Lattice* lattice,
                   const bool tim_rev,
                   ELPH_cmplx* restrict elph_mat_Sq)
{
    /*
    WARNING !! : Dmats must be for S symmetry and not for S^-1
    * qpt is unrotated q-point in crystal coordinates

    We need to compute gmn(S^-1k)*D_nn' * D^*_mm'

    D_nn'   = D^\dagger_nn'(S,S^-1*k) -- Eq.1
    D_mm'   = D^\dagger_mm'(S, S^-1*k + q) -Eq.2
    elph_mat_Sq and elph_mat_q : (( nmodes, nspin, nbands, nbands))
    Dmats = (nspin, nbands, nbands)
    */
    const ND_int nbnds = lattice->nbnds;
    const ND_int nmodes = lattice->nmodes;
    const ND_int nspin = lattice->nspin;

    ELPH_cmplx* tmp_buffer = calloc(nbnds * nbnds, sizeof(ELPH_cmplx));
    CHECK_ALLOC(tmp_buffer);

    char g_conj = 'T';
    if (tim_rev)
    {
        g_conj = 'C';
    }

    for (ND_int ispin = 0; ispin < nspin; ++ispin)
    {
        const ELPH_cmplx* D_r = Dmats_r + ispin * nbnds * nbnds;
        const ELPH_cmplx* D_l = Dmats_l + ispin * nbnds * nbnds;

        for (ND_int imode = 0; imode < nmodes; ++imode)
        {
            ND_int gshift = nbnds * nbnds * (ispin + nspin * imode);

            ELPH_cmplx* g_Sq = elph_mat_Sq + gshift;
            const ELPH_cmplx* g_q = elph_mat_q + gshift;

            // 1st zero out the buffer else U.B
            for (ND_int i = 0; i < (nbnds * nbnds); ++i)
            {
                g_Sq[i] = 0;
            }
            matmul_cmplx(g_conj, 'C', g_q, D_r, tmp_buffer, 1.0, 0.0,
                         nbnds, nbnds, nbnds, nbnds, nbnds, nbnds);

            matmul_cmplx('T', 'T', tmp_buffer, D_l, g_Sq, 1.0, 0.0,
                         nbnds, nbnds, nbnds, nbnds, nbnds, nbnds);
        }
    }

    free(tmp_buffer);
}
