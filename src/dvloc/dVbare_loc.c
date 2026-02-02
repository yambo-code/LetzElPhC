#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "common/ELPH_timers.h"
#include "common/constants.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/omp_pragma_def.h"
#include "common/parallel.h"
#include "dvloc.h"
#include "elphC.h"
#include "fft/fft.h"

void dVlocq(const ELPH_float* qpt, struct Lattice* lattice,
            struct Pseudo* pseudo, const ELPH_cmplx* eigVec, ELPH_cmplx* Vlocr,
            MPI_Comm commK)
{
    /**
     * @brief Computes the change in local bare potential in real space due to
     * atomic displacements.
     *
     * @details
     * This function calculates the derivative of the local pseudopotential
     * (dVloc/dtau) in real space for a given q-point, transforming from
     * reciprocal space (G-vectors) using an inverse FFT. The result is
     * projected onto phonon modes using the provided eigenvectors.
     *
     * ### Mathematical Formulation:
     * - The local potential derivative in reciprocal space is:
     *   \[
     *   \frac{dV_{loc}(q+G)}{d\tau} = -i (q+G) \cdot V_{loc}(|q+G|) \cdot
     * e^{-i(q+G)\cdot\tau}
     *   \]
     * - After interpolation and mode projection, an inverse FFT is applied to
     * obtain the real-space result.
     *
     * @param[in] qpt         q-point in crystal coordinates (range: [-1, 1]).
     * @param[in] lattice     Pointer to Lattice structure containing:
     *                        - alat_vec/blat_vec (lattice/reciprocal vectors),
     *                        - atomic_pos (Cartesian positions),
     *                        - atom_type (atomic species indices),
     *                        - fft_dims (FFT grid dimensions),
     *                        - dimension (cutoff type: '3', '2', or '1' for
     * 3D/2D/1D).
     * @param[in] pseudo      Pointer to Pseudo structure with pseudopotential
     * tables:
     *                        - vloc_table (interpolation data for Vloc(G)),
     *                        - loc_pseudo[itype].Zval (valence charge per
     * species).
     * @param[in] eigVec      Phonon eigenvectors (shape: [nmodes, 3*natom]),
     *                        used to project onto modes.
     * @param[out] Vlocr      Output array (shape: [nmodes, nffts_loc]) storing
     * the real-space dVloc/dtau for each mode.
     * @param[in] commK       MPI communicator for parallel FFT grid
     * distribution.
     *
     * @note
     * - The q-point must be in fractional coordinates (range [-1, 1]).
     * - Memory for `Vlocr` must be pre-allocated by the caller.
     * - Uses spline interpolation for Vloc(G) and includes a long-range Coulomb
     * correction.
     * - For 2D systems (`cutoff == '2'`), an analytic cutoff is applied to
     * remove z-periodicity.
     *
     * @warning
     * - Terminates with `error_msg` if:
     *   - No G-vectors are assigned to a process (`G_vecs_xy < 1`).
     *   - |q+G| exceeds the interpolation table bounds.
     */

    ELPH_start_clock("dV_bare");
    const ELPH_float* latvec = lattice->alat_vec;
    const ND_int ntype = pseudo->ntype;
    const ELPH_float* atom_pos = lattice->atomic_pos;
    const int* atom_type = lattice->atom_type;
    const char cutoff = lattice->dimension;
    const ND_int natom = lattice->natom;
    const ND_int nmodes = lattice->nmodes;
    const ND_int FFTx = lattice->fft_dims[0];
    const ND_int FFTy = lattice->fft_dims[1];
    const ND_int FFTz = lattice->fft_dims[2];
    const ELPH_float volume = lattice->volume;
    const ELPH_float* blat = lattice->blat_vec;
    const ELPH_float eta = 1.0;
    /* This is spread of nuclear charge */

    const ND_int nffts_loc = FFTx * FFTy * lattice->nfftz_loc;

    for (ND_int i = 0; i < (nffts_loc * nmodes); ++i)
    {
        Vlocr[i] = 0;
    }

    ND_int G_vecs_xy, Gxy_shift;

    G_vecs_xy = get_mpi_local_size_idx(FFTx * FFTy, &Gxy_shift, commK);
    if (G_vecs_xy < 1)
    {
        error_msg("Some process do not contain G vectors");
    }
    const ND_int size_G_vecs = G_vecs_xy * lattice->fft_dims[2];
    const ND_int size_VG = nmodes * size_G_vecs;

    // this must be calloc
    ELPH_cmplx* VlocG = calloc(size_VG, sizeof(ELPH_cmplx));
    CHECK_ALLOC(VlocG);

    int* gvecs = malloc(3 * size_G_vecs * sizeof(int));
    CHECK_ALLOC(gvecs);

    const ND_int npts_co = pseudo->vloc_table->npts_co;
    const ELPH_float* g_co = pseudo->vloc_table->g_co;
    const ELPH_float* vlocg = pseudo->vloc_table->vlocg;
    const ELPH_float* vploc_co = pseudo->vloc_table->vploc_co;

    const ELPH_float dg =
        pseudo->vloc_table->dg;  // spacing between two g points in g_co

    ELPH_float* VlocGtype = malloc(sizeof(ELPH_float) * (ntype));
    CHECK_ALLOC(VlocGtype);

    ELPH_cmplx* VlocGz_cart = malloc(FFTz * nmodes * sizeof(ELPH_cmplx));
    CHECK_ALLOC(VlocGz_cart);

    // Note. No openmp par for, not thread safe for this outer most loop.
    for (ND_int ig = 0; ig < G_vecs_xy; ++ig)
    {
        const ND_int fft_glob_idx = Gxy_shift + ig;
        const ND_int ix = fft_glob_idx / FFTy;
        const ND_int jy = fft_glob_idx % FFTy;

        for (ND_int kz = 0; kz < FFTz; ++kz)
        {
            int* g_temp_set = gvecs + FFTz * 3 * ig + kz * 3;
            g_temp_set[0] = ix;
            g_temp_set[1] = jy;
            g_temp_set[2] = kz;

            ELPH_float qGtemp[3] = {
                get_miller_idx(ix, FFTx) + qpt[0],
                get_miller_idx(jy, FFTy) + qpt[1],
                get_miller_idx(kz, FFTz) + qpt[2]};  // | q + G|

            ELPH_float qGtempCart[3];  // in cartisian coordinate

            MatVec3f(blat, qGtemp, false,
                     qGtempCart);  // 2*pi is included here //

            ELPH_float qGnorm = sqrt(qGtempCart[0] * qGtempCart[0] +
                                     qGtempCart[1] * qGtempCart[1] +
                                     qGtempCart[2] * qGtempCart[2]);

            if (qGnorm > g_co[npts_co - 1])
            {
                error_msg("Vlocg Interpolation failed due to out of bounds");
            }
            ND_int gidx = floor(qGnorm / dg);

            ELPH_cmplx* tmp_ptr = VlocGz_cart + kz * nmodes;

            ELPH_float cutoff_fac = 1;
            /* using analytic cutoff which works only when z periodicity is
             * broken */
            if (cutoff == '2')
            {
                ELPH_float Gp_norm = sqrt(qGtempCart[0] * qGtempCart[0] +
                                          qGtempCart[1] * qGtempCart[1]);
                ELPH_float qGp = latvec[8] * Gp_norm;
                ELPH_float cos_sin_fac = cos(qGtempCart[2] * latvec[8] * 0.5);
                if (Gp_norm > ELPH_EPS)
                {
                    cos_sin_fac =
                        cos_sin_fac - sin(qGtempCart[2] * latvec[8] * 0.5) *
                                          qGtempCart[2] / Gp_norm;
                }
                cutoff_fac = 1 - exp(-qGp * 0.5) * cos_sin_fac;
            }

            for (ND_int itype = 0; itype < ntype; ++itype)
            {
                VlocGtype[itype] = spline_interpolate(
                    qGnorm, gidx, g_co, vlocg + itype * npts_co,
                    vploc_co + itype * npts_co);
                // add long range part back
                if (qGnorm >= ELPH_EPS)
                {
                    VlocGtype[itype] -= (4 * ELPH_PI * ELPH_e2) *
                                        ((pseudo->loc_pseudo)[itype].Zval) *
                                        exp(-qGnorm * qGnorm * 0.25 / eta) *
                                        cutoff_fac / (qGnorm * qGnorm * volume);
                }
                // printf("%f \n",VlocGtype[itype]);
            }
            ELPH_OMP_PAR_SIMD
            for (ND_int ia = 0; ia < natom; ++ia)
            {
                const ND_int itype = atom_type[ia];
                const ELPH_float* pos_temp = atom_pos + 3 * ia;
                ELPH_float qdottau = qGtempCart[0] * pos_temp[0] +
                                     qGtempCart[1] * pos_temp[1] +
                                     qGtempCart[2] * pos_temp[2];
                ELPH_cmplx factor = -I * cexp(-I * qdottau);
                factor *= VlocGtype[itype];
                tmp_ptr[ia * 3] = factor * qGtempCart[0];
                tmp_ptr[ia * 3 + 1] = factor * qGtempCart[1];
                tmp_ptr[ia * 3 + 2] = factor * qGtempCart[2];
            }
        }
        // VlocG -> (nffs,atom,3),
        //  get it in mode basis @ (nu,atom,3) @(nffs,atom,3)^T
        matmul_cmplx('N', 'T', eigVec, VlocGz_cart, VlocG + FFTz * ig, 1.0, 0.0,
                     nmodes, nmodes, size_G_vecs, nmodes, FFTz, nmodes);
    }

    free(VlocGtype);
    free(VlocGz_cart);

    struct ELPH_fft_plan fft_plan;

    wfc_plan(&fft_plan, size_G_vecs, lattice->nfftz_loc, G_vecs_xy, gvecs,
             lattice->fft_dims, FFTW_MEASURE, commK);

    invfft3D(&fft_plan, nmodes, VlocG, Vlocr, false);

    wfc_destroy_plan(&fft_plan);
    free(gvecs);
    free(VlocG);
    ELPH_stop_clock("dV_bare");
}
