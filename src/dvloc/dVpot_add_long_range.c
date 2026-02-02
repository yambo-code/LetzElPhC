#include <complex.h>
#include <math.h>
#include <stdbool.h>
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
#include "wfc/wfc.h"

void dV_add_longrange(const ELPH_float* qpt, struct Lattice* lattice,
                      struct Phonon* phonon, const ELPH_float* Zvals,
                      const ELPH_cmplx* eigVec, ELPH_cmplx* dVscf,
                      const ND_int sign, const bool only_induced_part,
                      const ELPH_float EcutRy, const bool* nmags_add,
                      const ELPH_float eta_bare, const ELPH_float eta_induced,
                      MPI_Comm commK)
{
    /**
     * Add or subtracts the long range part of the deformation potential in real
    space

    Input
    qpt       : q-point in crystal coordinates
    if sign < 0: long range part is subtracted to the dvscf
    else       : long range part is added to the advscf
    dVscf : (nmodes, nmag, nfft_loc)
    Assumes dVscf is already initialized, else UB.
    only_induced_part : If only_induced_part is true, then only dipoles and
     quadrupoles are considered else Valance charge( Z_val/r),
    dipoles and quadrupoles are consider.
    nmags_add : is list of bools for which the potential is add or subtraced
    eta_bare is ewald summation parameter for bare
    eta_induced is ewald summation parameter for induced part
    if nmags_add[imag] == true, the potential is added/sub to dvscf[:,imag,:]
     *
     * @param[in] qpt         q-point in crystal coordinates (range: [-1, 1]).
     * @param[in] lattice     Pointer to Lattice structure containing:
     *                        - alat_vec/blat_vec (lattice/reciprocal vectors),
     *                        - atomic_pos (Cartesian positions),
     *                        - fft_dims (FFT grid dimensions),
     *                        - dimension (cutoff type: '3', '2', or '1' for
     * 3D/2D/1D).
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

    if (lattice->dimension == '2' && fabs(qpt[2]) > ELPH_EPS)
    {
        error_msg(
            "In 2D, only qz == 0 points are accepted when interpolating.");
    }

    ELPH_start_clock("dV_longrange");
    //
    const ELPH_float* latvec = lattice->alat_vec;
    const ELPH_float* atom_pos = lattice->atomic_pos;
    const char cutoff = lattice->dimension;
    const ND_int natom = lattice->natom;
    const ND_int nmodes = lattice->nmodes;
    const ELPH_float volume = lattice->volume;
    const ELPH_float* blat = lattice->blat_vec;

    const ND_int nffts_loc =
        lattice->fft_dims[0] * lattice->fft_dims[1] * lattice->nfftz_loc;

    // compute the small fft box for the Cutoff for which we compute
    // dV_longrage(G)
    //
    ND_int Gbox[3];
    get_fft_box(EcutRy, blat, Gbox, commK);
    // make sure we donot exceed fft_dims
    for (int i = 0; i < 3; ++i)
    {
        Gbox[i] = MIN(Gbox[i], lattice->fft_dims[i]);
    }

    ND_int G_vecs_xy, Gxy_shift;
    G_vecs_xy = get_mpi_local_size_idx(Gbox[0] * Gbox[1], &Gxy_shift, commK);
    if (G_vecs_xy < 1)
    {
        error_msg(
            "Some process do not contain G vectors. reduce PW "
            "parallelization.");
    }
    const ND_int size_G_vecs = G_vecs_xy * Gbox[2];
    const ND_int size_VG = nmodes * size_G_vecs;

    // this must be calloc
    ELPH_cmplx* VlocG = calloc(size_VG, sizeof(ELPH_cmplx));
    CHECK_ALLOC(VlocG);

    int* gvecs = malloc(3 * size_G_vecs * sizeof(*gvecs));
    CHECK_ALLOC(gvecs);

    ELPH_float* Gvecs_z = malloc(3 * Gbox[2] * sizeof(*Gvecs_z));
    CHECK_ALLOC(Gvecs_z);

    ELPH_cmplx* VlocGz_cart = calloc(Gbox[2] * nmodes, sizeof(ELPH_cmplx));
    CHECK_ALLOC(VlocGz_cart);

    ELPH_float qcart[3];
    MatVec3f(blat, qpt, false, qcart);
    for (int ii = 0; ii < 3; ++ii)
    {
        qcart[ii] /= (2.0 * ELPH_PI);
    }
    //
    const ELPH_float* Zvals_tmp = Zvals;
    if (only_induced_part)
    {
        Zvals_tmp = NULL;
    }
    const ELPH_float* epslion = phonon->epsilon;
    const ELPH_float* Zeu = phonon->Zborn;
    const ELPH_float* Qpole = phonon->Qpole;
    // Note. No openmp par for, not thread safe for this outer most loop.
    for (ND_int ig = 0; ig < G_vecs_xy; ++ig)
    {
        const ND_int G_glob_idx = Gxy_shift + ig;
        ND_int ix = G_glob_idx / Gbox[1];
        ND_int jy = G_glob_idx % Gbox[1];
        ix = get_miller_idx(ix, Gbox[0]);
        jy = get_miller_idx(jy, Gbox[1]);

        for (ND_int kz = 0; kz < Gbox[2]; ++kz)
        {
            ND_int Gkz = get_miller_idx(kz, Gbox[2]);
            int* g_temp_set = gvecs + Gbox[2] * 3 * ig + kz * 3;
            //
            g_temp_set[0] = ix;
            g_temp_set[1] = jy;
            g_temp_set[2] = Gkz;

            ELPH_float qGtemp[3] = {ix, jy, Gkz};  //|G|
            //
            MatVec3f(blat, qGtemp, false, Gvecs_z + 3 * kz);
            // 2*pi is included here //

            for (int ii = 0; ii < 3; ++ii)
            {
                Gvecs_z[kz * 3 + ii] /= (2.0 * ELPH_PI);
            }
        }
        dVlong_range_kernel(qcart, Gvecs_z, Gbox[2], Zvals_tmp, epslion, Zeu,
                            Qpole, natom, atom_pos, cutoff, volume, latvec[8],
                            EcutRy, eta_bare, eta_induced, VlocGz_cart);
        // VlocG -> (nffs,atom,3),
        //  get it in mode basis @ (nu,atom,3) @(nffs,atom,3)^T
        matmul_cmplx('N', 'T', eigVec, VlocGz_cart, VlocG + Gbox[2] * ig, 1.0,
                     0.0, nmodes, nmodes, size_G_vecs, nmodes, Gbox[2], nmodes);
    }

    free(Gvecs_z);
    free(VlocGz_cart);

    struct ELPH_fft_plan fft_plan;

    wfc_plan(&fft_plan, size_G_vecs, lattice->nfftz_loc, G_vecs_xy, gvecs,
             lattice->fft_dims, FFTW_MEASURE, commK);

    ELPH_cmplx* Vlocr = calloc(nffts_loc, sizeof(*Vlocr));
    CHECK_ALLOC(Vlocr);

    ELPH_float add_sub_factor = 1.0;
    if (sign < 0)
    {
        add_sub_factor = -1.0;
    }
    for (ND_int imode = 0; imode < nmodes; ++imode)
    {
        invfft3D(&fft_plan, 1, VlocG + imode * size_G_vecs, Vlocr, false);
        for (ND_int imag = 0; imag < lattice->nmag; ++imag)
        {
            if (nmags_add[imag])
            {
                ELPH_cmplx* dvscf_tmp =
                    dVscf + (imode * lattice->nmag + imag) * nffts_loc;
                for (ND_int ift = 0; ift < nffts_loc; ++ift)
                {
                    dvscf_tmp[ift] =
                        dvscf_tmp[ift] + add_sub_factor * Vlocr[ift];
                }
            }
        }
    }

    wfc_destroy_plan(&fft_plan);
    free(gvecs);
    free(VlocG);
    ELPH_stop_clock("dV_longrange");
}
