// adds or removed e^{ikr} phase to fft grid. Inplace operation
//
//
#include <complex.h>
#include <stdlib.h>

#include "common/constants.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/omp_pragma_def.h"
#include "dvloc.h"
#include "elphC.h"

void multiply_eikr(ELPH_cmplx* pot_grid, const ELPH_float* qpt_crys,
                   const struct Lattice* lattice, const ND_int nsets,
                   const ND_int sign)
{
    // sign < 0. e^{-ikr} is multiplied
    // else e^{ikr} is multiplied
    // qpt in crystal coordinates
    //
    const ND_int nffts =
        lattice->fft_dims[0] * lattice->fft_dims[1] * lattice->nfftz_loc;

    ELPH_cmplx factor = 2 * I * ELPH_PI;
    if (sign < 0)
    {
        factor = -factor;
    }

    ELPH_cmplx* eikr = malloc(sizeof(*eikr) * nffts);
    CHECK_ALLOC(eikr);

    ND_int FFTyz = lattice->fft_dims[1] * lattice->nfftz_loc;
    // (Nx, Ny, Nz_loc)  (Ny*Nz_loc, Nz_loc, 1)
    ELPH_OMP_PAR_FOR_SIMD
    for (ND_int ig = 0; ig < nffts; ++ig)
    {
        ND_int ix = ig / FFTyz;
        ND_int iy = (ig % FFTyz) / lattice->nfftz_loc;
        ND_int iz =
            lattice->nfftz_loc_shift + (ig % FFTyz) % lattice->nfftz_loc;

        ELPH_float ri = ((ELPH_float)ix) / lattice->fft_dims[0];
        ELPH_float rj = ((ELPH_float)iy) / lattice->fft_dims[1];
        ELPH_float rk = ((ELPH_float)iz) / lattice->fft_dims[2];

        ELPH_float qdotr =
            qpt_crys[0] * ri + qpt_crys[1] * rj + qpt_crys[2] * rk;
        eikr[ig] = cexp(factor * qdotr);
    }

    for (ND_int iset = 0; iset < nsets; ++iset)
    {
        ELPH_cmplx* pot_grid_iset = pot_grid + iset * nffts;
        ELPH_OMP_PAR_FOR_SIMD
        for (ND_int ig = 0; ig < nffts; ++ig)
        {
            pot_grid_iset[ig] *= eikr[ig];
        }
    }

    free(eikr);
}
