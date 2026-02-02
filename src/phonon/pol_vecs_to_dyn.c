// Convert polarization vectors to dynamical matrix (will be in row major)
//
//
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "common/error.h"
#include "common/numerical_func.h"
#include "elphC.h"
#include "phonon.h"

void pol_vecs_to_dyn(const ELPH_float* omega, const ND_int natom,
                     const ELPH_float* atomic_masses, ELPH_cmplx* pol_vecs)
{
    ND_int nmodes = 3 * natom;
    ELPH_cmplx* tmp_buf = malloc(2 * sizeof(*tmp_buf) * nmodes * nmodes);
    CHECK_ALLOC(tmp_buf);
    // mass normalize
    mass_normalize_pol_vecs(atomic_masses, nmodes, natom, 0.5, pol_vecs);
    //
    memcpy(tmp_buf, pol_vecs, sizeof(*tmp_buf) * nmodes * nmodes);
    memcpy(tmp_buf + nmodes * nmodes, pol_vecs,
           sizeof(*tmp_buf) * nmodes * nmodes);
    //
    for (ND_int i = 0; i < nmodes; ++i)
    {
        ELPH_float factor = omega[i] > 0 ? 1.0 : -1.0;
        factor = factor * omega[i] * omega[i];
        for (ND_int j = 0; j < nmodes; ++j)
        {
            tmp_buf[i * nmodes + j] = factor * conj(tmp_buf[i * nmodes + j]);
        }
    }

    matmul_cmplx('T', 'N', tmp_buf + nmodes * nmodes, tmp_buf, pol_vecs, 1.0,
                 0.0, nmodes, nmodes, nmodes, nmodes, nmodes, nmodes);

    free(tmp_buf);
}
