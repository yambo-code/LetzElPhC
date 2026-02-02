// does a mass normalization of phonon eigenvectors
//
//
#include <complex.h>
#include <math.h>

#include "elphC.h"
#include "phonon.h"

void mass_normalize_pol_vecs(const ELPH_float* atomic_masses,
                             const ND_int nsets, const ND_int natoms,
                             const ELPH_float power, ELPH_cmplx* pol_vecs)
{
    // pol_vecs (nsets, natom,3);
    // atomic_masses (natom)
    // does  e^v_ix = Mi^power * e^v_ix
    //
    for (ND_int imode = 0; imode < nsets; ++imode)
    {
        for (ND_int ia = 0; ia < natoms; ++ia)
        {
            ELPH_float mass_tmp = pow(atomic_masses[ia], power);
            for (ND_int i = 0; i < 3; ++i)
            {
                pol_vecs[imode * 3 * natoms + ia * 3 + i] *= mass_tmp;
            }
        }
    }
}
