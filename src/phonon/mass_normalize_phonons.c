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

void mass_normalize_force_constants(const ELPH_float* atomic_masses,
                                    const ND_int nsets, const ND_int natoms,
                                    const ELPH_float power, ELPH_cmplx* frc)
{
    // Force constants (nsets, natom, 3, natom,3);
    // atomic_masses (natom)
    // FC(iset, ia, i, jb, j) *= Ma^power * Mb^power
    ND_int nmodes = 3 * natoms;
    //
    for (ND_int ia = 0; ia < natoms; ++ia)
    {
        const ELPH_float Ma = pow(atomic_masses[ia], power);
        for (ND_int ib = 0; ib < natoms; ++ib)
        {
            const ELPH_float Mb = pow(atomic_masses[ib], power);
            for (ND_int iset = 0; iset < nsets; ++iset)
            {
                ELPH_cmplx* frc_tmp =
                    frc + iset * nmodes * nmodes + ia * 9 * natoms + 3 * ib;
                for (ND_int i = 0; i < 3; ++i)
                {
                    for (ND_int j = 0; j < 3; ++j)
                    {
                        frc_tmp[j + i * nmodes] *= (Ma * Mb);
                    }
                }
            }
        }
    }
}
