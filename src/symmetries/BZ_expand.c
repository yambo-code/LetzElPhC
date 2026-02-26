#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "elphC.h"
#include "symmetries.h"

ND_int bz_expand(const ND_int Nibz, const ND_int Nsym,
                 const ELPH_float* ibz_kpts, const struct symmetry* symms,
                 const ELPH_float* lat_vec, ELPH_float* kpoints, ND_int* kstar,
                 int* kmap)
{
    /*
    Expands k/q points from iBZ to full BZ
    ibz_kpts must be in cart units and
    outpur "kpoints" is in crystal units

    Note kpoints array must be allocated before
    */
    ND_int nkBZ_found = 0;

    for (ND_int ikpt = 0; ikpt < Nibz; ++ikpt)
    {
        const ELPH_float* kpt_tmp = ibz_kpts + 3 * ikpt;
        int nkstar = 0;

        for (ND_int isym = 0; isym < Nsym; ++isym)
        {
            ELPH_float kstar_cart[3], kstar_crys[3];

            MatVec3f(symms[isym].Rmat, kpt_tmp, false, kstar_cart);

            MatVec3f(lat_vec, kstar_cart, true, kstar_crys);

            if (find_kidx_in_list(nkBZ_found, kpoints, kstar_crys) < 0)
            {
                memcpy(kpoints + nkBZ_found * 3, kstar_crys,
                       3 * sizeof(ELPH_float));
                kmap[2 * nkBZ_found] = ikpt;
                kmap[2 * nkBZ_found + 1] = isym;
                ++nkstar;
                ++nkBZ_found;
            }
        }
        if (kstar)
        {
            kstar[ikpt] = nkstar;
        }
    }

    // In case the user has given full BZ points for iBZ, then we
    // simple ignore symmetry expansion and set those variables)
    if (Nsym > 1 && nkBZ_found == Nibz)
    {
        // find index of identity
        int iden_idx = -1;
        for (ND_int isym = 0; isym < Nsym; ++isym)
        {
            ELPH_float sum_diff = 3;
            for (ND_int ix = 0; ix < 9; ++ix)
            {
                sum_diff += symms[isym].Rmat[ix] * symms[isym].Rmat[ix];
            }
            for (ND_int ix = 0; ix < 3; ++ix)
            {
                sum_diff -= 2.0 * symms[isym].Rmat[ix * 4];
            }
            sum_diff = sqrt(fabs(sum_diff)) / 3;
            if (sum_diff < 1e-5)
            {
                iden_idx = isym;
                break;
            }
        }
        if (iden_idx < 0)
        {
            error_msg("Identity matrix not found.");
        }
        //
        for (ND_int ikpt = 0; ikpt < Nibz; ++ikpt)
        {
            kmap[2 * ikpt] = ikpt;
            kmap[2 * ikpt + 1] = iden_idx;
            if (kstar)
            {
                kstar[ikpt] = 1;
            }
            MatVec3f(lat_vec, ibz_kpts + 3 * ikpt, true, kpoints + 3 * ikpt);
        }
    }

    return nkBZ_found;
}
