#include <stdbool.h>
#include <string.h>

#include "../common/dtypes.h"
#include "../common/numerical_func.h"
#include "../elphC.h"
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

    return nkBZ_found;
}
