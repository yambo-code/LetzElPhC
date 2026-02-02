// give a uniform gamma centred, get generate kpoints in iBZ
//
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "common/constants.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "elphC.h"
#include "symmetries.h"

ND_int generate_iBZ_kpts(const ND_int* kgrid, const ND_int Nsym,
                         const struct symmetry* symms,
                         const ELPH_float* lat_vec, const ELPH_float* blat,
                         ELPH_float* ibz_kpts, const bool crystal)
{
    // the size of ibz_kpt must be atleat (product(kgrid),3)
    // if crystal is true, then
    // //
    for (ND_int i = 0; i < kgrid[0]; ++i)
    {
        for (ND_int j = 0; j < kgrid[1]; ++j)
        {
            for (ND_int k = 0; k < kgrid[2]; ++k)
            {
                ELPH_float* ik_ibz_tmp =
                    ibz_kpts + 3 * (k + j * kgrid[2] + i * kgrid[1] * kgrid[2]);
                ik_ibz_tmp[0] = ((ELPH_float)i) / kgrid[0];
                ik_ibz_tmp[1] = ((ELPH_float)j) / kgrid[1];
                ik_ibz_tmp[2] = ((ELPH_float)k) / kgrid[2];

                ik_ibz_tmp[0] -= rint(ik_ibz_tmp[0]);
                ik_ibz_tmp[1] -= rint(ik_ibz_tmp[1]);
                ik_ibz_tmp[2] -= rint(ik_ibz_tmp[2]);
            }
        }
    }
    //
    ND_int kgrid_product = kgrid[0] * kgrid[1] * kgrid[2];
    //
    if (!symms || Nsym == 1)
    {
        if (!crystal)
        {
            for (ND_int i = 0; i < kgrid_product; ++i)
            {
                ELPH_float qtmp[3];
                MatVec3f(blat, ibz_kpts + 3 * i, false, qtmp);
                for (ND_int xi = 0; xi < 3; ++xi)
                {
                    ibz_kpts[3 * i + xi] = qtmp[xi] / (2.0 * ELPH_PI);
                }
            }
        }
        return kgrid_product;
    }

    bool* iBZ_list = calloc(kgrid_product, sizeof(*iBZ_list));
    CHECK_ALLOC(iBZ_list);

    // Most compilers will remove this loop
    for (ND_int ikpt = 0; ikpt < kgrid_product; ++ikpt)
    {
        iBZ_list[ikpt] = false;
    }

    for (ND_int ikpt = 0; ikpt < kgrid_product; ++ikpt)
    {
        if (iBZ_list[ikpt])
        {
            continue;
        }
        const ELPH_float* kpt_tmp = ibz_kpts + 3 * ikpt;
        ELPH_float kcart_tmp[3];
        MatVec3f(blat, kpt_tmp, false, kcart_tmp);
        for (ND_int ix = 0; ix < 3; ++ix)
        {
            kcart_tmp[ix] /= (2 * ELPH_PI);
        }

        for (ND_int isym = 0; isym < Nsym; ++isym)
        {
            ELPH_float kstar_cart[3], kstar_crys[3];

            MatVec3f(symms[isym].Rmat, kcart_tmp, false, kstar_cart);

            MatVec3f(lat_vec, kstar_cart, true, kstar_crys);

            ND_int ikpt_rot_idx =
                find_kidx_in_list(kgrid_product - ikpt, kpt_tmp, kstar_crys);
            if (ikpt_rot_idx <= 0)
            {
                continue;
            }
            else
            {
                iBZ_list[ikpt + ikpt_rot_idx] = true;
            }
        }
    }
    ND_int nkiBZ_found = 0;

    for (ND_int ikpt = 0; ikpt < kgrid_product; ++ikpt)
    {
        if (!iBZ_list[ikpt])
        {
            if (nkiBZ_found != ikpt)
            {
                memcpy(ibz_kpts + 3 * nkiBZ_found, ibz_kpts + 3 * ikpt,
                       3 * sizeof(*ibz_kpts));
            }
            ++nkiBZ_found;
        }
    }
    free(iBZ_list);

    if (!crystal)
    {
        for (ND_int i = 0; i < nkiBZ_found; ++i)
        {
            ELPH_float qtmp[3];
            MatVec3f(blat, ibz_kpts + 3 * i, false, qtmp);
            for (ND_int xi = 0; xi < 3; ++xi)
            {
                ibz_kpts[3 * i + xi] = qtmp[xi] / (2.0 * ELPH_PI);
            }
        }
    }

    return nkiBZ_found;
}
