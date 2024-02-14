/*
This file contains function that finds the index of S^-1 in the group
*/

#include "symmetries.h"

ND_int find_inv_symm_idx(ND_int nsym, const ELPH_float* Smat,
                         const ELPH_float* point_group, const bool trans)
{
    // set trans = true if symmetries in point group array must be transpose i.e
    // if trans = true => Smat@point_group.T == 1 (compared to Identity)
    // else Smat@point_group is compared to Identity
    // return -1 if idx not found

    ND_int idx = -1;

    char Sym_trans = 'N';
    if (trans)
    {
        Sym_trans = 'T';
    }

    for (ND_int i = 0; i < nsym; ++i)
    {
        ELPH_float Mat_tmp[9];

        Gemm3x3f(Smat, 'N', point_group + i * 9, Sym_trans, Mat_tmp);

        Mat_tmp[0] -= 1;
        Mat_tmp[4] -= 1;
        Mat_tmp[8] -= 1;

        ELPH_float sum_tmp = 0;
        for (int xi = 0; xi < 9; ++xi)
        {
            sum_tmp += Mat_tmp[xi] * Mat_tmp[xi];
        }
        sum_tmp = sqrt(sum_tmp);
        if (sum_tmp < ELPH_EPS)
        {
            idx = i;
            break;
        }
    }

    return idx;
}
