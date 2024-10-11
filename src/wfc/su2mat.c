/*
This file contains SU2mat functions which computes spinor rotation matrix for a
given symmetry matrix
*/
#include <complex.h>
#include <math.h>
#include <stdbool.h>

#include "../common/constants.h"
#include "../common/error.h"
#include "../common/numerical_func.h"
#include "../elphC.h"
#include "wfc.h"

/* SU2 mat */
void SU2mat(const ELPH_float* sym_in, const ND_int nspinor,
            const bool invert_sym, const bool time_rev, ELPH_cmplx* su2mat)
{
    if (nspinor == 1)
    {
        su2mat[0] = 1.0;
        return;
    }

    /* Convert to proper rotations */
    ELPH_float det = det3x3(sym_in);

    // printf("det : %f \n",det);
    if (fabs(fabs(det) - 1) > 1e-4)
    {
        error_msg("Wrong symmetric matrices");
    }

    ELPH_float sym[9];

    for (int i = 0; i < 9; ++i)
    {
        sym[i] = sym_in[i] / det;
    }

    ELPH_cmplx pauli0[4] = {1, 0, 0, 1};
    ELPH_cmplx paulix[4] = {0, 1, 1, 0};
    ELPH_cmplx pauliy[4] = {0, -I, I, 0};
    ELPH_cmplx pauliz[4] = {1, 0, 0, -1};

    ELPH_float alpha, beta, delta;

    if (fabs(sym[6] + 1) < ELPH_EPS)
    {
        alpha = 0;
        beta = ELPH_PI / 2.0;
        delta = atan2(sym[1], sym[2]);
    }
    else if (fabs(sym[6] - 1) < ELPH_EPS)
    {
        alpha = 0;
        beta = -ELPH_PI / 2.0;
        delta = atan2(-sym[1], -sym[2]);
    }
    else
    {
        beta = -asin(sym[6]);
        alpha = atan2(sym[3] / cos(beta), sym[0] / cos(beta));
        delta = atan2(sym[7] / cos(beta), sym[8] / cos(beta));
    }
    alpha = alpha / 2;
    beta = beta / 2;
    delta = delta / 2;

    ELPH_cmplx su_delta[4], su_beta[4], su_alpha[4];

    for (int i = 0; i < 4; ++i)
    {
        su_delta[i] = cos(delta) * pauli0[i] - I * sin(delta) * paulix[i];
        su_beta[i] = cos(beta) * pauli0[i] - I * sin(beta) * pauliy[i];
        su_alpha[i] = cos(alpha) * pauli0[i] - I * sin(alpha) * pauliz[i];
    }

    ELPH_cmplx tempsu[4];
    matmul_Cmpl2x2(su_beta, su_delta, tempsu);
    matmul_Cmpl2x2(su_alpha, tempsu, su2mat);

    // In case of time reversel we have -I*sigma_y
    if (time_rev)
    { /* Note this is not commutative */
        ELPH_cmplx a1, b1, c1, d1;
        a1 = su2mat[0];
        b1 = su2mat[1];
        c1 = su2mat[2];
        d1 = su2mat[3];

        su2mat[0] = -c1;
        su2mat[1] = -d1;
        su2mat[2] = a1;
        su2mat[3] = b1;

/* The below if def is due to backward compatiblity. In Yambo ver <= 5.12, the
pauli matrix are not correct. so due to this it picks up an additional phase
factor of -1 for time_rev matrix */
#if defined(YAMBO_LT_5_1)
        for (int is11 = 0; is11 < 4; ++is11)
        {
            su2mat[is11] = -su2mat[is11];
        }
#endif
    }

    if (invert_sym)
    {
        ELPH_cmplx b1 = conj(su2mat[1]);

        su2mat[0] = conj(su2mat[0]);
        su2mat[3] = conj(su2mat[3]);
        su2mat[1] = conj(su2mat[2]);
        su2mat[2] = b1;
    }
}
