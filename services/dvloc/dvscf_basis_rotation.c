#include <stdlib.h>
#include <string.h>

#include "common/error.h"
#include "common/numerical_func.h"
#include "dvloc.h"
#include "elphC.h"
// inplace rotation of dvscf
//
//
//
void dVscf_change_basis(ELPH_cmplx* dvscf, const ELPH_cmplx* rot_vecs,
                        const ND_int nsets, const ND_int nmodes,
                        const ND_int nmag, const ND_int Nx, const ND_int Ny,
                        const ND_int Nz, const char blas_char)
{
    // dVscf, (nsets, nmodes, nmag, Nx, Ny, Nz_loc)
    // rot_vecs (modes, modes)
    // if blas_char == 'N" rot_vecs@dVscf
    // == 'C', rot_vecs^dagger @dvscf
    // === 'T' rot_vecs^T@dvscf

    ND_int tmp_buf_len = Ny * Nz * nmodes;
    ELPH_cmplx* tmp_buf = calloc(tmp_buf_len, sizeof(*tmp_buf));
    CHECK_ALLOC(tmp_buf);

    ND_int dvscf_stride = nmag * Nx * Ny * Nz;

    for (ND_int iset = 0; iset < nsets; ++iset)
    {
        for (ND_int img = 0; img < nmag; ++img)
        {
            for (ND_int ix = 0; ix < Nx; ++ix)
            {
                ELPH_cmplx* dvscf_tmp = dvscf + iset * nmodes * dvscf_stride +
                                        img * Nx * Ny * Nz + ix * Ny * Nz;
                matmul_cmplx(blas_char, 'N', rot_vecs, dvscf_tmp, tmp_buf, 1.0,
                             0.0, nmodes, dvscf_stride, Ny * Nz, nmodes,
                             Ny * Nz, nmodes);
                for (ND_int imode = 0; imode < nmodes; ++imode)
                {
                    memcpy(dvscf_tmp + imode * dvscf_stride,
                           tmp_buf + imode * Ny * Nz,
                           Ny * Nz * sizeof(*tmp_buf));
                }
            }
        }
    }

    free(tmp_buf);
}
