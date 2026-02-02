// compute non analytical term for phonon
//
// 1) Only the dipole-dipole term is include//
//

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/ELPH_timers.h"
#include "common/constants.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/omp_pragma_def.h"
#include "elphC.h"
#include "phonon.h"

static void add_ph_dyn_long_range_internal(
    const ELPH_float* qpt, struct Lattice* lattice, struct Phonon* phonon,
    const ND_int* Ggrid, const ND_int sign, const ELPH_float* atomic_masses,
    const ELPH_float eta, ELPH_cmplx* dyn_mat);

void add_ph_dyn_long_range(const ELPH_float* qpt, struct Lattice* lattice,
                           struct Phonon* phonon, const ND_int* Ggrid,
                           const ND_int sign, const ELPH_float* atomic_masses,
                           ELPH_cmplx* dyn_mat_asr, const ELPH_float eta,
                           ELPH_cmplx* dyn_mat)
{
    // compute the long range part
    add_ph_dyn_long_range_internal(qpt, lattice, phonon, Ggrid, sign,
                                   atomic_masses, eta, dyn_mat);
    // now subtract the asr correction
    //
    ELPH_float factor = -1;
    if (sign < 0)
    {
        factor = -factor;
    }

    ND_int nmodes = lattice->natom * 3;
    //
    for (ND_int ia = 0; ia < lattice->natom; ++ia)
    {
        ELPH_cmplx* tmp_buf_ia = dyn_mat_asr + ia * 9;
        for (ND_int i = 0; i < 3; ++i)
        {
            for (ND_int j = 0; j < 3; ++j)
            {
                dyn_mat[(ia * 3 + i) * nmodes + ia * 3 + j] +=
                    factor * tmp_buf_ia[3 * i + j];
            }
        }
    }
}

void compute_dyn_lr_asr_correction(struct Lattice* lattice,
                                   struct Phonon* phonon, const ND_int* Ggrid,
                                   const ELPH_float* atomic_masses,
                                   const ELPH_float eta,
                                   ELPH_cmplx* dyn_mat_asr)
{
    // ASR for the lr non-analytical dynamical matrix.
    // (dyn_mat_asr) : (natom, 3, 3)
    // See X. Gonze  PhysRevB.50.13035
    //
    for (ND_int i = 0; i < lattice->natom * 9; ++i)
    {
        dyn_mat_asr[i] = 0.0;
    }

    if (!phonon->epsilon || !phonon->Zborn)
    {
        return;
    }
    ND_int nmodes = 3 * lattice->natom;
    ELPH_cmplx* tmp_dyn_mat = calloc(nmodes * nmodes, sizeof(*tmp_dyn_mat));
    CHECK_ALLOC(tmp_dyn_mat);
    // let the comiplers remove this loop.
    for (ND_int i = 0; i < nmodes * nmodes; ++i)
    {
        tmp_dyn_mat[i] = 0.0;
    }

    ELPH_float qpt_zero[3] = {0.0, 0.0, 0.0};
    add_ph_dyn_long_range_internal(qpt_zero, lattice, phonon, Ggrid, 1,
                                   atomic_masses, eta, tmp_dyn_mat);

    //
    for (ND_int ia = 0; ia < lattice->natom; ++ia)
    {
        ELPH_cmplx* tmp_buf_ia = dyn_mat_asr + ia * 9;
        for (ND_int ja = 0; ja < lattice->natom; ++ja)
        {
            for (ND_int i = 0; i < 3; ++i)
            {
                for (ND_int j = 0; j < 3; ++j)
                {
                    tmp_buf_ia[3 * i + j] +=
                        tmp_dyn_mat[(ia * 3 + i) * nmodes + ja * 3 + j];
                }
            }
        }
    }
    free(tmp_dyn_mat);
}

static void add_ph_dyn_long_range_internal(
    const ELPH_float* qpt, struct Lattice* lattice, struct Phonon* phonon,
    const ND_int* Ggrid, const ND_int sign, const ELPH_float* atomic_masses,
    const ELPH_float eta, ELPH_cmplx* dyn_mat)
{
    // adds or subtracts non-analytical term to the dynamical matrix.
    // sign < 0 : subtract else add
    // qpt in crystal coordinater
    //
    // (3D) X. Gonze et al  Phys. Rev. B 50, 13035(R)
    // (2D ) T. Sohier et al Nano Lett. 2017, 17, 6, 3758–3763
    // (2D) S. Ponce et al PHYSICAL REVIEW B 107, 155424 (2023)
    //
    if (!phonon->epsilon || (!phonon->Zborn && !phonon->Qpole))
    {
        return;
    }
    //
    ELPH_start_clock("dyn lr part");

    if (lattice->dimension == '2' && fabs(qpt[2]) > ELPH_EPS)
    {
        error_msg(
            "In 2D, only qz == 0 points are accepted when interpolating.");
    }

    ELPH_cmplx factor = 4.0 * ELPH_PI * ELPH_e2 / lattice->volume;  // prefactor
    if (sign < 0)
    {
        factor = -factor;
    }
    //
    ELPH_float eps_alpha[9];

    memcpy(eps_alpha, phonon->epsilon, sizeof(eps_alpha));

    ELPH_float zlat = lattice->alat_vec[8];
    //
    if (lattice->dimension == '2')
    {
        factor = factor * zlat / 2.0;
        // compute alpha = c/2 * (eps-1)
        for (int i = 0; i < 9; ++i)
        {
            eps_alpha[i] = 0.0;
        }
        // inplane
        eps_alpha[0] = 0.5 * zlat * (phonon->epsilon[0] - 1);
        eps_alpha[1] = 0.5 * zlat * (phonon->epsilon[1]);
        eps_alpha[3] = 0.5 * zlat * (phonon->epsilon[3]);
        eps_alpha[4] = 0.5 * zlat * (phonon->epsilon[4] - 1);
    }

    ND_int GridZ = Ggrid[2];
    if (lattice->dimension == '2')
    {
        GridZ = 1;
    }
    //
    ND_int Gvec_size = Ggrid[0] * Ggrid[1] * GridZ;

    ND_int natom = lattice->natom;

    ELPH_cmplx* Zdotq_tau =
        calloc(2 * Gvec_size * 3 * natom, sizeof(*Zdotq_tau));
    CHECK_ALLOC(Zdotq_tau);
    // leave it to the compilers(most will very likely remove it )
    for (ND_int i = 0; i < (2 * Gvec_size * 3 * natom); ++i)
    {
        Zdotq_tau[i] = 0.0;
    }
    //
    ELPH_OMP_PAR_FOR_SIMD
    for (ND_int ig = 0; ig < Gvec_size; ++ig)
    {
        ELPH_cmplx* out_tmp1 = Zdotq_tau + ig * 3 * natom;
        ELPH_cmplx* out_tmp2 =
            Zdotq_tau + ig * 3 * natom + Gvec_size * 3 * natom;
        const ND_int Gx = ig / (Ggrid[1] * GridZ);
        const ND_int Gy = (ig % (Ggrid[1] * GridZ)) / GridZ;
        const ND_int Gz = (ig % (Ggrid[1] * GridZ)) % GridZ;
        //
        ELPH_float qplusG[3], tmp_buf[3];
        tmp_buf[0] = (qpt[0] + get_miller_idx(Gx, Ggrid[0]));
        tmp_buf[1] = (qpt[1] + get_miller_idx(Gy, Ggrid[1]));
        tmp_buf[2] = (qpt[2] + get_miller_idx(Gz, GridZ));
        // COnvert to cart units (2*pi) is included
        MatVec3f(lattice->blat_vec, tmp_buf, false, qplusG);
        //
        ELPH_float qplusG_norm = sqrt(dot3_macro(qplusG, qplusG));
        if (qplusG_norm < ELPH_EPS)
        {
            continue;
        }
        //
        MatVec3f(eps_alpha, qplusG, false, tmp_buf);
        ELPH_float q_eps_q = dot3_macro(tmp_buf, qplusG);

        ELPH_float decay_fac = exp(-q_eps_q * q_eps_q * 0.125 / eta);

        if (lattice->dimension == '2')
        {
            q_eps_q += qplusG_norm;
            decay_fac = exp(-qplusG_norm * qplusG_norm * 0.125 / eta);
        }
        //
        for (ND_int ia = 0; ia < natom; ++ia)
        {
            const ELPH_float* Z_k =
                phonon->Zborn ? (phonon->Zborn + 9 * ia) : NULL;
            const ELPH_float* Q_k =
                phonon->Qpole ? (phonon->Qpole + 27 * ia) : NULL;
            const ELPH_float* tau_k = lattice->atomic_pos + 3 * ia;

            ELPH_float Qpole_buf[3] = {0.0, 0.0, 0.0};
            //
            if (Z_k)
            {
                MatVec3f(Z_k, qplusG, true, tmp_buf);
            }
            else
            {
                memcpy(tmp_buf, Qpole_buf, sizeof(tmp_buf));
            }
            //
            ELPH_cmplx* out_tmp_buf = out_tmp1 + 3 * ia;
            ELPH_cmplx* out_tmp_buf2 = out_tmp2 + 3 * ia;

            // compute (q+G)_x Q_xyz * (q+G)_y
            if (Q_k)
            {
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        for (int k = 0; k < 3; ++k)
                        {
                            Qpole_buf[k] =
                                Qpole_buf[k] +
                                qplusG[i] * qplusG[j] * Q_k[k + 3 * j + 9 * i];
                        }
                    }
                }
                if (lattice->dimension == '2')
                {
                    // subtract Qzz|q+G|^2
                    for (int i = 0; i < 3; ++i)
                    {
                        Qpole_buf[i] -=
                            (qplusG_norm * qplusG_norm * Q_k[i + 24]);
                    }
                }
            }
            // e^{-iq.tau}
            ELPH_cmplx qdot_tau = cexp(-I * (dot3_macro(qplusG, tau_k)));
            qdot_tau /= q_eps_q;
            qdot_tau *= decay_fac;
            qdot_tau /= sqrt(atomic_masses[ia]);
            //
            for (ND_int i = 0; i < 3; ++i)
            {
                out_tmp_buf[i] =
                    (tmp_buf[i] - I * 0.5 * Qpole_buf[i]) * qdot_tau;
                out_tmp_buf2[i] = out_tmp_buf[i] * q_eps_q;
            }
        }
    }
    //
    ND_int nmodes = 3 * natom;
    matmul_cmplx('C', 'N', Zdotq_tau + Gvec_size * 3 * natom, Zdotq_tau,
                 dyn_mat, factor, 1.0, nmodes, nmodes, nmodes, nmodes, nmodes,
                 Gvec_size);
    //
    free(Zdotq_tau);
    //
    //
    ELPH_stop_clock("dyn lr part");
    //
    return;
}
