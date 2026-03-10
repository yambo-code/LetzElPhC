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
    // If atomic_masses as not NULL, this implies, we are passing mass
    // normalized dynamical matrices, else we are passing dynamical matrices
    // without 1/sqrt(Ma Mb) factor

    // compute the long range part
    add_ph_dyn_long_range_internal(qpt, lattice, phonon, Ggrid, sign,
                                   atomic_masses, eta, dyn_mat);
    //
    if (!dyn_mat_asr)
    {
        return;
    }
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
        ELPH_float ia_mass = 1.0;
        if (atomic_masses)
        {
            ia_mass = 1.0 / atomic_masses[ia];
        }
        //
        for (ND_int i = 0; i < 3; ++i)
        {
            for (ND_int j = 0; j < 3; ++j)
            {
                dyn_mat[(ia * 3 + i) * nmodes + ia * 3 + j] +=
                    (factor * tmp_buf_ia[3 * i + j] * ia_mass);
            }
        }
    }
}

void compute_dyn_lr_asr_correction(struct Lattice* lattice,
                                   struct Phonon* phonon, const ND_int* Ggrid,
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

    if (!phonon->epsilon || (!phonon->Zborn && !phonon->Qpole))
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
    // We donot want mass normalized.
    add_ph_dyn_long_range_internal(qpt_zero, lattice, phonon, Ggrid, 1, NULL,
                                   eta, tmp_dyn_mat);

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
    // (2D ) M Royo et al PHYSICAL REVIEW X 11, 041027 (2021)
    //
    // if atomic_masses is NULL, then the dynamical matrices are not mass
    // normalized.
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

    // In the case of 2D, we need to seperate out the perperndicular
    // contribution
    ELPH_cmplx* Zdotq_tau_per = NULL;
    if (lattice->dimension == '2')
    {
        Zdotq_tau_per =
            calloc(2 * Gvec_size * 3 * natom, sizeof(*Zdotq_tau_per));
        CHECK_ALLOC(Zdotq_tau_per);
    }
    // leave it to the compilers(most will very likely remove it )
    for (ND_int i = 0; i < (2 * Gvec_size * 3 * natom); ++i)
    {
        Zdotq_tau[i] = 0.0;
    }
    if (Zdotq_tau_per)
    {
        for (ND_int i = 0; i < (2 * Gvec_size * 3 * natom); ++i)
        {
            Zdotq_tau_per[i] = 0.0;
        }
    }
    //
    ELPH_OMP_PAR_FOR_SIMD
    for (ND_int ig = 0; ig < Gvec_size; ++ig)
    {
        ELPH_cmplx* out_tmp1 = Zdotq_tau + ig * 3 * natom;
        ELPH_cmplx* out_tmp2 =
            Zdotq_tau + ig * 3 * natom + Gvec_size * 3 * natom;
        //
        ELPH_cmplx* out_tmp1_per = NULL;
        ELPH_cmplx* out_tmp2_per = NULL;
        if (Zdotq_tau_per)
        {
            out_tmp1_per = Zdotq_tau_per + ig * 3 * natom;
            out_tmp2_per =
                Zdotq_tau_per + ig * 3 * natom + Gvec_size * 3 * natom;
        }
        //
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
        ELPH_float q_eps_q_per = 1.0;

        ELPH_float decay_fac = exp(-fabs(q_eps_q) * 0.125 / eta);

        if (lattice->dimension == '2')
        {
            // In 2D f(q) = 1-tanh(K*eta/2) is the decay factor
            // Evaluates f(|q|) using L = eta_induced*4*pi*alpha_per*1.001
            // Note eta >= 1 to satisfy stablity condition
            ELPH_float Leff =
                eta * zlat * 1.001 * (1.0 - 1.0 / phonon->epsilon[8]);
            decay_fac = 1 - tanh(qplusG_norm * 0.5 * Leff);
            // L = eta * zlat/2
            if (fabs(decay_fac) < ELPH_EPS)
            {
                continue;
            }

            q_eps_q *= decay_fac;
            q_eps_q += qplusG_norm;
            /* q_eps_q_per = 1 - 0.5 * zlat * qplusG_norm * decay_fac * */
            /*                       (1.0 - 1.0 / phonon->epsilon[8]); */
            q_eps_q_per = 1.0;
            // Note: The exact Royo 1/eps_perp causes unphysical divergence far
            // from Gamma. Because x = 2*PI*|q|*f(q)*alpha_perp < 1, we Taylor
            // expand 1/(1-x) ~ 1 + x. The out-of-plane dipole already has a |q|
            // prefactor. Keeping the x ~ |q| term creates unphysical O(|q|^2)
            // growth that breaks Wannier interpolation. Truncating to 0th-order
            // (1.0) keeps the exact linear physics at Gamma and stabilizes the
            // boundary.
            if (fabs(q_eps_q_per) < ELPH_EPS)
            {
                q_eps_q_per = 0.0;
            }
            else
            {
                q_eps_q_per = -qplusG_norm / q_eps_q_per;
            }
            decay_fac = sqrt(decay_fac);
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
            if (atomic_masses)
            {
                qdot_tau /= sqrt(atomic_masses[ia]);
            }
            //
            for (ND_int i = 0; i < 3; ++i)
            {
                out_tmp_buf[i] =
                    (tmp_buf[i] - I * 0.5 * Qpole_buf[i]) * qdot_tau;
                out_tmp_buf2[i] = out_tmp_buf[i] * q_eps_q;
            }
            // Out of plane term for 2D
            if (Zdotq_tau_per)
            {
                qdot_tau *= q_eps_q;
                for (ND_int i = 0; i < 3; ++i)
                {
                    if (Z_k)
                    {
                        out_tmp1_per[3 * ia + i] += Z_k[i + 6] * qdot_tau *
                                                    q_eps_q_per /
                                                    phonon->epsilon[8];
                        out_tmp2_per[3 * ia + i] +=
                            Z_k[i + 6] * qdot_tau / phonon->epsilon[8];
                    }
                    if (Q_k)
                    {
                        for (ND_int ixx = 0; ixx < 3; ++ixx)
                        {
                            out_tmp1_per[3 * ia + i] -= I * qplusG[ixx] *
                                                        Q_k[i + 18 + 3 * ixx] *
                                                        qdot_tau * q_eps_q_per;
                            out_tmp2_per[3 * ia + i] -= qplusG[ixx] * I *
                                                        Q_k[i + 18 + 3 * ixx] *
                                                        qdot_tau;
                        }
                    }
                }
            }
        }
    }
    //
    ND_int nmodes = 3 * natom;
    matmul_cmplx('C', 'N', Zdotq_tau + Gvec_size * 3 * natom, Zdotq_tau,
                 dyn_mat, factor, 1.0, nmodes, nmodes, nmodes, nmodes, nmodes,
                 Gvec_size);
    //
    if (Zdotq_tau_per)
    {
        matmul_cmplx('C', 'N', Zdotq_tau_per + Gvec_size * 3 * natom,
                     Zdotq_tau_per, dyn_mat, factor, 1.0, nmodes, nmodes,
                     nmodes, nmodes, nmodes, Gvec_size);
    }
    free(Zdotq_tau);
    free(Zdotq_tau_per);
    //
    //
    ELPH_stop_clock("dyn lr part");
    //
    return;
}
