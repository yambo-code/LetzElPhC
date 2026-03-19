#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/dtypes.h"
#include "common/error.h"
#include "elphC.h"
#include "phonon.h"

static void quadrupole_3d_to_2d_internal(const ND_int natom,
                                         const ELPH_float* epsilon,
                                         const ELPH_float* tau,
                                         const ELPH_float* Zborn,
                                         const ELPH_float* Qpoles_3D,
                                         ELPH_float* Qpoles_2D);

void quadrupole_3d_to_2d(struct Lattice* lattice, struct Phonon* phonon)
{
    if (lattice->dimension != '2')
    {
        return;
    }
    if (phonon->epsilon && (phonon->Zborn || phonon->Qpole))
    {
        ELPH_float* Qpoles_2D = calloc(27 * lattice->natom, sizeof(*Qpoles_2D));
        CHECK_ALLOC(Qpoles_2D);

        // let the compiler remove this remove.
        for (ND_int i = 0; i < (27 * lattice->natom); ++i)
        {
            Qpoles_2D[i] = 0.0;
        }
        //
        quadrupole_3d_to_2d_internal(lattice->natom, phonon->epsilon,
                                     lattice->atomic_pos, phonon->Zborn,
                                     phonon->Qpole, Qpoles_2D);

        if (phonon->Qpole)
        {
            memcpy(phonon->Qpole, Qpoles_2D,
                   27 * lattice->natom * sizeof(*Qpoles_2D));
            free(Qpoles_2D);
        }
        else
        {
            phonon->Qpole = Qpoles_2D;
            // Note we donot need to free here, as when we free phonon struct.
            // it gets automatically freed.
        }
    }
    return;
}

static void quadrupole_3d_to_2d_internal(
    const ND_int natom, const ELPH_float* epsilon, const ELPH_float* tau,
    const ELPH_float* Zborn, const ELPH_float* Qpoles_3D, ELPH_float* Qpoles_2D)
{
    // Convert 3D quadrupoles to 2D
    // See Appendix B of PHYS. REV. X 11, 041027 (2021)
    //
    if (!tau || !Qpoles_2D || (!Zborn && !Qpoles_3D))
    {
        return;
    }
    //
    ELPH_float inv_eps = 1.0 / epsilon[8];

    for (ND_int kappa = 0; kappa < natom; kappa++)
    {
        ELPH_float tau_kz = tau[kappa * 3 + 2];
        for (ND_int beta = 0; beta < 3; beta++)
        {
            /* ---------------------------------------------------------------
             * Step 1: Out-of-plane Q^{zz}_{kappa beta} (Equation B5b)
             * ---------------------------------------------------------------
             */
            ELPH_float Z_z = Zborn ? Zborn[kappa * 9 + 2 * 3 + beta] : 0.0;
            ELPH_float Q_zz =
                Qpoles_3D ? Qpoles_3D[kappa * 27 + 2 * 9 + 2 * 3 + beta] : 0.0;
            ELPH_float Qhat_zz = (Q_zz + 2.0 * tau_kz * Z_z) * inv_eps;
            //
            Qpoles_2D[kappa * 27 + 2 * 9 + 2 * 3 + beta] = Qhat_zz;
            /* ---------------------------------------------------------------
             * Step 2: Mixed components Q^{z alpha} and Q^{alpha z} (Equation
             * B5a)
             * ---------------------------------------------------------------
             */
            for (ND_int alpha = 0; alpha < 2; alpha++)
            {
                ELPH_float Z_a =
                    Zborn ? Zborn[kappa * 9 + alpha * 3 + beta] : 0.0;
                ELPH_float Q_za =
                    Qpoles_3D ? Qpoles_3D[kappa * 27 + 2 * 9 + alpha * 3 + beta]
                              : 0.0;
                ELPH_float Q_az =
                    Qpoles_3D ? Qpoles_3D[kappa * 27 + alpha * 9 + 2 * 3 + beta]
                              : 0.0;
                //
                Qpoles_2D[kappa * 27 + 2 * 9 + alpha * 3 + beta] =
                    (Q_za + tau_kz * Z_a) * inv_eps;

                Qpoles_2D[kappa * 27 + alpha * 9 + 2 * 3 + beta] =
                    (Q_az + tau_kz * Z_a) * inv_eps;
            }

            /* ---------------------------------------------------------------
             * Step 3: In-plane Q^{alpha gamma} (Equation B10)
             * ---------------------------------------------------------------
             */
            for (ND_int alpha = 0; alpha < 2; alpha++)
            {
                for (ND_int gamma = 0; gamma < 2; gamma++)
                {
                    ELPH_float Q_ag = Qpoles_3D
                                          ? Qpoles_3D[kappa * 27 + alpha * 9 +
                                                      gamma * 3 + beta]
                                          : 0.0;
                    ELPH_float delta = (alpha == gamma) ? 1.0 : 0.0;
                    //
                    Qpoles_2D[kappa * 27 + alpha * 9 + gamma * 3 + beta] =
                        Q_ag + (delta - epsilon[alpha * 3 + gamma]) * Qhat_zz;
                }
            }
        }
    }
}
