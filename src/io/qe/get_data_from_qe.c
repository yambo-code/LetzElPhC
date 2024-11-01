#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../common/constants.h"
#include "../../common/dtypes.h"
#include "../../common/error.h"
#include "../../common/numerical_func.h"
#include "../../common/parallel.h"
#include "../../elphC.h"
#include "../../symmetries/symmetries.h"
#include "../mpi_bcast.h"
#include "qe_io.h"

void get_data_from_qe(struct Lattice* lattice, struct Phonon* phonon,
                      const char* ph_save_dir, char*** pseudo_pots,
                      const struct ELPH_MPI_Comms* Comm)
{
    /*
    This functions gets basic ground state data and phonon data from q.e
    that are not saved by YAMBO.
    */
    int mpi_error;

    char* tmp_buffer = NULL;
    if (Comm->commW_rank == 0)
    {
        tmp_buffer = calloc(1024 + strlen(ph_save_dir), 1);
        CHECK_ALLOC(tmp_buffer);
    }
    if (Comm->commW_rank == 0)
    {
        strcpy(tmp_buffer, ph_save_dir);
        strcat(tmp_buffer, "/dyn0");
        read_qpts_qe(tmp_buffer, &phonon->nq_iBZ, &phonon->nq_BZ,
                     &phonon->qpts_iBZ);
    }
    // Bcast data
    mpi_error = MPI_Bcast(&phonon->nq_iBZ, 1, ELPH_MPI_ND_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&phonon->nq_BZ, 1, ELPH_MPI_ND_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    // divide the qpoints in iBZ over q pools
    phonon->nq_iBZ_loc = distribute_to_grps(phonon->nq_iBZ, Comm->nqpools,
                                            Comm->commW_rank / Comm->commQ_size,
                                            &phonon->nq_shift);

    if (phonon->nq_iBZ_loc < 1)
    {
        error_msg(
            "There are no qpoints in some qpools, Make sure nqpool < # of "
            "qpoints in iBZ.");
    }

    // Bcast phonon->qpts_iBZ
    if (Comm->commW_rank != 0)
    {
        phonon->qpts_iBZ = malloc(sizeof(ELPH_float) * 3 * phonon->nq_iBZ);
        CHECK_ALLOC(phonon->qpts_iBZ);
    }
    mpi_error = MPI_Bcast(phonon->qpts_iBZ, 3 * phonon->nq_iBZ, ELPH_MPI_float,
                          0, Comm->commW);
    MPI_error_msg(mpi_error);
    // qpts must be divided to get then cart units
    // now get the basic info from
    ELPH_float alat[3];
    char* PSEUDO_DIR = NULL;
    bool ph_tim_rev;

    ELPH_float* ph_sym_mats = NULL;
    ELPH_float* ph_sym_tau = NULL;

    ELPH_float lat_vec[9];  // a[:,i] is ith lattice vector
    if (Comm->commW_rank == 0)
    {
        strcpy(tmp_buffer, ph_save_dir);
        strcat(tmp_buffer, "/data-file-schema.xml");
        parse_qexml(tmp_buffer, lat_vec, alat, &lattice->dimension,
                    &lattice->is_soc_present, &lattice->nmag, lattice->fft_dims,
                    &phonon->nph_sym, &ph_sym_mats, &ph_sym_tau, &ph_tim_rev,
                    &PSEUDO_DIR, pseudo_pots);
        // free pseudo pots, no longer need
        free(PSEUDO_DIR);
        PSEUDO_DIR = NULL;
    }

    // Bcast all the variables

    mpi_error = MPI_Bcast(lat_vec, 9, ELPH_MPI_float, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error =
        MPI_Bcast(lattice->fft_dims, 3, ELPH_MPI_ND_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(alat, 3, ELPH_MPI_float, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&lattice->dimension, 1, MPI_CHAR, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&lattice->nmag, 1, ELPH_MPI_ND_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&phonon->nph_sym, 1, ELPH_MPI_ND_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error =
        MPI_Bcast(&lattice->is_soc_present, 1, MPI_C_BOOL, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&ph_tim_rev, 1, MPI_C_BOOL, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    ELPH_float blat[9];
    reciprocal_vecs(lat_vec, blat);
    for (int ix = 0; ix < 9; ++ix)
    {
        blat[ix] /= (2.0f * ELPH_PI);
    }

    // allocate memory for phonon symmetric matrices on rest of the cpus
    phonon->ph_syms = malloc(sizeof(struct symmetry) * 2 * phonon->nph_sym);
    CHECK_ALLOC(phonon->ph_syms);
    // we create 2*phonon->nph_sym sets of symmetries (factor 2 to store the
    // time rev case) Note: the second half([phonon->nph_sym:]) are only
    // used when time reversal is present
    // --
    // --
    // now convert phonon symmetric matrices and phonon frac. trans. vec to cart
    // units.
    if (Comm->commW_rank == 0)
    {
        for (ND_int isym = 0; isym < phonon->nph_sym; ++isym)
        {
            // Note we also fill the second half but are only used when tim_rev
            // is present
            ELPH_float* sym_tmp = ph_sym_mats + isym * 9;
            ELPH_float* sym_tmp_trev =
                ph_sym_mats + (isym + phonon->nph_sym) * 9;

            // It should be noted that we use Sx + v convention, but qe uses
            // S(x+v) so our v = S*tau_qe
            ELPH_float* vec_tmp = ph_sym_tau + isym * 3;
            ELPH_float* vec_tmp_trev =
                ph_sym_tau + (isym + phonon->nph_sym) * 3;

            // compute lat_vec@S@b^T
            // first  S and b^T
            Gemm3x3f(sym_tmp, 'N', blat, 'T',
                     sym_tmp_trev);  // sym_tmp_trev is used as tmp buffer
            Gemm3x3f(lat_vec, 'N', sym_tmp_trev, 'N', sym_tmp);

            for (int ix = 0; ix < 9; ++ix)
            {
                if (fabs(sym_tmp[ix]) < ELPH_EPS)
                {
                    sym_tmp[ix] = ELPH_EPS;
                    // This is to make sure that we donot observe any
                    // platform/compiler specific things
                }
                sym_tmp_trev[ix] = -sym_tmp[ix];
            }
            // first to cartisian units
            MatVec3f(lat_vec, vec_tmp, false, vec_tmp_trev);
            // compute S*tau
            MatVec3f(sym_tmp, vec_tmp_trev, false, vec_tmp);
            // we also negate the frac .tras. vec (just a convention used in
            // this code)
            for (int ix = 0; ix < 3; ++ix)
            {
                vec_tmp_trev[ix] = -vec_tmp[ix];
            }

            memcpy(phonon->ph_syms[isym].Rmat, sym_tmp, sizeof(ELPH_float) * 9);
            memcpy(phonon->ph_syms[isym].tau, vec_tmp, sizeof(ELPH_float) * 3);
            phonon->ph_syms[isym].time_rev = false;

            memcpy(phonon->ph_syms[isym + phonon->nph_sym].Rmat, sym_tmp_trev,
                   sizeof(ELPH_float) * 9);
            memcpy(phonon->ph_syms[isym + phonon->nph_sym].tau, vec_tmp_trev,
                   sizeof(ELPH_float) * 3);
            phonon->ph_syms[isym + phonon->nph_sym].time_rev = true;
        }
    }

    /* In case of time reversal symmetry :
     we double the number of symmetries and also use the second half of the
     symmetries */
    if (ph_tim_rev)
    {
        phonon->nph_sym *= 2;
    }

    // find the inverse index for symmetries
    // Note this is only true for the spacial symmetries.
    // Warning: The inverse of time reversal symmetries is not correct. (just
    // donot use it !)
    if (Comm->commW_rank == 0)
    {
        for (ND_int isym = 0; isym < phonon->nph_sym; ++isym)
        {
            phonon->ph_syms[isym].inv_idx =
                find_inv_symm_idx(phonon->nph_sym, phonon->ph_syms[isym].Rmat,
                                  ph_sym_mats, false);

            if (phonon->ph_syms[isym].inv_idx < 0)
            {
                error_msg("Phonon symmetries do not form a group");
            }
        }
    }

    // bcast symetries
    Bcast_symmetries(phonon->nph_sym, phonon->ph_syms, 0, Comm->commW);

    free(ph_sym_mats);
    free(ph_sym_tau);
    free(tmp_buffer);

    phonon->qpts_BZ = malloc(phonon->nq_BZ * 3 * sizeof(ELPH_float));
    CHECK_ALLOC(phonon->qpts_BZ);

    phonon->qmap = malloc(phonon->nq_BZ * 2 * sizeof(int));
    CHECK_ALLOC(phonon->qmap);

    phonon->nqstar = malloc(phonon->nq_iBZ * sizeof(ND_int));
    CHECK_ALLOC(phonon->nqstar);

    // convert to iBZ qpts to cart. coordinates.
    for (ND_int iqpt = 0; iqpt < phonon->nq_iBZ; ++iqpt)
    {
        ELPH_float* qpt_tmp = phonon->qpts_iBZ + 3 * iqpt;
        for (int ix = 0; ix < 3; ++ix)
        {
            qpt_tmp[ix] /= alat[ix];
        }
    }

    // expand to full BZ
    ND_int expand_nBZ = bz_expand(
        phonon->nq_iBZ, phonon->nph_sym, phonon->qpts_iBZ, phonon->ph_syms,
        lat_vec, phonon->qpts_BZ, phonon->nqstar, phonon->qmap);

    if (expand_nBZ != phonon->nq_BZ)
    {
        error_msg("Expansion of q-points in full BZ failed");
    }
    //// convert iBZ qpts to crystal coordinates
    for (ND_int iqpt = 0; iqpt < phonon->nq_iBZ; ++iqpt)
    {
        ELPH_float* qpt_tmp = phonon->qpts_iBZ + 3 * iqpt;
        ELPH_float qcart_tmp[3];
        memcpy(qcart_tmp, qpt_tmp, sizeof(ELPH_float) * 3);
        // convert to crystal coordinates
        MatVec3f(lat_vec, qcart_tmp, true, qpt_tmp);
    }
}
