#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/constants.h"
#include "common/cwalk/cwalk.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/parallel.h"
#include "elphC.h"
#include "io/mpi_bcast.h"
#include "parse_upf.h"
#include "qe_io.h"
#include "symmetries/symmetries.h"

void get_data_from_qe(struct Lattice* lattice, struct Phonon* phonon,
                      struct Pseudo* pseudo, const char* ph_save_dir,
                      ELPH_float* alat_param, const struct ELPH_MPI_Comms* Comm)
{
    /*
    This functions gets basic ground state data and phonon data from q.e
    that are not saved by YAMBO.
    */
    int mpi_error;

    char* tmp_buffer = NULL;
    size_t temp_str_len = 1;
    if (Comm->commW_rank == 0)
    {
        temp_str_len = 1024 + strlen(ph_save_dir);
        tmp_buffer = calloc(temp_str_len, 1);
        CHECK_ALLOC(tmp_buffer);
    }
    if (Comm->commW_rank == 0)
    {
        cwk_path_join(ph_save_dir, "dyn0", tmp_buffer, temp_str_len);
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
    char** pseudo_pots = NULL;
    char* PSEUDO_DIR = NULL;
    bool ph_tim_rev;

    ND_int natoms = 0;
    pseudo->ntype = 0;

    ELPH_float* ph_sym_mats = NULL;
    ELPH_float* ph_sym_tau = NULL;
    bool* ph_trevs = NULL;
    bool ph_mag_symm = false;

    phonon->Zborn = NULL;
    phonon->epsilon = NULL;
    phonon->Qpole = NULL;
    bool phTensors_present[3] = {false, false, false};
    // {epsilon, Zborn, Qpole}

    ELPH_float* lat_vec = lattice->alat_vec;  // a[:,i] is ith lattice vector
    if (Comm->commW_rank == 0)
    {
        cwk_path_join(ph_save_dir, "data-file-schema.xml", tmp_buffer,
                      temp_str_len);
        parse_qexml(tmp_buffer, &natoms, lat_vec, alat, &lattice->dimension,
                    &lattice->is_soc_present, &lattice->nmag, lattice->fft_dims,
                    &phonon->nph_sym, &ph_sym_mats, &ph_sym_tau, &ph_trevs,
                    &ph_mag_symm, &ph_tim_rev, &PSEUDO_DIR, &pseudo_pots,
                    &lattice->nspinor, &pseudo->ntype, &lattice->atom_type,
                    &lattice->atomic_pos);

        // free pseudo pots, no longer need
        free(PSEUDO_DIR);
        PSEUDO_DIR = NULL;

        // read dielectric and other relevent tensors from tensors.xml (if
        // exists)
        cwk_path_join(ph_save_dir, "tensors.xml", tmp_buffer, temp_str_len);
        read_ph_tensors_qe(tmp_buffer, natoms, phonon);
        // read quadrupole
        cwk_path_join(ph_save_dir, "quadrupole.fmt", tmp_buffer, temp_str_len);
        read_quadrupole_fmt(tmp_buffer, &phonon->Qpole, natoms);

        if (phonon->epsilon)
        {
            phTensors_present[0] = true;
        }
        if (phonon->Zborn)
        {
            phTensors_present[1] = true;
        }
        if (phonon->Qpole)
        {
            phTensors_present[2] = true;
        }
    }

    // Bcast all the variables
    //
    mpi_error = MPI_Bcast(&natoms, 1, ELPH_MPI_ND_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);
    //
    lattice->natom = natoms;
    lattice->nmodes = 3 * lattice->natom;

    mpi_error = MPI_Bcast(&lattice->nspinor, 1, MPI_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&pseudo->ntype, 1, ELPH_MPI_ND_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    if (Comm->commW_rank)
    {
        lattice->atomic_pos = malloc(3 * natoms * sizeof(*lattice->atomic_pos));
        CHECK_ALLOC(lattice->atomic_pos);

        lattice->atom_type = malloc(natoms * sizeof(*lattice->atom_type));
        CHECK_ALLOC(lattice->atom_type);
    }

    mpi_error = MPI_Bcast(lattice->atomic_pos, 3 * natoms, ELPH_MPI_float, 0,
                          Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(lattice->atom_type, natoms, MPI_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    //
    mpi_error = MPI_Bcast(phTensors_present, 3, MPI_C_BOOL, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(lat_vec, 9, ELPH_MPI_float, 0, Comm->commW);
    MPI_error_msg(mpi_error);
    //
    lattice->volume = fabs(det3x3(lat_vec));
    reciprocal_vecs(lat_vec, lattice->blat_vec);

    mpi_error =
        MPI_Bcast(lattice->fft_dims, 3, ELPH_MPI_ND_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(alat, 3, ELPH_MPI_float, 0, Comm->commW);
    MPI_error_msg(mpi_error);
    if (alat_param)
    {
        memcpy(alat_param, alat, sizeof(alat));
    }

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

    mpi_error = MPI_Bcast(&ph_mag_symm, 1, MPI_C_BOOL, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    ELPH_float blat[9];
    for (int ix = 0; ix < 9; ++ix)
    {
        blat[ix] = lattice->blat_vec[ix] / (2.0f * ELPH_PI);
    }

    // allocate and bcast dielectric, born charges, Quadrapoles (if present)
    // Quadrapoles not supported yet
    if (phTensors_present[0])
    {
        if (Comm->commW_rank)
        {
            phonon->epsilon = malloc(9 * sizeof(*phonon->epsilon));
            CHECK_ALLOC(phonon->epsilon);
        }
        mpi_error =
            MPI_Bcast(phonon->epsilon, 9, ELPH_MPI_float, 0, Comm->commW);
        MPI_error_msg(mpi_error);
    }
    if (phTensors_present[1])
    {
        if (Comm->commW_rank)
        {
            phonon->Zborn = malloc(9 * natoms * sizeof(*phonon->Zborn));
            CHECK_ALLOC(phonon->Zborn);
        }
        mpi_error = MPI_Bcast(phonon->Zborn, natoms * 9, ELPH_MPI_float, 0,
                              Comm->commW);
        MPI_error_msg(mpi_error);
    }
    if (phTensors_present[2])
    {
        if (Comm->commW_rank)
        {
            phonon->Qpole = malloc(27 * natoms * sizeof(*phonon->Qpole));
            CHECK_ALLOC(phonon->Qpole);
        }
        mpi_error = MPI_Bcast(phonon->Qpole, natoms * 27, ELPH_MPI_float, 0,
                              Comm->commW);
        MPI_error_msg(mpi_error);
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
            // Note we also fill the second half but are only used when no
            // magnetic symmetries and inversion are present
            ELPH_float* sym_tmp = ph_sym_mats + isym * 9;
            ELPH_float* sym_tmp_trev =
                ph_sym_mats + (isym + phonon->nph_sym) * 9;

            bool isym_trev = ph_trevs[isym];

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
                if (isym_trev)
                {
                    sym_tmp[ix] = -sym_tmp[ix];
                }
                if (fabs(sym_tmp[ix]) < ELPH_EPS)
                {
                    sym_tmp[ix] = fabs(sym_tmp[ix]);
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
            // this code) for time reversal symmetry. THis is aleady done when
            // we do S*tau if S is timerev. as S already contained -negation.
            // The trev counterpart case must be reversed
            for (int ix = 0; ix < 3; ++ix)
            {
                vec_tmp_trev[ix] = -vec_tmp[ix];
            }

            memcpy(phonon->ph_syms[isym].Rmat, sym_tmp, sizeof(ELPH_float) * 9);
            memcpy(phonon->ph_syms[isym].tau, vec_tmp, sizeof(ELPH_float) * 3);
            phonon->ph_syms[isym].time_rev = isym_trev;

            memcpy(phonon->ph_syms[isym + phonon->nph_sym].Rmat, sym_tmp_trev,
                   sizeof(ELPH_float) * 9);
            memcpy(phonon->ph_syms[isym + phonon->nph_sym].tau, vec_tmp_trev,
                   sizeof(ELPH_float) * 3);
            phonon->ph_syms[isym + phonon->nph_sym].time_rev = !isym_trev;
        }
    }

    /* In case of time reversal symmetry :
     we double the number of symmetries and also use the second half of the
     symmetries */
    if (ph_tim_rev && !ph_mag_symm)
    {
        phonon->nph_sym *= 2;
    }

    // bcast symetries
    Bcast_symmetries(phonon->nph_sym, phonon->ph_syms, 0, Comm->commW);

    free(ph_trevs);
    free(ph_sym_mats);
    free(ph_sym_tau);

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
    // Read local part of pseudo information
    pseudo->loc_pseudo = malloc(pseudo->ntype * sizeof(struct local_pseudo));
    CHECK_ALLOC(pseudo->loc_pseudo);
    if (Comm->commW_rank == 0)
    {
        for (ND_int itype = 0; itype < pseudo->ntype; ++itype)
        {
            cwk_path_join(ph_save_dir, pseudo_pots[itype], tmp_buffer,
                          temp_str_len);

            parse_upf(tmp_buffer, pseudo->loc_pseudo + itype);
            free(pseudo_pots[itype]);
        }
        free(pseudo_pots);
        pseudo_pots = NULL;
    }
    // Bcast all the pseudo information
    pseudo->ngrid_max = 0;
    for (ND_int itype = 0; itype < pseudo->ntype; ++itype)
    {
        Bcast_local_pseudo(pseudo->loc_pseudo + itype, true, 0, Comm->commW);
        ND_int ngrid_pot = pseudo->loc_pseudo[itype].ngrid;
        pseudo->ngrid_max = MAX(pseudo->ngrid_max, ngrid_pot);
    }
    free(tmp_buffer);
}
