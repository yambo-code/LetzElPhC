#include <complex.h>
#include <math.h>
#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "common/ELPH_timers.h"
#include "common/constants.h"
#include "common/cwalk/cwalk.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/free_dtypes.h"
#include "common/init_dtypes.h"
#include "common/numerical_func.h"
#include "common/parallel.h"
#include "common/print_info.h"
#include "dvloc/dvloc.h"
#include "elphC.h"
#include "interpolation.h"
#include "interpolation_utilities.h"
#include "io/io.h"
#include "io/qe/qe_io.h"
#include "parser/parser.h"
#include "phonon/phonon.h"
#include "symmetries/symmetries.h"
#include "wfc/wfc.h"
#include "wigner_seitz.h"

void interpolation_driver(const char* ELPH_input_file,
                          enum ELPH_dft_code dft_code, MPI_Comm comm_world)
{
    //
    struct interpolation_usr_input* input_data;
    read_interpolation_input_file(ELPH_input_file, &input_data, comm_world);
    //
    const char* ph_save = input_data->ph_save_dir;
    const char* ph_save_interpolated = input_data->ph_save_interpolation_dir;
    const ND_int* qgrid_new = input_data->qgrid_fine;

    init_ELPH_clocks();

    struct ELPH_MPI_Comms* mpi_comms = malloc(sizeof(struct ELPH_MPI_Comms));
    CHECK_ALLOC(mpi_comms);

    int mpi_error = MPI_SUCCESS;
    // Only plain wave plarallization is supported.
    create_parallel_comms(1, 1, comm_world, mpi_comms);

    print_ELPH_logo(mpi_comms->commW_rank, stdout);
    print_info_msg(mpi_comms->commW_rank,
                   "********** Interpolation Program started **********");

    // first create ph_interpolation directory and copy necessary file
    if (0 == mpi_comms->commW_rank)
    {
        if (dft_code == DFT_CODE_QE)
        {
            copy_ph_save_to_ph_interpolation_qe(ph_save, ph_save_interpolated);
        }
    }
    // Barrier is important to ensure we only move forward,
    // once copying is done.
    mpi_error = MPI_Barrier(mpi_comms->commW);
    MPI_error_msg(mpi_error);

    struct Lattice* lattice = malloc(sizeof(struct Lattice));
    CHECK_ALLOC(lattice);
    init_lattice_type(lattice);

    struct Phonon* phonon = malloc(sizeof(struct Phonon));
    CHECK_ALLOC(phonon);
    init_phonon_type(phonon);

    struct Pseudo* pseudo = malloc(sizeof(*pseudo));
    CHECK_ALLOC(pseudo);
    init_Pseudo_type(pseudo);

    const bool interpolate_dvscf = input_data->interpolate_dvscf;
    const bool write_dVbare = input_data->write_dVbare;

    const enum asr_kind asr_fc =
        asr_kind_from_string(input_data->asr, !mpi_comms->commW_rank);
    const enum asr_kind asr_born =
        asr_kind_from_string(input_data->zasr, !mpi_comms->commW_rank);
    //
    // this the ewald summation parameter and is always set to 1
    const ELPH_float eta_bare = 1.0;
    //
    ELPH_float* Zvals = NULL;
    ELPH_float alat_scale[3];
    if (dft_code == DFT_CODE_QE)
    {
        get_interpolation_data_from_qe(lattice, phonon, pseudo, ph_save, &Zvals,
                                       alat_scale, mpi_comms);
    }
    else
    {
        error_msg("Only qe is supported currently.");
    }
    //
    // Apply acoustic sum rule for born charges
    apply_acoustic_sum_rule_born_charges(asr_born, phonon->Zborn,
                                         lattice->natom);
    //
    // In case of 2D, we need to convert 3D quadrupoles to 2D
    if (lattice->dimension == '2')
    {
        quadrupole_3d_to_2d(lattice, phonon);
    }
    // We need atomic masses
    ELPH_float* atomic_masses = malloc(sizeof(*atomic_masses) * lattice->natom);
    CHECK_ALLOC(atomic_masses);

    ELPH_float* dummy1 = malloc(sizeof(*dummy1) * lattice->nmodes);
    CHECK_ALLOC(dummy1);

    if (dft_code == DFT_CODE_QE)
    {
        if (0 == mpi_comms->commW_rank)
        {
            ELPH_cmplx* ref_pat_basis = malloc(
                sizeof(*ref_pat_basis) * lattice->nmodes * lattice->nmodes);
            CHECK_ALLOC(ref_pat_basis);

            char read_buf[1024];
            cwk_path_join(ph_save, "dyn1", read_buf, sizeof(read_buf));
            ELPH_float qpt_tmp[3];
            ND_int iq_read = read_dyn_qe(read_buf, lattice, qpt_tmp, dummy1,
                                         ref_pat_basis, atomic_masses);
            if (iq_read != 1)
            {
                error_msg("More than 1 dynmat read.");
            }
            free(ref_pat_basis);
        }

        //
        mpi_error = MPI_Bcast(atomic_masses, lattice->natom, ELPH_MPI_float, 0,
                              mpi_comms->commW);
        MPI_error_msg(mpi_error);
    }

    // *******************************************************************
    // ******************  setup wigner seitz ****************************
    // *******************************************************************
    //
    // get the coarse q-grid
    ND_int q_grid_co[3];
    find_qpt_grid(phonon->nq_BZ, phonon->qpts_BZ, q_grid_co);
    // first setup wigner seitz vectors.
    ND_int* ws_vecs_dyn = NULL;
    ND_int* ws_degen_dyn = NULL;
    ND_int n_ws_vecs_dyn = 0;

    n_ws_vecs_dyn = build_wigner_seitz_vectors(
        q_grid_co, lattice->alat_vec, ELPH_EPS, lattice->atomic_pos,
        lattice->natom, lattice->atomic_pos, lattice->natom, &ws_vecs_dyn,
        &ws_degen_dyn);
    //
    ND_int* ws_vecs_dvscf = NULL;
    ND_int* ws_degen_dvscf = NULL;
    ND_int n_ws_vecs_dvscf = 0;
    if (interpolate_dvscf)
    {
        n_ws_vecs_dvscf = build_wigner_seitz_vectors(
            q_grid_co, lattice->alat_vec, ELPH_EPS, lattice->atomic_pos,
            lattice->natom, NULL, 0, &ws_vecs_dvscf, &ws_degen_dvscf);
    }
    // find qBZ to fft grid indices
    ND_int* indices_q2fft = malloc(2 * phonon->nq_BZ * sizeof(*indices_q2fft));
    CHECK_ALLOC(indices_q2fft);
    //
    Sorted_qpts_idxs(phonon->nq_BZ, phonon->qpts_BZ,
                     indices_q2fft + phonon->nq_BZ);
    for (ND_int iq = 0; iq < phonon->nq_BZ; ++iq)
    {
        indices_q2fft[indices_q2fft[phonon->nq_BZ + iq]] = iq;
    }
    //
    // read dvscf and eigen_vectors
    // allocate large buffers
    ELPH_cmplx* dVscfs_co = NULL;
    ELPH_cmplx* dyns_co = NULL;

    ND_int dvscf_loc_len = lattice->nmodes * lattice->nmag *
                           lattice->nfftz_loc * lattice->fft_dims[0] *
                           lattice->fft_dims[1];

    ND_int nfft_loc =
        lattice->nfftz_loc * lattice->fft_dims[0] * lattice->fft_dims[1];
    //
    if (interpolate_dvscf)
    {
        dVscfs_co = malloc(phonon->nq_BZ * dvscf_loc_len * sizeof(*dVscfs_co));
        CHECK_ALLOC(dVscfs_co);
    }

    dyns_co = malloc(phonon->nq_BZ * lattice->nmodes * lattice->nmodes *
                     sizeof(*dyns_co));
    CHECK_ALLOC(dyns_co);

    ELPH_float* omega_ph_co =
        malloc(phonon->nq_BZ * lattice->nmodes * sizeof(*omega_ph_co));
    CHECK_ALLOC(omega_ph_co);

    ELPH_float EcutRy = 30;
    bool dvscf_composite_form = false;
    // If true, then we use dvscf is V*s0 + Bx*s1 + By*s2 + Bz*s3
    // where si are pauli matrices (2x2)
    // if false, then dvscf is [V,Bx, By, Bz]
    // Always initiate to false
    bool nmags_add_long_range[4] = {false, false, false, false};
    bool only_induced_part_long_range = false;
    // In case dvscf potential contains only part from Hartree +
    // exchange-correlation then we need to remove long_range couloumb only due
    // to change density + induced potential due to dipoles and higher order
    // multipole terms. So  only_induced_part_long_range = true; means we remove
    // only long-range due to higher order multipole (dipole + quadrupole ...)
    // if false, we also need to remove longrange part from valance change i.e (
    // monopole + dipole + quadrupole...)
    if (dft_code == DFT_CODE_QE)
    {
        // q.e stores dvscf in [V,Bx,By,Bz]
        dvscf_composite_form = false;
        only_induced_part_long_range = false;
        nmags_add_long_range[0] = true;
        if (lattice->nmag == 2)
        {
            nmags_add_long_range[1] = true;
        }
    }

    // *******************************************************************
    // dvscf and dyn IO and expand to ful BZ and remove long-range parts *
    // *******************************************************************
    //
    ND_int iqpt_tmp = 0;
    for (ND_int iqco = 0; iqco < phonon->nq_iBZ; ++iqco)
    {
        // read eigen vectors and dvscf
        // we will for now read eigenvectors in dyns_co buffer
        // latter we convert it to dynamical matrices
        ND_int iq_fft_idx = indices_q2fft[iqpt_tmp];
        //
        ELPH_cmplx* dV_co_tmp =
            dVscfs_co ? dVscfs_co + iq_fft_idx * dvscf_loc_len : NULL;
        ELPH_cmplx* eigs_co =
            dyns_co + iq_fft_idx * lattice->nmodes * lattice->nmodes;
        ELPH_float* omega_co_tmp = omega_ph_co + iq_fft_idx * lattice->nmodes;
        //
        if (dft_code == DFT_CODE_QE)
        {
            get_dvscf_dyn_qe(ph_save, lattice, iqco, eigs_co, dV_co_tmp,
                             omega_co_tmp, mpi_comms);
        }
        if (dV_co_tmp)
        {
            // remore long range
            dV_add_longrange(phonon->qpts_iBZ + iqco * 3, lattice, phonon,
                             Zvals, eigs_co, dV_co_tmp, -1,
                             only_induced_part_long_range, EcutRy,
                             nmags_add_long_range, eta_bare,
                             input_data->eta_induced, mpi_comms->commK);
            // THe potential here is lattice periodic and not q peridoic.
            // so multiply with e^iqr factor to make it q periodic.
            // Before we make it q periodic, we must rotate the dvscfs
        }
        ++iqpt_tmp;
        // we will remove the long range part of dynmats later
        //
        for (ND_int istar = 1; istar < phonon->nqstar[iqco]; ++istar)
        {
            iq_fft_idx = indices_q2fft[iqpt_tmp];
            //
            ELPH_cmplx* dV_co_star =
                dVscfs_co ? dVscfs_co + iq_fft_idx * dvscf_loc_len : NULL;
            ELPH_cmplx* eigs_co_star =
                dyns_co + iq_fft_idx * lattice->nmodes * lattice->nmodes;
            ELPH_float* omega_co_tmp_star =
                omega_ph_co + iq_fft_idx * lattice->nmodes;

            //
            ND_int iq_iBZ = phonon->qmap[2 * iqpt_tmp];
            ND_int idx_qsym = phonon->qmap[2 * iqpt_tmp + 1];
            if (iq_iBZ != iqco)
            {
                error_msg("Wrong qstar.");
            }
            // rotate dvscf

            memcpy(omega_co_tmp_star, omega_co_tmp,
                   lattice->nmodes * sizeof(*omega_ph_co));
            //
            struct symmetry* sym_star = phonon->ph_syms + idx_qsym;
            rotate_eig_vecs(sym_star, lattice, phonon->qpts_iBZ + iqco * 3,
                            eigs_co, eigs_co_star);
            if (dVscfs_co)
            {
                rotate_dvscf(dV_co_tmp, sym_star, lattice, dvscf_composite_form,
                             dV_co_star, mpi_comms->commK);
            }
            ++iqpt_tmp;
        }
    }

    //
    // *******************************************************************
    // ******************  q->R (dvscf) **********************************
    // *******************************************************************
    if (dVscfs_co)
    {
        for (ND_int i = 0; i < phonon->nq_BZ; ++i)
        {
            ND_int iq = indices_q2fft[i];
            ELPH_cmplx* rot_vecs =
                dyns_co + iq * lattice->nmodes * lattice->nmodes;
            // remove mass normalization in the dynmats
            mass_normalize_pol_vecs(atomic_masses, lattice->nmodes,
                                    lattice->natom, 1.0, rot_vecs);

            dVscf_change_basis(dVscfs_co + iq * dvscf_loc_len, rot_vecs, 1,
                               lattice->nmodes, lattice->nmag,
                               lattice->fft_dims[0], lattice->fft_dims[1],
                               lattice->nfftz_loc, 'C');
            // get back to previous normalization
            mass_normalize_pol_vecs(atomic_masses, lattice->nmodes,
                                    lattice->natom, -1.0, rot_vecs);
            //
            // Now remove the e^-iqr factor to make the dvscf q periodic
            // 1) get the rotated q point
            ND_int iq_iBZ = phonon->qmap[2 * i];
            ND_int idx_qsym = phonon->qmap[2 * i + 1];
            //
            ELPH_float tmp_qpt[3], qpt_cart_iq[3];
            MatVec3f(lattice->blat_vec, phonon->qpts_iBZ + iq_iBZ * 3, false,
                     tmp_qpt);
            MatVec3f(phonon->ph_syms[idx_qsym].Rmat, tmp_qpt, false,
                     qpt_cart_iq);
            // remove 2*pi
            qpt_cart_iq[0] /= (2 * ELPH_PI);
            qpt_cart_iq[1] /= (2 * ELPH_PI);
            qpt_cart_iq[2] /= (2 * ELPH_PI);
            // convert to crystal units
            MatVec3f(lattice->alat_vec, qpt_cart_iq, true, tmp_qpt);
            //
            // 2) multiply with e^iqr
            multiply_eikr(dVscfs_co + iq * dvscf_loc_len, tmp_qpt, lattice,
                          lattice->nmodes * lattice->nmag, 1);
        }
        //
        // Now perform fourier transform
        fft_q2R(dVscfs_co, q_grid_co, dvscf_loc_len);
    }

    // Get a grid to perform summation
    ND_int Ggrid_phonon[3];
    get_fft_box(EcutRy, lattice->blat_vec, Ggrid_phonon, mpi_comms->commK);

    // first compute the long_range asr term
    ELPH_cmplx* dyn_mat_asr_lr =
        calloc(lattice->natom * 9, sizeof(*dyn_mat_asr_lr));
    CHECK_ALLOC(dyn_mat_asr_lr);

    compute_dyn_lr_asr_correction(lattice, phonon, Ggrid_phonon,
                                  input_data->eta_ph, dyn_mat_asr_lr);

    ELPH_cmplx* fc_lr = NULL;
    if (asr_fc == ASR_ALL || asr_fc == ASR_ALL_HUANG)
    {
        // In case, we apply rotational sum rules, then we also need long range
        // force constants.
        // See C. Lin et al.: npj Comput. Mater. 8, 236 (2022)
        fc_lr = malloc(phonon->nq_BZ * lattice->nmodes * lattice->nmodes *
                       sizeof(*fc_lr));
        CHECK_ALLOC(fc_lr);
    }

    // *******************************************************************
    // ******************  Polarization vectors to dynmats ***************
    // *******************************************************************
    //
    // FIX ME need to parallize
    for (ND_int i = 0; i < phonon->nq_BZ; ++i)
    {
        ND_int iq = indices_q2fft[i];
        // convert polarization vectors to dynamical matrix and remove long
        // range part
        ELPH_cmplx* pol_vecs_iq =
            dyns_co + iq * lattice->nmodes * lattice->nmodes;

        ELPH_cmplx* fc_lr_iq = NULL;

        pol_vecs_to_dyn(omega_ph_co + iq * lattice->nmodes, lattice->natom,
                        atomic_masses, pol_vecs_iq);
        // Here dynamical matrices have 1/sqrt(Ma*Mb) factor
        if (fc_lr)
        {
            fc_lr_iq = fc_lr + iq * lattice->nmodes * lattice->nmodes;
            memcpy(fc_lr_iq, pol_vecs_iq,
                   lattice->nmodes * lattice->nmodes * sizeof(*fc_lr));
        }
        const ELPH_float* qpt_iq_tmp = phonon->qpts_BZ + 3 * i;
        // remove long range part
        // Note the long-range part also contains 1/sqrt(Ma*Mb)
        add_ph_dyn_long_range(qpt_iq_tmp, lattice, phonon, Ggrid_phonon, -1,
                              atomic_masses, dyn_mat_asr_lr, input_data->eta_ph,
                              pol_vecs_iq);
        //
        // for long range force constants
        if (fc_lr)
        {
            for (ND_int ii = 0; ii < (lattice->nmodes * lattice->nmodes); ++ii)
            {
                fc_lr_iq[ii] -= pol_vecs_iq[ii];
            }
        }
    }

    // *******************************************************************
    // ******************  q->R (dyn to fc) ******************************
    // *******************************************************************
    //
    // fourier transform phonons
    fft_q2R(dyns_co, q_grid_co, lattice->nmodes * lattice->nmodes);
    //
    // IN case of rotational sum rules, we also need to add long range force
    // constants
    if (fc_lr)
    {
        fft_q2R(fc_lr, q_grid_co, lattice->nmodes * lattice->nmodes);
        //
        for (ND_int ii = 0;
             ii < (phonon->nq_BZ * lattice->nmodes * lattice->nmodes); ++ii)
        {
            dyns_co[ii] += fc_lr[ii];
        }
    }

    // *******************************************************************
    // ****************** Acoustic sum rule ******************************
    // *******************************************************************
    //
    // Here the force constant are mass normalized.
    // i.e \tilde{C}(R)_{ab} = 1/sqrt(Ma*Mb) * C(R)_{ab}, so we remove
    // normalization
    mass_normalize_force_constants(atomic_masses, phonon->nq_BZ, lattice->natom,
                                   0.5, dyns_co);

    // Apply Acoustic sum rule for force constants
    apply_acoustic_sum_rule_fc(asr_fc, q_grid_co, lattice->natom, dyns_co,
                               lattice->atomic_pos, lattice->alat_vec,
                               ws_vecs_dyn, n_ws_vecs_dyn, ws_degen_dyn);

    // Normalize the force constants, so that we donot need to normalize it
    // for every dynamical matrix
    mass_normalize_force_constants(atomic_masses, phonon->nq_BZ, lattice->natom,
                                   -0.5, dyns_co);
    //
    // now remove the long range force constants incase they are added
    if (fc_lr)
    {
        for (ND_int ii = 0;
             ii < (phonon->nq_BZ * lattice->nmodes * lattice->nmodes); ++ii)
        {
            dyns_co[ii] -= fc_lr[ii];
        }
    }
    // free fc_lr, no longer needed.
    free(fc_lr);
    fc_lr = NULL;
    //
    //
    ND_int nqpts_to_interpolate = qgrid_new[0] * qgrid_new[1] * qgrid_new[2];
    // this will be over written lattern with number of qpts in iBZ

    ELPH_float* qpts_interpolation = NULL;
    // in crystal coordinates

    // *******************************************************************
    // ****************** Generate/read q-points *************************
    // *******************************************************************
    //
    bool qpts_usr_provide = false;
    if (0 == strlen(input_data->qlist_file))
    {
        qpts_interpolation =
            malloc(sizeof(*qpts_interpolation) * 3 * nqpts_to_interpolate);
        CHECK_ALLOC(qpts_interpolation);
        // in crystal coordinates
        struct symmetry* symms_iBZexpand =
            input_data->nosym ? NULL : phonon->ph_syms;
        nqpts_to_interpolate = generate_iBZ_kpts(
            qgrid_new, phonon->nph_sym, symms_iBZexpand, lattice->alat_vec,
            lattice->blat_vec, qpts_interpolation, true);
    }
    else
    {
        qpts_usr_provide = true;
        if (0 == mpi_comms->commW_rank)
        {
            qpts_interpolation = parse_qpt_entries(input_data->qlist_file,
                                                   &nqpts_to_interpolate);
            if (!qpts_interpolation)
            {
                error_msg("Reading qpoint file failed.");
            }
        }

        mpi_error = MPI_Bcast(&nqpts_to_interpolate, 1, ELPH_MPI_ND_INT, 0,
                              mpi_comms->commW);
        MPI_error_msg(mpi_error);

        if (0 != mpi_comms->commW_rank)
        {
            qpts_interpolation =
                malloc(sizeof(*qpts_interpolation) * 3 * nqpts_to_interpolate);
            CHECK_ALLOC(qpts_interpolation);
        }

        mpi_error = MPI_Bcast(qpts_interpolation, 3 * nqpts_to_interpolate,
                              ELPH_MPI_float, 0, mpi_comms->commW);
        MPI_error_msg(mpi_error);
    }
    // local part
    ELPH_cmplx* Vlocr = NULL;

    // netcdf variable for writing dVbare (only used when the user asks to
    // write)
    int ncid_dVbare = 0, nc_err = 0;
    int ncvar_dVbare = 0, ncvar_ph_freq = 0, ncvar_ph_eig = 0;
    //
    // *******************************************************************
    // ******************** dVbare IO ************************************
    // *******************************************************************
    //
    if (write_dVbare)
    {
        // In case the user wants to dum dVbare, we contruct it. Note that
        // This is not actually added. Instead, in the long_range_term, the
        // long_range monopole term is added to completely avoid reconstruting
        // the full bare
        Vlocr = malloc(sizeof(*Vlocr) * lattice->nmodes * nfft_loc);
        CHECK_ALLOC(Vlocr);
        // buffer to store local part of the pseudo potential
        //  compute the Vlocg table
        //// first find the qmax for vloc table
        ELPH_float qmax_val = fabs(qpts_interpolation[0]);
        for (ND_int imax = 0; imax < (nqpts_to_interpolate * 3); ++imax)
        {
            if (fabs(qpts_interpolation[imax]) > qmax_val)
            {
                qmax_val = fabs(qpts_interpolation[imax]);
            }
        }
        // Note this needs to be set before compute the Vlocg table else U.B
        pseudo->vloc_table->qmax_abs = ceil(fabs(qmax_val)) + 1;
        create_vlocg_table(lattice, pseudo, mpi_comms);

        // Create a netcdf and define netcdf dimensions.
        // World comm must open it (even though for now it is same as CommK and
        // commQ)
        if ((nc_err =
                 nc_create_par("ndb.dVbare", NC_NETCDF4 | NC_CLOBBER,
                               mpi_comms->commW, MPI_INFO_NULL, &ncid_dVbare)))
        {
            ERR(nc_err);
        }
        // Donot do prefilling
        if ((nc_err = ncsetfill(ncid_dVbare, NC_NOFILL)))
        {
            fprintf(stderr, "Error setting nc_fill to ndb.dVbare file.");
            ERR(nc_err);
        }
        // Define variables (dVbare, freq, eigs, qpts_reduced)
        // dVbare
        ND_int dims[6] = {nqpts_to_interpolate, lattice->nmodes,
                          lattice->fft_dims[0], lattice->fft_dims[1],
                          lattice->fft_dims[2], 2};
        //
        size_t chunksize[6] = {1, 1, lattice->fft_dims[0], lattice->fft_dims[1],
                               1, 2};
        def_ncVar(ncid_dVbare, &ncvar_dVbare, 6, ELPH_NC4_IO_FLOAT, dims,
                  "dVbare_local",
                  (char*[]){"nq", "nmodes", "Nx", "Ny", "Nz", "re_im"},
                  chunksize);
        // Collective IO
        if ((nc_err =
                 nc_var_par_access(ncid_dVbare, ncvar_dVbare, NC_COLLECTIVE)))
        {
            ERR(nc_err);
        }
        // freq
        def_ncVar(ncid_dVbare, &ncvar_ph_freq, 2, ELPH_NC4_IO_FLOAT, dims,
                  "FREQ", (char*[]){"nq", "nmodes"}, NULL);
        // eigs
        dims[2] = lattice->natom;
        dims[3] = 3;
        dims[4] = 2;
        def_ncVar(ncid_dVbare, &ncvar_ph_eig, 5, ELPH_NC4_IO_FLOAT, dims,
                  "POLARIZATION_VECTORS",
                  (char*[]){"nq", "nmodes", "atom", "pol", "re_im"}, NULL);
        // qpts_reduced
        int ncvar_qpt_tmp = 0;
        dims[1] = 3;
        def_ncVar(ncid_dVbare, &ncvar_qpt_tmp, 2, ELPH_NC4_IO_FLOAT, dims,
                  "qpoints", (char*[]){"nq", "pol"}, NULL);
        // write qpoints now
        if (0 == mpi_comms->commW_rank)
        {
            // qpts in reduced units
            if ((nc_err = nc_put_var(ncid_dVbare, ncvar_qpt_tmp,
                                     qpts_interpolation)))
            {
                ERR(nc_err);
            }
        }
    }

    // allocations
    ELPH_cmplx* dvscf_interpolated = NULL;
    if (dVscfs_co)
    {
        dvscf_interpolated =
            malloc(dvscf_loc_len * sizeof(*dvscf_interpolated));
        CHECK_ALLOC(dvscf_interpolated);
    }

    ELPH_cmplx* dyn_interpolated =
        malloc(sizeof(*dyn_interpolated) * lattice->nmodes * lattice->nmodes);
    CHECK_ALLOC(dyn_interpolated);

    ELPH_float* ph_freq_iq_interp = NULL;
    //
    if (write_dVbare)
    {
        ph_freq_iq_interp =
            malloc(lattice->nmodes * sizeof(*ph_freq_iq_interp));
        CHECK_ALLOC(ph_freq_iq_interp);
    }

    // *******************************************************************
    // ************************** LO-TO setup ****************************
    // *******************************************************************
    const ELPH_float loto_mag = 1e-5;
    // eps for q->0
    //
    if (input_data->loto)
    {
        ELPH_float loto_dir_norm =
            sqrt(dot3_macro(input_data->loto_dir, input_data->loto_dir));
        // No LO-TO splitting for user provided qpts. The user must instead
        // provide add a small shift in the direction
        if (qpts_usr_provide)
        {
            if (mpi_comms->commW_rank == 0)
            {
                fprintf(stdout,
                        "Warning: When qlist_file is provided, the loto and "
                        "loto_dir variables are ignored. "
                        "Please refer to qlist_file variable documentation for "
                        "more details.\n");
            }
            input_data->loto = false;
        }
        else if (loto_dir_norm < ELPH_EPS)
        {
            if (mpi_comms->commW_rank == 0)
            {
                fprintf(stdout,
                        "Warning: Too small LO-TO direction vector, Turning "
                        "off LO-TO splitting.\n");
            }
            input_data->loto = false;
        }
        else
        {
            ELPH_float tmp_loto_dir[3];
            for (ND_int ix = 0; ix < 3; ++ix)
            {
                tmp_loto_dir[ix] =
                    loto_mag * input_data->loto_dir[ix] / loto_dir_norm;
            }
            MatVec3f(lattice->alat_vec, tmp_loto_dir, true,
                     input_data->loto_dir);
        }
    }

    // *******************************************************************
    // ************************** Interpolation loop  ********************
    // *******************************************************************
    //
    // now interpolate
    for (ND_int iq = 0; iq < nqpts_to_interpolate; ++iq)
    {
        char read_buf[1024];
        char dvscf_dyn_name[32];

        ELPH_float* qpt_interpolate = qpts_interpolation + 3 * iq;

        ELPH_float qpt_interpolate_cart[3];
        MatVec3f(lattice->blat_vec, qpt_interpolate, false,
                 qpt_interpolate_cart);
        //
        // *******************************************************************
        // ******************** dVscf interpolation **************************
        // *******************************************************************
        if (dVscfs_co)
        {
            snprintf(dvscf_dyn_name, sizeof(dvscf_dyn_name), "dvscf%lld",
                     (long long)(iq + 1));
            cwk_path_join(ph_save_interpolated, dvscf_dyn_name, read_buf,
                          sizeof(read_buf));

            fft_R2q_dvscf(dVscfs_co, qpt_interpolate, q_grid_co, lattice->natom,
                          dvscf_loc_len / lattice->nmodes, ws_vecs_dvscf,
                          n_ws_vecs_dvscf, ws_degen_dvscf, dvscf_interpolated);
            //
            // remove e^iqr phase to make it lattice periodic
            multiply_eikr(dvscf_interpolated, qpt_interpolate, lattice,
                          lattice->nmodes * lattice->nmag, -1);
            // add long range back
            //
            dV_add_longrange(qpt_interpolate, lattice, phonon, Zvals, NULL,
                             dvscf_interpolated, 1,
                             only_induced_part_long_range, EcutRy,
                             nmags_add_long_range, eta_bare,
                             input_data->eta_induced, mpi_comms->commK);

            // write to file
            if (dft_code == DFT_CODE_QE)
            {
                // we have dvscf in cart basis and we want to store them
                // in cart basis
                write_dvscf_qe(read_buf, lattice, dvscf_interpolated,
                               mpi_comms->commK);
                // write patterns
                // Write placeholder pattern file indicating that it is an
                // identity matrix.
                if (0 == mpi_comms->commQ_rank)
                {
                    snprintf(dvscf_dyn_name, sizeof(dvscf_dyn_name),
                             "patterns.%lld.xml", (long long)(iq + 1));
                    cwk_path_join(ph_save_interpolated, dvscf_dyn_name,
                                  read_buf, sizeof(read_buf));
                    write_identity_patterns_xml(read_buf);
                }
            }
        }

        // *******************************************************************
        // ******************** Phonon interpolation *************************
        // *******************************************************************
        //
        fft_R2q_dyn(dyns_co, qpt_interpolate, q_grid_co, lattice->natom,
                    ws_vecs_dyn, n_ws_vecs_dyn, ws_degen_dyn, dyn_interpolated);

        // in case of LO-TO, add small vector to gamma
        bool loto_eps_added = false;
        if (input_data->loto)
        {
            ELPH_float sum = 0.0;
            for (ND_int ix = 0; ix < 3; ++ix)
            {
                ELPH_float iq_diff_tmp =
                    qpt_interpolate[ix] - rint(qpt_interpolate[ix]);
                sum += iq_diff_tmp * iq_diff_tmp;
            }
            sum = sqrt(sum);
            if (sum < ELPH_EPS)
            {
                loto_eps_added = true;
                for (ND_int ix = 0; ix < 3; ++ix)
                {
                    qpt_interpolate[ix] += input_data->loto_dir[ix];
                }
            }
        }
        // add back the long range part
        add_ph_dyn_long_range(qpt_interpolate, lattice, phonon, Ggrid_phonon, 1,
                              atomic_masses, dyn_mat_asr_lr, input_data->eta_ph,
                              dyn_interpolated);
        //
        if (loto_eps_added)
        {
            for (ND_int ix = 0; ix < 3; ++ix)
            {
                qpt_interpolate[ix] -= input_data->loto_dir[ix];
            }
        }

        if ((!write_dVbare || dVscfs_co) && 0 == mpi_comms->commW_rank)
        {
            // interpolate dyn file
            snprintf(dvscf_dyn_name, sizeof(dvscf_dyn_name), "dyn%lld",
                     (long long)(iq + 1));
            cwk_path_join(ph_save_interpolated, dvscf_dyn_name, read_buf,
                          sizeof(read_buf));
            // write dyn file
            if (dft_code == DFT_CODE_QE)
            {
                ELPH_float qpt_tmp_iq[3];
                for (ND_int ix = 0; ix < 3; ++ix)
                {
                    qpt_tmp_iq[ix] = alat_scale[ix] * qpt_interpolate_cart[ix] /
                                     (2 * ELPH_PI);
                }
                write_dyn_qe(read_buf, lattice->natom, qpt_tmp_iq,
                             dyn_interpolated, atomic_masses);
            }
        }
        // *******************************************************************
        // ******************** dVbare interpolation *************************
        // *******************************************************************
        //
        if (write_dVbare)
        {
            // Symmetrize the matrix
            for (ND_int idim1 = 0; idim1 < lattice->nmodes; idim1++)
            {
                for (ND_int jdim1 = 0; jdim1 <= idim1; jdim1++)
                {
                    dyn_interpolated[idim1 * lattice->nmodes + jdim1] =
                        0.5 *
                        (dyn_interpolated[idim1 * lattice->nmodes + jdim1] +
                         conj(dyn_interpolated[jdim1 * lattice->nmodes +
                                               idim1]));
                }
            }
            // Note:
            // We should diagonalize on on single cpu and BCast eig_vecs. This
            // is because, in case of small numerical diff, we might get
            // different eig values on different cpus which will be diaster.
            if (0 == mpi_comms->commQ_rank)
            {
                int lpack_info = diagonalize_hermitian(
                    'V', 'U', lattice->nmodes, lattice->nmodes,
                    dyn_interpolated, ph_freq_iq_interp);
                if (lpack_info)
                {
                    error_msg("Error diagonalizing dynamical matrix");
                }
            }
            mpi_error = MPI_Bcast(ph_freq_iq_interp, lattice->nmodes,
                                  ELPH_MPI_float, 0, mpi_comms->commQ);
            MPI_error_msg(mpi_error);

            mpi_error =
                MPI_Bcast(dyn_interpolated, lattice->nmodes * lattice->nmodes,
                          ELPH_MPI_cmplx, 0, mpi_comms->commQ);
            MPI_error_msg(mpi_error);

            //
            for (ND_int imode = 0; imode < lattice->nmodes; imode++)
            {
                ELPH_float tmp_fq = sqrt(fabs(ph_freq_iq_interp[imode]));
                if (ph_freq_iq_interp[imode] < 0)
                {
                    ph_freq_iq_interp[imode] = -tmp_fq;
                }
                else
                {
                    ph_freq_iq_interp[imode] = tmp_fq;
                }

                ELPH_cmplx* eig_tmp_ptr =
                    dyn_interpolated + imode * lattice->nmodes;
                for (ND_int jmode = 0; jmode < lattice->nmodes; jmode++)
                {
                    ND_int ia = jmode / 3;
                    eig_tmp_ptr[jmode] /= sqrt(atomic_masses[ia]);
                }
            }
            dVlocq(qpt_interpolate, lattice, pseudo, dyn_interpolated, Vlocr,
                   mpi_comms->commW);
            // now write

            size_t startp[6] = {iq, 0, 0, 0, lattice->nfftz_loc_shift, 0};
            size_t countp[6] = {1,
                                lattice->nmodes,
                                lattice->fft_dims[0],
                                lattice->fft_dims[1],
                                lattice->nfftz_loc,
                                2};
            //
            if ((nc_err = nc_put_vara(ncid_dVbare, ncvar_dVbare, startp, countp,
                                      Vlocr)))
            {
                ERR(nc_err);
            }
            if (0 == mpi_comms->commW_rank)
            {
                startp[4] = 0;
                countp[2] = lattice->natom;
                countp[3] = 3;
                countp[4] = 2;
                if ((nc_err = nc_put_vara(ncid_dVbare, ncvar_ph_eig, startp,
                                          countp, dyn_interpolated)))
                {
                    ERR(nc_err);
                }
                if ((nc_err = nc_put_vara(ncid_dVbare, ncvar_ph_freq, startp,
                                          countp, ph_freq_iq_interp)))
                {
                    ERR(nc_err);
                }
            }
        }
        if (dft_code == DFT_CODE_QE)
        {
            // convert the interpolated qpoints in weird q.e units
            // ELPH_float* qpt_interpolate = qpts_interpolation + 3 * iq;
            for (ND_int ix = 0; ix < 3; ++ix)
            {
                qpt_interpolate[ix] =
                    alat_scale[ix] * qpt_interpolate_cart[ix] / (2 * ELPH_PI);
            }
        }
    }
    //
    //
    // Write dyn0 file in case required
    if ((!write_dVbare || dVscfs_co) && 0 == mpi_comms->commW_rank &&
        !qpts_usr_provide)
    {
        char read_buf[1024];
        cwk_path_join(ph_save_interpolated, "dyn0", read_buf, sizeof(read_buf));

        write_qpts_qe(read_buf, nqpts_to_interpolate, qpts_interpolation,
                      qgrid_new);
    }

    // close the netcdf file incase opened
    if (write_dVbare)
    {
        if ((nc_err = nc_close(ncid_dVbare)))
        {
            ERR(nc_err);
        }
    }

    // *******************************************************************
    // ************************** Cleanup ********************************
    // *******************************************************************

    free(ph_freq_iq_interp);

    free(ws_vecs_dvscf);
    free(ws_degen_dvscf);

    free(ws_vecs_dyn);
    free(ws_degen_dyn);

    free(dyn_mat_asr_lr);
    free(indices_q2fft);
    free(dvscf_interpolated);
    free(dyn_interpolated);
    free(qpts_interpolation);
    int World_rank_tmp = mpi_comms->commW_rank;

    free(atomic_masses);
    free(dummy1);
    free(Vlocr);

    free(omega_ph_co);
    free(dVscfs_co);
    free(dyns_co);

    free(Zvals);
    // free user input
    free_interpolation_usr_input(input_data);
    //
    free_Pseudo_type(pseudo);
    free_lattice_type(lattice);
    free_phonon_type(phonon);

    free(lattice);
    free(phonon);
    free(pseudo);
    //
    free_parallel_comms(mpi_comms);
    free(mpi_comms);
    // From here nothing should happen apart from last minite things such
    // as printing clocks
    //

    if (0 == World_rank_tmp)
    {
        print_ELPH_clock_summary();
    }
    // cleanup the clocks
    cleanup_ELPH_clocks();
    // done with the calculation
    print_info_msg(World_rank_tmp,
                   "********** Interpolation Program ended **********");

    return;
}
