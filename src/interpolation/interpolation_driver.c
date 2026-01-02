#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "common/ELPH_timers.h"
#include "common/constants.h"
#include "common/cwalk/cwalk.h"
#include "common/dtypes.h"
#include "common/error.h"
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

    struct Lattice* lattice = malloc(sizeof(struct Lattice));
    CHECK_ALLOC(lattice);
    init_lattice_type(lattice);

    struct Phonon* phonon = malloc(sizeof(struct Phonon));
    CHECK_ALLOC(phonon);
    init_phonon_type(phonon);

    struct Pseudo* pseudo = malloc(sizeof(*pseudo));
    CHECK_ALLOC(pseudo);
    init_Pseudo_type(pseudo);

    bool interpolate_dvscf = input_data->interpolate_dvscf;
    bool write_dVbare = false;

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
    //
    // We need atomic masses
    ELPH_float* atomic_masses = malloc(sizeof(*atomic_masses) * lattice->natom);
    CHECK_ALLOC(atomic_masses);

    ELPH_float* dummy1 = malloc(sizeof(*dummy1) * lattice->nmodes);
    CHECK_ALLOC(dummy1);

    // These are reference patern basis.
    ELPH_cmplx* ref_pat_basis =
        malloc(sizeof(*ref_pat_basis) * lattice->nmodes * lattice->nmodes);
    CHECK_ALLOC(ref_pat_basis);

    if (dft_code == DFT_CODE_QE)
    {
        if (0 == mpi_comms->commW_rank)
        {
            char read_buf[1024];
            cwk_path_join(ph_save, "dyn1", read_buf, sizeof(read_buf));
            ELPH_float qpt_tmp[3];
            ND_int iq_read = read_dyn_qe(read_buf, lattice, qpt_tmp, dummy1,
                                         ref_pat_basis, atomic_masses);
            if (iq_read != 1)
            {
                error_msg("More than 1 dynmat read.");
            }

            if (interpolate_dvscf)
            {
                // read the first pattern file
                cwk_path_join(ph_save, "patterns.1.xml", read_buf,
                              sizeof(read_buf));
                read_pattern_qe(read_buf, lattice, ref_pat_basis);
            }
        }

        //
        mpi_error = MPI_Bcast(atomic_masses, lattice->natom, ELPH_MPI_float, 0,
                              mpi_comms->commW);
        MPI_error_msg(mpi_error);

        mpi_error = MPI_Bcast(ref_pat_basis, lattice->nmodes * lattice->nmodes,
                              ELPH_MPI_cmplx, 0, mpi_comms->commW);
        MPI_error_msg(mpi_error);

        MPI_error_msg(mpi_error);
    }

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
    // local part
    ELPH_cmplx* Vlocr = NULL;

    if (interpolate_dvscf || write_dVbare)
    {
        Vlocr = malloc(sizeof(*Vlocr) * lattice->nmodes * nfft_loc);
        // buffer to store local part of the pseudo potential
        CHECK_ALLOC(Vlocr);
        //  compute the Vlocg table
        //// first find the qmax for vloc table
        ELPH_float qmax_val = fabs(phonon->qpts_iBZ[0]);
        for (ND_int imax = 0; imax < (phonon->nq_iBZ * 3); ++imax)
        {
            if (fabs(phonon->qpts_iBZ[imax]) > qmax_val)
            {
                qmax_val = fabs(phonon->qpts_iBZ[imax]);
            }
        }
        // Note this needs to be set before compute the Vlocg table else U.B
        pseudo->vloc_table->qmax_abs = ceil(fabs(qmax_val)) + 1;
        create_vlocg_table(lattice, pseudo, mpi_comms);
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
    bool only_induced_part_long_range = true;
    bool add_dVbare = false;
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
        only_induced_part_long_range = true;
        nmags_add_long_range[0] = true;
        if (lattice->nmag == 2)
        {
            nmags_add_long_range[1] = true;
        }
        add_dVbare = true;
    }

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
            // add the local part incase required
            if (add_dVbare)
            {
                dVlocq(phonon->qpts_iBZ + iqco * 3, lattice, pseudo, eigs_co,
                       Vlocr, mpi_comms->commK);
                // add local part from nuclei
                // lattice->nmodes  * nfft_loc
                // lattice->nmodes * lattice->nmag
                for (ND_int imode = 0; imode < lattice->nmodes; ++imode)
                {
                    ELPH_cmplx* Vlocr_tmp = Vlocr + imode * nfft_loc;
                    for (ND_int imag = 0; imag < lattice->nmag; ++imag)
                    {
                        if (nmags_add_long_range[imag])
                        {
                            ELPH_cmplx* dvscf_tmp =
                                dV_co_tmp +
                                (imode * lattice->nmag + imag) * nfft_loc;
                            for (ND_int ift = 0; ift < nfft_loc; ++ift)
                            {
                                dvscf_tmp[ift] += Vlocr_tmp[ift];
                            }
                        }
                    }
                }
            }
            // remore long range
            dV_add_longrange(phonon->qpts_iBZ + iqco * 3, lattice, phonon,
                             Zvals, eigs_co, dV_co_tmp, -1,
                             only_induced_part_long_range, EcutRy,
                             nmags_add_long_range, input_data->eta_bare,
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

    compute_dyn_lr_asr_correction(lattice, phonon, Ggrid_phonon, atomic_masses,
                                  input_data->eta_ph, dyn_mat_asr_lr);

    // FIX ME need to parallize
    for (ND_int i = 0; i < phonon->nq_BZ; ++i)
    {
        ND_int iq = indices_q2fft[i];
        // convert polarization vectors to dynamical matrix and remove long
        // range part
        ELPH_cmplx* pol_vecs_iq =
            dyns_co + iq * lattice->nmodes * lattice->nmodes;

        pol_vecs_to_dyn(omega_ph_co + iq * lattice->nmodes, lattice->natom,
                        atomic_masses, pol_vecs_iq);

        const ELPH_float* qpt_iq_tmp = phonon->qpts_BZ + 3 * i;
        // remove long range part
        add_ph_dyn_long_range(qpt_iq_tmp, lattice, phonon, Ggrid_phonon, -1,
                              atomic_masses, dyn_mat_asr_lr, input_data->eta_ph,
                              pol_vecs_iq);
        //
    }

    // fourier transform phonons
    fft_q2R(dyns_co, q_grid_co, lattice->nmodes * lattice->nmodes);
    //
    ND_int nqpts_to_interpolate = qgrid_new[0] * qgrid_new[1] * qgrid_new[2];
    // this will be over written lattern with number of qpts in iBZ

    ELPH_float* qpts_interpolation =
        malloc(sizeof(*qpts_interpolation) * 3 * nqpts_to_interpolate);
    // in crystal coordinates
    nqpts_to_interpolate = generate_iBZ_kpts(
        qgrid_new, phonon->nph_sym, phonon->ph_syms, lattice->alat_vec,
        lattice->blat_vec, qpts_interpolation, true);

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
            // change to pattern basis
            // In principle, we should first construct pattern basis based on
            // the litle group of q but since we only use dvscfs internally, we
            // multiply with the dvscfs in cart basis to pattern basis with the
            // first pattern basis as when computing electron-phonon we again
            // unstrip the pattern basis, so it does not effect what we
            // use. This allows to skip a function to write pattern.xml
            // file as both are not needed for our purposes.
            dVscf_change_basis(dvscf_interpolated, ref_pat_basis, 1,
                               lattice->nmodes, lattice->nmag,
                               lattice->fft_dims[0], lattice->fft_dims[1],
                               lattice->nfftz_loc, 'N');
            //
            // remove e^iqr phase to make it lattice periodic
            multiply_eikr(dvscf_interpolated, qpt_interpolate, lattice,
                          lattice->nmodes * lattice->nmag, -1);
            // add long range back
            //
            dV_add_longrange(qpt_interpolate, lattice, phonon, Zvals,
                             ref_pat_basis, dvscf_interpolated, 1,
                             only_induced_part_long_range, EcutRy,
                             nmags_add_long_range, input_data->eta_bare,
                             input_data->eta_induced, mpi_comms->commK);

            // write to file
            if (dft_code == DFT_CODE_QE)
            {
                write_dvscf_qe(read_buf, lattice, dvscf_interpolated,
                               mpi_comms->commK);
            }
        }

        if (0 == mpi_comms->commW_rank)
        {
            // interpolate dyn file
            snprintf(dvscf_dyn_name, sizeof(dvscf_dyn_name), "dyn%lld",
                     (long long)(iq + 1));
            cwk_path_join(ph_save_interpolated, dvscf_dyn_name, read_buf,
                          sizeof(read_buf));

            fft_R2q_dyn(dyns_co, qpt_interpolate, q_grid_co, lattice->natom,
                        ws_vecs_dyn, n_ws_vecs_dyn, ws_degen_dyn,
                        dyn_interpolated);

            // add back the long range part
            add_ph_dyn_long_range(qpt_interpolate, lattice, phonon,
                                  Ggrid_phonon, 1, atomic_masses,
                                  dyn_mat_asr_lr, input_data->eta_ph,
                                  dyn_interpolated);
            // write dyn file
            // FIX me write qpoint in alat units q.e
            // Convert qpoints to
            if (dft_code == DFT_CODE_QE)
            {
                write_dyn_qe(read_buf, lattice->natom, qpt_interpolate,
                             dyn_interpolated, atomic_masses);
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

    if (0 == mpi_comms->commW_rank)
    {
        char read_buf[1024];
        cwk_path_join(ph_save_interpolated, "dyn0", read_buf, sizeof(read_buf));

        write_qpts_qe(read_buf, nqpts_to_interpolate, qpts_interpolation,
                      qgrid_new);
    }

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
    free(ref_pat_basis);
    free(Vlocr);

    free(omega_ph_co);
    free(dVscfs_co);
    free(dyns_co);

    free(Zvals);
    // free user input
    free_interpolation_usr_input(input_data);
    //
    free_save_data(NULL, lattice, pseudo, phonon);
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
