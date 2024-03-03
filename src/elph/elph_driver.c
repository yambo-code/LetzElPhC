/*
THe starting point for the entire code
*/
#include "elph.h"

void elph_driver(const char* ELPH_input_file,
                 enum ELPH_dft_code dft_code,
                 MPI_Comm comm_world)
{

    struct usr_input* input_data;
    // read the input file
    read_input_file(ELPH_input_file, &input_data, comm_world);
    // Note input parameters are broadcasted internally
    // All the parameters in input_data must be available for all cpus in
    // comm_world

    char* kernel = input_data->kernel;

    struct ELPH_MPI_Comms* mpi_comms = malloc(sizeof(struct ELPH_MPI_Comms));
    CHECK_ALLOC(mpi_comms);

    create_parallel_comms(input_data->nqpool, input_data->nkpool,
                          comm_world, mpi_comms);

    struct Lattice* lattice = malloc(sizeof(struct Lattice));
    CHECK_ALLOC(lattice);

    struct Pseudo* pseudo = malloc(sizeof(struct Pseudo));
    CHECK_ALLOC(pseudo);

    struct Phonon* phonon = malloc(sizeof(struct Phonon));
    CHECK_ALLOC(phonon);

    struct WFC* wfcs;

    // read the SAVE data and phonon related data.
    read_and_alloc_save_data(input_data->save_dir, mpi_comms,
                             input_data->start_bnd, input_data->end_bnd, &wfcs,
                             input_data->ph_save_dir, lattice, pseudo, phonon,
                             dft_code);

    //======= Now start the real computation =========
    // a) COmpute the D_mats and store them in the netcdf file
    // ============= Dmats =====================
    bool dmat_file_found = false; /// FIX ME
    if (!dmat_file_found)
    {
        compute_and_write_dmats("ndb.Dmats", wfcs, lattice, phonon->nph_sym,
                                phonon->ph_syms, mpi_comms);
    }
    // b) Compute elph
    // ============= ELPH iBZ computation =============
    ND_int nmodes = lattice->nmodes;
    ND_int nfft_loc = lattice->fft_dims[0] * lattice->fft_dims[1] * lattice->nfftz_loc;

    ELPH_cmplx* eigVec = malloc(sizeof(ELPH_cmplx) * nmodes * nmodes);
    //(nmodes,nmodes) // buffer to store eigen vectors
    CHECK_ALLOC(eigVec);

    ELPH_cmplx* dVscf = malloc(sizeof(ELPH_cmplx) * nmodes * lattice->nmag * nfft_loc);
    // (nmodes,nmag,Nx,Ny,Nz) // bufffer to store dVscf
    CHECK_ALLOC(dVscf);

    ELPH_float* omega_ph = malloc(sizeof(ELPH_float) * nmodes);
    // buffer for storing phonon freq
    CHECK_ALLOC(omega_ph);

    int ncid_elph, ncid_dmat, nc_err;
    int varid_eig, varid_elph, varid_omega, varid_dmat;
    // Define netcdf variables
    if (mpi_comms->commK_rank == 0)
    {
        // open Dmat file
        if ((nc_err = nc_open_par("ndb.Dmats", NC_NOWRITE, mpi_comms->commR,
                                  MPI_INFO_NULL, &ncid_dmat)))
        {
            ERR(nc_err);
        }

        // get dmat var id for dmats
        if ((nc_err = nc_inq_varid(ncid_dmat, "Dmats", &varid_dmat)))
        {
            ERR(nc_err);
        }

        // create elph file
        if ((nc_err = nc_create_par("ndb.elph", NC_NETCDF4, mpi_comms->commR,
                                    MPI_INFO_NULL, &ncid_elph)))
        {
            fprintf(stderr, "Error creating ndb.elph file.");
            ERR(nc_err);
        }

        // set no fill mode (to avoid writting twice)
        if ((nc_err = ncsetfill(ncid_elph, NC_NOFILL)))
        {
            fprintf(stderr, "Error setting nc_fill to ndb.elph file.");
            ERR(nc_err);
        }

        def_ncVar(
            ncid_elph, &varid_eig, 5, ELPH_NC4_IO_FLOAT,
            (ND_int[]) { phonon->nq_BZ, nmodes, nmodes / 3, 3, 2 },
            "POLARIZATION_VECTORS", (char*[]) { "nq", "nmodes", "atom", "pol", "re_im" },
            NULL);

        def_ncVar(ncid_elph, &varid_omega, 2, ELPH_NC4_IO_FLOAT,
                  (ND_int[]) { phonon->nq_BZ, nmodes }, "FREQ",
                  (char*[]) { "nq", "nmodes" }, NULL);

        def_ncVar(
            ncid_elph, &varid_elph, 7, ELPH_NC4_IO_FLOAT,
            (ND_int[]) { phonon->nq_BZ, lattice->nkpts_BZ, nmodes,
                         lattice->nspin, lattice->nbnds, lattice->nbnds, 2 },
            "elph_mat",
            (char*[]) { "nq", "nk", "nmodes", "nspin", "initial_band", "final_band_PH_abs", "re_im" },
            (size_t[]) { 1, 1, nmodes, lattice->nspin, lattice->nbnds,
                         lattice->nbnds, 2 });
    }

    ELPH_cmplx* eig_Sq = NULL;

    if (mpi_comms->commQ_rank == 0)
    {
        eig_Sq = calloc(nmodes * nmodes, sizeof(ELPH_cmplx));
        CHECK_ALLOC(eig_Sq);
    }

    for (ND_int iqpt = 0; iqpt < phonon->nq_iBZ_loc; ++iqpt)
    {
        ND_int iqpt_iBZg = iqpt + phonon->nq_shift;
        // read dynamical matrix and dvscf for the iBZ qpt
        if (dft_code == DFT_CODE_QE)
        {
            get_dvscf_dyn_qe(input_data->ph_save_dir, lattice, iqpt_iBZg,
                             eigVec, dVscf, omega_ph, mpi_comms);
            // qe dvscf only contains dV_Ha + dV_xc, we need to add the local
            // part of pseudo
        }
        else
        {
            error_msg("Currently only quantum espresso supported");
        }
        ELPH_cmplx* Vlocr = malloc(sizeof(ELPH_cmplx) * nmodes * nfft_loc);
        // buffer to store local part of the pseudo potential
        CHECK_ALLOC(Vlocr);
        //
        // compute the local part of the bare
        dVlocq(phonon->qpts_iBZ + iqpt_iBZg * 3, lattice, pseudo, eigVec, Vlocr,
               mpi_comms->commK);

        if (dft_code == DFT_CODE_QE)
        {
            add_dvscf_qe(dVscf, Vlocr, lattice); // add bare local to induce part
        }
        else
        {
            error_msg("Currently only quantum espresso supported");
        }

        free(Vlocr);

        ND_int qpos = 0; // positon of this iBZ qpoint in full q point list
        for (ND_int i = 0; i < iqpt_iBZg; ++i)
        {
            qpos += phonon->nqstar[i];
        }
        // write eigen vectors and frequencies
        if (mpi_comms->commQ_rank == 0)
        {
            size_t startp[5] = { qpos, 0, 0, 0, 0 };
            size_t countp[5] = { 1, nmodes, nmodes / 3, 3, 2 };
            if ((nc_err = nc_put_vara(ncid_elph, varid_eig, startp, countp, eigVec)))
            {
                ERR(nc_err);
            }

            if ((nc_err = nc_put_vara(ncid_elph, varid_omega, startp, countp, omega_ph)))
            {
                ERR(nc_err);
            }
            // write down the rotate eigen vectors;
            for (ND_int istar = 1; istar < phonon->nqstar[iqpt_iBZg]; ++istar)
            {
                ND_int qpos_star = qpos + istar;

                struct symmetry* sym_rot = phonon->ph_syms + phonon->qmap[2 * qpos_star + 1];

                rotate_eig_vecs(sym_rot, lattice, phonon->qpts_iBZ + iqpt_iBZg * 3,
                                eigVec, eig_Sq);
                //
                startp[0] = qpos_star;
                if ((nc_err = nc_put_vara(ncid_elph, varid_eig, startp, countp, eig_Sq)))
                {
                    ERR(nc_err);
                }

                if ((nc_err = nc_put_vara(ncid_elph, varid_omega, startp, countp, omega_ph)))
                {
                    ERR(nc_err);
                }
                //
            }
        }
        // Now compute and write the electron-phonon matrix elements
        compute_and_write_elphq(wfcs, lattice, pseudo, phonon, iqpt_iBZg,
                                eigVec, dVscf, ncid_elph, varid_elph,
                                ncid_dmat, varid_dmat, true,
                                input_data->kminusq, mpi_comms);
    }

    free(eig_Sq);

    if (mpi_comms->commK_rank == 0)
    {
        // close files
        if ((nc_err = nc_close(ncid_elph)))
        {
            ERR(nc_err);
        }

        if ((nc_err = nc_close(ncid_dmat)))
        {
            ERR(nc_err);
        }
    }

    // ELPH_cmplx Ry2Ha = pow(2,-1.5);
    free(omega_ph);
    free(eigVec);
    free(dVscf);

    // cleanup
    free_usr_input(input_data);
    free_save_data(wfcs, lattice, pseudo, phonon);
    free(lattice);
    free(pseudo);
    free(phonon);
    free_parallel_comms(mpi_comms);
    free(mpi_comms);
    fftw_fun(cleanup)();
}
