/*
 Yambo calling point
*/
#include <complex.h>
#include <fftw3.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/ELPH_timers.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/init_dtypes.h"
#include "common/parallel.h"
#include "common/print_info.h"
#include "dvloc/dvloc.h"
#include "dvG_utils.h"
#include "elph.h"
#include "elphC.h"
#include "fft/fft.h"
#include "io/io.h"
#include "io/qe/qe_io.h"
#include "parser/parser.h"
#include "symmetries/symmetries.h"
#include "common/string_func.h"

/*
 * Callback-enabled driver: identical flow to elph_driver but writes gkkp
 * data via fill_fn(iq_BZ, ik_BZ, data, ...) instead of ndb.elph.
 * ndb.Dmats is still created and used internally for BZ symmetry expansion.
 * Eigenvectors / frequencies are NOT written (yambo reads them elsewhere).
 */
void elph_driver_cb(const char* ELPH_input_file, enum ELPH_dft_code dft_code,
                    MPI_Comm comm_world, elph_fill_fn fill_fn)
{
    struct elph_usr_input* input_data;
    init_ELPH_clocks();

    read_elph_input_file(ELPH_input_file, &input_data, comm_world);

    struct kernel_info* kernel = malloc(sizeof(struct kernel_info));
    init_kernel(kernel);
    set_kernel(input_data->kernel_str, kernel);

    struct ELPH_MPI_Comms* mpi_comms = malloc(sizeof(struct ELPH_MPI_Comms));
    CHECK_ALLOC(mpi_comms);

    create_parallel_comms(input_data->nqpool, input_data->nkpool, comm_world,
                          mpi_comms);

    print_ELPH_logo(mpi_comms->commW_rank, elph_get_log_file());
    print_info_msg(mpi_comms->commW_rank,
                   "********** Program started **********");
    print_input_info(input_data->save_dir, input_data->ph_save_dir,
                     input_data->kernel_str, input_data->kminusq, dft_code,
                     mpi_comms);

    struct Lattice* lattice = malloc(sizeof(struct Lattice));
    CHECK_ALLOC(lattice);
    init_lattice_type(lattice);

    struct Pseudo* pseudo = malloc(sizeof(struct Pseudo));
    CHECK_ALLOC(pseudo);
    init_Pseudo_type(pseudo);

    struct Phonon* phonon = malloc(sizeof(struct Phonon));
    CHECK_ALLOC(phonon);
    init_phonon_type(phonon);

    struct WFC* wfcs;

    if (dft_code == DFT_CODE_QE)
    {
        get_data_from_qe(lattice, phonon, pseudo, input_data->ph_save_dir, NULL,
                         mpi_comms);
    }
    else
    {
        error_msg("Only QE supported");
    }
    read_and_alloc_save_data(input_data->save_dir, mpi_comms,
                             input_data->start_bnd, input_data->end_bnd, &wfcs,
                             input_data->ph_save_dir, lattice, pseudo, phonon);
    char DM_name[100];
    strlcpy_custom(DM_name, input_data->ph_save_dir,100);
    strlcat(DM_name,"/ndb.Dmats", 100);

    print_lattice_info(mpi_comms, lattice);
    print_phonon_info(mpi_comms, phonon);

    print_info_msg(mpi_comms->commW_rank, "");
    print_info_msg(mpi_comms->commW_rank,
                   "=== Computing Dmats for phonon symmetries ===");
    compute_and_write_dmats(DM_name, wfcs, lattice, phonon->nph_sym,
                            phonon->ph_syms, mpi_comms);

    ND_int nmodes = lattice->nmodes;
    ND_int nfft_loc =
        lattice->fft_dims[0] * lattice->fft_dims[1] * lattice->nfftz_loc;

    ELPH_cmplx* eigVec = malloc(sizeof(ELPH_cmplx) * nmodes * nmodes);
    CHECK_ALLOC(eigVec);

    ELPH_cmplx* dVscf =
        malloc(sizeof(ELPH_cmplx) * nmodes * lattice->nmag * nfft_loc);
    CHECK_ALLOC(dVscf);

    ELPH_float* omega_ph = malloc(sizeof(ELPH_float) * nmodes);
    CHECK_ALLOC(omega_ph);

    int ncid_dmat, nc_err;
    int varid_dmat;

    if (mpi_comms->commK_rank == 0)
    {
        if ((nc_err = nc_open_par(DM_name, NC_NOWRITE, mpi_comms->commR,
                                  MPI_INFO_NULL, &ncid_dmat)))
        {
            ERR(nc_err);
        }

        if ((nc_err = nc_inq_varid(ncid_dmat, "Dmats", &varid_dmat)))
        {
            ERR(nc_err);
        }

        size_t Dmat_counts[6] = {0, 0, 0, 0, 0, 0};
        if ((nc_err = nc_get_vara(ncid_dmat, varid_dmat, Dmat_counts,
                                  Dmat_counts, NULL)))
        {
            ERR(nc_err);
        }
    }

    print_info_msg(mpi_comms->commW_rank,
                   "=== Computing Electron-phonon matrix elements ===");
    print_info_msg(mpi_comms->commW_rank, "");

    for (ND_int iqpt = 0; iqpt < phonon->nq_iBZ_loc; ++iqpt)
    {
        print_info_msg(mpi_comms->commW_rank, "### q-point : %d/%d",
                       (int)(iqpt + 1), (int)phonon->nq_iBZ_loc);

        ND_int iqpt_iBZg = iqpt + phonon->nq_shift;

        if (dft_code == DFT_CODE_QE)
        {
            ELPH_cmplx* dvscf_read = NULL;
            if (kernel->screening == ELPH_DFPT_SCREENING)
            {
                dvscf_read = dVscf;
            }
            else
            {
                ND_int dvscf_num = nmodes * lattice->nmag * nfft_loc;
                for (ND_int ix = 0; ix < dvscf_num; ++ix)
                {
                    dVscf[ix] = 0.0;
                }
            }
            get_dvscf_dyn_qe(input_data->ph_save_dir, lattice, iqpt_iBZg,
                             eigVec, dvscf_read, omega_ph, mpi_comms);
        }
        else
        {
            error_msg("Currently only quantum espresso supported");
        }

        ELPH_cmplx* Vlocr = malloc(sizeof(ELPH_cmplx) * nmodes * nfft_loc);
        CHECK_ALLOC(Vlocr);
        dVlocq(phonon->qpts_iBZ + iqpt_iBZg * 3, lattice, pseudo, eigVec, Vlocr,
               mpi_comms->commK);

        if (dft_code == DFT_CODE_QE)
        {
            add_dvscf_qe(dVscf, Vlocr, lattice);
        }
        else
        {
            error_msg("Currently only quantum espresso supported");
        }
        free(Vlocr);

        compute_and_write_elphq(wfcs, lattice, pseudo, phonon, iqpt_iBZg,
                                eigVec, dVscf, 0, 0, ncid_dmat,
                                varid_dmat, kernel->non_loc,
                                input_data->kminusq, mpi_comms, fill_fn,
                                iqpt_iBZg);
    }

    if (mpi_comms->commK_rank == 0)
    {
        if ((nc_err = nc_close(ncid_dmat)))
        {
            ERR(nc_err);
        }
    }

    free(omega_ph);
    free(eigVec);
    free(dVscf);

    int World_rank_tmp = mpi_comms->commW_rank;

    free(kernel);
    free_elph_usr_input(input_data);
    free_save_data(wfcs, lattice, pseudo, phonon);
    free(lattice);
    free(pseudo);
    free(phonon);
    free_parallel_comms(mpi_comms);
    free(mpi_comms);
    fftw_fun(cleanup)();

    if (0 == World_rank_tmp)
    {
        print_ELPH_clock_summary();
    }
    cleanup_ELPH_clocks();
    print_info_msg(World_rank_tmp, "********** Program ended **********");
}


