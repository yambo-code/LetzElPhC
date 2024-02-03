/*
THe starting point for the entire code
*/
#include "elph.h"

int main(int argc, char* argv[])
{   
    
    #if defined(ELPH_OMP_PARALLEL_BUILD)
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    #else
    MPI_Init(&argc, &argv);
    #endif
    
    int mpi_error;

    struct usr_input * input_data;
    // read the input file
    read_input_file(argv[1], &input_data, MPI_COMM_WORLD);
    // Note input parameters are broadcasted internally
    // All the parameters in input_data must be available for all cpus in MPI_COMM_WORLD
    
    ND_int NK = input_data->nkpool;
    ND_int NQ  = input_data->nqpool;
    char * SAVEDIR  = input_data->save_dir ;
    ND_int FIRST_BAND   = input_data->start_bnd ;
    ND_int LAST_BAND  = input_data->end_bnd ;
    char * PH_SAVE = input_data->ph_save_dir;
    char * kernel = input_data->kernel;
    
    struct ELPH_MPI_Comms * mpi_comms = malloc(sizeof(struct ELPH_MPI_Comms));

    create_parallel_comms(NQ, NK, MPI_COMM_WORLD, mpi_comms);

    struct Lattice * lattice    = malloc(sizeof(struct Lattice));
    struct Pseudo  * pseudo     = malloc(sizeof(struct Pseudo));
    struct Phonon  * phonon     = malloc(sizeof(struct Phonon));
    struct WFC * wfcs;

    char ** pseudo_pots = NULL;
    // read phonon stuff from qe (for abinit a similar routine needs to be written)
    get_data_from_qe(lattice, phonon, PH_SAVE, &pseudo_pots, mpi_comms);

    // read the YAMBO SAVE data
    read_and_alloc_save_data(SAVEDIR, mpi_comms, FIRST_BAND, LAST_BAND, \
            &wfcs, PH_SAVE, pseudo_pots, lattice, pseudo, phonon, "qe");
    
    // free pseudo pots, no longer need
    if (mpi_comms->commW_rank ==0)
    {
        for (ND_int ipot = 0; ipot < pseudo->ntype; ++ipot) free(pseudo_pots[ipot]);
        free(pseudo_pots);
    }

    ND_int nk_totalBZ = lattice->kmap->dims[0];
    if (nk_totalBZ/mpi_comms->nkpools < 1)
        error_msg("There are no kpoints in some cpus, Make sure nkpool < # of kpoints in full BZ.");
    
    
    //======= Now we got all we need. start the real computation =========
    // a) COmpute the D_mats and store them in the netcdf file
    // ============= Dmats =====================
    bool dmat_file_found = false; /// FIX ME
    if (!dmat_file_found)
    {   
        compute_and_write_dmats("ndb.Dmats", wfcs, lattice, phonon->nph_sym, phonon->ph_sym_mats, \
                                    phonon->ph_sym_tau, phonon->time_rev_array, mpi_comms);
    }       
    // b) Compute elph 
    // ============= ELPH iBZ computation =============
    ND_array(Nd_cmplxS) eigVec[1], dVscf[1];

    ND_int nmodes = lattice->atomic_pos->dims[0]*3;
    
    ND_function(init, Nd_cmplxS)(eigVec, 2, nd_idx{nmodes,nmodes}); 
    //(nmodes,nmodes)
    ND_function(init, Nd_cmplxS)(dVscf, 5, nd_idx{nmodes,lattice->nmag, \
        lattice->fft_dims[0],lattice->fft_dims[1],lattice->nfftz_loc}); 
    // (nmodes,nmag,Nx,Ny,Nz)

    // allocate the memory
    ND_function(malloc, Nd_cmplxS)(eigVec);
    ND_function(malloc, Nd_cmplxS)(dVscf);
    ELPH_float * omega_ph = malloc(sizeof(ELPH_float)*nmodes);

    for (ND_int iqpt = 0; iqpt < phonon->nq_iBZ_loc; ++iqpt)
    {   
        ND_int iqpt_iBZg = iqpt + phonon->nq_shift;
        // read dynamica matrix and dvscf for the iBZ qpt
        get_dvscf_dyn_qe(PH_SAVE, lattice, iqpt_iBZg, eigVec->data, dVscf->data, omega_ph, mpi_comms);
        // Now compute the electron-phonon matrix elements
        //compute_elphq(wfcs, &lattice, &pseudo, qpt,  &eigVec, &dVscf, elph_kq, mpi_comms);
    }
    // ELPH_cmplx Ry2Ha = pow(2,-1.5);

    free(omega_ph);
    ND_function(destroy,Nd_cmplxS)(eigVec);
    ND_function(destroy,Nd_cmplxS)(dVscf);

    // cleanup
    free_usr_input(input_data);
    free_save_data(wfcs, lattice, pseudo);
    free_phonon_data(phonon);
    free(lattice); free(pseudo); free(phonon);
    free_parallel_comms(mpi_comms);
    free(mpi_comms);
    fftw_fun(cleanup)();
    MPI_Finalize();

    return 0;

}



