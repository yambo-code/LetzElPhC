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

    // read the input file
    read_input_file(argv[1], &input_data, MPI_COMM_WORLD);
    // Note input parameters are broadcasted internally

    struct usr_input * input_data;
    
    ND_int NK = input_data->nkpool;
    ND_int NQ  = input_data->nqpool;
    char * SAVEDIR  = input_data->save_dir ;
    ND_int FIRST_BAND   = input_data->start_bnd ;
    ND_int LAST_BAND  = input_data->end_bnd ;
    char * PH_SAVE = input_data->ph_save_dir;
    char * kernel = input_data->kernel;
    
    struct ELPH_MPI_Comms * mpi_comms = malloc(sizeof(struct ELPH_MPI_Comms));

    create_parallel_comms(NQ, NK, MPI_COMM_WORLD, mpi_comms);

    /*======= Section from here is q.e specific ===========*/
    // Now get the basic parameters from q.e or any dft code.
    // For now only qe supported
    ND_int nqpt_iBZ, nqpt_fullBZ;
    
    ELPH_float * qpts_iBZ;

    if (mpi_comms->commW_rank ==0)
    {
        read_qpts_qe(const char * dyn0_file, &nqpt_iBZ, &nqpt_fullBZ, &qpts_iBZ);
    }
    // Bcast data
    mpi_error = MPI_Bcast(&nqpt_iBZ,    1, ELPH_MPI_ND_INT, 0, mpi_comms->commW);
    mpi_error = MPI_Bcast(&nqpt_fullBZ, 1, ELPH_MPI_ND_INT, 0, mpi_comms->commW);
    
    // divide the qpoints in iBZ over q pools
    ND_int qshift ;
    ND_int nq_iBZ_this_pool = distribute_to_grps(nqpt_iBZ, NQ, \
                            mpi_comms->commW_rank/mpi_comms->commQ_size, &qshift);

    if (nq_iBZ_this_pool < 1) 
        error_msg("There are no qpoints in some qpools, Make sure nqpool < # of qpoints in iBZ.");
    
    // Bcast qpts_iBZ
    if (mpi_comms->commW_rank !=0)
    {
        ELPH_float * qpts_iBZ = malloc(sizeof(ELPH_float)*3*nqpt_iBZ);
    }
    mpi_error = MPI_Bcast(qpts_iBZ, 3*nqpt_iBZ, ELPH_MPI_float, 0, mpi_comms->commW);
    // qpts must be divided to get then cart units

    // now get the basic info from  
    ND_int FFT_dims[3];
    ELPH_float alat[3];

    // PSEUDO_DIR and pseudo_pots are required only on main root node
    char * PSEUDO_DIR   = NULL;
    char ** pseudo_pots = NULL;

    /* Phonon symmetries not to be confused with the ones used in yambo for wfcs*/
    bool ph_tim_rev;
    ELPH_float * ph_sym_tau;
    ELPH_float * ph_sym_mats;
    ND_int nph_sym;
    
    struct Lattice * lattice = malloc(sizeof(struct Lattice));
    struct Pseudo  * pseudo  = malloc(sizeof(struct Pseudo));
    struct WFC * wfcs;
    
    if (mpi_comms->commW_rank ==0)
    {
        parse_qexml(const char * xml_file, alat, &lattice->dimension, \
        &lattice->is_soc_present, &lattice->nmag, FFT_dims, &nph_sym, \
        &ph_sym_mats, &ph_sym_tau, &ph_tim_rev, &PSEUDO_DIR, &pseudo_pots);
    }
    
    // Bcast all the variables
    mpi_error = MPI_Bcast(FFT_dims, 3, ELPH_MPI_ND_INT, 0, mpi_comms->commW); 
    mpi_error = MPI_Bcast(alat, 3, ELPH_MPI_float, 0, mpi_comms->commW );
    
    mpi_error = MPI_Bcast(&lattice->nmag, 1, ELPH_MPI_ND_INT, 0, mpi_comms->commW); 
    mpi_error = MPI_Bcast(&nph_sym, 1, ELPH_MPI_ND_INT, 0, mpi_comms->commW); 

    mpi_error = MPI_Bcast(&lattice->is_soc_present, 1, MPI_C_BOOL, 0, mpi_comms->commW); 
    mpi_error = MPI_Bcast(&ph_tim_rev,     1, MPI_C_BOOL, 0, mpi_comms->commW); 
    
    // allocate memory for phonon symmetric matrices on rest of the cpus
    if (mpi_comms->commW_rank !=0)
    {   
        // we create 2*nph_sym sets of symmetries (factor 2 to store the time rev case)
        // Note: the second half([nph_sym:]) are only used when time reversal is present
        ph_sym_mats = malloc(sizeof(ELPH_float)*3*3*2*nph_sym);
        ph_sym_tau  = malloc(sizeof(ELPH_float)*3*2*nph_sym); 
    }

    if (mpi_comms->commW_rank ==0 ) printf("Reading SAVE data \n");
    
    read_and_alloc_save_data(SAVEDIR, mpi_comms, FIRST_BAND, LAST_BAND, \
                            &wfcs, PSEUDO_DIR, pseudo_pots, lattice, pseudo,FFT_dims);
    

    // free pseudo pots, no longer need
    if (mpi_comms->commW_rank ==0)
    {
        free(PSEUDO_DIR);
        for (ND_int ipot = 0; ipot < pseudo->ntype; ++ipot) free(pseudo_pots[ipot]);
        free(pseudo_pots);
        PSEUDO_DIR = NULL;
        pseudo_pots = NULL;
    }

    // convert qpts to crystal coordinates
    //qpts_iBZ, 3*nqpt_iBZ
    for (ND_int iqpt=0; iqpt<nqpt_iBZ; ++iqpt)
    {   
        ELPH_float * qpt_tmp = qpts_iBZ + 3*iqpt;
        ELPH_float qcart_tmp[3];
        for (int ix =0; ix < 3; ++ix) qcart_tmp[ix] = qpt_tmp[ix]/alat[ix] ;
        // convert to crystal coordinates
        MatVec3f(lattice->alat_vec->data, qcart_tmp, true, qpt_tmp); 
    }
    


    // now convert phonon symmetric matrices and phonon frac. trans. vec to cart units.
    if (mpi_comms->commW_rank ==0)
    {   
        ELPH_float blat[9];
        reciprocal_vecs(lattice->alat_vec->data,blat);
        for (int ix = 0; ix < 9; ++ix) blat[ix] /= (2.0f*ELPH_PI);

        for (ND_int isym=0; isym < nph_sym; ++isym)
        {   
            
            // Note we also fill the second half but are only used when tim_rev is present
            ELPH_float * sym_tmp      = ph_sym_mats + isym*9;
            ELPH_float * sym_tmp_trev = ph_sym_mats + (isym + nph_sym)*9;

            ELPH_float * vec_tmp      = ph_sym_tau +  isym*3;
            ELPH_float * vec_tmp_trev = ph_sym_tau +  (isym + nph_sym)*3;
            
            // compute lat_vec@S^T@b^T
            // first  S^T and b^T
            Gemm3x3f(sym_tmp, 'T', blat,'T', sym_tmp_trev); //sym_tmp_trev is used as tmp buffer
            Gemm3x3f(lattice->alat_vec->data,'N', sym_tmp_trev,'N', sym_tmp);
            
            for (int ix = 0; ix < 9; ++ix) sym_tmp_trev[ix] = -sym_tmp[ix];
            MatVec3f(lattice->alat_vec->data, vec_tmp, false, vec_tmp_trev);

            // we also negate the frac .tras. vec (just a convention used in this code)
            for (int ix = 0; ix < 3; ++ix)
            {
                vec_tmp[ix]      =  vec_tmp_trev[ix];
                vec_tmp_trev[ix] = -vec_tmp_trev[ix];
            }

        }
    }
    
    // create a array of bool which tells whether the symmetry is time_rev
    bool * ph_tim_rev_array = malloc(sizeof(bool)*2*nph_sym);
    for (ND_int isym=0; isym < nph_sym; ++isym)
    {   
        ph_tim_rev_array[isym]          = false; // 1st half are normal
        ph_tim_rev_array[isym+nph_sym]  = true; // second hald are time rev
    }

    /* In case of time reversal symmetry : 
     we double the number of symmetries and also use the second half of the symmetries */
    if (ph_tim_rev) nph_sym *= 2; 
    // bcast 
    mpi_error = MPI_Bcast(ph_sym_mats, 3*3*nph_sym, ELPH_MPI_float, 0, mpi_comms->commW );
    mpi_error = MPI_Bcast(ph_sym_tau,  3*nph_sym,   ELPH_MPI_float, 0, mpi_comms->commW );
    
    //======= Now we got all we need. start the real computation =========
    // a) COmpute the D_mats and store them in the netcdf file
    
    ND_int nk_totalBZ = lattice->kmap->dims[0];
    if (nk_totalBZ/mpi_comms->nkpools < 1)
        error_msg("There are no kpoints in some cpus, Make sure nkpool < # of kpoints in full BZ.");
    
    // ============= Dmats =====================
    bool dmat_file_found = false; /// FIX ME
    if (!dmat_file_found)
    {   
        compute_and_write_dmats("ndb.Dmats", wfcs, lattice, nph_sym, ph_sym_mats, \
                                    ph_sym_tau, ph_tim_rev_array, mpi_comms);
    }       

    // b) Compute elph in iBZ
    // ============= ELPH iBZ computation =============
    
    for (ND_int iqpt = 0; iqpt < nq_iBZ_this_pool; ++iqpt)
    {   
        ND_int iqpt_iBZg = iqpt + qshift;

        compute_elphq(wfcs, &lattice, &pseudo, qpt,  &eigVec, &dVscf, elph_kq, mpi_comms);
    }
    



    ELPH_cmplx Ry2Ha = pow(2,-1.5);
    
    
    // cleanup
    free(ph_tim_rev_array);
    free(ph_sym_tau);
    free(ph_sym_mats);
    free(qpts_iBZ);
    free_usr_input(input_data);
    free_save_data(wfcs, &lattice, &pseudo);
    free(lattice); free(pseudo);
    free_parallel_comms(mpi_comms);
    free(mpi_comms);
    fftw_fun(cleanup)();
    MPI_Finalize();

    return 0;

}
