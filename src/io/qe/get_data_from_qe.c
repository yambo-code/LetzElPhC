#include "qe_io.h"

void get_data_from_qe(struct Lattice * lattice, \
    struct Phonon * phonon, const char * ph_save_dir, \
    char *** pseudo_pots, const struct ELPH_MPI_Comms * Comm)
{
    /*
    This functions gets basic ground state data and phonon data from q.e 
    that are not saved by YAMBO.
    */
    int mpi_error;
    
    char * tmp_buffer = NULL;
    if (Comm->commW_rank ==0)
    {
        tmp_buffer = calloc(1024+strlen(ph_save_dir),1);
    }
    if (Comm->commW_rank ==0)
    {   
        strcpy(tmp_buffer,ph_save_dir);
        strcat(tmp_buffer, "/dyn0");
        read_qpts_qe(tmp_buffer, &phonon->nq_iBZ, &phonon->nq_BZ, &phonon->qpts_iBZ);
    }
    // Bcast data
    mpi_error = MPI_Bcast(&phonon->nq_iBZ,1, ELPH_MPI_ND_INT, 0, Comm->commW);
    mpi_error = MPI_Bcast(&phonon->nq_BZ, 1, ELPH_MPI_ND_INT, 0, Comm->commW);
    
    // divide the qpoints in iBZ over q pools
    phonon->nq_iBZ_loc = distribute_to_grps(phonon->nq_iBZ, Comm->nqpools, \
                    Comm->commW_rank/Comm->commQ_size, &phonon->nq_shift);

    if (phonon->nq_iBZ_loc < 1) 
        error_msg("There are no qpoints in some qpools, Make sure nqpool < # of qpoints in iBZ.");
    
    // Bcast phonon->qpts_iBZ
    if (Comm->commW_rank !=0)
    {
        phonon->qpts_iBZ = malloc(sizeof(ELPH_float)*3*phonon->nq_iBZ);
    }
    mpi_error = MPI_Bcast(phonon->qpts_iBZ, 3*phonon->nq_iBZ, ELPH_MPI_float, 0, Comm->commW);
    // qpts must be divided to get then cart units
    // now get the basic info from  
    ELPH_float alat[3];
    char * PSEUDO_DIR   = NULL;
    bool ph_tim_rev;
    
    ELPH_float lat_vec[9]; // a[:,i] is ith lattice vector
    if (Comm->commW_rank ==0)
    {   
        strcpy(tmp_buffer,ph_save_dir);
        strcat(tmp_buffer, "/data-file-schema.xml");
        parse_qexml(tmp_buffer, lat_vec, alat, &lattice->dimension, \
        &lattice->is_soc_present, &lattice->nmag, lattice->fft_dims, &phonon->nph_sym, \
        &phonon->ph_sym_mats, &phonon->ph_sym_tau, &ph_tim_rev, &PSEUDO_DIR, pseudo_pots);
        // free pseudo pots, no longer need
        free(PSEUDO_DIR);
        PSEUDO_DIR = NULL;
    }
    
    // Bcast all the variables
    mpi_error = MPI_Bcast(lat_vec, 9, ELPH_MPI_float, 0, Comm->commW );
    mpi_error = MPI_Bcast(lattice->fft_dims, 3, ELPH_MPI_ND_INT, 0, Comm->commW );
    mpi_error = MPI_Bcast(alat, 3, ELPH_MPI_float, 0, Comm->commW );
    
    mpi_error = MPI_Bcast(&lattice->nmag, 1, ELPH_MPI_ND_INT, 0, Comm->commW); 
    mpi_error = MPI_Bcast(&phonon->nph_sym, 1, ELPH_MPI_ND_INT, 0, Comm->commW); 

    mpi_error = MPI_Bcast(&lattice->is_soc_present, 1, MPI_C_BOOL, 0, Comm->commW); 
    mpi_error = MPI_Bcast(&ph_tim_rev,     1, MPI_C_BOOL, 0, Comm->commW); 

    ELPH_float blat[9];
    reciprocal_vecs(lat_vec,blat);
    for (int ix = 0; ix < 9; ++ix) blat[ix] /= (2.0f*ELPH_PI);

    // allocate memory for phonon symmetric matrices on rest of the cpus
    if (Comm->commW_rank !=0)
    {   
        // we create 2*phonon->nph_sym sets of symmetries (factor 2 to store the time rev case)
        // Note: the second half([phonon->nph_sym:]) are only used when time reversal is present
        phonon->ph_sym_mats = malloc(sizeof(ELPH_float)*3*3*2*phonon->nph_sym);
        phonon->ph_sym_tau  = malloc(sizeof(ELPH_float)*3*2*phonon->nph_sym); 
    }

    // convert qpts to crystal coordinates
    //phonon->qpts_iBZ, 3*phonon->nq_iBZ
    for (ND_int iqpt=0; iqpt<phonon->nq_iBZ; ++iqpt)
    {   
        ELPH_float * qpt_tmp = phonon->qpts_iBZ + 3*iqpt;
        ELPH_float qcart_tmp[3];
        for (int ix =0; ix < 3; ++ix) qcart_tmp[ix] = qpt_tmp[ix]/alat[ix] ;
        // convert to crystal coordinates
        MatVec3f(lat_vec, qcart_tmp, true, qpt_tmp); 
    }
    
    // now convert phonon symmetric matrices and phonon frac. trans. vec to cart units.
    if (Comm->commW_rank ==0)
    {   
        for (ND_int isym=0; isym < phonon->nph_sym; ++isym)
        {   
            // Note we also fill the second half but are only used when tim_rev is present
            ELPH_float * sym_tmp      = phonon->ph_sym_mats + isym*9;
            ELPH_float * sym_tmp_trev = phonon->ph_sym_mats + (isym + phonon->nph_sym)*9;

            ELPH_float * vec_tmp      = phonon->ph_sym_tau +  isym*3;
            ELPH_float * vec_tmp_trev = phonon->ph_sym_tau +  (isym + phonon->nph_sym)*3;
            
            // compute lat_vec@S^T@b^T
            // first  S^T and b^T
            Gemm3x3f(sym_tmp, 'T', blat,'T', sym_tmp_trev); //sym_tmp_trev is used as tmp buffer
            Gemm3x3f(lat_vec,'N', sym_tmp_trev,'N', sym_tmp);
            
            for (int ix = 0; ix < 9; ++ix) sym_tmp_trev[ix] = -sym_tmp[ix];
            MatVec3f(lat_vec, vec_tmp, false, vec_tmp_trev);

            // we also negate the frac .tras. vec (just a convention used in this code)
            for (int ix = 0; ix < 3; ++ix)
            {
                vec_tmp[ix]      =  vec_tmp_trev[ix];
                vec_tmp_trev[ix] = -vec_tmp_trev[ix];
            }

        }
    }
    
    // create a array of bool which tells whether the symmetry is time_rev
    phonon->time_rev_array = malloc(sizeof(bool)*2*phonon->nph_sym);
    for (ND_int isym=0; isym < phonon->nph_sym; ++isym)
    {   
        phonon->time_rev_array[isym]                  = false; // 1st half are normal
        phonon->time_rev_array[isym+phonon->nph_sym]  = true; // second hald are time rev
    }

    /* In case of time reversal symmetry : 
     we double the number of symmetries and also use the second half of the symmetries */
    if (ph_tim_rev) phonon->nph_sym *= 2; 
    // bcast 
    mpi_error = MPI_Bcast(phonon->ph_sym_mats, 3*3*phonon->nph_sym, ELPH_MPI_float, 0, Comm->commW );
    mpi_error = MPI_Bcast(phonon->ph_sym_tau,  3*phonon->nph_sym,   ELPH_MPI_float, 0, Comm->commW );
    
    if (Comm->commW_rank ==0) free(tmp_buffer);

    phonon->qpts_BZ = malloc(phonon->nq_BZ*3*sizeof(ELPH_float));
    phonon->qmap    = malloc(phonon->nq_BZ*2*sizeof(int));
    phonon->nqstar  = malloc(phonon->nq_iBZ*sizeof(ND_int));
    
    // expand the qpoints into full BZ
    ND_int nqBZ_found = 0;
    for (ND_int iqpt = 0 ; iqpt < phonon->nq_iBZ; ++iqpt)
    {   
        ELPH_float qpt_iBZ_tmp[3];
        MatVec3f(blat,phonon->qpts_iBZ + 3*iqpt, false, qpt_iBZ_tmp);
        phonon->nqstar[iqpt] = 0;
        for (ND_int isym = 0; isym < phonon->nph_sym; ++isym)
        {
            ELPH_float qstar_cart[3], qstar_crys[3];
            MatVec3f(phonon->ph_sym_mats + isym*9, qpt_iBZ_tmp, false, qstar_cart);
            MatVec3f(lat_vec, qstar_cart, true, qstar_crys); 

            bool qstar_present = false;
            for (ND_int ifound = 0; ifound<nqBZ_found; ++ifound)
            {
                ELPH_float * qtmp = phonon->qpts_BZ + ifound*3;
                ELPH_float sum_norm=0;
                for (int ix = 0; ix < 3; ++ix)
                {
                    ELPH_float diff_tmp = qtmp[ix]-qstar_crys[ix];
                    diff_tmp = diff_tmp- rint(diff_tmp);
                    sum_norm += diff_tmp*diff_tmp;
                }
                sum_norm = sqrt(sum_norm);
                if (sum_norm < ELPH_EPS) 
                {
                    qstar_present = true;
                    break;
                }
            }
        
            if(!qstar_present)
            {   
                memcpy(phonon->qpts_BZ + nqBZ_found*3, qstar_crys,3*sizeof(ELPH_float));
                phonon->qmap[2*nqBZ_found] = iqpt;
                phonon->qmap[2*nqBZ_found+1] = isym;
                ++phonon->nqstar[iqpt];
                ++nqBZ_found;
            }
        }
    }

    // some sanity checks
    ND_int sum_nq_star = 0;
    for (ND_int i = 0; i<phonon->nq_iBZ; ++i) sum_nq_star += phonon->nqstar[i];
    if (sum_nq_star != nqBZ_found || nqBZ_found != phonon->nq_BZ)
        error_msg("Expansion of q points to full BZ failed.");

}




