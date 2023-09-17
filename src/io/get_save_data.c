/* This function reads all the lattice, pseudo and wfcs data from SAVE DIR */

#include "io.h"


/*static functions */
static void quick_read(const int ncid, char* var_name, void * data_out);

static void alloc_and_set_Gvec(ND_array(Nd_floatS) * gvec, ND_int ik, ND_array(Nd_floatS) * totalGvecs, \
                    ND_array(Nd_floatS) * Gvecidxs, ELPH_float * lat_param, ND_int nG, ND_int nG_shift);

static void quick_read_sub(const int ncid, char* var_name, const size_t * startp, \
                            const size_t * countp, void * data_out);

static void get_wfc_from_save(ND_int spin_stride_len, ND_int ik, ND_int nkiBZ, \
            ND_int nspin, ND_int nspinor, ND_int start_band, ND_int nbnds, \
            ND_int nG, ND_int G_shift, const char * save_dir, char * work_array, \
            ELPH_cmplx * out_wfc, MPI_Comm comm);

/* Function body */
void read_and_alloc_save_data(char * SAVEdir, MPI_Comm commQ, MPI_Comm commK,  \
                ND_int start_band, ND_int end_band, struct WFC ** wfcs,char * pseudo_dir, \
                char ** pseudo_pots, struct Lattice * lattice, struct Pseudo * pseudo, \
                const ND_int * FFT_dims)
{
    /* This function allocates and reads data from SAVE dir.
    The following data is read : wfcs(in iBZ), lattice and pseudo
    start_band, end_band are give in fortran indices i.e 1st band starts 
    from 1 instead of 0;
    // pseudo_pots  : list of pseudopotential files (for now only upf2 is supported)
    These variables are not allocated here 
    ---
    Lattice :
    dimesnion : (read from input )
    ---
    pseudo : 
    */
    /*
    Expect wfcs, all are read by the single IO and broadcasted 

    FFT_dims // only rank 0 has to pass this. 
    */

    int my_rank, comm_size, mpi_error;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int cpus_qpool, qcolor, qrank;

    MPI_Comm_rank(commQ, &qrank);
    MPI_Comm_size(commQ, &cpus_qpool);
    qcolor = my_rank/cpus_qpool ;
    int nqpools = comm_size/cpus_qpool;

    int npw_cpus, krank;
    MPI_Comm_rank(commK, &krank);
    MPI_Comm_size(commK, &npw_cpus);
    int kcolor = qrank/npw_cpus;

    if (my_rank == 0) memcpy(lattice->fft_dims,FFT_dims,sizeof(ND_int)*3);
    // Bcast fft_dims
    mpi_error = MPI_Bcast(lattice->fft_dims, 3, MPI_INT, 0, MPI_COMM_WORLD );

    ND_int nffts = lattice->fft_dims[0]*lattice->fft_dims[1]*lattice->fft_dims[2];

    ND_int nffts_per_core = nffts/npw_cpus;
    ND_int nffts_rem = nffts%npw_cpus;
    if (nffts_per_core <1) error_msg("Some cpus do not contain plane waves ");
    
    
    ND_int nffts_in_this_cpu = nffts_per_core;
    if (krank < nffts_rem) ++nffts_in_this_cpu;

    lattice->nffts_loc = nffts_in_this_cpu;

    lattice->nfft_shift_loc = krank*nffts_per_core;
    if (krank < nffts_rem) lattice->nfft_shift_loc += krank;
    else                   lattice->nfft_shift_loc += nffts_rem;

    int dbid, ppid, tempid, retval; // file ids for ns.db1 , pp_pwscf*

    char * temp_str = malloc(sizeof(char) * (strlen(pseudo_dir) + strlen(SAVEdir) + 100));

    int nkBZ ; // total kpoints in BZ

    char * elements = malloc(sizeof(char)*3*104); // coded 104 elements 
    {   
        char * temp = "NA\0H \0He\0Li\0Be\0B \0C \0N \0O " \
        "\0F \0Ne\0Na\0Mg\0Al\0Si\0P \0S \0Cl\0Ar\0K \0Ca" \
        "\0Sc\0Ti\0V \0Cr\0Mn\0Fe\0Co\0Ni\0Cu\0Zn\0Ga\0Ge" \
        "\0As\0Se\0Br\0Kr\0Rb\0Sr\0Y \0Zr\0Nb\0Mo\0Tc\0Ru" \
        "\0Rh\0Pd\0Ag\0Cd\0In\0Sn\0Sb\0Te\0I \0Xe\0Cs\0Ba" \
        "\0La\0Ce\0Pr\0Nd\0Pm\0Sm\0Eu\0Gd\0Tb\0Dy\0Ho\0Er" \
        "\0Tm\0Yb\0Lu\0Hf\0Ta\0W \0Re\0Os\0Ir\0Pt\0Au\0Hg" \
        "\0Tl\0Pb\0Bi\0Po\0At\0Rn\0Fr\0Ra\0Ac\0Th\0Pa\0U " \
        "\0Np\0Pu\0Am\0Cm\0Bk\0Cf\0Es\0Fm\0Md\0No\0Lr\0" ;
        memcpy(elements,temp,sizeof(char)*3*104); // printf(elements+3*Z) will give symbol for Z
    }
    /*****/
    if (my_rank == 0)
    {
        sprintf(temp_str, "%s/ndb.kindx", SAVEdir) ; 
        NC_open_file(temp_str, 'r', &tempid);

        ELPH_float kindx_pars[7];
        quick_read(tempid, "PARS", kindx_pars);
        nkBZ = (int)rint(kindx_pars[0]); // FIX ME !! or kindx_pars[5] ?

        NC_close_file(tempid);
    } 
    /* broad cast (int nkBZ) */
    mpi_error = MPI_Bcast(&nkBZ, 1, MPI_INT, 0, MPI_COMM_WORLD );
    /*******/

    //printf("Debug-%d \n",1);
    ELPH_float dimensions[18];
    if (my_rank == 0)
    {
        sprintf(temp_str, "%s/ns.db1", SAVEdir) ; 
        NC_open_file(temp_str, 'r', &dbid);
        quick_read(dbid, "DIMENSIONS", dimensions);
    }
    /* bcast ELPH_float dimensions[18] */
    mpi_error = MPI_Bcast(dimensions, 18, ELPH_MPI_float, 0, MPI_COMM_WORLD );

    lattice->nspinor = (int)rint(dimensions[11]);
    lattice->nspin   = (int)rint(dimensions[12]);
    lattice->timerev = (int)rint(dimensions[9]);
    lattice->total_bands = (int) rint(dimensions[5]);

    if (start_band < 1 || end_band<1)
    {
        if (my_rank == 0) printf("Warning : invalid bands used in calculation. computing matrix elements for all bands, \
                Bands index belong to [1,nbnds] \n");
        start_band = 1;
        end_band = lattice->total_bands;
    }
    if (start_band>lattice->total_bands  || end_band>lattice->total_bands  || start_band>=end_band)
    {
        if (my_rank == 0) printf("Warning : invalid bands used in calculation. computing matrix elements for all bands \n");
        start_band = 1;
        end_band = lattice->total_bands;
    }

    lattice->start_band = start_band;
    lattice->end_band = end_band;
    lattice->nbnds = end_band-start_band+1;

    int nibz = (int) rint(dimensions[6]);

    /* Free me in free function */
    ND_array(Nd_floatS) * lattice_data = malloc( 9*sizeof(ND_array(Nd_floatS)) );
    ND_array(i) * kmap = malloc( 1* sizeof(ND_array(i)) );

    lattice->alat_vec           = lattice_data;
    lattice->atomic_pos         = lattice_data+1;
    lattice->kpt_iredBZ         = lattice_data+2;
    lattice->kpt_fullBZ         = lattice_data+3;
    lattice->kpt_fullBZ_crys    = lattice_data+4;
    lattice->sym_mat            = lattice_data+5;
    lattice->frac_trans         = lattice_data+6;
    lattice->kmap               = kmap;


    ND_array(Nd_floatS) * sym_temp  = lattice_data+7;
    ND_array(Nd_floatS) * kibz_temp = lattice_data+8;


    if (my_rank == 0)
    {
        ND_function(init,Nd_floatS) (lattice->alat_vec,          0, NULL); // alat_vec 'r'
    }

    ND_function(init,Nd_floatS) (lattice->kpt_iredBZ,        2, nd_idx{nibz,3} ); // ibZ kpts 'c'
    ND_function(init,Nd_floatS) (lattice->kpt_fullBZ,        2, nd_idx{nkBZ,3} ); // full bZ kpts cart 'c'
    ND_function(init,Nd_floatS) (lattice->kpt_fullBZ_crys,   2, nd_idx{nkBZ,3} ); // full bZ kpts crystal  'c'
    ND_function(calloc,Nd_floatS) (lattice->kpt_iredBZ);
    ND_function(calloc,Nd_floatS) (lattice->kpt_fullBZ);
    ND_function(calloc,Nd_floatS) (lattice->kpt_fullBZ_crys);

    ND_function(set_all,Nd_floatS) (lattice->kpt_fullBZ,0.0); // zero initialize buffer

    nd_init_i(lattice->kmap, 2, nd_idx{nkBZ,2}); // 'c'
    nd_calloc_i(lattice->kmap);


    ELPH_float lat_param[3];
    if (my_rank == 0)
    {
        quick_read(dbid, "LATTICE_PARAMETER", lat_param);
    }
    /*Bcast ELPH_float lat_param[3] */
    mpi_error = MPI_Bcast(lat_param, 3, ELPH_MPI_float, 0, MPI_COMM_WORLD );


    if (my_rank == 0)
    {
        ND_function(init,Nd_floatS) (sym_temp,  0, NULL); // symtemp 'r'
        ND_function(init,Nd_floatS) (kibz_temp, 0, NULL); // kpt_iredBZ.T 'r'
    }
    /* read */
    if (my_rank == 0)
    {
        ND_function(readVar, Nd_floatS) (dbid, "LATTICE_VECTORS", lattice->alat_vec  );
        ND_function(readVar, Nd_floatS) (dbid, "SYMMETRY", sym_temp);
        ND_function(readVar, Nd_floatS) (dbid, "K-POINTS", kibz_temp);
    }
    /* Bcast */
    Bcast_ND_arrayFloat(lattice->alat_vec, true, 0, MPI_COMM_WORLD);
    Bcast_ND_arrayFloat(sym_temp, true, 0, MPI_COMM_WORLD);
    Bcast_ND_arrayFloat(kibz_temp, true, 0, MPI_COMM_WORLD);

    /* allocate symmetric matrices*/
    ND_function(init,Nd_floatS) (lattice->sym_mat, *(sym_temp->rank), sym_temp->dims ); // ibZ kpts 'c'
    ND_function(calloc,Nd_floatS) (lattice->sym_mat);

    /* tau . Yambo do not accept non-symmorphic, so tau is 0 */
    ND_function(init,Nd_floatS) (lattice->frac_trans,   2, nd_idx{sym_temp->dims[0],3} ); // (nsym,3)
    ND_function(malloc,Nd_floatS) (lattice->frac_trans);
    ND_function(set_all,Nd_floatS) (lattice->frac_trans,0.0);


    /* Get kpoints to cart coordinates */
    for (ND_int i = 0 ; i<nibz ; ++i )
    {
        kibz_temp->data[i + 0*nibz] = kibz_temp->data[i + 0*nibz]/lat_param[0];
        kibz_temp->data[i + 1*nibz] = kibz_temp->data[i + 1*nibz]/lat_param[1];
        kibz_temp->data[i + 2*nibz] = kibz_temp->data[i + 2*nibz]/lat_param[2];
        /* Transpose to (nibz,3)*/
        lattice->kpt_iredBZ->data[0 + i*3] = kibz_temp->data[i + 0*nibz];
        lattice->kpt_iredBZ->data[1 + i*3] = kibz_temp->data[i + 1*nibz];
        lattice->kpt_iredBZ->data[2 + i*3] = kibz_temp->data[i + 2*nibz];
    }

    /* Expand kpoints to full BZ */
    bz_expand(kibz_temp, sym_temp, lattice->alat_vec, lattice->kpt_fullBZ_crys , lattice->kmap);
    
    /* setup time reversal array */
    lattice->time_rev_array = malloc(sizeof(bool)*sym_temp->dims[0]);

    /* Transpose symmetries as yambo stores in reverse order and compute time_rev_array */
    for (ND_int i =0 ; i<sym_temp->dims[0]; ++i)
    {
        transpose3x3f(sym_temp->data + 9*i, lattice->sym_mat->data + 9*i);
        
        /* fill the time-rev array*/
        if (i >= (sym_temp->dims[0]/(lattice->timerev+1)) ) lattice->time_rev_array[i] =true;
        else lattice->time_rev_array[i] =false;

    }

    {   /* create a scope to free the stack after usage */
        ELPH_float blat[9] ;
        reciprocal_vecs(lattice->alat_vec->data,blat); // b[:,i]  are blat. blat comes with 2*pi factor
    
        ND_function(matmulX, Nd_floatS) ('N', 'T', lattice->kpt_fullBZ_crys->data, blat,\
                lattice->kpt_fullBZ->data, 1.0f/(2.0f*ELPH_PI) , 0.0, 3, 3, 3, nkBZ, 3, 3);
    
        ELPH_float ntype;
        if (my_rank == 0)
        {
            quick_read(dbid, "number_of_atom_species", &ntype);
        }
        /* Bcast ELPH_float ntype */
        mpi_error = MPI_Bcast(&ntype, 1, ELPH_MPI_float, 0, MPI_COMM_WORLD );

        pseudo->ntype = (ND_int)rint(ntype) ;
    }

    ND_function(destroy, Nd_floatS)(sym_temp);
    ND_function(destroy, Nd_floatS)(kibz_temp);

    /* Read atomic positions */
    ND_array(Nd_floatS) atomic_map, atom_pos_temp ;

    
    if (my_rank == 0)
    {
        ND_function(init,Nd_floatS) (&atomic_map,  0, NULL); // 'r'
        ND_function(init,Nd_floatS) (&atom_pos_temp, 0, NULL); //  'r'
    }

    ELPH_float * natom_per_type = malloc(sizeof(ELPH_float) * 2 * pseudo->ntype );
    ELPH_float * atomic_numbers = natom_per_type + pseudo->ntype ;
    
    
    if (my_rank == 0)
    {
        quick_read(dbid, "N_ATOMS", natom_per_type);
        quick_read(dbid, "atomic_numbers", atomic_numbers);
        ND_function(readVar, Nd_floatS) (dbid, "ATOM_MAP",  &atomic_map);
        ND_function(readVar, Nd_floatS) (dbid, "ATOM_POS",  &atom_pos_temp);
    }
    /* Bcast atomic_map and atom_pos_temp */
    Bcast_ND_arrayFloat(&atomic_map, true, 0, MPI_COMM_WORLD);
    Bcast_ND_arrayFloat(&atom_pos_temp, true, 0, MPI_COMM_WORLD);
    /* Bcast ELPH_float natom_per_type,  atomic_numbers */
    mpi_error = MPI_Bcast(natom_per_type, pseudo->ntype, ELPH_MPI_float, 0, MPI_COMM_WORLD );
    mpi_error = MPI_Bcast(atomic_numbers, pseudo->ntype, ELPH_MPI_float, 0, MPI_COMM_WORLD );


    ND_int total_atoms =0 ;
    for (ND_int ia = 0; ia<pseudo->ntype; ++ia) total_atoms += (ND_int)rint(natom_per_type[ia]) ;

    lattice->atom_type = malloc(sizeof(int)*total_atoms);
    char * atom_symbols = malloc(sizeof(char)*3*pseudo->ntype);
    //
    ND_function(init,Nd_floatS) (lattice->atomic_pos,   2, nd_idx{total_atoms,3} ); // full bZ kpts crystal  'c'
    ND_function(calloc,Nd_floatS) (lattice->atomic_pos);
    for (ND_int it = 0 ; it<pseudo->ntype ; ++it)
    {   
        ND_int ia_no = rint(atomic_numbers[it]); 
        memcpy(atom_symbols + 3*it, elements + 3*ia_no, 3*sizeof(char));

        ND_int nspec = rint(natom_per_type[it]);
        for (ND_int ispec=0 ; ispec <nspec; ++ispec )
        {
            ND_int iatom = rint( (*ND_function(ele,Nd_floatS)(&atomic_map,nd_idx{it,ispec}) ) -1) ;

            ELPH_float * get_ptr = ND_function(ele,Nd_floatS)(&atom_pos_temp,nd_idx{it,ispec,0}) ;
            ELPH_float * set_ptr = ND_function(ele,Nd_floatS)(lattice->atomic_pos,nd_idx{iatom,0}) ;
            memcpy(set_ptr, get_ptr, 3*sizeof(ELPH_float));
            (lattice->atom_type)[iatom] = it;
        }
    }

    //
    free(natom_per_type);
    ND_function(destroy, Nd_floatS)(&atomic_map);
    ND_function(destroy, Nd_floatS)(&atom_pos_temp);

    ELPH_float * nGmax = malloc(sizeof(ELPH_float) * nibz ); // max number of gvectors for each wfc in iBZ
    *wfcs = malloc(sizeof(struct WFC) * nibz); // wfcs in iBZ
    struct WFC * wfc_temp = *wfcs ;

    /* allocate arrays of arrays for wfc, gvsc, Fk */
    ND_array(Nd_cmplxS) * wfc_alloc_arrays  = malloc(nibz*sizeof(ND_array(Nd_cmplxS))); // free me 
    ND_array(Nd_floatS) * gvec_alloc_arrays = malloc(nibz*sizeof(ND_array(Nd_floatS))); // free me
    ND_array(Nd_floatS) * Fk_alloc_arrays   = malloc(nibz*sizeof(ND_array(Nd_floatS))); // free me

    if (my_rank == 0) quick_read(dbid, "WFC_NG", nGmax);
    /* Bcast ELPH_float * nGmax */
    mpi_error = MPI_Bcast(nGmax, nibz, ELPH_MPI_float, 0, MPI_COMM_WORLD );

    ND_array(Nd_floatS) totalGvecs, Gvecidxs ;
    
    if (my_rank == 0) 
    {
        ND_function(init,Nd_floatS) (&totalGvecs,  0, NULL); // all gvectors
        ND_function(init,Nd_floatS) (&Gvecidxs,    0, NULL); // gvector indices in the above gvectors

        ND_function(readVar, Nd_floatS) (dbid, "G-VECTORS", &totalGvecs);
        ND_function(readVar, Nd_floatS) (dbid, "WFC_GRID",  &Gvecidxs);
        NC_close_file(dbid); // close ns.db1 file
    }
    Bcast_ND_arrayFloat(&totalGvecs, true, 0, MPI_COMM_WORLD);
    Bcast_ND_arrayFloat(&Gvecidxs, true, 0, MPI_COMM_WORLD);

    // ! Warning, Only read only mode for opening files
    for (ND_int ik = 0 ; ik <nibz ; ++ik)
    {   
        /*set total pws */
        (wfc_temp+ik)->npw_total = rint(nGmax[ik]);

        ND_int pw_per_core = ((wfc_temp+ik)->npw_total)/npw_cpus;
        ND_int pw_rem      = ((wfc_temp+ik)->npw_total)%npw_cpus;
        ND_int pw_this_cpu = pw_per_core;
        if (krank < pw_rem) ++pw_this_cpu ;
            
        /* set the local number of pw's */
        (wfc_temp+ik)->npw_loc = pw_this_cpu;

        ND_int G_shift = pw_per_core*krank;
            
        if (pw_rem !=0)
        {
            if (krank < pw_rem) G_shift += krank ;
            else G_shift += pw_rem;
        }

        (wfc_temp+ik)->gvec = gvec_alloc_arrays+ik ;
            
        alloc_and_set_Gvec((wfc_temp+ik)->gvec, ik, &totalGvecs, &Gvecidxs, \
                                lat_param, pw_this_cpu,G_shift);
            
        /* initiate, allocate and load wfcs*/
        (wfc_temp+ik)->wfc = wfc_alloc_arrays+ik;
        ND_function(init,Nd_cmplxS)((wfc_temp+ik)->wfc, 4, \
        nd_idx{lattice->nspin, lattice->nbnds, lattice->nspinor, pw_this_cpu}); 
        // (nspin, bands, nspinor, npw)
        ND_function(malloc,Nd_cmplxS)((wfc_temp+ik)->wfc);

        get_wfc_from_save((wfc_temp+ik)->wfc->strides[0], ik, nibz, \
        lattice->nspin, lattice->nspinor, lattice->start_band, \
        lattice->nbnds, pw_this_cpu,G_shift, SAVEdir, temp_str, \
        (wfc_temp+ik)->wfc->data, commK);

        /* initiate, allocate and load Fk (Kleinbylander Coefficients)*/
        (wfc_temp+ik)->Fk = Fk_alloc_arrays+ik ;
        ND_function(init,Nd_floatS)((wfc_temp+ik)->Fk, 0, NULL); 
        sprintf(temp_str, "%s/ns.kb_pp_pwscf_fragment_%d", SAVEdir, (int)(ik+1) ) ;  // fix it for abinit 
        /* Abinit has a aditional spin dimension instead of 2*n projectors */
        if ((retval = nc_open_par(temp_str, NC_NOWRITE, commK, MPI_INFO_NULL, &ppid))) ERR(retval);
        sprintf(temp_str, "PP_KB_K%d", (int)(ik+1)) ;  // fix be for abinit
        int varid_temp;
        if ((retval = nc_inq_varid(ppid, temp_str, &varid_temp))) ERR(retval); // get the varible id of the file
        // collective IO
        if ((retval = nc_var_par_access(ppid, varid_temp, NC_COLLECTIVE))) ERR(retval); // NC_COLLECTIVE or NC_INDEPENDENT
        ND_function(readVar_sub, Nd_floatS)(ppid, temp_str, (wfc_temp+ik)->Fk, ND_ALL,ND_ALL,nd_idx{G_shift,pw_this_cpu,1} );
        
        NC_close_file(ppid);
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    /* Free temp gvec arrays */
    ND_function(destroy, Nd_floatS)(&totalGvecs);
    ND_function(destroy, Nd_floatS)(&Gvecidxs);

    /* pseudo data */

    ELPH_float * Zval = malloc(sizeof(ELPH_float)*pseudo->ntype); // this needs to read from pseudo pots // data not set here
    pseudo->Zval = Zval ;
    ND_array(Nd_floatS)* kb_arrays = malloc(sizeof(ND_array(Nd_floatS)) *( (3*pseudo->ntype) + 2 ) );

    pseudo->Vloc_atomic = kb_arrays;
    pseudo->r_grid      = kb_arrays + pseudo->ntype ;
    pseudo->rab_grid    = kb_arrays + (2*pseudo->ntype);
    pseudo->PP_table    = kb_arrays + (3*pseudo->ntype);
    pseudo->Fsign       = kb_arrays + (3*pseudo->ntype) + 1 ;
    
    if (my_rank == 0)
    {
        ND_function(init,Nd_floatS) (pseudo->PP_table,  0, NULL); 
        ND_function(init,Nd_floatS) (pseudo->Fsign,    0, NULL); 

        sprintf(temp_str, "%s/ns.kb_pp_pwscf", SAVEdir) ;  // fix be for abinit 
        NC_open_file(temp_str, 'r', &ppid);
    
        {
            ELPH_float kb_pars[4];
            quick_read(ppid, "PARS", kb_pars);
            pseudo->lmax = rint(kb_pars[3]); // yambo stores lmax+1 ! FIX ME
        }

        ND_function(readVar, Nd_floatS) (ppid, "PP_TABLE", pseudo->PP_table);
        ND_function(readVar, Nd_floatS) (ppid, "PP_KBS",   pseudo->Fsign);

        NC_close_file(ppid); 
    }
    /* Bcast PP_table and Fsign */
    Bcast_ND_arrayFloat(pseudo->PP_table, true, 0, MPI_COMM_WORLD);
    Bcast_ND_arrayFloat(pseudo->Fsign, true, 0, MPI_COMM_WORLD);
    // Bcast int pseudo->lmax
    mpi_error = MPI_Bcast(&(pseudo->lmax), 1, MPI_INT, 0, MPI_COMM_WORLD );

    // compute f-coeffs
    alloc_and_Compute_f_Coeff(lattice, pseudo); 

    /* Read upf data */
    /* First get the pseudo pots order */
    ND_int * pseudo_order = malloc(sizeof(ND_int)*pseudo->ntype);

    if (my_rank == 0)
    {
        for (ND_int ipot1 = 0 ; ipot1<pseudo->ntype ; ++ipot1)
        {
            char temp_ele[3];
        
            sprintf(temp_str, "%s/%s", pseudo_dir,pseudo_pots[ipot1]) ; 
            /* read elements from pseudo pots */
            get_upf_element(temp_str, temp_ele); // only single process !
            bool found = false ;

            for (ND_int ipot2 = 0 ; ipot2<pseudo->ntype ; ++ipot2)
            {
                if(strcmp(temp_ele,atom_symbols + 3*ipot2)  == 0)
                {
                    found = true;
                    pseudo_order[ipot1] = ipot2;
                    break;
                }
            }
            if (!found)
            {
                printf("Pseudo for element %s not found \n",temp_ele);
                exit(EXIT_FAILURE);
            }
        }
    }
    // Bcast pseudo_order[pseudo->ntype]
    mpi_error = MPI_Bcast(pseudo_order, pseudo->ntype, ELPH_MPI_ND_INT, 0, MPI_COMM_WORLD);
    
    /* Get data from upfs */
    if (my_rank == 0)
    {
        for (ND_int ipot = 0 ; ipot<pseudo->ntype ; ++ipot)
        {   
            ND_int iorder = pseudo_order[ipot];
            sprintf(temp_str, "%s/%s", pseudo_dir,pseudo_pots[ipot]) ; 
            parse_upf2(temp_str, Zval + iorder, pseudo->Vloc_atomic + iorder, \
                pseudo->r_grid + iorder, pseudo->rab_grid + iorder);
        }
    }

    // Bcast all the pseudo information
    mpi_error = MPI_Bcast(Zval, pseudo->ntype, ELPH_MPI_float, 0, MPI_COMM_WORLD );
    for (ND_int itype = 0; itype<pseudo->ntype; ++itype)
    {
        Bcast_ND_arrayFloat(pseudo->Vloc_atomic, true, 0, MPI_COMM_WORLD);
        Bcast_ND_arrayFloat(pseudo->r_grid, true, 0, MPI_COMM_WORLD);
        Bcast_ND_arrayFloat(pseudo->rab_grid, true, 0, MPI_COMM_WORLD);
    }

    // reuse pseudo_order buffer to find the ngrid max
    for (ND_int ipot = 0 ; ipot<pseudo->ntype ; ++ipot)
    {
        pseudo_order[ipot] = (pseudo->Vloc_atomic+ipot)->dims[0];
    }

    pseudo->ngrid_max = find_maxint(pseudo_order,pseudo->ntype );
    lattice->npw_max  = find_maxfloat(nGmax, nibz) ; // find the max number of pws i.e max(nGmax)

    // free all buffers
    free(pseudo_order);
    free(atom_symbols);
    free(elements);
    free(nGmax);
    free(temp_str);
}


void free_save_data(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo)
{   
    /* free fCoeff*/
    free_f_Coeff(lattice, pseudo); // NOTE this function must be called before freeing PP_table

    ND_int nkiBZ = lattice->kpt_iredBZ->dims[0];
    /* Free wavefunctions */
    for (ND_int ik =0 ; ik<nkiBZ; ++ik)
    {
        ND_function(destroy, Nd_cmplxS)((wfcs+ik)->wfc);
        ND_function(destroy, Nd_floatS)((wfcs+ik)->gvec);
        ND_function(destroy, Nd_floatS)((wfcs+ik)->Fk);
    }

    free(wfcs->wfc);
    free(wfcs->gvec);
    free(wfcs->Fk);
    free(wfcs);

    // free local and non local pseudo data
    for (ND_int i = 0 ; i< ((3*pseudo->ntype) + 2); ++i)
    {
        ND_function(destroy, Nd_floatS)(pseudo->Vloc_atomic + i);
    }
    free(pseudo->Vloc_atomic);
    free(pseudo->Zval);

    /* free lattice data */
    for (int i = 0; i<7; ++i) ND_function(destroy, Nd_floatS)(lattice->alat_vec + i); // 7 and 8 are already freed
    free(lattice->alat_vec);
    free(lattice->time_rev_array);
    ND_function(destroy, i)(lattice->kmap);
    free(lattice->kmap);
    free(lattice->atom_type);

}






// ----


static void alloc_and_set_Gvec(ND_array(Nd_floatS) * gvec, ND_int ik, ND_array(Nd_floatS) * totalGvecs, \
                    ND_array(Nd_floatS) * Gvecidxs, ELPH_float * lat_param, ND_int nG, ND_int nG_shift)
{   
    ND_function(init,Nd_floatS)(gvec, 2, nd_idx{nG,3} );
    ND_function(malloc,Nd_floatS)(gvec);

    ELPH_float * gidx_temp = ND_function(ele,Nd_floatS)(Gvecidxs,nd_idx{ik,0});
    // #pragma openmp parallel for ? FIX ME
    for (ND_int ig = 0 ; ig <nG ; ++ig )
    {
        ND_int gidx = rint(gidx_temp[ig+nG_shift]-1) ;
        if (gidx < 0)
        {
            printf("Wrong g indices \n");
            exit(EXIT_FAILURE);
        }
        ELPH_float * gvec_temp = ND_function(ele,Nd_floatS)(gvec, nd_idx{ig,0}) ;

        gvec_temp[0] = (*ND_function(ele,Nd_floatS)(totalGvecs,nd_idx{0,gidx}))/lat_param[0] ;
        gvec_temp[1] = (*ND_function(ele,Nd_floatS)(totalGvecs,nd_idx{1,gidx}))/lat_param[1] ;
        gvec_temp[2] = (*ND_function(ele,Nd_floatS)(totalGvecs,nd_idx{2,gidx}))/lat_param[2] ;

    }
}



static void get_wfc_from_save(ND_int spin_stride_len, ND_int ik, ND_int nkiBZ, \
            ND_int nspin, ND_int nspinor, ND_int start_band, ND_int nbnds, \
            ND_int nG, ND_int G_shift, const char * save_dir, char * work_array, \
            ELPH_cmplx * out_wfc, MPI_Comm comm)
{   
    int wfID, retval ;
    // NO OPENMP !! , Not thread safe
    for (ND_int is = 0; is<nspin; ++is)
    {   
        sprintf(work_array, "%s/ns.wf_fragments_%d_1", save_dir, (int)( is*nkiBZ + (ik+1) ) ) ;

        //NC_open_file(work_array, 'r', &wfID);
        if ((retval = nc_open_par(work_array, NC_NOWRITE, comm, MPI_INFO_NULL, &wfID))) ERR(retval);

        sprintf(work_array, "WF_COMPONENTS_@_SP_POL%d_K%d_BAND_GRP_1", (int)(is+1) , (int)(ik+1) ) ;

        //// (nspin, bands, nspinor, npw)
        size_t startp[4] = {start_band-1, 0      , G_shift , 0} ;
        size_t countp[4] = {nbnds       , nspinor, nG      , 2};
        quick_read_sub(wfID, work_array, startp, countp, out_wfc + is*spin_stride_len);

        NC_close_file(wfID); 
    }
}


static void quick_read(const int ncid, char* var_name, void * data_out)
{   /*  Serial read
        load the entire varible data to data_out pointer 
    */
    int varid, retval;

    if ((retval = nc_inq_varid(ncid, var_name, &varid))) ERR(retval); // get the varible id of the file
    
    if ((retval = nc_get_var(ncid, varid, data_out))) ERR(retval); //get data in floats

}


static void quick_read_sub(const int ncid, char* var_name, const size_t * startp, \
                            const size_t * countp, void * data_out)
{   /*  Serial read
        load the slice of varible data to data_out pointer 
    */
    int varid, retval;

    if ((retval = nc_inq_varid(ncid, var_name, &varid))) ERR(retval); // get the varible id of the file

    // collective IO
    if ((retval = nc_var_par_access(ncid, varid, NC_COLLECTIVE))) ERR(retval); // NC_COLLECTIVE or NC_INDEPENDENT

    if ((retval = nc_get_vara(ncid, varid, startp, countp,data_out))) ERR(retval); //get data in floats

}


