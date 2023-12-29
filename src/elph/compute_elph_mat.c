/*
This function computes electron-matrix elements for all q's
*/
#include "elph.h"

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    struct usr_input * input_data;

    read_input_file(argv[1], &input_data);

    ELPH_float qpt[3] = {0.0,0.0,0.0};

    ND_int NK = input_data->nkpool;
    ND_int NQ  = input_data->nqpool;
    char * SAVEDIR  = input_data->save_dir ;
    ND_int FIRST_BAND   = input_data->start_bnd ;
    ND_int LAST_BAND  = input_data->end_bnd ;
    char * PSEUDO_DIR  = input_data->pseudo_dir ;
    char * DVSCF_NC  = input_data->dvscf_file ;

    MPI_Comm commQ, commK ;

    create_parallel_comms(NQ, NK, &commQ, &commK);

    struct Lattice lattice ;
    struct Pseudo pseudo ;
    struct WFC * wfcs;

    lattice.dimension = input_data->dimension;
    char ** pseudo_pots = input_data->pseudos ;

    ND_int FFT_dims[3];
    ND_int nq;
    
    get_FFT_dims(DVSCF_NC, &nq, FFT_dims);


    int npw_cpus, krank;
    MPI_Comm_rank(commK, &krank);

    if (krank ==0 ) printf("Reading \n");
    
    read_and_alloc_save_data(SAVEDIR, commQ, commK, FIRST_BAND, LAST_BAND, \
                &wfcs, PSEUDO_DIR, pseudo_pots, &lattice, &pseudo,FFT_dims);

    ND_array(Nd_cmplxS) eigVec,  dVscf;

    // ND_function(init, Nd_cmplxS) (&eigVec, 0, NULL);
    // ND_function(init, Nd_cmplxS) (&dVscf, 0, NULL);
    read_dvscfq(DVSCF_NC, &eigVec, &lattice, &dVscf,0, commK);

    if (krank ==0 ) printf("Computing el-ph \n");
    
    /*compute el-ph*/

    ELPH_cmplx * elph_kq = NULL;


    if (krank ==0 )
    {   
        ND_int size_elph = lattice.kmap->dims[0]*eigVec.dims[0]*lattice.nspin*lattice.nbnds*lattice.nbnds;
        elph_kq = malloc(sizeof(ELPH_cmplx)*size_elph);
        for (ND_int isize = 0 ; isize<size_elph ; ++isize) elph_kq[isize] = 0;
    }
    //((k, nmodes, nspin, nbands, nbands))

    compute_elph(wfcs, &lattice, &pseudo, qpt,  &eigVec, &dVscf, elph_kq, commQ, commK);


    ELPH_cmplx Ry2Ha = pow(2,-1.5);
    
    int nk11 = 1 ;
    int nb11 = 26 ;
    int nb22 = 36 ; 

    if (krank ==0 )
    {   
        ND_array(Nd_cmplxS) elph_w;
        ND_function(init, Nd_cmplxS)(&elph_w, 5, nd_idx{lattice.kmap->dims[0], eigVec.dims[0], lattice.nspin, lattice.nbnds, lattice.nbnds  } );
        elph_w.data = elph_kq;
        ND_function(write, Nd_cmplxS)("nc.elph", "elph_mat", &elph_w, (char * [5]) {"nk", "modes", "nspin", "nbndK", "nbndKq"},NULL);
        ND_function(uninit, Nd_cmplxS)(&elph_w);

        for (int imode = 0; imode < 9 ; ++imode)
        {   
            ND_int shi_idx = nb22 + lattice.nbnds*nb11 + (lattice.nbnds*lattice.nbnds)*(0 + imode*lattice.nspin + nk11*lattice.nspin*eigVec.dims[0]);
            ELPH_cmplx val = Ry2Ha * elph_kq[shi_idx];  // (nq,nk,nv,is,nbk,nbk+q)
            printf("%10e + %10e I \n",creal(val),cimag(val));
        }
    }

    ND_array(Nd_cmplxS) Dmat_w;
    if (krank ==0 )
    {   
        printf("Computing electronic representation matrices...\n");
        ND_function(init, Nd_cmplxS)(&Dmat_w, 5, nd_idx{lattice.sym_mat->dims[0],lattice.kmap->dims[0], lattice.nspin, lattice.nbnds, lattice.nbnds});
        ND_function(malloc, Nd_cmplxS) (&Dmat_w);
    }

    for (ND_int isym=0; isym< lattice.sym_mat->dims[0]; ++isym)
    {
        for (ND_int ikBZ=0; ikBZ< lattice.kmap->dims[0]; ++ikBZ)
        {   
            ELPH_cmplx * Dkmn_rep_ptr = NULL;
            if (krank ==0 ) Dkmn_rep_ptr = ND_function(ele, Nd_cmplxS)(&Dmat_w,nd_idx{isym,ikBZ,0,0,0});

            electronic_reps(wfcs, &lattice, lattice.sym_mat->data + 9*isym,  (ELPH_float []){0,0,0}, \
            lattice.time_rev_array[isym], ikBZ, Dkmn_rep_ptr, commK);
        }

    }

    if (krank ==0 ) 
    {
        ND_function(write, Nd_cmplxS)("nc.Dmats", "Dmats", &Dmat_w, (char * [5]) {"nsym", "nk", "nspin", "nbndb", "nbnda"},NULL);
        ND_function(destroy, Nd_cmplxS) (&Dmat_w);
    }


    // cleanup
    free_usr_input(input_data);
    free_save_data(wfcs, &lattice, &pseudo);
    free(elph_kq);
    ND_function(destroy,Nd_cmplxS)(&eigVec);
    ND_function(destroy,Nd_cmplxS)(&dVscf);

    free_parallel_comms(&commQ, &commK);

    fftw_fun(cleanup)();

    MPI_Finalize();

    return 0;

}
