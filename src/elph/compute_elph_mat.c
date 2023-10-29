/*
This function computes electron-matrix elements for all q's
*/
#include "elph.h"

#define NK 1
#define NQ 1
#define SAVEDIR "/Users/murali/phd/one_phonon_raman/wse2/SAVE" //"/Users/murali/phd/one_phonon_raman/si/bse/si_data/raman/silicon.save/SAVE" 
#define FIRST_BAND  1
#define LAST_BAND 40
#define PSEUDO_DIR "/Users/murali/phd/one_phonon_raman/wse2"
#define DVSCF_NC "/Users/murali/phd/one_phonon_raman/wse2/nc.dVscf"

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    MPI_Comm commQ, commK ;

    create_parallel_comms(NQ, NK, &commQ, &commK);

    struct Lattice lattice ;
    struct Pseudo pseudo ;
    struct WFC * wfcs;

    lattice.dimension = '2';

    char ** pseudo_pots = malloc(sizeof(char *)*10);
    //pseudo_pots[0] = "Si.upf" ;//"Si_PBE_SR.SG15v1.2.UPF" ;//"Mo.upf";
    pseudo_pots[0] = "W_PBE_nof_FR.SG15v1.2.UPF";
    pseudo_pots[1] = "Se_PBE_FR.SG15v1.2.UPF";


    ND_int FFT_dims[3];
    ND_int nq;
    
    
    get_FFT_dims(DVSCF_NC, &nq, FFT_dims);

    printf("Reading \n");
    
    read_and_alloc_save_data(SAVEDIR, commQ, commK, FIRST_BAND, LAST_BAND, \
                &wfcs, PSEUDO_DIR, pseudo_pots, &lattice, &pseudo,FFT_dims);

    ND_array(Nd_cmplxS) eigVec,  dVscf;

    // ND_function(init, Nd_cmplxS) (&eigVec, 0, NULL);
    // ND_function(init, Nd_cmplxS) (&dVscf, 0, NULL);
    read_dvscfq(DVSCF_NC, &eigVec, &lattice, &dVscf,0, commK);

    printf("Computing el-ph \n");
    ELPH_float qpt[3] = {0.0,0.0,0.0};
    

    /*compute el-ph*/

    ELPH_cmplx * elph_kq = NULL;
    int npw_cpus, krank;
    MPI_Comm_rank(commK, &krank);

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
        for (int imode = 0; imode < 9 ; ++imode)
        {   
            ND_int shi_idx = nb22 + lattice.nbnds*nb11 + (lattice.nbnds*lattice.nbnds)*(0 + imode*lattice.nspin + nk11*lattice.nspin*eigVec.dims[0]);
            ELPH_cmplx val = Ry2Ha * elph_kq[shi_idx];  // (nq,nk,nv,is,nbk,nbk+q)
            printf("%10e + %10e I \n",creal(val),cimag(val));
        }
    }

    free(elph_kq);
    ND_function(destroy,Nd_cmplxS)(&eigVec);
    ND_function(destroy,Nd_cmplxS)(&dVscf);

    free_parallel_comms(&commQ, &commK);

    fftw_fun(cleanup)();

    MPI_Finalize();

    return 0;

}
