/*
This function computes electron-matrix elements for all q's
*/
#include "elph.h"

#define NK 2
#define NQ 1
#define SAVEDIR "/Users/murali/phd/one_phonon_raman/wse2/SAVE" //"/Users/murali/phd/one_phonon_raman/si/bse/si_data/raman/silicon.save/SAVE" 
#define FIRST_BAND  1
#define LAST_BAND 40
#define PSEUDO_DIR "/Users/murali/phd/one_phonon_raman/wse2"
#define DVSCF_NC "/Users/murali/phd/one_phonon_raman/wse2/nc.dVscf_new"

// int main(int argc, char* argv[])
// {
//     MPI_Init(&argc, &argv);

//     MPI_Comm commQ, commK ;

//     create_parallel_comms(NQ, NK, &commQ, &commK);


//     struct Lattice lattice ;
//     struct Pseudo pseudo ;
//     struct WFC * wfcs;

//     lattice.dimension = '2';

//     char ** pseudo_pots = malloc(sizeof(char *)*10);
//     //pseudo_pots[0] = "Si.upf" ;//"Si_PBE_SR.SG15v1.2.UPF" ;//"Mo.upf";
//     pseudo_pots[0] = "W_PBE_nof_FR.SG15v1.2.UPF";
//     pseudo_pots[1] = "Se_PBE_FR.SG15v1.2.UPF";
//     const ND_int FFT_dims[3] = {45,45,270};

//     printf("Reading \n");
    
//     read_and_alloc_save_data(SAVEDIR, commQ, commK, FIRST_BAND, LAST_BAND, \
//                 &wfcs, PSEUDO_DIR, pseudo_pots, &lattice, &pseudo,FFT_dims);

//     ND_array(Nd_cmplxS) eigVec,  dVscf , eigVecq,  dVscfq, elphkq;

//     ND_function(init, Nd_cmplxS) (&eigVec, 0, NULL);
//     ND_function(init, Nd_cmplxS) (&dVscf, 0, NULL);

//     int ncid_dvscf ;
//     NC_open_file(DVSCF_NC, 'r', &ncid_dvscf);

//     ND_function(readVar, Nd_cmplxS) (ncid_dvscf, "dVscfs", &dVscf);
//     ND_function(readVar, Nd_cmplxS) (ncid_dvscf, "ph_pol_vec", &eigVec);

//     ND_function(init_strip_dims, Nd_cmplxS) (&dVscf, 1, &dVscfq);
//     ND_function(init_strip_dims, Nd_cmplxS) (&eigVec, 1, &eigVecq);

//     ND_function(strip_dims, Nd_cmplxS) (&dVscf,  1, nd_idx{0},  &dVscfq);
//     ND_function(strip_dims, Nd_cmplxS) (&eigVec, 1, nd_idx{0}, &eigVecq);


//     printf("dims %lld , %lld \n",dVscf.dims[1],eigVec.dims[1]);
//     NC_close_file(ncid_dvscf); 

//     printf("Computing el-ph \n");
//     ELPH_float qpt[3] = {0.0,0.0,0.0};
    
//     int nq = 1;

//     ND_function(init, Nd_cmplxS) (&elphkq, 6, nd_idx{nq,  lattice.kmap->dims[0],  eigVecq.dims[0],  lattice.nspin, lattice.nbnds, lattice.nbnds});
//     ND_function(malloc, Nd_cmplxS) (&elphkq);
//     ND_function(set_all, Nd_cmplxS) (&elphkq,0.0);
//     /*compute el-ph*/
    
//     free_parallel_comms(&commQ, &commK);
//     MPI_Finalize()

//     return 0;

// }