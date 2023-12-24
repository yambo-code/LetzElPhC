// #include "symmetries.h"


// void electronic_reps(struct WFC * wfcs, struct Lattice * lattice, \
//     ELPH_float * restrict sym_mat,  ELPH_float * restrict sym_v,  \
//     const ND_int ikBZ, ELPH_cmplx * restrict Dkmn_rep, MPI_Comm commK)
// {
//     /*
//     This is a function to compute representation matrices for the given symmetry and k point
//     */
    
//     /*
//     First, we perform a binary search and get the indices of rotate indices
//     */

    
//     int ikibz    = *(lattice->kmap->data + ikBZ*2)      ;
//     int ksym     = *(lattice->kmap->data + ikBZ*2 + 1)  ;

//     // compute the Rk vector and find it in the list of k-points

//     ELPH_cmplx Rk_vec[3] = {0,0,0};

//     // start of small scope
//     { 
//     ELPH_cmplx Rk_tmp[3] = {0,0,0};
//     MatVec3f(sym_mat, lattice->kpt_fullBZ->data + 3*ikBZ, false, Rk_tmp);
//     // convert to crystal coordinates
//     MatVec3f(lattice->alat_vec->data,Rk_tmp, true, Rk_vec);
//     } 
//     // end of scope



// }





