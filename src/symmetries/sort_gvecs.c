// #include "symmetries.h"

// /*
// This contains 3 function
// 1) function to build kd_tree for gvecs
// 2) function to sorted indices for rotated g 
// 3) function to destroy the tree

// !!! Note all these functions are must run on
// only 1 cpu
// */


// void build_gtree(void ** gtree, const ND_int ngvecs, \
//     const ELPH_float * restrict lat_vec, \
//     const ELPH_float * restrict gvecs)
// {   
//     /*
//     Builds kd_tree for given gvecs
//     lat_vec --> a[:,i]

//     */

//     void * gvec_tree = *gtree;
//     gvec_tree = kd_create(3);

//     // build the tree
//     for (ND_int i = 0 ; i<ngvecs ;  ++i)
//     {   
//         // create a single precision tree
//         ELPH_float gcrys[3] = {0,0,0};

//         assert(kd_insert3f(gvec_tree, gvecs[3*i], gvecs[3*i+1], gvecs[3*i+2], gvecs + 3*i) == 0 );
//     }


// }


// void sort_gvec()