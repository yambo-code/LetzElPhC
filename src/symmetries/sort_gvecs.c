#include "symmetries.h"
#include "../common/kdtree/kdtree.h"

/*
This contains 3 function
1) function to build kd_tree for gvecs
2) function to get indices for rotated g 
3) free kd_tree
*/

void build_gtree(void ** gtree, const ND_int ngvecs, \
    const ELPH_float * restrict lat_vec, \
    ELPH_float * restrict gvecs)
{   
    /*
    Builds kd_tree for given gvecs
    lat_vec --> a[:,i]
    gvecs : cart units

    Note: the output gtree must be freed with kd_free(gtree);
    */
    if (gvecs == NULL) return;

    *gtree = kd_create(3);
    void * gvec_tree = *gtree;

    // build the tree
    for (ND_int i = 0 ; i<ngvecs ;  ++i)
    {   
        // create a single precision tree
        ELPH_float gcrys[3] = {0,0,0};
        MatVec3f(lat_vec, gvecs + 3*i, true, gcrys);
        int kd_err = kd_insert3(gvec_tree, gcrys[0], gcrys[1], gcrys[2], gvecs + 3*i);
        if (kd_err != 0) error_msg("Building kdtree for gvecs failed");
    }
}

void free_gtree(void * gtree)
{   
    /*
    Frees the  allocated gvector tree
    */
    kd_free(gtree);
}

void sort_gvec(void *gtree, const ELPH_float * restrict tree_gvecs,
    const ND_int ngvecs, const ELPH_float * restrict sym_mat, \
    const ELPH_float * restrict G0, const ELPH_float * restrict lat_vec, \
    const ELPH_float * restrict gvecs_out, ND_int * restrict out_idx)
{
    /*
    This functions returns indices of the Sym_mat@gvecs_out + G0 in the 
    tree_gvecs(kdtree is build for these tree_gvecs).
    if some gvecs in gvecs_out are not found in kdtree, then index 
    is set to -1
    tree_gvecs : starting gvec ptr used to create gtree
    ngvecs : number of gvecs_out
    G0 and gvecs_out are in cart

    Complexity : ~O(nlogn)
    */
    if (out_idx == NULL) return;
    // first compute lat^T@sym_mat
    ELPH_float trans_mat[9];
    ELPH_float ulm_vec[3]={0,0,0};

    if (sym_mat != NULL) Gemm3x3f(lat_vec, 'T', sym_mat,  'N', trans_mat);
    else transpose3x3f(lat_vec, trans_mat);

    if (G0 != NULL) MatVec3f(lat_vec, G0, true, ulm_vec);

    for (ND_int i = 0 ; i<ngvecs ;  ++i)
    {   
        // create a single precision tree
        ELPH_float gcrys[3] = {0,0,0};
        MatVec3f(trans_mat, gvecs_out + 3*i, false, gcrys);
        void * kd_out = kd_nearest3(gtree, gcrys[0]+ulm_vec[0], gcrys[1]+ulm_vec[1], gcrys[2]+ulm_vec[2]);
        if (kd_res_size(kd_out) == 0) out_idx[i] = -1;
        else
        {   
            double nearest_pt[3];
            ELPH_float * tmp_ptr = kd_res_item(kd_out,nearest_pt);
            // check distance
            double dist = 0;
            for (int i =0; i<3; ++i) dist += (gcrys[i]-nearest_pt[i])*(gcrys[i]-nearest_pt[i]);
            dist = sqrt(dist);
            if (dist < ELPH_EPS) out_idx[i] = (tmp_ptr-tree_gvecs)/3 ;
            else out_idx[i] = -1;
        }
        kd_res_free(kd_out);
    }
}

