#include "gsort.h"
/*
This file contain g vectors sorting functions
*/

#define COMP_AB(A,B) ( (A > B)-(A < B) )

static int gsort_cmp(const void *a, const void *b);



void Sorted_gvecs_idxs(const ND_int npw, ELPH_float * gvecs, ND_int * indices )
{
    /*
    Sort gvecs and group all gvecs with same (x,y).

    // gvecs in crystal coordinates
    */
    if (gvecs == NULL) return ; 

    ELPH_float ** gvec_ptrs = malloc(sizeof(ELPH_float *)*npw);

    /* fill the struct */
    for (ND_int i = 0 ; i< npw; ++i)
    {
        gvec_ptrs[i] = gvecs + 3*i; 
    }

    qsort(gvec_ptrs, npw, sizeof(ELPH_float *), gsort_cmp);

    // store the sorted indices
    for (ND_int i = 0 ; i< npw; ++i) indices[i] = (gvec_ptrs[i]-gvecs)/3;
    
    free(gvec_ptrs);

}



static int gsort_cmp(const void *a, const void *b)
{   
    /*
    First sorts along x, then y and finally z
    */
    ELPH_float * v1 = *(ELPH_float **)a;
    ELPH_float * v2 = *(ELPH_float **)b;

    ND_int Gax, Gay, Gaz;
    ND_int Gbx, Gby, Gbz;
    Gax = rint(v1[0]);
    Gay = rint(v1[1]);
    Gaz = rint(v1[2]);

    Gbx = rint(v2[0]);
    Gby = rint(v2[1]);
    Gbz = rint(v2[2]);

    if (Gax == Gbx)  
    {
        if (Gay == Gby) return COMP_AB(Gaz,Gbz);
        else            return COMP_AB(Gay,Gby);
    }
    else return  COMP_AB(Gax,Gbx); 
}
