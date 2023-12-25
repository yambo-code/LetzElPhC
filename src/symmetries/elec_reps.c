#include "symmetries.h"


void electronic_reps(struct WFC * wfcs, struct Lattice * lattice, \
    ELPH_float * restrict sym_mat,  ELPH_float * restrict sym_v,  \
    const ND_int ikBZ, ELPH_cmplx * restrict Dkmn_rep, MPI_Comm commK)
{
    /*
    This is a function to compute representation matrices for the given symmetry and k point
    */
    
    // compute the Rk vector and find it in the list of k-points
    ELPH_float Rk_vec[3] = {0,0,0};

    // start of small scope
    { 
    ELPH_float Rk_tmp[3] = {0,0,0};
    MatVec3f(sym_mat, lattice->kpt_fullBZ->data + 3*ikBZ, false, Rk_tmp);
    // convert to crystal coordinates
    MatVec3f(lattice->alat_vec->data,Rk_tmp, true, Rk_vec);
    } 
    // end of scope
    
    // now find the index of the rotated k point
    ND_int iRkBZ = -1;
    // find the index
    for (ND_int ik = 0; ik < lattice->kpt_fullBZ_crys->dims[0], ++ik)
    {
        ELPH_float * restrict ik_vec_tmp = lattice->kpt_fullBZ_crys->data + 3*ik ;
        ELPH_float sum = 0;
        for (int i = 0; i < 3 ; ++i)
        {   
            ELPH_float diff_tmp = ik_vec_tmp[i]-Rk_vec[i];
            diff_tmp = diff_tmp-rint(diff_tmp);
            sum += diff_tmp*diff_tmp;
        }
        sum = sqrt(sum);
        if (sum < ELPH_EPS)
        {
            iRkBZ = ik;
            break;
        }
    }

    if (iRkBZ < 0) error_msg("Rotated k point not found. Either Wrong Phonon symmetry or using non-uniform kgrid");


    /* 
    Now we compute the Dmats i.e $ \langle Rk | U(R) | k \rangle $
    if k = Sym1 * ik1 and Rk = Sym2 * ik2 then 
    Dmats = $ \langle Sym2 * ik2 | U(R) | Sym1 * ik1  \rangle $
    */

    // now get the corresponding iBZ point and symmetry for k and Rk
    int ik1      = *(lattice->kmap->data + ikBZ*2)      ;
    int Sym1     = *(lattice->kmap->data + ikBZ*2 + 1)  ;

    int ik2      = *(lattice->kmap->data + iRkBZ*2)      ;
    int Sym2     = *(lattice->kmap->data + iRkBZ*2 + 1)  ;

    /* 
    note to avoid any copying of memory, we apply all symmetry operations 
    to only one wave-function. here we apply on ik2
    

    => Dmats = $ \langle ik2 | {U^\dagger(Sym2) U(R) U(Sym1)} | ik1  \rangle $
    where operator in brackets is applied on left
    if Sym1 is time reversal then we must conjugate i.e Dmats = conj(Dmats)
    
    (U^\dagger(Sym2) U(R) U(Sym1))^\dagger  is applied to ik2
    */

    


}







