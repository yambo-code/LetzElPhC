#include "symmetries.h"
#include "../wfc/wfc.h"

#define bool2int(bval) (bval ? 1 : 0)
/*
Note that the bool2int() is not needed in practice, the standard 
anyways mandate that true = 1 and false = 0. this is 
just for readability. Any decent compiler will remove this.
*/


void electronic_reps(const struct WFC * wfcs, const struct Lattice * lattice, \
    const ELPH_float * Rsym_mat,  const ELPH_float * Rsym_v, \
    const bool tim_revR, const ND_int ikBZ, ELPH_cmplx * restrict Dkmn_rep, MPI_Comm commK)
{
    /*
    This is a function to compute representation matrices for the given symmetry and k point
    */
    
    // compute the Rk vector and find it in the list of k-points
    ELPH_float Rk_vec[3] = {0,0,0};

    // start of small scope
    { 
    ELPH_float Rk_tmp[3] = {0,0,0};
    MatVec3f(Rsym_mat, lattice->kpt_fullBZ->data + 3*ikBZ, false, Rk_tmp);
    // convert to crystal coordinates
    MatVec3f(lattice->alat_vec->data,Rk_tmp, true, Rk_vec);
    } 
    // end of scope
    
    // now find the index of the rotated k point
    ND_int iRkBZ = -1;
    // find the index
    for (ND_int ik = 0; ik < lattice->kpt_fullBZ_crys->dims[0]; ++ik)
    {
        ELPH_float * ik_vec_tmp = lattice->kpt_fullBZ_crys->data + 3*ik ;
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
    Dmats = $ \langle Sym2 * k2 | U(R) | Sym1 * k1  \rangle $
    */

    // now get the corresponding iBZ point and symmetry for k and Rk
    const int ik1      = *(lattice->kmap->data + ikBZ*2)      ;
    const int iSym1     = *(lattice->kmap->data + ikBZ*2 + 1)  ;

    const int ik2      = *(lattice->kmap->data + iRkBZ*2)      ;
    const int iSym2     = *(lattice->kmap->data + iRkBZ*2 + 1)  ;

    const ELPH_float * Sym1  = lattice->sym_mat->data + 9*iSym1;
    const ELPH_float * Sym2  = lattice->sym_mat->data + 9*iSym2;

    const bool tr1 = lattice->time_rev_array[iSym1];
    const bool tr2 = lattice->time_rev_array[iSym2];

    const ELPH_float * tau1  = ND_function(ele,Nd_floatS)(lattice->frac_trans, nd_idx{iSym1,0});
    const ELPH_float * tau2  = ND_function(ele,Nd_floatS)(lattice->frac_trans, nd_idx{iSym2,0});
    
    /* 
    note to reduce copying of memory, we apply all symmetry operations 
    to only one wave-function. here we apply on ik2
    

    => Dmats = $ \langle ik2 | {U^\dagger(Sym2) U(R) U(Sym1)} | ik1  \rangle $
    where operator in brackets is applied on left
    if R or Sym1 is time reversal then we must conjugate i.e Dmats = conj(Dmats)

    (U^\dagger(Sym2) U(R) U(Sym1))^\dagger  is applied to ik2
    */

    // FIX ME R or Sym1 is time reversal then we must conjugate i.e Dmats = conj(Dmats)

    // Get the wfcs in iBZ
    //
    const ELPH_cmplx * wfc_k1 = (wfcs+ik1)->wfc->data ;
    const ELPH_cmplx * wfc_k2 = (wfcs+ik2)->wfc->data ;

    const ELPH_float * gvecs_k1 = (wfcs+ik1)->gvec->data ;
    const ELPH_float * gvecs_k2 = (wfcs+ik2)->gvec->data ;

    const ND_int npw_k1_loc = (wfcs+ik1)->npw_loc ;
    const ND_int npw_k2_loc = (wfcs+ik2)->npw_loc ;

    const ND_int npw_k1_total = (wfcs+ik1)->npw_total ;
    const ND_int npw_k2_total = (wfcs+ik2)->npw_total ;

    // 
    

}







/*
    // In general for three symmetric operations : (A3,v3)@(A2,v2)@(A1,v1) = (A3@A2@A1, v3 + A3*v2 + A3@A2*v1)

    // In our case : A1 = Sym2, A2 = R^{-1}, A3 = Sym1^{-1}

    // If Ax + v is symmetry then its inverse operation is A^{-1}x - A^{-1}v
    
    // In case of time reversal symmetries, 

    // v1 -> v1*(-1)^{T(A1) + T(A2) + T(A3)}
    // v2 -> v2*(-1)^{T(A2) + T(A3)}
    // v3 -> v3*(-1)^{T(A3)}

    // where T(A) = 1 if A is time reversal symmetry. else 0

    // Wavefunction co-efficients Cg = conj(Cg) if (T(A1) + T(A2) + T(A3))%2 = 1


    // */
    
    // //{T(Sym1) + T(R) + T(Sym2)} 
    // bool tr_123 =  (bool2int(tr1) + bool2int(tim_revR) + bool2int(tr2) + 1)%2 ;
    // // {T(Sym1) + T(R)}
    // bool tr_23 =  (bool2int(tr1) + bool2int(tim_revR) + 1)%2 ;
    // // {T(Sym1)}
    // bool tr_3 =  (bool2int(tr1) + 1)%2 ;

    // /* Note : In the above 3 bools, we added +1 on rhs to account for 
    // complex conjugation of the wave-function when computing the sandwich */
    
    // ELPH_cmplx SU2_mat123[4]={1,0,0,1};  // SU^\dagger(S1)@SU^\dagger(R)@SU(S2)

    // // initiate to I_2x2. if nspinor == 1 then only first element in considered
    // // compute SU(2) mats
    // if (lattice->nspinor == 2)
    // { 
    //     ELPH_cmplx SU2_tempS1[4]={1,0,0,1}; //SU^dagger(S1)
    //     ELPH_cmplx SU2_tempR[4]={1,0,0,1}; //SU^dagger(R)
    //     ELPH_cmplx SU2_tempS2[4]={1,0,0,1}; //SU(S2)

    //     ELPH_cmplx SU2_temp[4]={1,0,0,1}; // temp
        
    //     SU2mat(Sym1, lattice->nspinor, true, tr1, SU2_tempS1); //SU^dagger(S1)
    //     SU2mat(Rsym_mat, lattice->nspinor, true, tim_revR, SU2_tempR); //SU^dagger(R)
    //     SU2mat(Sym2, lattice->nspinor, false, tr2, SU2_tempS2); //SU(S2)

    //     // conjugate if required (due to time reversal symmetries)
    //     for (ND_int i = 0; i < (lattice->nspinor*lattice->nspinor); ++i)
    //     {
    //         if (tr_3) SU2_tempS1[i] = conj(SU2_tempS1[i]);
    //         if (tr_23) SU2_tempR[i] = conj(SU2_tempR[i]);
    //         if (tr_123) SU2_tempS2[i] = conj(SU2_tempS2[i]);
    //     }

    //     // compute the product of the 3 su(2) mats i.e 
    //     // SU^\dagger(S1)@SU^\dagger(R)@SU(S2)
        
    //     matmul_Cmpl2x2(SU2_tempS1, SU2_tempR, SU2_temp); // SU^\dagger(S1)@SU^\dagger(R)
    //     matmul_Cmpl2x2(SU2_temp, SU2_tempS2, SU2_mat123); // SU^\dagger(S1)@SU^\dagger(R)@SU(S2)
    // } 


    // // Compute the total rotation matrix and frac trans

    // ELPH_float Sym123[9] = {0,0,0,0,0,0,0,0,0};
    // ELPH_float tau_123[3] = {0,0,0};
    // // !!! FIX ME frac tran fix

    // { // create a small scope
    //     // 1st compute the rotational matrix
    //     ELPH_float Sym_1R_temp[9] = {0,0,0,0,0,0,0,0,0}; //Sym1^{-1}@R^{-1} i.e Sym1^{T}@R^{T}
    //     Gemm3x3f(Sym1, 'T', Rsym_mat,'T',Sym_1R_temp);
    //     Gemm3x3f(Sym_1R_temp, 'N', Sym2,'N',Sym123); //Sym1^{-1}@R^{-1}@Sym2 i.e,  Sym1^{T}@R^{T}@Sym2
    //     // fractional translation
    //     // v1 = tau_2, v2 = -R^{-1}Rsym_v, v3 = -Sym1^{-1}tau_1

    //     ELPH_float v1_tmp[3] = {0,0,0};
    //     ELPH_float v2_tmp[3] = {0,0,0};
        
    //     ELPH_float v1_fac = 1; 
    //     ELPH_float v2_fac = 1; 
    //     ELPH_float v3_fac = 1; 

    //     if (tr_123) v1_fac = -1;
    //     if(tr_23) v2_fac = -1;
    //     if(tr_3) v3_fac = -1;

    //     // v1_tmp = Sym1^{T}@R^{T}*(tau_2-Rsym_v)
    //     for (int i = 0 ; i<3 ; ++i) v1_tmp[i] = v1_fac*tau2[i]-Rsym_v[i]*v2_fac;
    //     MatVec3f(Sym_1R_temp, tau_2, false, v1_tmp);

    // } // end of scope
