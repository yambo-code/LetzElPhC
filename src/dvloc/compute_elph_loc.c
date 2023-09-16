#include "dvloc.h"


void compute_elphLocal_q(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo, \
            ELPH_float * qpt, ND_array(Nd_cmplxS) * eigVec, ND_array(Nd_cmplxS) * dVscf_full, 
            ELPH_cmplx * elph_kq)
{
    /*
    dVscf -> (nmodes,nmag,Nx,Ny,Nz)
    ((k, nmodes, nspin, nbands, nbands))
    */

    /*
    First compute the 
    */
    
    int * kmap  = lattice->kmap->data ; 
    int * KplusQidxs = malloc((lattice->kmap->dims[0])*sizeof(int));
    get_KplusQ_idxs(lattice->kpt_fullBZ_crys, KplusQidxs , qpt, lattice->alat_vec, true);

    
    /* Compute elph-matrix elements*/
    for (ND_int i =0 ; i <lattice->kmap->dims[0]; ++i )
    {   
        int ik    = *(kmap + i*2)      ;
        int ksym  = *(kmap + i*2 + 1)  ;
        int ikq   = *(kmap + KplusQidxs[i]*2)      ;
        int kqsym = *(kmap + KplusQidxs[i]*2 + 1)  ;
        ELPH_cmplx * elph_kq_mn = elph_kq + i*elph_kstride_k + iv*elph_kstride_mode ;

        elphLocal(qpt, wfcs, lattice, ikq, ik, kqsym, ksym, true, &dVscf, elph_kq_mn);
    }
    
    free(KplusQidxs);

}

