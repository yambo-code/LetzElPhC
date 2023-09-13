#include "Vnonloc.h"


void elphNL_q(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo, \
            ELPH_float * qpt, ND_array(Nd_cmplxS) * eigVec, ELPH_cmplx * elph_kq,
            MPI_Comm commK, MPI_Comm commQ)
{
    /*
    COmpute non local part for all kpoints and modes for particular q point 

    eigVec : eigen vectors , (nu,atom,3)
    qpt is in crystal
    out : elph_kq ((k, nmodes, nspin, nbands, nbands))

    !! Warning . elph_kq must be initialized else Undefined behaviour !
    */

    int my_rank, comm_size, mpi_error;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int cpus_qpool, qcolor, qrank;

    MPI_Comm_rank(commQ, &qrank);
    MPI_Comm_size(commQ, &cpus_qpool);

    qcolor = my_rank/cpus_qpool ;
    int nqpools = comm_size/cpus_qpool;

    int npw_cpus, krank;
    MPI_Comm_rank(commK, &krank);
    MPI_Comm_size(commK, &npw_cpus);

    int kcolor = qrank/npw_cpus;
    int nkpools = cpus_qpool/npw_cpus ;


    ND_int nk_totalBZ = lattice->kmap->dims[0];

    int nk_per_color = nk_totalBZ/nkpools ;
    int k_rem      = nk_totalBZ%nkpools ;

    int k_in_this_color = nk_per_color;

    if (kcolor < k_rem) ++k_in_this_color;

    
    ND_int nbnd = wfcs->wfc->dims[1] ; 
    ND_int elph_kstride = (eigVec->dims[0]) * (lattice->nspin) *nbnd * nbnd; // 1st stride value of elph_kq

    int * kmap  = lattice->kmap->data ; 
    /* Compute elph-matrix elements*/
    int * KplusQidxs = malloc((lattice->kmap->dims[0])*sizeof(int));

    get_KplusQ_idxs(lattice->kpt_fullBZ_crys, KplusQidxs , qpt, lattice->alat_vec, true);

    for (ND_int ii =0 ; ii <k_in_this_color; ++ii )
    {   
        ND_int i  = kcolor*nk_per_color
        if (kcolor < k_rem) i += kcolor;
        else                i += k_rem;

        int ik    = *(kmap + i*2)      ;
        int ksym  = *(kmap + i*2 + 1)  ;
        int ikq   = *(kmap + KplusQidxs[i]*2)      ;
        int kqsym = *(kmap + KplusQidxs[i]*2 + 1)  ;
        ELPH_cmplx * elph_kq_mn = elph_kq + i*elph_kstride ;
        
        elphNonLocal(wfcs, lattice, pseudo, ikq, ik, kqsym, ksym, true, eigVec, elph_kq_mn, commK);
    }

    free(KplusQidxs);
}