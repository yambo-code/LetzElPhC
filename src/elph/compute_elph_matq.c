#include "elph.h"


void compute_elphq(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo, \
            ELPH_float * qpt, ND_array(Nd_cmplxS) * eigVec, ND_array(Nd_cmplxS) * dVscfq, 
            ELPH_cmplx * elph_kq, const struct ELPH_MPI_Comms * Comm)
{   

    /*
    dVscf -> (nmodes,nmag,Nx,Ny,Nz)
    ((k, nmodes, nspin, nbands, nbands))
    */

    /* distribute k points */
    int mpi_error;
    ND_int nk_totalBZ = lattice->kmap->dims[0];

    ND_int kshift, nk_this_pool;
    
    nk_this_pool = distribute_to_grps(nk_totalBZ, Comm->nkpools, \
                    Comm->commQ_rank/Comm->commK_size, &kshift);
    
    if (nk_this_pool < 1)
        error_msg("There are no kpoints in some cpus, Make sure nkpool < # of kpoints in full BZ.");

    /* Computing the  change in local potential */
    ND_array(Nd_cmplxS) Vlocr[1];
    // allocate memory for Vlocr
    ND_int dim_Buffer[4] = {eigVec->dims[0],dVscfq->dims[2], dVscfq->dims[3], dVscfq->dims[4] }; //(nmodes, Nx, Ny, Nz_loc)

    ND_function(init, Nd_cmplxS) (Vlocr, 4, dim_Buffer);
    ND_function(malloc, Nd_cmplxS)  (Vlocr);

    /* compute the local part of the bare */
    dVlocq(qpt, lattice, pseudo, eigVec, Vlocr, Comm->commK);
    /* add bare local to induce part*/
    add_dvscf(dVscfq, Vlocr); 
    /* now we can destroy Vlocr */
    ND_function(destroy, Nd_cmplxS) (Vlocr);

    /* get the k point indices and symmetries */
    int * kmap  = lattice->kmap->data ; 
    int * KplusQidxs = malloc((lattice->kmap->dims[0])*sizeof(int));
    
    get_KplusQ_idxs(lattice->kpt_fullBZ_crys, KplusQidxs , qpt, lattice->alat_vec, true);
    
    ND_int nbnd = wfcs->wfc->dims[1] ; 
    ND_int elph_kstride = (eigVec->dims[0]) * (lattice->nspin) *nbnd * nbnd; // 1st stride value of elph_kq
    
    /* Now Compute elph-matrix elements for each kpoint */
    for (ND_int ii =0 ; ii <nk_this_pool; ++ii )
    {   
        /* compute the global k index */
        ND_int i  = kshift + ii;
        int ik    = *(kmap + i*2)      ;
        int ksym  = *(kmap + i*2 + 1)  ;
        int ikq   = *(kmap + KplusQidxs[i]*2)      ;
        int kqsym = *(kmap + KplusQidxs[i]*2 + 1)  ;
        
        ELPH_cmplx * elph_kq_mn = elph_kq + ii*elph_kstride ;
        
        elphLocal(qpt, wfcs, lattice, ikq, ik, kqsym, ksym, dVscfq, Comm, elph_kq_mn);
        /* add the non local part elements */
        // WARNING : non local must be added after elphLocal else U.B
        add_elphNonLocal(wfcs, lattice, pseudo, ikq, ik, kqsym, ksym, eigVec, elph_kq_mn, Comm);
    }
    // free wfc buffers
    free(KplusQidxs);

}


