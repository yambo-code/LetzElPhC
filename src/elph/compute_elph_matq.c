#include "elph.h"


void compute_elph(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo, \
            ELPH_float * qpt, ND_array(Nd_cmplxS) * eigVec, ND_array(Nd_cmplxS) * dVscfq, 
            ELPH_cmplx * elph_kq, MPI_Comm commQ, MPI_Comm commK)
{
    /*
    dVscf -> (nmodes,nmag,Nx,Ny,Nz)
    ((k, nmodes, nspin, nbands, nbands))
    */

    /* distribute k points */
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

    /* Computing the  change in local potential */
    ND_array(Nd_cmplxS) Vlocr[1];
    // allocate memory for Vlocr
    ND_int dim_Buffer[2] = {eigVec->dims[0],lattice->nffts_loc}; //(nmodes, nffts_loc)

    ND_function(init, Nd_cmplxS) (Vlocr, 2, dim_Buffer);
    ND_function(malloc, Nd_cmplxS)  (Vlocr);

    /* compute the local part of the bare */
    dVlocq(qpt, lattice, pseudo, eigVec, Vlocr, commK);
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


    struct wfcBox wfcRspace[1];
    // allocate wfc buffers
    {
        const ND_int * fft_dims = lattice->fft_dims;
        const ND_int dimensions[6] = {lattice->nspin,nbnd,lattice->nspinor,fft_dims[0],fft_dims[1],fft_dims[2]};
        alloc_wfcBox(wfcRspace, dimensions, lattice->npw_max, FFTW_MEASURE, commK);
    }

    /* Now Compute elph-matrix elements for each kpoint */
    for (ND_int ii =0 ; ii <k_in_this_color; ++ii )
    {   
        /* compute the global k index */
        ND_int i  = kcolor*nk_per_color;
        if (kcolor < k_rem) i += kcolor;
        else                i += k_rem;
        i += ii;
        int ik    = *(kmap + i*2)      ;
        int ksym  = *(kmap + i*2 + 1)  ;
        int ikq   = *(kmap + KplusQidxs[i]*2)      ;
        int kqsym = *(kmap + KplusQidxs[i]*2 + 1)  ;
        
        ELPH_cmplx * elph_kq_mn = elph_kq + i*elph_kstride ;
        
        elphLocal(qpt, wfcs, lattice, ikq, ik, kqsym, ksym, dVscfq, commK, wfcRspace, elph_kq_mn);
        /* add the non local part elements */
        add_elphNonLocal(wfcs, lattice, pseudo, ikq, ik, kqsym, ksym, eigVec, elph_kq_mn, commK);
    }
    // free wfc buffers
    free_wfcBox(wfcRspace);
    free(KplusQidxs);

}

