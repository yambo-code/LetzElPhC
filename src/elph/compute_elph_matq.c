#include "elph.h"

void compute_and_write_elphq(struct WFC* wfcs, struct Lattice* lattice,
                             struct Pseudo* pseudo, struct Phonon* phonon,
                             const ND_int iqpt, ELPH_cmplx* eigVec,
                             ELPH_cmplx* dVscfq, const int ncid,
                             const int varid, const bool non_loc,
                             const bool kminusq,
                             const struct ELPH_MPI_Comms* Comm)
{
    /*
    dVscf -> (nmodes,nmag,Nx,Ny,Nz)
    ((k, nmodes, nspin, nbands, nbands))
    */
    /* distribute k points */
    int mpi_error;
    const ND_int nk_totalBZ = lattice->nkpts_BZ;

    const ELPH_float* qpt = phonon->qpts_iBZ + iqpt * 3;

    ND_int qpos = 0; // positon of this iBZ qpoint in full q point list
    for (ND_int i = 0; i < iqpt; ++i)
    {
        qpos += phonon->nqstar[i];
    }

    ND_int kshift;
    ND_int nk_this_pool = distribute_to_grps(nk_totalBZ, Comm->nkpools,
                                             Comm->commQ_rank / Comm->commK_size, &kshift);

    if (nk_this_pool < 1)
    {
        error_msg(
            "There are no kpoints in some cpus, Make sure nkpool < # of "
            "kpoints in full BZ.");
    }
    /* get the k point indices and symmetries */
    int* kmap = lattice->kmap;
    int* KplusQidxs = malloc(nk_totalBZ * sizeof(int));

    get_KplusQ_idxs(nk_totalBZ, lattice->kpt_fullBZ_crys, qpt, KplusQidxs);

    ND_int nbnds = lattice->nbnds;
    ND_int nmodes = lattice->natom * 3;

    ELPH_cmplx* elph_kq_mn = NULL;
    if (Comm->commK_rank == 0)
    {
        elph_kq_mn = calloc(nbnds * nbnds * lattice->nspin * nmodes, sizeof(ELPH_cmplx));
    }
    //// (nu, nspin, mk, nk+q)
    /* Now Compute elph-matrix elements for each kpoint */

    size_t startp[7] = { 0, 0, 0, 0, 0, 0, 0 };
    size_t countp[7] = { 1, 1, nmodes, lattice->nspin, nbnds, nbnds, 2 };

    for (ND_int ii = 0; ii < nk_this_pool; ++ii)
    {
        /* compute the global k index */
        ND_int i = kshift + ii;
        int ik = *(kmap + i * 2);
        int ksym = *(kmap + i * 2 + 1);
        int ikq = *(kmap + KplusQidxs[i] * 2);
        int kqsym = *(kmap + KplusQidxs[i] * 2 + 1);

        elphLocal(qpt, wfcs, lattice, ikq, ik, kqsym, ksym, dVscfq, Comm,
                  elph_kq_mn);
        /* add the non local part elements */
        // WARNING : non local must be added after elphLocal else U.B
        if (non_loc)
        {
            add_elphNonLocal(wfcs, lattice, pseudo, ikq, ik, kqsym, ksym,
                             eigVec, elph_kq_mn, Comm);
        }
        startp[0] = qpos;
        startp[1] = i;
        int nc_err;
        if (Comm->commK_rank == 0)
        {
            if ((nc_err = nc_put_vara(ncid, varid, startp, countp, elph_kq_mn)))
            {
                ERR(nc_err);
            }
        }
    }
    // free wfc buffers
    if (Comm->commK_rank == 0)
    {
        free(elph_kq_mn);
    }

    free(KplusQidxs);
}
