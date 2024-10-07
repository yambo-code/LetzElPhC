#include "../common/dtypes.h"
#include "../common/error.h"
#include "../common/numerical_func.h"
#include "../common/parallel.h"
#include "../common/progess_bar.h"
#include "../dvloc/dvloc.h"
#include "../elphC.h"
#include "../nonloc/Vnonloc.h"
#include "../symmetries/symmetries.h"
#include "elph.h"
#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void compute_and_write_elphq(struct WFC* wfcs, struct Lattice* lattice,
                             struct Pseudo* pseudo, struct Phonon* phonon,
                             const ND_int iqpt, ELPH_cmplx* eigVec,
                             ELPH_cmplx* dVscfq,
                             const int ncid_elph, const int varid_elph,
                             const int ncid_dmat, const int varid_dmat,
                             const bool non_loc,
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

    // sanity check
    if (phonon->qmap[2 * qpos] != iqpt || phonon->qmap[2 * qpos + 1] != 0)
    {
        error_msg("Qpoint in iBZ cannot be traced in phonon qmap");
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
    CHECK_ALLOC(KplusQidxs);

    ELPH_float qpt_tmp[3];
    memcpy(qpt_tmp, qpt, sizeof(ELPH_float) * 3);
    // incase if we want yambo convention
    if (kminusq)
    {
        for (int xi = 0; xi < 3; ++xi)
        {
            qpt_tmp[xi] = -qpt_tmp[xi];
        }
    }
    get_KplusQ_idxs(nk_totalBZ, lattice->kpt_fullBZ_crys, qpt_tmp, KplusQidxs);

    ND_int nbnds = lattice->nbnds;
    ND_int nmodes = lattice->nmodes;

    ELPH_cmplx* elph_kq_mn = NULL;
    if (Comm->commK_rank == 0)
    {
        elph_kq_mn = calloc(nbnds * nbnds * lattice->nspin * nmodes, sizeof(ELPH_cmplx));
        CHECK_ALLOC(elph_kq_mn);
    }
    //// (nu, nspin, mk, nk+q)
    /* Now Compute elph-matrix elements for each kpoint */

    size_t startp[7] = { 0, 0, 0, 0, 0, 0, 0 };
    size_t countp[7] = { 1, 1, nmodes, lattice->nspin, nbnds, nbnds, 2 };

    ELPH_cmplx* gSq_buff = NULL;
    ELPH_cmplx* D_mat_l = NULL;
    ELPH_cmplx* D_mat_r = NULL;

    if (Comm->commK_rank == 0)
    {
        gSq_buff = calloc(nbnds * nbnds * lattice->nspin * nmodes, sizeof(ELPH_cmplx));
        CHECK_ALLOC(gSq_buff);

        D_mat_l = calloc(nbnds * nbnds * lattice->nspin, sizeof(ELPH_cmplx));
        CHECK_ALLOC(D_mat_l);

        D_mat_r = calloc(nbnds * nbnds * lattice->nspin, sizeof(ELPH_cmplx));
        CHECK_ALLOC(D_mat_r);
    }

    // start the progress bar
    struct progress_bar pbar[1];
    start_progressbar(pbar, Comm->commW_rank, nk_this_pool);
    // compute electron-phonon matrix elements for each kpoint
    for (ND_int ii = 0; ii < nk_this_pool; ++ii)
    {
        /* compute the global k index */
        ND_int i = kshift + ii;

        int idx_k = i;
        int idx_kq = KplusQidxs[i];

        if (kminusq)
        {
            swap_ints(&idx_k, &idx_kq);
        }

        int ik = kmap[2 * idx_k];
        int ksym = kmap[2 * idx_k + 1];

        int ikq = kmap[2 * idx_kq];
        int kqsym = kmap[2 * idx_kq + 1];

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
            if ((nc_err = nc_put_vara(ncid_elph, varid_elph, startp, countp, elph_kq_mn)))
            {
                ERR(nc_err);
            }
        }

        // expand the el-ph matrix elements in full BZ
        const ELPH_float* kpt_BZ = lattice->kpt_fullBZ + 3 * i;

        // only master node does this
        if (Comm->commK_rank == 0)
        {
            for (ND_int istar = 1; istar < phonon->nqstar[iqpt]; ++istar)
            {
                ND_int qpos_star = qpos + istar;
                // get the symmetry
                int istar_symm = phonon->qmap[2 * qpos_star + 1];

                size_t D_mat_sp[6] = { istar_symm, idx_kq, 0, 0, 0, 0 };
                size_t D_mat_cp[6] = { 1, 1, lattice->nspin, lattice->nbnds, lattice->nbnds, 2 };
                // read D_mats //const int ncid_dmat, const int varid_dmat,
                if ((nc_err = nc_get_vara(ncid_dmat, varid_dmat, D_mat_sp, D_mat_cp, D_mat_l)))
                { // k + q
                    ERR(nc_err);
                }

                D_mat_sp[1] = idx_k;
                if ((nc_err = nc_get_vara(ncid_dmat, varid_dmat, D_mat_sp, D_mat_cp, D_mat_r)))
                { // k
                    ERR(nc_err);
                }

                // gSq_buff
                elph_q_rotate(D_mat_l, elph_kq_mn, D_mat_r, lattice,
                              phonon->ph_syms[istar_symm].time_rev, gSq_buff);

                ELPH_float Sk_crys[3];
                ELPH_float Sk_cart[3];

                MatVec3f(phonon->ph_syms[istar_symm].Rmat, kpt_BZ,
                         false, Sk_cart);
                // convert to crsytal units
                MatVec3f(lattice->alat_vec, Sk_cart, true, Sk_crys);

                // now get the indices of S*k

                ND_int idx_Sk = find_kidx_in_list(lattice->nkpts_BZ,
                                                  lattice->kpt_fullBZ_crys, Sk_crys);

                if (idx_Sk < 0)
                {
                    error_msg("Unable to to find S*k index "
                              "This is due to incommenserate k and q grids.");
                }

                startp[0] = qpos_star;
                startp[1] = idx_Sk;
                // Write it for Sq and Sk point
                if ((nc_err = nc_put_vara(ncid_elph, varid_elph, startp, countp, gSq_buff)))
                {
                    ERR(nc_err);
                }
            }
        }
        mpi_error = MPI_Barrier(Comm->commK);
        MPI_error_msg(mpi_error);
        // update the progress bar
        print_progressbar(pbar);
    }
    // free wfc buffers
    if (Comm->commK_rank == 0)
    {
        free(elph_kq_mn);
        free(gSq_buff);
        free(D_mat_l);
        free(D_mat_r);
    }

    free(KplusQidxs);
}
