/*
This function gets the dvscf and eig vectors for q point from qe
*/
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../common/cwalk/cwalk.h"
#include "../../common/dtypes.h"
#include "../../common/error.h"
#include "../../elphC.h"
#include "qe_io.h"

void get_dvscf_dyn_qe(const char* ph_save_dir, struct Lattice* lattice,
                      ND_int iq_BZ, ELPH_cmplx* eig, ELPH_cmplx* dvscf,
                      ELPH_float* omega_ph, const struct ELPH_MPI_Comms* Comm)
{
    int mpi_error;
    // ph_save_dir must be available on all cpus else U.B
    // if dvscf == NULL, then only eig vecs are returned.
    ND_int nmodes = lattice->nmodes;

    char small_buf[32];

    size_t tmp_char_buf_len = 1024 + strlen(ph_save_dir);

    char* tmp_char_buf = malloc(tmp_char_buf_len);
    CHECK_ALLOC(tmp_char_buf);

    ELPH_cmplx* pat_vecs = NULL;
    if (dvscf != NULL)
    {
        pat_vecs = malloc(sizeof(ELPH_cmplx) * nmodes * nmodes);
        CHECK_ALLOC(pat_vecs);
    }

    if (Comm->commQ_rank == 0)
    {
        if (dvscf != NULL)
        {
            snprintf(small_buf, 32, "patterns.%d.xml", (int)(iq_BZ + 1));
            cwk_path_join(ph_save_dir, small_buf, tmp_char_buf,
                          tmp_char_buf_len);

            read_pattern_qe(tmp_char_buf, lattice, pat_vecs);
        }
        ELPH_float qpts[3];

        snprintf(small_buf, 32, "dyn%d", (int)(iq_BZ + 1));
        cwk_path_join(ph_save_dir, small_buf, tmp_char_buf, tmp_char_buf_len);

        ND_int ndyn_read =
            read_dyn_qe(tmp_char_buf, lattice, qpts, omega_ph, eig);
        if (ndyn_read != 1)
        {
            error_msg(
                "Wrong number of dynmats read from dynamical matrix file");
        }
    }

    // Bcast variables

    mpi_error = MPI_Bcast(eig, nmodes * nmodes, ELPH_MPI_cmplx, 0, Comm->commQ);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(omega_ph, nmodes, ELPH_MPI_float, 0, Comm->commQ);
    MPI_error_msg(mpi_error);

    if (dvscf != NULL)
    {
        mpi_error = MPI_Bcast(pat_vecs, nmodes * nmodes, ELPH_MPI_cmplx, 0,
                              Comm->commQ);
        MPI_error_msg(mpi_error);

        snprintf(small_buf, 32, "dvscf%d", (int)(iq_BZ + 1));
        cwk_path_join(ph_save_dir, small_buf, tmp_char_buf, tmp_char_buf_len);

        if (Comm->commRq_rank == 0)
        {
            read_dvscf_qe(tmp_char_buf, lattice, eig, pat_vecs, dvscf,
                          Comm->commK);
        }
        // broad cast dvscf
        ND_int dvscf_len = nmodes * lattice->nmag * lattice->fft_dims[0] *
                           lattice->fft_dims[1] * lattice->nfftz_loc;
        ND_int max_int_val = ((ND_int)INT_MAX) - 10;

        // check for buffer overflow
        if (dvscf_len > max_int_val)
        {
            error_msg("Buffer Overflow in dvscf_len. Contact developer");
        }

        mpi_error =
            MPI_Bcast(dvscf, dvscf_len, ELPH_MPI_cmplx, 0, Comm->commRq);
        MPI_error_msg(mpi_error);

        free(pat_vecs);
    }

    free(tmp_char_buf);
}
