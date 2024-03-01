#include <limits.h>

#include "qe_io.h"
/*
This file contains routines to read dvscf/drho from qe output files
*/

void read_dvscf_qe(const char* dvscf_file, struct Lattice* lattice,
                   const ELPH_cmplx* eig, const ELPH_cmplx* pats,
                   ELPH_cmplx* restrict dvscf_out, MPI_Comm commK)
{
    /*
    takes input dvscf file name, patterns and phonon polarization vectors
    output : dvscf in mode basis (nmodes, nmag, FFTx, FFTy, FFTz)
    *
    eig : (nmodes, natom, pol)
    pattern : (nmodes, natom, pol)
    */
    // double complex // nmag

    // allocate a temporary buffer to read
    // qe store stuff in doubles as an array of shape (nmodes, nmag, FFTz, FFTy,
    // FFTx)
    int my_rank, Comm_size, mpi_error;

    mpi_error = MPI_Comm_size(commK, &Comm_size);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_rank(commK, &my_rank);
    MPI_error_msg(mpi_error);

    const ND_int* fft_dims = lattice->fft_dims;
    const ND_int nmodes = lattice->nmodes;

    ND_int read_buffer_count = fft_dims[0] * fft_dims[1] * lattice->nfftz_loc;
    // check for buffer overflow of read_buffer_count. This is because mpi takes
    // parameter as int
    if (read_buffer_count >= ((ND_int)INT_MAX))
    {
        error_msg("Buffer overflow in mpi count argument");
    }

    double complex* read_buf = malloc(sizeof(double complex) * read_buffer_count);
    CHECK_ALLOC(read_buf);

    MPI_File handle;

    if (MPI_File_open(commK, dvscf_file, MPI_MODE_RDONLY, MPI_INFO_NULL,
                      &handle)
        != MPI_SUCCESS)
    {
        if (my_rank == 0)
        {
            fprintf(stderr, "Unable to find the %s file \n", dvscf_file);
        }
        error_msg("unable to open dvscf file");
    }

    // set the file pointer to 0
    // mpi_error = MPI_File_seek(handle,0,MPI_SEEK_SET);

    // now read stuff
    ND_int nsets = nmodes * lattice->nmag;
    for (ND_int iset = 0; iset < nsets; ++iset)
    {
        ELPH_cmplx* dvscf_tmp = dvscf_out + iset * read_buffer_count;

        MPI_Offset offset = sizeof(double complex) * fft_dims[0] * fft_dims[1] * (iset * fft_dims[2] + lattice->nfftz_loc_shift);

        mpi_error = MPI_File_read_at_all(handle, offset, read_buf, read_buffer_count,
                             MPI_C_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        MPI_error_msg(mpi_error);
        // transpose from (FFTz, FFTy, FFTx)->(FFTx, FFTy, FFTz)
        ND_int ld_read_buf = fft_dims[0] * fft_dims[1];
        ND_int ld_dvscf_buf = lattice->nfftz_loc * fft_dims[1];

        for (ND_int iy = 0; iy < fft_dims[1]; ++iy)
        {
            double complex* read_buf_y_ptr = read_buf + iy * fft_dims[0];
            ELPH_cmplx* dvscf_y_ptr = dvscf_tmp + iy * lattice->nfftz_loc;

            for (ND_int ix = 0; ix < fft_dims[0]; ++ix)
            {
                for (ND_int iz = 0; iz < lattice->nfftz_loc; ++iz)
                {
                    dvscf_y_ptr[ix * ld_dvscf_buf + iz] = read_buf_y_ptr[iz * ld_read_buf + ix];
                }
            }
        }
    }

    // free the read buffer
    free(read_buf);

    // the read buffer is in pattern basis. convert to cart basis then to mode
    // mode basis

    /*
    conversion to cart einsum('vsijk, vax-> sijkax', dVscf,np.conj(pattern_vec))
    i.e pattern_vec^\dagger@dVscf
    // cart to mode basis
    => Output = eig@pattern_vec^\dagger@dVscf

    So, we first compute eig@pattern_vec^\dagger and multiply to the read array
    */

    // allocate buffer for eig@pattern_vec^\dagger

    ELPH_cmplx* eig_pat = calloc(nmodes * nmodes, sizeof(ELPH_cmplx));
    // const ELPH_cmplx * eig, const ELPH_cmplx * pats
    CHECK_ALLOC(eig_pat);

    // compute eig@pattern_vec^\dagger
    matmul_cmplx('N', 'C', eig, pats, eig_pat, 1.0, 0.0, nmodes, nmodes, nmodes,
                 nmodes, nmodes, nmodes);

    // Now convert to mode basis

    // we do not want to create a another large buffer, so we create a
    // relatively smaller buffer and do multiple matmuls (nmodes, FFTy,
    // FFTz_loc)
    ELPH_cmplx* dvscf_mode = calloc(nmodes * fft_dims[1] * lattice->nfftz_loc, sizeof(ELPH_cmplx));
    CHECK_ALLOC(dvscf_mode);

    ND_int niter = lattice->nmag * fft_dims[0]; // FFTx*nmag
    ND_int ld_mode = niter * fft_dims[1] * lattice->nfftz_loc; // nmag*FFTx*FFTy*FFTz_loc
    for (ND_int iter = 0; iter < niter; ++iter)
    {
        ELPH_cmplx* dvscf_pat_tmp = dvscf_out + iter * fft_dims[1] * lattice->nfftz_loc;

        matmul_cmplx('N', 'N', eig_pat, dvscf_pat_tmp, dvscf_mode, 1.0, 0.0,
                     nmodes, ld_mode, fft_dims[1] * lattice->nfftz_loc, nmodes,
                     fft_dims[1] * lattice->nfftz_loc, nmodes);

        // copy back the output to input buf
        for (ND_int imode = 0; imode < nmodes; ++imode)
        {
            ELPH_cmplx* tmp_copy = dvscf_pat_tmp + imode * ld_mode;
            memcpy(tmp_copy,
                   dvscf_mode + imode * fft_dims[1] * lattice->nfftz_loc,
                   fft_dims[1] * lattice->nfftz_loc * sizeof(ELPH_cmplx));
        }
    }
    // free space
    free(dvscf_mode);
    free(eig_pat);

    // close the file
    if (MPI_File_close(&handle) != MPI_SUCCESS)
    {
        if (my_rank == 0)
        {
            fprintf(stderr, "Unable to close %s file \n", dvscf_file);
        }
        error_msg("unable to close dvscf file");
    }
}
