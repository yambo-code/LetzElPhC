// THis file contains write functions
// to read fake dyn files and dvscfs.
// The file can only be read by Letzelphc
//
#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/error.h"
#include "elphC.h"
#include "qe_io.h"

void write_dyn_qe(const char* file_name, ND_int natom, const ELPH_float* qpts,
                  const ELPH_cmplx* dyn_mat, const ELPH_float* atomic_masses)
{
    //
    FILE* fp = fopen(file_name, "w");
    if (!fp)
    {
        error_msg("Error creating dyn file");
    }

    // Write comment lines
    fprintf(fp, "Dynamical matrix file\n");
    fprintf(fp, "%s\n",
            "This is fake dyn file created by LetzElphC, only to be read by "
            "LetzElphC.");

    ND_int ntype = natom;  // set it same as natoms
    ND_int ibrav = -1;     // set to negative (should not be 0)
    //
    //  "%3lld %5lld %4lld %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f",
    fprintf(fp, "%3lld%5lld%4lld%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n",
            (long long)ntype, (long long)natom, (long long)ibrav, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0);

    // Write atom types and masses
    for (ND_int ia = 0; ia < natom; ia++)
    {
        fprintf(fp, " %4lld 'X%lld' %14.6f\n", (long long)(ia + 1),
                (long long)(ia + 1), atomic_masses[ia]);
    }

    // Write atom positions
    for (ND_int ia = 0; ia < natom; ia++)
    {
        fprintf(fp, "%5lld%5lld%18.10f%18.10f%18.10f\n", (long long)(ia + 1),
                (long long)(ia + 1), 0.0, 0.0, 0.0);
    }

    // For each q-point
    ND_int nq = 1;
    for (ND_int iq = 0; iq < nq; iq++)
    {
        const ELPH_float* qpt = qpts + iq * 3;
        const ELPH_cmplx* dyn_mat_q = dyn_mat + iq * 9 * natom * natom;

        // Write q-poND_int header
        fprintf(fp, "\n     Dynamical  Matrix in cartesian axes\n\n");
        fprintf(fp, "     q = ( %14.9f%14.9f%14.9f ) \n\n", qpt[0], qpt[1],
                qpt[2]);

        // Write dynamical matrix blocks
        for (ND_int ia = 0; ia < natom; ia++)
        {
            for (ND_int ib = 0; ib < natom; ib++)
            {
                fprintf(fp, "%5lld%5lld\n", (long long)(ia + 1),
                        (long long)(ib + 1));

                ELPH_float mass_sqrt =
                    sqrt(atomic_masses[ia] * atomic_masses[ib]);
                for (ND_int ix = 0; ix < 3; ix++)
                {
                    for (ND_int iy = 0; iy < 3; iy++)
                    {
                        // our dynamical matrix is in row major
                        // but in dynma we store in coloumn major.
                        ND_int idx =
                            iy + ib * 3 + ix * 3 * natom + ia * 3 * 3 * natom;
                        fprintf(fp, "%12.8f %12.8f   ",
                                creal(dyn_mat_q[idx]) * mass_sqrt,
                                cimag(dyn_mat_q[idx]) * mass_sqrt);
                    }
                    fprintf(fp, "\n");
                }
            }
        }
    }

    fprintf(fp, "\n");
    fclose(fp);
}

void write_dvscf_qe(const char* dvscf_file, struct Lattice* lattice,
                    const ELPH_cmplx* dvscf_in, MPI_Comm commK)
{
    /*
    Takes dvscf in mode basis (nmodes, nmag, FFTx, FFTy, FFTz_loc)
    and writes it to file in QE format (nmodes, nmag, FFTz, FFTy, FFTx)
    where FFTz is the global dimension
    */
    int my_rank, Comm_size, mpi_error;

    mpi_error = MPI_Comm_size(commK, &Comm_size);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_rank(commK, &my_rank);
    MPI_error_msg(mpi_error);

    const ND_int* fft_dims = lattice->fft_dims;
    const ND_int nmodes = lattice->nmodes;

    ND_int write_buffer_count = fft_dims[0] * fft_dims[1] * lattice->nfftz_loc;
    // check for buffer overflow of write_buffer_count
    if (write_buffer_count >= ((ND_int)INT_MAX))
    {
        error_msg("Buffer overflow in mpi count argument");
    }

    double complex* write_buf =
        malloc(sizeof(double complex) * write_buffer_count);
    CHECK_ALLOC(write_buf);

    MPI_File handle;

    if (MPI_File_open(commK, dvscf_file, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                      MPI_INFO_NULL, &handle) != MPI_SUCCESS)
    {
        error_msg("unable to open dvscf file for writing");
    }

    // Set the file size (important for MPI-IO)
    MPI_Offset file_size = sizeof(double complex) * nmodes * lattice->nmag *
                           fft_dims[0] * fft_dims[1] * fft_dims[2];
    mpi_error = MPI_File_set_size(handle, file_size);
    MPI_error_msg(mpi_error);

    // Now write stuff
    ND_int nsets = nmodes * lattice->nmag;
    for (ND_int iset = 0; iset < nsets; ++iset)
    {
        const ELPH_cmplx* dvscf_tmp = dvscf_in + iset * write_buffer_count;

        // Transpose from (FFTx, FFTy, FFTz_loc)->(FFTz, FFTy, FFTx) in
        // write_buf
        ND_int ld_dvscf_buf = lattice->nfftz_loc * fft_dims[1];
        ND_int ld_write_buf = fft_dims[0] * fft_dims[1];

        for (ND_int iy = 0; iy < fft_dims[1]; ++iy)
        {
            double complex* write_buf_y_ptr = write_buf + iy * fft_dims[0];
            const ELPH_cmplx* dvscf_y_ptr = dvscf_tmp + iy * lattice->nfftz_loc;

            for (ND_int ix = 0; ix < fft_dims[0]; ++ix)
            {
                for (ND_int iz = 0; iz < lattice->nfftz_loc; ++iz)
                {
                    write_buf_y_ptr[iz * ld_write_buf + ix] =
                        dvscf_y_ptr[ix * ld_dvscf_buf + iz];
                }
            }
        }

        MPI_Offset offset = sizeof(double complex) * fft_dims[0] * fft_dims[1] *
                            (iset * fft_dims[2] + lattice->nfftz_loc_shift);

        mpi_error =
            MPI_File_write_at_all(handle, offset, write_buf, write_buffer_count,
                                  MPI_C_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
        MPI_error_msg(mpi_error);
    }

    // Free the write buffer
    free(write_buf);
    // close the file
    if (MPI_File_close(&handle) != MPI_SUCCESS)
    {
        error_msg("unable to close dvscf file");
    }
}

void write_qpts_qe(const char* dyn_file, ND_int nqpt_iBZ,
                   const ELPH_float* qpts, const ND_int* qgrid)
{
    /*
    Writes q-points to dyn0 file in the format readable by read_qpts_qe
    Input:
        dyn_file - output file name
        nqpt_iBZ - number of q-points in irreducible BZ
        qpts - array of q-points in iBZ (size 3*nqpt_iBZ)
        qgrid - [3] array specifying the q-point mesh
    */
    FILE* fp = fopen(dyn_file, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "Opening file %s for writing failed \n", dyn_file);
        error_msg("Unable to open the dyn0 file for writing");
    }

    // Write qgrid
    fprintf(fp, "%lld %lld %lld\n", (long long)qgrid[0], (long long)qgrid[1],
            (long long)qgrid[2]);

    // Write number of q-points in iBZ
    fprintf(fp, "%lld\n", (long long)nqpt_iBZ);

    // Write q-points
    for (ND_int i = 0; i < nqpt_iBZ; ++i)
    {
        fprintf(fp, "%.10f %.10f %.10f\n", qpts[3 * i], qpts[3 * i + 1],
                qpts[3 * i + 2]);
    }

    fclose(fp);
}
