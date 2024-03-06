#include "qe_io.h"

/*
File contains function to read dynamical matrix file (in the old format)
and outputs phonon polarization vectors
*/

void zheev_(char* jobz, char* uplo, int* n, double complex* a, int* lda,
            double* w, double complex* work, int* lwork, double* rwork,
            int* info);

// only one cpu calls the routine

#define DYN_READ_BUF_SIZE 10000
#define DYN_FLOAT_BUF_SIZE 200

ND_int read_dyn_qe(const char* dyn_file, struct Lattice* lattice,
                   ELPH_float* restrict qpts, ELPH_float* restrict omega,
                   ELPH_cmplx* restrict pol_vecs)
{
    /*
    // reads all the dynamical matrices in the file
    The function return value is number of dynmats found and read
    the qpoints and pol_vecs are updated according
    */

    ND_int nq_found = 0; // function return value

    char* fgets_err; // error code for fgets

    // First, open the file
    FILE* fp = fopen(dyn_file, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Opening file %s failed \n", dyn_file);
        error_msg("Unable to open the dynmat file");
    }

    char* read_buf = malloc(sizeof(char) * DYN_READ_BUF_SIZE);
    CHECK_ALLOC(read_buf);

    ELPH_float* read_fbuf = malloc(sizeof(ELPH_float) * DYN_FLOAT_BUF_SIZE);
    CHECK_ALLOC(read_fbuf);

    // two comment lines
    fgets(read_buf, DYN_READ_BUF_SIZE, fp);
    fgets(read_buf, DYN_READ_BUF_SIZE, fp);
    fgets(read_buf, DYN_READ_BUF_SIZE, fp); // this line has 9 floats
    if (parser_doubles_from_string(read_buf, read_fbuf) != 9)
    {
        error_msg("Error reading line 3 in dyn file");
    }
    ND_int ntype = rint(read_fbuf[0]);
    ND_int natom = rint(read_fbuf[1]);
    int nmodes = natom * 3;
    if (natom != lattice->natom)
    {
        error_msg("Wrong number of atoms in dyn file");
    }
    ND_int ibrav = rint(read_fbuf[2]);
    if (ibrav == 0)
    {
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // scratch
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // lat vec 1
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // lat vec 2
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // lat vec 3
    }

    ELPH_float* atm_mass = malloc(sizeof(ELPH_float) * (natom + ntype));
    CHECK_ALLOC(atm_mass);

    // read atomic mass for each type
    ELPH_float* atm_mass_type = atm_mass + natom;
    for (ND_int i = 0; i < ntype; ++i)
    {
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);
        int tmpi_read;
        double tmpf_read;
        char tmps_read[128];
        if (sscanf(read_buf, "%d '%[^']' %lf", &tmpi_read, tmps_read, &tmpf_read) != 3)
        {
            error_msg("Failed to read atomic masses from dyn file");
        }
        atm_mass_type[i] = tmpf_read;
        if (fabs(atm_mass_type[i]) < ELPH_EPS)
        {
            error_msg("Zero masses in dynamical file");
        }
    }
    // now set masses for each atom
    for (ND_int i = 0; i < natom; ++i)
    {
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);
        if (parser_doubles_from_string(read_buf, read_fbuf) != 5)
        {
            error_msg("Failed to read atomic masses from dyn file");
        }

        ND_int itype = rint(read_fbuf[1]);
        atm_mass[i] = atm_mass_type[itype - 1];
    }

    // allocate a tmp buffer for diagonalization
    double complex* dyn_mat_tmp = malloc(sizeof(double complex) * nmodes * nmodes);
    CHECK_ALLOC(dyn_mat_tmp);

    // allocate workspace for zheev
    double* omega2 = malloc(sizeof(double) * 4 * nmodes); // // (3N-2 + N) = 4N-2
    CHECK_ALLOC(omega2);

    double* rwork = omega2 + nmodes;

    int info_z, lwork;
    double complex tmp_work_var;
    lwork = -1; // set up a query request
    zheev_("V", "U", &nmodes, dyn_mat_tmp, &nmodes, omega2, &tmp_work_var,
           &lwork, rwork, &info_z);
    if (info_z != 0)
    {
        error_msg(
            "Error in query request for zheev, when diagonalizing dyn mat");
    }
    lwork = (int)rint(creal(tmp_work_var));

    double complex* work_array = malloc(sizeof(double complex) * lwork);
    CHECK_ALLOC(work_array);

    // start reading the dynamical matrices
    while (!nq_found) // while(true) // put while(true) to read all the
                      // dynamical matrices from file.
    {
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // empty
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // content strong
        if (!string_start_with(read_buf, "Dynamical", true))
        {
            break; // break the loop if dynmaical not found
        }
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // empty
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // qpoint

        ELPH_float* restrict qpt_tmp = qpts + nq_found * 3;

        if (parser_doubles_from_string(read_buf, qpt_tmp) != 3)
        {
            error_msg("error reading qpoint from dynamat files");
        }

        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // empty

        ELPH_cmplx* restrict eig = pol_vecs + nq_found * nmodes * nmodes;
        ELPH_float* restrict omega_q = omega + nq_found * nmodes;

        for (ND_int ia = 0; ia < natom; ++ia)
        {
            for (ND_int ib = 0; ib < natom; ++ib)
            {
                fgets(read_buf, DYN_READ_BUF_SIZE, fp); // read ia and ib
                int itmp, jtmp;
                sscanf(read_buf, "%d  %d", &itmp, &jtmp);
                --itmp;
                --jtmp;
                if (itmp != ia || jtmp != ib)
                {
                    error_msg(
                        "error reading dynamical matrix from dynamat files");
                }

                // we must divide by 1/sqrt(Ma*Mb)
                ELPH_float inv_mass_sqtr = sqrt(atm_mass[ia] * atm_mass[ib]);
                inv_mass_sqtr = 1.0 / inv_mass_sqtr;

                for (int ix = 0; ix < 3; ++ix)
                {
                    fgets(read_buf, DYN_READ_BUF_SIZE, fp); // read dynmat
                    // get the values of dynamical matrix
                    if (parser_doubles_from_string(read_buf, read_fbuf) != 6)
                    {
                        error_msg(
                            "error reading dynamical matrix from dynamat "
                            "files");
                    }
                    for (int i = 0; i < 6; ++i)
                    {
                        read_fbuf[i] *= inv_mass_sqtr;
                    }
                    for (int iy = 0; iy < 3; ++iy) // // (ib,iy,ia,ix)
                    {
                        // we store in column major format (this is because
                        // lapack routines are column major)
                        dyn_mat_tmp[ix + ia * 3 + iy * nmodes + ib * 3 * nmodes] = read_fbuf[2 * iy] + I * read_fbuf[2 * iy + 1];
                    }
                }
            }
        }

        // symmetrize the matrix
        for (ND_int idim1 = 0; idim1 < nmodes; ++idim1)
        {
            for (ND_int jdim1 = 0; jdim1 <= idim1; ++jdim1)
            {
                dyn_mat_tmp[idim1 * nmodes + jdim1] = 0.5 * (dyn_mat_tmp[idim1 * nmodes + jdim1] + conj(dyn_mat_tmp[jdim1 * nmodes + idim1]));
            }
        }

        // diagonalize the dynamical matrix and divide with sqrt of masses
        zheev_("V", "U", &nmodes, dyn_mat_tmp, &nmodes, omega2, work_array,
               &lwork, rwork, &info_z);
        if (info_z != 0)
        {
            error_msg("Error when diagonalizing dyn mat");
        }
        // now store them in eig and omega, (neig,natom,pol)
        for (ND_int imode = 0; imode < nmodes; ++imode)
        {
            omega_q[imode] = sqrt(fabs(omega2[imode]));
            if (omega2[imode] < 0)
            {
                omega_q[imode] = -omega_q[imode];
            }

            ELPH_cmplx* eig_tmp_ptr = eig + imode * nmodes;
            double complex* dyn_tmp_ptr = dyn_mat_tmp + imode * nmodes;
            for (ND_int jmode = 0; jmode < nmodes; ++jmode)
            {
                ND_int ia = jmode / 3;
                eig_tmp_ptr[jmode] = dyn_tmp_ptr[jmode] / sqrt(atm_mass[ia]);
            }
        }
        // update the counter
        ++nq_found;
    }

    free(work_array);
    free(omega2);
    free(dyn_mat_tmp);
    free(atm_mass);
    free(read_fbuf);
    free(read_buf);
    // close the file
    fclose(fp);

    if (nq_found == 0)
    {
        error_msg("No dynamical matrices found in the dyn file");
    }

    return nq_found;
}

void read_qpts_qe(const char* dyn0_file, ND_int* nqpt_iBZ, ND_int* nqpt_fullBZ,
                  ELPH_float** qpts)
{
    /*
    Reads q-points from dyn0 file. (in 2*pi/alat units)
    also gives number of qpoints in iBZ and full BZ
    memory for qpts is allocated inside and must be freed outside
    */
    FILE* fp = fopen(dyn0_file, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Opening file %s failed \n", dyn0_file);
        error_msg("Unable to open the dyn0 file");
    }

    char* read_buf = malloc(DYN_READ_BUF_SIZE);
    CHECK_ALLOC(read_buf);

    int qgrid[3];
    fgets(read_buf, DYN_READ_BUF_SIZE, fp);
    if (sscanf(read_buf, "%d %d %d", qgrid, qgrid + 1, qgrid + 2) != 3)
    {
        error_msg("Error reading qgrid from dyn0 file");
    }

    int nq_iBZ_tmp;
    fgets(read_buf, DYN_READ_BUF_SIZE, fp);
    if (sscanf(read_buf, "%d", &nq_iBZ_tmp) != 1)
    {
        error_msg("Error reading 2nd line from dyn0 file");
    }

    *nqpt_iBZ = nq_iBZ_tmp;
    *nqpt_fullBZ = qgrid[0] * qgrid[1] * qgrid[2];

    ELPH_float* iBZ_qpts = malloc(sizeof(ELPH_float) * 3 * nq_iBZ_tmp);
    CHECK_ALLOC(iBZ_qpts);

    *qpts = iBZ_qpts;

    // now read list of qpoints
    for (ND_int i = 0; i < *nqpt_iBZ; ++i)
    {
        float qpt_tmp[3];
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);
        if (sscanf(read_buf, "%f %f %f", qpt_tmp, qpt_tmp + 1, qpt_tmp + 2) != 3)
        {
            error_msg("Error reading qpoint from dyn0 file");
        }

        iBZ_qpts[3 * i] = qpt_tmp[0];
        iBZ_qpts[3 * i + 1] = qpt_tmp[1];
        iBZ_qpts[3 * i + 2] = qpt_tmp[2];
    }

    free(read_buf);
    fclose(fp);
}
