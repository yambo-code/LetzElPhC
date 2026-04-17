#include <complex.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/cblas.h"
#include "common/constants.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/string_func.h"
#include "elphC.h"
#include "io/ezxml/ezxml.h"
#include "qe_io.h"

/*
File contains function to read dynamical matrix file (in the old format)
and outputs phonon polarization vectors
*/

// only one cpu calls the routine

#define DYN_READ_BUF_SIZE 10000
#define DYN_FLOAT_BUF_SIZE 200

static ND_int read_dyn_qe_old(FILE* fp, struct Lattice* lattice,
                              ELPH_float* qpts, ELPH_float* omega,
                              ELPH_cmplx* pol_vecs, ELPH_float* amass);

static ND_int read_dyn_xml(FILE* fp, struct Lattice* lattice, ELPH_float* qpts,
                           ELPH_float* omega, ELPH_cmplx* pol_vecs,
                           ELPH_float* amass);

static bool is_dyn_xml(FILE* fp);

ND_int read_dyn_qe(const char* dyn_file, struct Lattice* lattice,
                   ELPH_float* qpts, ELPH_float* omega, ELPH_cmplx* pol_vecs,
                   ELPH_float* amass)
{
    // if amass (atomic masses in au) is NULL, will be ignored
    //
    // First, open the file
    FILE* fp = fopen(dyn_file, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Opening file %s failed \n", dyn_file);
        error_msg("Unable to open the dynmat file");
    }

    ND_int nq_found = 0;

    bool is_xml_format = is_dyn_xml(fp);
    // First check if it is xml file

    if (is_xml_format)
    {
        nq_found = read_dyn_xml(fp, lattice, qpts, omega, pol_vecs, amass);
    }
    else
    {
        nq_found = read_dyn_qe_old(fp, lattice, qpts, omega, pol_vecs, amass);
    }
    fclose(fp);
    return nq_found;
}

static ND_int read_dyn_qe_old(FILE* fp, struct Lattice* lattice,
                              ELPH_float* qpts, ELPH_float* omega,
                              ELPH_cmplx* pol_vecs, ELPH_float* amass)
{
    /*
    // reads all the dynamical matrices in the file
    The function return value is number of dynmats found and read
    the qpoints and pol_vecs are updated according
    */

    ND_int nq_found = 0;  // function return value

    char* read_buf = malloc(sizeof(char) * DYN_READ_BUF_SIZE);
    CHECK_ALLOC(read_buf);

    ELPH_float* read_fbuf = malloc(sizeof(ELPH_float) * DYN_FLOAT_BUF_SIZE);
    CHECK_ALLOC(read_fbuf);

    // two comment lines
    fgets(read_buf, DYN_READ_BUF_SIZE, fp);
    fgets(read_buf, DYN_READ_BUF_SIZE, fp);
    fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // this line has 9 floats
    if (parse_floats_from_string(read_buf, read_fbuf, DYN_FLOAT_BUF_SIZE) != 9)
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
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // scratch
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // lat vec 1
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // lat vec 2
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // lat vec 3
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
        if (sscanf(read_buf, "%d '%[^']' %lf", &tmpi_read, tmps_read,
                   &tmpf_read) != 3)
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
        if (parse_floats_from_string(read_buf, read_fbuf, DYN_FLOAT_BUF_SIZE) !=
            5)
        {
            error_msg("Failed to read atomic masses from dyn file");
        }

        ND_int itype = rint(read_fbuf[1]);
        atm_mass[i] = atm_mass_type[itype - 1];
    }

    // allocate a tmp buffer for diagonalization
    double _Complex* dyn_mat_tmp =
        malloc(sizeof(double _Complex) * nmodes * nmodes);
    CHECK_ALLOC(dyn_mat_tmp);

    // allocate workspace for zheev
    double* omega2 =
        malloc(sizeof(double) * 4 * nmodes);  // // (3N-2 + N) = 4N-2
    CHECK_ALLOC(omega2);

    double* rwork = omega2 + nmodes;

    int info_z, lwork;
    double _Complex tmp_work_var;
    lwork = -1;  // set up a query request
    zheev_("V", "U", &nmodes, dyn_mat_tmp, &nmodes, omega2, &tmp_work_var,
           &lwork, rwork, &info_z);
    if (info_z != 0)
    {
        error_msg(
            "Error in query request for zheev, when diagonalizing dyn mat");
    }
    lwork = (int)rint(creal(tmp_work_var * 1.005));

    double _Complex* work_array = malloc(sizeof(double _Complex) * lwork);
    CHECK_ALLOC(work_array);

    // start reading the dynamical matrices
    while (!nq_found)  // while(true) // put while(true) to read all the
                       // dynamical matrices from file.
    {
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // empty
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // content strong
        if (!string_start_with(read_buf, "Dynamical", true))
        {
            break;  // break the loop if dynmaical not found
        }
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // empty
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // qpoint

        ELPH_float* qpt_tmp = qpts + nq_found * 3;

        if (parse_floats_from_string(read_buf, qpt_tmp, 3) != 3)
        {
            error_msg("error reading qpoint from dynamat files");
        }

        fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // empty

        ELPH_cmplx* eig = pol_vecs + nq_found * nmodes * nmodes;
        ELPH_float* omega_q = omega + nq_found * nmodes;

        for (ND_int ia = 0; ia < natom; ++ia)
        {
            for (ND_int ib = 0; ib < natom; ++ib)
            {
                fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // read ia and ib
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
                    fgets(read_buf, DYN_READ_BUF_SIZE, fp);  // read dynmat
                    // get the values of dynamical matrix
                    if (parse_floats_from_string(read_buf, read_fbuf,
                                                 DYN_FLOAT_BUF_SIZE) != 6)
                    {
                        error_msg(
                            "error reading dynamical matrix from dynamat "
                            "files");
                    }
                    for (int i = 0; i < 6; ++i)
                    {
                        read_fbuf[i] *= inv_mass_sqtr;
                    }
                    for (int iy = 0; iy < 3; ++iy)  // // (ib,iy,ia,ix)
                    {
                        // we store in column major format (this is because
                        // lapack routines are column major)
                        dyn_mat_tmp[ix + ia * 3 + iy * nmodes +
                                    ib * 3 * nmodes] =
                            read_fbuf[2 * iy] + I * read_fbuf[2 * iy + 1];
                    }
                }
            }
        }

        // symmetrize the matrix
        for (ND_int idim1 = 0; idim1 < nmodes; ++idim1)
        {
            for (ND_int jdim1 = 0; jdim1 <= idim1; ++jdim1)
            {
                dyn_mat_tmp[idim1 * nmodes + jdim1] =
                    0.5 * (dyn_mat_tmp[idim1 * nmodes + jdim1] +
                           conj(dyn_mat_tmp[jdim1 * nmodes + idim1]));
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
            double _Complex* dyn_tmp_ptr = dyn_mat_tmp + imode * nmodes;
            for (ND_int jmode = 0; jmode < nmodes; ++jmode)
            {
                ND_int ia = jmode / 3;
                eig_tmp_ptr[jmode] = dyn_tmp_ptr[jmode] / sqrt(atm_mass[ia]);
            }
        }
        // update the counter
        ++nq_found;
    }

    if (amass)
    {
        memcpy(amass, atm_mass, natom * sizeof(*amass));
    }

    free(work_array);
    free(omega2);
    free(dyn_mat_tmp);
    free(atm_mass);
    free(read_fbuf);
    free(read_buf);
    // close the file

    if (nq_found == 0)
    {
        error_msg("No dynamical matrices found in the dyn file");
    }

    return nq_found;
}

static ND_int read_dyn_xml(FILE* fp, struct Lattice* lattice, ELPH_float* qpts,
                           ELPH_float* omega, ELPH_cmplx* pol_vecs,
                           ELPH_float* amass)
{
    /*
    Reads dynamical matrices from QE XML format file.
    Returns number of q-points found and read.
    Updates qpts, omega, and pol_vecs arrays.
    */

    ND_int nq_found = 0;

    ezxml_t xml = ezxml_parse_fp(fp);
    if (xml == NULL)
    {
        error_msg("Error parsing dynamical xml file");
    }

    ezxml_t geom_info = ezxml_get(xml, "GEOMETRY_INFO", -1);
    if (geom_info == NULL)
    {
        error_msg("No GEOMETRY_INFO tag found in XML file");
    }

    // Read basic parameters
    ezxml_t xml_ntype = ezxml_get(geom_info, "NUMBER_OF_TYPES", -1);
    if (!xml_ntype)
    {
        error_msg("Error parsing number of types");
    }
    int ntypes = atoi(xml_ntype->txt);

    ezxml_t xml_natoms = ezxml_get(geom_info, "NUMBER_OF_ATOMS", -1);
    if (!xml_natoms)
    {
        error_msg("Error parsing number of atoms");
    }
    int natoms = atoi(xml_natoms->txt);

    int nmodes = natoms * 3;

    if (natoms != lattice->natom)
    {
        error_msg("Wrong number of atoms in dyn file");
    }

    // Read atomic masses
    ELPH_float* atm_mass = malloc(sizeof(ELPH_float) * (natoms + ntypes));
    CHECK_ALLOC(atm_mass);

    ELPH_float* atm_mass_type = atm_mass + natoms;
    for (int it = 0; it < ntypes; it++)
    {
        char tag[50];
        snprintf(tag, sizeof(tag), "MASS.%d", (int)(it + 1));

        ezxml_t mass_tag = ezxml_get(geom_info, tag, -1);
        if (!mass_tag)
        {
            error_msg("Cannot parse atomic mass for dyn.xml files");
        }
        // in xml files. it is in amu
        atm_mass_type[it] = atof(mass_tag->txt) * 911.444243096;
    }
    // Find all atoms of this type and set their masses
    for (int ia = 0; ia < natoms; ia++)
    {
        char atom_tag[50];
        snprintf(atom_tag, sizeof(atom_tag), "ATOM.%d", (int)(ia + 1));
        ezxml_t atom = ezxml_get(geom_info, atom_tag, -1);
        if (atom)
        {
            const char* index = ezxml_attr(atom, "INDEX");
            if (index)
            {
                int ia_idx = atoi(index) - 1;
                atm_mass[ia] = atm_mass_type[ia_idx];
            }
            else
            {
                error_msg(
                    "Cannot parse atomic index from atom for dyn.xml files");
            }
        }
        else
        {
            error_msg("Cannot parse atomic mass for dyn.xml files");
        }
    }

    // Allocate temporary buffers
    double _Complex* dyn_mat_tmp =
        malloc(sizeof(double _Complex) * nmodes * nmodes);
    CHECK_ALLOC(dyn_mat_tmp);

    double* omega2 = malloc(sizeof(double) * 4 * nmodes);  // (3N-2 + N) = 4N-2
    CHECK_ALLOC(omega2);

    double* rwork = omega2 + nmodes;

    // Prepare for diagonalization
    int info_z, lwork;
    double _Complex tmp_work_var;
    lwork = -1;  // set up a query request
    zheev_("V", "U", &nmodes, dyn_mat_tmp, &nmodes, omega2, &tmp_work_var,
           &lwork, rwork, &info_z);
    if (info_z != 0)
    {
        error_msg("Error in query request for zheev");
    }
    lwork = (int)rint(creal(tmp_work_var * 1.005));

    double _Complex* work_array = malloc(sizeof(double _Complex) * lwork);
    CHECK_ALLOC(work_array);

    // Read dynamical matrices for each q-point
    while (!nq_found)
    {
        // we only need one dynamical matrix.
        char dyn_tag[50];
        snprintf(dyn_tag, sizeof(dyn_tag), "DYNAMICAL_MAT_.%d",
                 (int)(nq_found + 1));
        ezxml_t dyn_mat = ezxml_get(xml, dyn_tag, -1);
        if (dyn_mat == NULL)
        {
            break;  // No more q-points
        }

        ELPH_float* qpt_tmp = qpts + nq_found * 3;
        ELPH_cmplx* eig = pol_vecs + nq_found * nmodes * nmodes;
        ELPH_float* omega_q = omega + nq_found * nmodes;

        // Read q-point
        ezxml_t qpt_xml = ezxml_get(dyn_mat, "Q_POINT", -1);
        if (!qpt_xml)
        {
            error_msg("Error reading q-point from XML");
        }
        const char* q_str = qpt_xml->txt;
        if (parse_floats_from_string(q_str, qpt_tmp, 3) != 3)
        {
            error_msg("Error reading q-point from XML");
        }

        // Read dynamical matrix
        for (int ia = 0; ia < natoms; ++ia)
        {
            for (int ib = 0; ib < natoms; ++ib)
            {
                char phi_tag[50];
                snprintf(phi_tag, sizeof(phi_tag), "PHI.%d.%d", (int)(ia + 1),
                         (int)(ib + 1));
                ezxml_t dynr_xml = ezxml_get(dyn_mat, phi_tag, -1);
                if (!dynr_xml)
                {
                    error_msg("Error parsing dynamical matrix from XML");
                }
                const char* phi_str = dynr_xml->txt;

                ELPH_float phi_vals[18];
                if (parse_floats_from_string(phi_str, phi_vals,
                                             ARRAY_LEN(phi_vals)) != 18)
                {
                    error_msg("Error reading dynamical matrix from XML");
                }

                // Fill dynamical matrix (include mass factors)
                ELPH_float inv_mass_sqtr = sqrt(atm_mass[ia] * atm_mass[ib]);
                inv_mass_sqtr = 1.0 / inv_mass_sqtr;

                for (int ix = 0; ix < 3; ix++)
                {
                    for (int iy = 0; iy < 3; iy++)
                    {
                        // Column-major storage for LAPACK
                        dyn_mat_tmp[(ix + ia * 3) + (iy + ib * 3) * nmodes] =
                            (phi_vals[2 * ix + 6 * iy] +
                             I * phi_vals[2 * ix + 6 * iy + 1]) *
                            inv_mass_sqtr;
                    }
                }
            }
        }

        // Symmetrize the matrix
        for (ND_int idim1 = 0; idim1 < nmodes; idim1++)
        {
            for (ND_int jdim1 = 0; jdim1 <= idim1; jdim1++)
            {
                dyn_mat_tmp[idim1 * nmodes + jdim1] =
                    0.5 * (dyn_mat_tmp[idim1 * nmodes + jdim1] +
                           conj(dyn_mat_tmp[jdim1 * nmodes + idim1]));
            }
        }

        // Diagonalize the dynamical matrix
        zheev_("V", "U", &nmodes, dyn_mat_tmp, &nmodes, omega2, work_array,
               &lwork, rwork, &info_z);
        if (info_z != 0)
        {
            error_msg("Error diagonalizing dynamical matrix");
        }

        // Store eigenvalues and eigenvectors
        for (ND_int imode = 0; imode < nmodes; imode++)
        {
            omega_q[imode] = sqrt(fabs(omega2[imode]));
            if (omega2[imode] < 0)
            {
                omega_q[imode] = -omega_q[imode];
            }

            ELPH_cmplx* eig_tmp_ptr = eig + imode * nmodes;
            double _Complex* dyn_tmp_ptr = dyn_mat_tmp + imode * nmodes;
            for (ND_int jmode = 0; jmode < nmodes; jmode++)
            {
                ND_int ia = jmode / 3;
                eig_tmp_ptr[jmode] = dyn_tmp_ptr[jmode] / sqrt(atm_mass[ia]);
            }
        }

        nq_found++;
    }

    // copy atomic masses
    if (amass)
    {
        memcpy(amass, atm_mass, natoms * sizeof(*amass));
    }

    // Clean up
    free(atm_mass);
    free(dyn_mat_tmp);
    free(omega2);
    free(work_array);
    ezxml_free(xml);

    if (nq_found == 0)
    {
        error_msg("No dynamical matrices found in the XML file");
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

    long long int qgrid[3];
    fgets(read_buf, DYN_READ_BUF_SIZE, fp);
    if (sscanf(read_buf, "%lld %lld %lld", qgrid, qgrid + 1, qgrid + 2) != 3)
    {
        error_msg("Error reading qgrid from dyn0 file");
    }

    long long int nq_iBZ_tmp;
    fgets(read_buf, DYN_READ_BUF_SIZE, fp);
    if (sscanf(read_buf, "%lld", &nq_iBZ_tmp) != 1)
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
        double qpt_tmp[3];
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);
        if (sscanf(read_buf, "%lf %lf %lf", qpt_tmp, qpt_tmp + 1,
                   qpt_tmp + 2) != 3)
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

static bool is_dyn_xml(FILE* fp)
{
    // must be rewind(fp) at the end.
    // the if(true) is to create a small scope
    bool is_xml_format = false;
    char line[1024];
    while (fgets(line, sizeof(line), fp))
    {
        // Skip leading whitespace
        char* p = line;
        while (isspace((unsigned char)(*p)))
        {
            p++;
        }

        // Skip empty lines
        if (*p == '\0')
        {
            continue;
        }

        // Check for '<?xml'
        if (strstr(p, "<?xml") != NULL)
        {
            is_xml_format = true;
        }
        break;
    }
    rewind(fp);
    return is_xml_format;
}
