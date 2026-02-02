/*
This file contains function that parses data-file-schema.xml file
*/

#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "common/constants.h"
#include "common/error.h"
#include "common/string_func.h"
#include "elphC.h"
#include "io/ezxml/ezxml.h"
#include "qe_io.h"

#define ELPH_XML_READ_LINE_SIZE 1000

void parse_qexml(const char* xml_file, ND_int* natoms, ELPH_float* lat_vec,
                 ELPH_float* alat, char* dim, bool* is_soc_present,
                 ND_int* nmag, ND_int* fft_dims, ND_int* nph_sym,
                 ELPH_float** ph_sym_mats, ELPH_float** ph_sym_tau,
                 bool** ph_trevs, bool* ph_mag_symm, bool* ph_tim_rev,
                 char** pseudo_dir, char*** pseudo_pots, int* nspinor,
                 ND_int* ntype, int** atomic_type, ELPH_float** atomic_positons)
{
    /*
    get pseudo potential information, fft information from xml file

    need to free pseudo_dir and pseudo_pots pointers
    */

    // set defaults
    *is_soc_present = false;
    *dim = '3';
    *nmag = 1;
    if (nspinor)
    {
        *nspinor = 1;
    }

    // check if there is any assume_isolated tag in xml file
    FILE* fp = fopen(xml_file, "r");
    if (!fp)
    {
        error_msg("Error opening data-file-schema.xml file");
    }

    ezxml_t qexml = ezxml_parse_fp(fp);
    if (qexml == NULL)
    {
        error_msg("Error parsing data-file-schema.xml file");
    }

    char* tmp_read = malloc(ELPH_XML_READ_LINE_SIZE);
    CHECK_ALLOC(tmp_read);

    const char* tmp_str;  // tmp variable for ezxml char pointer
    ezxml_t xml_tmp;

    xml_tmp = ezxml_get(qexml, "output", 0, "boundary_conditions", 0,
                        "assume_isolated", -1);
    if (xml_tmp)
    {
        char* assume_iso = xml_tmp->txt;
        if (string_start_with(assume_iso, "2D", true))
        {
            *dim = '2';
        }
    }
    // next get the pseudo pot directory
    xml_tmp =
        ezxml_get(qexml, "input", 0, "control_variables", 0, "pseudo_dir", -1);
    if (!xml_tmp)
    {
        error_msg("Parsing pseudo_dir from data-file-schema.xml file");
    }
    tmp_str = xml_tmp->txt;

    *pseudo_dir = malloc(strlen(tmp_str) +
                         1);  // we need to free this outside of this function
    CHECK_ALLOC(*pseudo_dir);

    strcpy(*pseudo_dir, tmp_str);
    // printf("pseudo dir : %s",*pseudo_dir);
    //  get ntypes
    ezxml_t atom_specs = ezxml_get(qexml, "output", 0, "atomic_species", -1);
    if (atom_specs == NULL)
    {
        error_msg("error reading atomic spices from data-file-schema.xml file");
    }

    tmp_str = ezxml_attr(atom_specs, "ntyp");
    if (!tmp_str)
    {
        error_msg("error ntyp attribute from data-file-schema.xml file");
    }
    ND_int atomic_ntypes = atoll(tmp_str);
    if (ntype)
    {
        *ntype = atomic_ntypes;
    }
    //
    *pseudo_pots = malloc(sizeof(char*) * (atomic_ntypes));
    CHECK_ALLOC(*pseudo_pots);

    char** pot_tmp = *pseudo_pots;

    char* atom_type_symbol =
        calloc((atomic_ntypes) * 16, sizeof(*atom_type_symbol));
    CHECK_ALLOC(atom_type_symbol);

    for (ND_int itype = 0; itype < atomic_ntypes; ++itype)
    {
        xml_tmp =
            ezxml_get(atom_specs, "species", (int)itype, "pseudo_file", -1);
        if (!xml_tmp)
        {
            error_msg("Parsing species from data-file-schema.xml file");
        }
        tmp_str = xml_tmp->txt;
        // printf("%d : %s \n",(int)itype, tmp_str);
        pot_tmp[itype] = malloc(1 + strlen(tmp_str));
        CHECK_ALLOC(pot_tmp[itype]);

        strcpy(pot_tmp[itype], tmp_str);

        // read atomic type
        xml_tmp = ezxml_get(atom_specs, "species", (int)itype, "");
        if (!xml_tmp)
        {
            error_msg(
                "error reading atomic spices from data-file-schema.xml "
                "file");
        }
        tmp_str = ezxml_attr(xml_tmp, "name");
        if (!tmp_str)
        {
            error_msg("error name attribute from data-file-schema.xml file");
        }
        // remove white spaces
        while (isspace((unsigned char)(*tmp_str)))
        {
            ++tmp_str;
        }
        if (strlen(tmp_str) > 15)
        {
            error_msg("Atomic Name is very long.");
        }
        strncpy_custom(atom_type_symbol + 16 * itype, tmp_str, 16);
    }

    // get number of atoms
    xml_tmp = ezxml_get(qexml, "output", 0, "atomic_structure", -1);
    if (!xml_tmp)
    {
        error_msg("Parsing atomic_structure from data-file-schema.xml file");
    }

    tmp_str = ezxml_attr(xml_tmp, "nat");
    if (!tmp_str)
    {
        error_msg("error nat attribute from data-file-schema.xml file");
    }
    *natoms = atoll(tmp_str);
    // get alat
    tmp_str = ezxml_attr(xml_tmp, "alat");
    if (!tmp_str)
    {
        error_msg("error alat attribute from data-file-schema.xml file");
    }
    alat[0] = atof(tmp_str);
    alat[1] = alat[0];
    alat[2] = alat[0];
    //
    // read atomic positons and set its type.
    int* atom_type = malloc((*natoms) * sizeof(*atom_type));
    CHECK_ALLOC(atom_type);
    if (atomic_type)
    {
        *atomic_type = atom_type;
    }
    //
    ELPH_float* atomic_pos = NULL;
    if (atomic_positons)
    {
        atomic_pos = malloc(3 * (*natoms) * sizeof(*atomic_pos));
        CHECK_ALLOC(atomic_pos);
        *atomic_positons = atomic_pos;
    }
    //
    for (ND_int ia = 0; ia < *natoms; ++ia)
    {
        xml_tmp = ezxml_get(qexml, "output", 0, "atomic_structure", 0,
                            "atomic_positions", 0, "atom", (int)ia, "");
        if (!xml_tmp)
        {
            error_msg(
                "Parsing atomic_positions from data-file-schema.xml file");
        }
        tmp_str = ezxml_attr(xml_tmp, "name");
        if (!tmp_str)
        {
            error_msg("error name attribute from data-file-schema.xml file");
        }
        while (isspace((unsigned char)(*tmp_str)))
        {
            ++tmp_str;
        }
        //
        int itype = -1;
        for (int it = 0; it < atomic_ntypes; ++it)
        {
            if (!strcmp(atom_type_symbol + 16 * it, tmp_str))
            {
                itype = it;
                break;
            }
        }
        if (itype < 0)
        {
            error_msg("Cannot find atomic type.");
        }
        //
        atom_type[ia] = itype;
        //
        if (atomic_positons)
        {
            tmp_str = xml_tmp->txt;
            if (parse_floats_from_string(tmp_str, atomic_pos + 3 * ia, 3) != 3)
            {
                error_msg("Error parsing atomic positions");
            }
        }
    }
    //
    //
    // get fft dims
    xml_tmp = ezxml_get(qexml, "output", 0, "basis_set", 0, "fft_grid", -1);
    if (!xml_tmp)
    {
        error_msg("Parsing fft_grid from data-file-schema.xml file");
    }

    tmp_str = ezxml_attr(xml_tmp, "nr1");
    if (!tmp_str)
    {
        error_msg("error nr1 attribute from data-file-schema.xml file");
    }
    fft_dims[0] = atoll(tmp_str);

    tmp_str = ezxml_attr(xml_tmp, "nr2");
    if (!tmp_str)
    {
        error_msg("error nr2 attribute from data-file-schema.xml file");
    }
    fft_dims[1] = atoll(tmp_str);

    tmp_str = ezxml_attr(xml_tmp, "nr3");
    if (!tmp_str)
    {
        error_msg("error nr3 attribute from data-file-schema.xml file");
    }
    fft_dims[2] = atoll(tmp_str);

    // get lattice vectors
    ELPH_float a_tmp_read[3];

    xml_tmp = ezxml_get(qexml, "output", 0, "atomic_structure", 0, "cell", 0,
                        "a1", -1);
    if (!xml_tmp)
    {
        error_msg("Parsing a1 from data-file-schema.xml file");
    }
    tmp_str = xml_tmp->txt;
    if (parse_floats_from_string(tmp_str, a_tmp_read, 3) != 3)
    {
        error_msg("Error parsing a1 vec from data-file-schema.xml");
    }
    for (int ix = 0; ix < 3; ++ix)
    {
        lat_vec[3 * ix + 0] = a_tmp_read[ix];  // a[:,i]
    }

    xml_tmp = ezxml_get(qexml, "output", 0, "atomic_structure", 0, "cell", 0,
                        "a2", -1);
    if (!xml_tmp)
    {
        error_msg("Parsing a2 from data-file-schema.xml file");
    }
    tmp_str = xml_tmp->txt;
    if (parse_floats_from_string(tmp_str, a_tmp_read, 3) != 3)
    {
        error_msg("Error parsing a2 vec from data-file-schema.xml");
    }
    for (int ix = 0; ix < 3; ++ix)
    {
        lat_vec[3 * ix + 1] = a_tmp_read[ix];
    }

    xml_tmp = ezxml_get(qexml, "output", 0, "atomic_structure", 0, "cell", 0,
                        "a3", -1);
    if (!xml_tmp)
    {
        error_msg("Parsing a3 from data-file-schema.xml file");
    }
    tmp_str = xml_tmp->txt;
    if (parse_floats_from_string(tmp_str, a_tmp_read, 3) != 3)
    {
        error_msg("Error parsing a3 vec from data-file-schema.xml");
    }
    for (int ix = 0; ix < 3; ++ix)
    {
        lat_vec[3 * ix + 2] = a_tmp_read[ix];
    }

    // check if soc is present
    xml_tmp =
        ezxml_get(qexml, "output", 0, "magnetization", 0, "spinorbit", -1);
    //
    if (!xml_tmp)
    {
        error_msg(
            "Parsing magnetization,spinorbit from data-file-schema.xml file");
    }
    tmp_str = xml_tmp->txt;
    //
    strncpy_custom(tmp_read, tmp_str, ELPH_XML_READ_LINE_SIZE);
    lowercase_str(tmp_read);

    if (strstr(tmp_read, "true"))
    {
        *is_soc_present = true;
    }

    // nmag
    // first check if this is a lsda calc
    bool lsda = false;
    xml_tmp = ezxml_get(qexml, "output", 0, "magnetization", 0, "lsda", -1);
    if (!xml_tmp)
    {
        error_msg("Parsing magnetization, lsda from data-file-schema.xml file");
    }
    tmp_str = xml_tmp->txt;
    strncpy_custom(tmp_read, tmp_str, ELPH_XML_READ_LINE_SIZE);
    lowercase_str(tmp_read);

    if (strstr(tmp_read, "true"))
    {
        lsda = true;
    }

    *ph_tim_rev = true;

    // check if there is noinv flag
    bool no_inv = false;
    xml_tmp = ezxml_get(qexml, "input", 0, "symmetry_flags", 0, "noinv", -1);
    if (!xml_tmp)
    {
        error_msg("Parsing noinv from data-file-schema.xml file");
    }
    tmp_str = xml_tmp->txt;
    strncpy_custom(tmp_read, tmp_str, ELPH_XML_READ_LINE_SIZE);
    lowercase_str(tmp_read);

    if (strstr(tmp_read, "true"))
    {
        no_inv = true;
    }

    if (no_inv)
    {
        *ph_tim_rev = false;
    }

    if (lsda)
    {
        *nmag = 2;
    }
    else
    {
        bool is_non_collinear = false;
        bool mag_system = false;

        xml_tmp =
            ezxml_get(qexml, "output", 0, "magnetization", 0, "noncolin", -1);
        if (!xml_tmp)
        {
            error_msg("Parsing noncolin from data-file-schema.xml file");
        }
        tmp_str = xml_tmp->txt;
        strncpy_custom(tmp_read, tmp_str, ELPH_XML_READ_LINE_SIZE);
        lowercase_str(tmp_read);

        if (strstr(tmp_read, "true"))
        {
            is_non_collinear = true;
            if (nspinor)
            {
                *nspinor = 2;
            }
        }

        if (is_non_collinear)
        {
            xml_tmp = ezxml_get(qexml, "output", 0, "magnetization", 0,
                                "do_magnetization", -1);
            if (!xml_tmp)
            {
                error_msg(
                    "Parsing do_magnetization from data-file-schema.xml file");
            }
            tmp_str = xml_tmp->txt;

            strncpy_custom(tmp_read, tmp_str, ELPH_XML_READ_LINE_SIZE);
            lowercase_str(tmp_read);

            if (strstr(tmp_read, "true"))
            {
                mag_system = true;
            }
        }
        //
        if (mag_system)
        {
            *nmag = 4;
        }
        else
        {
            *nmag = 1;
        }

        if (mag_system)
        {
            // NM : This is can be changed later in this function
            *ph_tim_rev = false;
        }
    }

    // Finally read the phonon symmetries from the xml file
    // first get number of symmetries
    xml_tmp = ezxml_get(qexml, "output", 0, "symmetries", 0, "nsym", -1);
    if (!xml_tmp)
    {
        error_msg("Parsing nsym from data-file-schema.xml file");
    }
    tmp_str = xml_tmp->txt;

    *nph_sym = atoi(tmp_str);

    // create a temp buffers for reading
    // we create 2*nph_sym sets of symmetries (factor 2 to store the time rev
    // case) The second half is used as temporary buffer
    // and also to store time_rev symmetries in nmag == 1 case
    // Note: the second half([nph_sym:]) are only used when nmag == 1
    // and these are not computed in this function. but we allocate
    // the storage
    //
    // in case of nmag == 2 or nmag == 4, we donot have a pure time reversal
    // symmetry but we can have a rotation + time rev
    //
    *ph_sym_mats = malloc(sizeof(ELPH_float) * 3 * 3 * 2 * (*nph_sym));
    CHECK_ALLOC(*ph_sym_mats);

    *ph_sym_tau = malloc(sizeof(ELPH_float) * 3 * 2 * (*nph_sym));
    CHECK_ALLOC(*ph_sym_tau);

    *ph_trevs = calloc(2 * (*nph_sym), sizeof(bool));
    CHECK_ALLOC(*ph_trevs);
    bool* trev_ptr = *ph_trevs;

    for (ND_int isym = 0; isym < 2 * (*nph_sym); ++isym)
    {
        // a  smart compiller will remove this loop.
        // but lets stick to standard and leave these to compilers
        trev_ptr[isym] = false;
    }

    bool inversion_sym = false;
    bool mag_sym_found = false;

    ELPH_float I3x3[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    for (ND_int isym = 0; isym < *nph_sym; ++isym)
    {
        ezxml_t sym_xml_tmp = ezxml_get(qexml, "output", 0, "symmetries", -1);
        if (sym_xml_tmp == NULL)
        {
            error_msg("Error reading in phonon symmetries");
        }
        //
        // rotation matrix (stored in transposed order) and frac .translation
        // are store in crystals coordinates parse rotation (in crystal units)
        // S_cart = (alat@S_crys^T@blat^T)
        xml_tmp = ezxml_get(sym_xml_tmp, "symmetry", (int)isym, "rotation", -1);
        if (!xml_tmp)
        {
            error_msg("Parsing rotation from data-file-schema.xml file");
        }
        tmp_str = xml_tmp->txt;
        if (parse_floats_from_string(tmp_str, (*ph_sym_mats) + 9 * isym, 9) !=
            9)
        {
            error_msg(
                "Error parsing symmetry matrices from data-file-schema.xml");
        }

        // parse translation (in crystal units) tau_cart = alat@tau_crys
        // It should be noted that we use Sx + v convention, but qe uses S(x+v)
        // so our v = S*tau_qe
        xml_tmp = ezxml_get(sym_xml_tmp, "symmetry", (int)isym,
                            "fractional_translation", -1);
        if (!xml_tmp)
        {
            error_msg(
                "Parsing fractional_translation from data-file-schema.xml "
                "file");
        }
        tmp_str = xml_tmp->txt;
        if (parse_floats_from_string(tmp_str, (*ph_sym_tau) + 3 * isym, 3) != 3)
        {
            error_msg(
                "Error parsing frac. trans. vecs from data-file-schema.xml");
        }
        // check if this symmetry corresponds to time_reversal symmetry for nmag
        // == 4
        xml_tmp = ezxml_get(sym_xml_tmp, "symmetry", (int)isym, "info", -1);
        if (!xml_tmp)
        {
            error_msg("Parsing symmetry info from data-file-schema.xml file");
        }
        const char* trev_tmp_str = ezxml_attr(xml_tmp, "time_reversal");
        if (trev_tmp_str)
        {
            strncpy_custom(tmp_read, trev_tmp_str, ELPH_XML_READ_LINE_SIZE);
            lowercase_str(tmp_read);

            if (strstr(tmp_read, "true"))
            {
                mag_sym_found = true;
                *ph_tim_rev = true;
                trev_ptr[isym] = true;
            }
            else
            {
                trev_ptr[isym] = false;
            }
        }
        //

        // check if inversion is present
        ELPH_float sum = 0;
        ELPH_float* sym_mat_tmpp = (*ph_sym_mats) + 9 * isym;
        for (int ix = 0; ix < 9; ++ix)
        {
            sum += fabs((I3x3[ix] + sym_mat_tmpp[ix]) *
                        (I3x3[ix] + sym_mat_tmpp[ix]));
        }
        sum = sqrt(sum);
        if (sum < ELPH_EPS)
        {
            inversion_sym = true;
        }
    }

    // incase magnetic symmetries found, we need to rearrange the symmetries and
    // push back
    if (mag_sym_found)
    {
        // send all timerev symmetries to end
        // first make a temporary copy
        ELPH_float* sym_src_ptr = *ph_sym_mats + (*nph_sym) * 9;
        ELPH_float* sym_des_ptr = *ph_sym_mats;
        memcpy(sym_src_ptr, sym_des_ptr, sizeof(ELPH_float) * 9 * (*nph_sym));

        ELPH_float* tau_src_ptr = *ph_sym_tau + (*nph_sym) * 3;
        ELPH_float* tau_des_ptr = *ph_sym_tau;
        memcpy(tau_src_ptr, tau_des_ptr, sizeof(ELPH_float) * 3 * (*nph_sym));

        bool* trev_src_ptr = *ph_trevs + (*nph_sym);
        bool* trev_des_ptr = *ph_trevs;
        memcpy(trev_src_ptr, trev_des_ptr, sizeof(bool) * (*nph_sym));

        ND_int t_sym_copiled = 0;
        ND_int n_sym_copiled = 0;

        if (*nph_sym % 2)
        {
            error_msg("Wrong number of symmetries");
        }
        for (ND_int isym = 0; isym < *nph_sym; ++isym)
        {
            ND_int dest_shift = 0;
            if (trev_src_ptr[isym])
            {
                dest_shift = (*nph_sym) / 2 + t_sym_copiled;
                ++t_sym_copiled;
            }
            else
            {
                dest_shift = n_sym_copiled;
                ++n_sym_copiled;
            }
            ELPH_float* sym_src_ptr_tmp = sym_src_ptr + 9 * isym;
            ELPH_float* sym_des_ptr_tmp = sym_des_ptr + 9 * dest_shift;
            memcpy(sym_des_ptr_tmp, sym_src_ptr_tmp, sizeof(ELPH_float) * 9);

            ELPH_float* tau_src_ptr_tmp = tau_src_ptr + 3 * isym;
            ELPH_float* tau_des_ptr_tmp = tau_des_ptr + 3 * dest_shift;
            memcpy(tau_des_ptr_tmp, tau_src_ptr_tmp, sizeof(ELPH_float) * 3);

            trev_des_ptr[dest_shift] = trev_src_ptr[isym];
        }
        if (t_sym_copiled != n_sym_copiled)
        {
            error_msg("Wrong number of symmetries");
        }
    }

    *ph_mag_symm = mag_sym_found;
    // set time reversal to false if inversion is present (no longer needed) for
    // non mangetic materials.
    if (inversion_sym && !mag_sym_found)
    {
        *ph_tim_rev = false;
    }

    if (!atomic_type)
    {
        free(atom_type);
    }
    if (!atomic_positons)
    {
        free(atomic_pos);
    }
    free(atom_type_symbol);
    free(tmp_read);
    ezxml_free(qexml);
    fclose(fp);
}
