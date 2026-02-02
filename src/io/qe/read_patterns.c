/*
This file contains function that parses pattern.xml files
*/
#include <complex.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/dtypes.h"
#include "common/error.h"
#include "common/string_func.h"
#include "elphC.h"
#include "io/ezxml/ezxml.h"
#include "qe_io.h"

void read_pattern_qe(const char* pat_file, struct Lattice* lattice,
                     ELPH_cmplx* pat_vecs)
{
    /*
    patvec dim (nmodes, natom, 3)
    */
    // open the xml file
    FILE* fp = fopen(pat_file, "r");
    if (fp == NULL)
    {
        error_msg("Opening pattern.xml file failed");
    }

    ezxml_t patxml = ezxml_parse_fp(fp);
    if (patxml == NULL)
    {
        error_msg("parsing pattern file failed \n");
    }

    // first get the number of irrep

    ezxml_t xml_get = ezxml_get(patxml, "IRREPS_INFO", 0, "NUMBER_IRR_REP", -1);
    if (!xml_get)
    {
        error_msg("parsing irreps info in pattern file failed \n");
    }
    int nirrep = atoi(xml_get->txt);

    ND_int natom = lattice->natom;
    ND_int nmodes = natom * 3;

    ELPH_float* pat_tmp_read = malloc(sizeof(ELPH_float) * 2 * nmodes);
    CHECK_ALLOC(pat_tmp_read);

    char rep_str[100];
    char pert_str[100];

    ND_int imode = 0;
    for (int irep = 0; irep < nirrep; ++irep)
    {
        snprintf(rep_str, 100, "REPRESENTION.%d", irep + 1);

        ezxml_t xml_get_pert = ezxml_get(patxml, "IRREPS_INFO", 0, rep_str, 0,
                                         "NUMBER_OF_PERTURBATIONS", -1);
        if (!xml_get_pert)
        {
            error_msg(
                "parsing NUMBER_OF_PERTURBATIONS in pattern file failed \n");
        }

        int npert = atoi(xml_get_pert->txt);

        for (int ipert = 0; ipert < npert; ++ipert)
        {
            snprintf(pert_str, 100, "PERTURBATION.%d", ipert + 1);

            ezxml_t xml_get_ipert =
                ezxml_get(patxml, "IRREPS_INFO", 0, rep_str, 0, pert_str, 0,
                          "DISPLACEMENT_PATTERN", -1);
            if (!xml_get_ipert)
            {
                error_msg(
                    "parsing DISPLACEMENT_PATTERN in pattern file failed \n");
            }
            char* per_vec_str = xml_get_ipert->txt;

            if (parse_floats_from_string(per_vec_str, pat_tmp_read,
                                         2 * nmodes) != (2 * nmodes))
            {
                error_msg("Reading patterns.xml file failed");
            }

            for (ND_int i = 0; i < nmodes; ++i)
            {
                pat_vecs[imode * nmodes + i] =
                    pat_tmp_read[2 * i] + I * pat_tmp_read[2 * i + 1];
            }

            ++imode;
        }
    }
    if (imode != nmodes)
    {
        error_msg("Parsing patterns.xml failed");
    }

    free(pat_tmp_read);
    ezxml_free(patxml);
    fclose(fp);
}
