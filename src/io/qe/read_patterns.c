/*
This file contains function that parses pattern.xml files
*/
#include "qe_io.h"

void read_pattern_qe(const char* pat_file, struct Lattice* lattice,
                     ELPH_cmplx* restrict pat_vecs)
{
    /*
    patvec dim (nmodes, natom, 3)
    */
    // open the xml file
    ezxml_t patxml = ezxml_parse_file(pat_file);
    if (patxml == NULL)
    {
        error_msg("Opening pattern file failed \n");
    }

    // first get the number of irrep

    int nirrep =
        atoll(ezxml_get(patxml, "IRREPS_INFO", 0, "NUMBER_IRR_REP", -1)->txt);

    ND_int natom = lattice->atomic_pos->dims[0];
    ND_int nmodes = natom * 3;

    ELPH_float* pat_tmp_read = malloc(sizeof(ELPH_float) * 2 * nmodes);
    char rep_str[100];
    char pert_str[100];

    ND_int imode = 0;
    for (int irep = 0; irep < nirrep; ++irep)
    {
        sprintf(rep_str, "REPRESENTION.%d", irep + 1);
        int npert = atoll(ezxml_get(patxml, "IRREPS_INFO", 0, rep_str, 0,
                                    "NUMBER_OF_PERTURBATIONS", -1)
                              ->txt);

        for (int ipert = 0; ipert < npert; ++ipert)
        {
            sprintf(pert_str, "PERTURBATION.%d", ipert + 1);

            char* per_vec_str =
                ezxml_get(patxml, "IRREPS_INFO", 0, rep_str, 0, pert_str, 0,
                          "DISPLACEMENT_PATTERN", -1)
                    ->txt;

            if (parser_doubles_from_string(per_vec_str, pat_tmp_read) !=
                (2 * nmodes))
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
}
