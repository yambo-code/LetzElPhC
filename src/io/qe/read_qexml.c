/*
This file contains function that parses data-file-schema.xml file
*/

#include "qe_io.h"
#include <string.h>

void parse_qexml(const char * xml_file, ELPH_float * alat, \
        ND_int * fft_dims, char ** pseudo_dir, char *** pseudo_pots)
{
    /*
    get pseudo potential information, fft information from xml file
    
    need to free pseudo_dir and pseudo_pots pointers
    */
    char * tmp_str; // tmp variable for ezxml char pointer

    ezxml_t qexml = ezxml_parse_file(xml_file);
    if (qexml == NULL) error_msg("Error parsing data-file-schema.xml file");

    // first get the pseudo pot directory
    tmp_str = ezxml_get(qexml,"input", 0, "control_variables",0, \
                        "pseudo_dir",-1)->txt;
    *pseudo_dir = malloc(strlen(tmp_str)+1); // we need to free this outside of this function
    strcpy(*pseudo_dir,tmp_str);
    printf("pseudo dir : %s",*pseudo_dir);
    // get ntypes
    ezxml_t atom_specs = ezxml_get(qexml, "input", 0, "atomic_species", -1);
    if (atom_specs == NULL) error_msg("error reading atomic spices from data-file-schema.xml file");
    
    ND_int ntype = atoll(ezxml_attr(atom_specs, "ntyp")) ;
    *pseudo_pots = malloc(sizeof(char *)*ntype);

    char ** pot_tmp = *pseudo_pots;

    for (ND_int itype = 0; itype < ntype; ++itype)
    {
        tmp_str = ezxml_get(atom_specs, "species", itype, "pseudo_file", -1)->txt;
        printf("%d : %s \n",(int)itype, tmp_str);
        pot_tmp[itype] = malloc(1+strlen(tmp_str));
        strcpy(pot_tmp[itype],tmp_str);
    }

    // get alat
    *alat = atof(ezxml_attr(ezxml_get(qexml, "input", 0, "atomic_structure", -1), "alat"));

    // get fft dims
    fft_dims[0] = atoll(ezxml_attr(ezxml_get(qexml, "output", 0, "basis_set", 0, "fft_grid", -1), "nr1"));
    fft_dims[1] = atoll(ezxml_attr(ezxml_get(qexml, "output", 0, "basis_set", 0, "fft_grid", -1), "nr2"));
    fft_dims[2] = atoll(ezxml_attr(ezxml_get(qexml, "output", 0, "basis_set", 0, "fft_grid", -1), "nr3"));

    ezxml_free(qexml);

}






