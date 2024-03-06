/*
This file contains functions that parses upf data from pseudo pots .
Only Local part, grids , valance electron info are read rest are available in
yambo
*/

#include "ezxml/ezxml.h"
#include "io.h"

static void parse_upf2(FILE* fp, struct local_pseudo* loc_pseudo);

static void parse_upf1(FILE* fp, struct local_pseudo* loc_pseudo);

static void get_upf2_element(FILE* fp, char* atomic_sym);

static void get_upf1_element(FILE* fp, char* atomic_sym);

// ===================================================================================

void parse_upf(const char* filename, struct local_pseudo* loc_pseudo)
{
    FILE* fp = fopen(filename, "r");
    if (fp == NULL)
    {
        error_msg("Unable to open the upf file");
    }

    int upf_version = 0;
    char start_str[100];
    fgets(start_str, 100, fp);
    if (strstr(start_str, "PP_INFO"))
    {
        upf_version = 1;
    }
    else if (strstr(start_str, "UPF version"))
    {
        upf_version = 2;
    }
    else
    {
        error_msg("Only UPF-1 or UPF-2 file formats supported");
    }

    if (fseek(fp, 0, SEEK_SET) != 0)
    {
        error_msg("error setting the start seek for upf file");
    }

    if (upf_version == 1)
    {
        parse_upf1(fp, loc_pseudo);
    }
    else if (upf_version == 2)
    {
        parse_upf2(fp, loc_pseudo);
    }
    else
    {
        error_msg("Only UPF-1 or UPF-2 file formats supported");
    }
    fclose(fp);
}

void get_upf_element(const char* filename, char* atomic_sym)
{
    FILE* fp = fopen(filename, "r");
    if (fp == NULL)
    {
        error_msg("Unable to open the upf file");
    }

    int upf_version = 0;
    char start_str[100];
    fgets(start_str, 100, fp);
    if (strstr(start_str, "<PP_INFO>"))
    {
        upf_version = 1;
    }
    else if (strstr(start_str, "<UPF version"))
    {
        upf_version = 2;
    }
    else
    {
        error_msg("Only UPF-1 or UPF-2 file formats supported");
    }

    if (fseek(fp, 0, SEEK_SET) != 0)
    {
        error_msg("error setting the start seek for upf file");
    }

    if (upf_version == 1)
    {
        get_upf1_element(fp, atomic_sym);
    }
    else if (upf_version == 2)
    {
        get_upf2_element(fp, atomic_sym);
    }
    else
    {
        error_msg("Only UPF-1 or UPF-2 file formats supported");
    }

    fclose(fp);
}

// ===================================================================================

static void parse_upf2(FILE* fp, struct local_pseudo* loc_pseudo)
{
    // parses upf format 2
    /*
    This function initiates and allocates memory for Vloc, r_grid, rab_grid
    */
    if (fseek(fp, 0, SEEK_SET) != 0)
    {
        error_msg("error setting the start seek for upf file");
    }

    ezxml_t upfFP = ezxml_parse_fp(fp);
    if (upfFP == NULL)
    {
        error_msg("Reading from pseudo potential file failed");
    }

    /* Read header details */
    ezxml_t header = ezxml_get(upfFP, "PP_HEADER", -1);
    if (header == NULL)
    {
        error_msg("Reading header from pseudo potential file failed");
    }

    ELPH_float Zval = atof(ezxml_attr(header, "z_valence"));
    ND_int ngrid = atoi(ezxml_attr(header, "mesh_size"));

    if (strcmp(ezxml_attr(header, "pseudo_type"), "NC"))
    {
        error_msg("Pseudo potential is not norm conserving");
    }

    /* Parse mesh grid*/
    ezxml_t mesh = ezxml_get(upfFP, "PP_MESH", 0, "PP_R", -1);
    if (mesh == NULL)
    {
        error_msg("Reading radial mesh from pseudo potential file failed");
    }

    /* parse rab */
    ezxml_t dmesh = ezxml_get(upfFP, "PP_MESH", 0, "PP_RAB", -1);
    if (dmesh == NULL)
    {
        error_msg("Reading radial rab from pseudo potential file failed");
    }

    /* parse local potential */
    ezxml_t loc_pot = ezxml_get(upfFP, "PP_LOCAL", -1);
    if (loc_pot == NULL)
    {
        error_msg("Reading local part from pseudo potential file failed");
    }

    loc_pseudo->ngrid = ngrid;
    loc_pseudo->Zval = Zval;

    loc_pseudo->Vloc_atomic = malloc(sizeof(ELPH_float) * ngrid);
    CHECK_ALLOC(loc_pseudo->Vloc_atomic);

    loc_pseudo->r_grid = malloc(sizeof(ELPH_float) * ngrid);
    CHECK_ALLOC(loc_pseudo->r_grid);

    loc_pseudo->rab_grid = malloc(sizeof(ELPH_float) * ngrid);
    CHECK_ALLOC(loc_pseudo->rab_grid);

    if (parser_doubles_from_string(loc_pot->txt, loc_pseudo->Vloc_atomic) != ngrid)
    {
        error_msg("Parsing local potential from upf-2 failed");
    }
    if (parser_doubles_from_string(mesh->txt, loc_pseudo->r_grid) != ngrid)
    {
        error_msg("Parsing radial mesh from upf-2 failed");
    }
    if (parser_doubles_from_string(dmesh->txt, loc_pseudo->rab_grid) != ngrid)
    {
        error_msg("Parsing Rab from upf-2 failed");
    }

    ezxml_free(upfFP);
}

static void parse_upf1(FILE* fp, struct local_pseudo* loc_pseudo)
{
    // parses upf 1
    /*
    This function initiates and allocates memory for Vloc, r_grid, rab_grid
    */
    // upf1 does not confine to XML format, so we make some kind of xml format
    // FILE* fp = fopen(filename, "r");
    // if (fp == NULL)
    //     error_msg("Unable to open the upf file");

    // header details cannot be read like XML
    // parse line by line and read the required header details
    if (fseek(fp, 0, SEEK_SET) != 0)
    {
        error_msg("error setting the start seek for upf file");
    }

    char* xml_buf = malloc(1000);
    CHECK_ALLOC(xml_buf);

    float Zval = -1;
    int ngrid = -1;
    bool read_header = false;
    while (fgets(xml_buf, 1000, fp))
    {
        if (strstr(xml_buf, "<PP_HEADER>"))
        {
            read_header = true;
            break;
        }
    }
    if (!read_header)
    {
        error_msg("Unable to parse header info from upf file");
    }

    // Now read line by line
    fgets(xml_buf, 1000, fp); // version number
    fgets(xml_buf, 1000, fp); // element label
    fgets(xml_buf, 1000, fp); // pseudo type
    if (!strstr(xml_buf, "NC"))
    {
        error_msg("Pseudo potential is not norm conserving");
    }

    fgets(xml_buf, 1000, fp); // nlcc
    fgets(xml_buf, 1000, fp); // XC info
    fgets(xml_buf, 1000, fp); // Valence electrons
    char* tmp_token = strtok(xml_buf, " ");
    sscanf(tmp_token, "%f", &Zval);
    fgets(xml_buf, 1000, fp); // total energy
    fgets(xml_buf, 1000, fp); // suggested cutoff
    fgets(xml_buf, 1000, fp); // max l
    fgets(xml_buf, 1000, fp); // max nmesh
    tmp_token = strtok(xml_buf, " ");

    sscanf(tmp_token, "%d", &ngrid);

    free(xml_buf);
    xml_buf = NULL;

    if (fseek(fp, 0, SEEK_END) != 0)
    {
        error_msg("error setting the end seek for upf file");
    }

    long length = ftell(fp);

    if (length < 0)
    {
        error_msg("error getting the lenght of total upf file");
    }

    if (fseek(fp, 0, SEEK_SET) != 0)
    {
        error_msg("error setting the start seek for upf file");
    }

    length += 15; // to add root elements and one for null terminator

    xml_buf = malloc(length+1);
    CHECK_ALLOC(xml_buf);

    xml_buf[length] = '\0';
    strcpy(xml_buf, "<root>");
    xml_buf[6] = '\n';

    if (fread(xml_buf + 7, 1, length - 15, fp) != (length - 15))
    {
        error_msg("Reading loading data from pseudo potential file failed");
    }

    strcpy(xml_buf + length - 8, "</root>");

    ezxml_t upfFP = ezxml_parse_str(xml_buf, length - 15);
    if (upfFP == NULL)
    {
        error_msg("Reading from pseudo potential file failed");
    }

    /* Parse mesh grid*/
    ezxml_t mesh = ezxml_get(upfFP, "PP_MESH", 0, "PP_R", -1);
    if (mesh == NULL)
    {
        error_msg("Reading radial mesh from pseudo potential file failed");
    }

    /* parse rab */
    ezxml_t dmesh = ezxml_get(upfFP, "PP_MESH", 0, "PP_RAB", -1);
    if (dmesh == NULL)
    {
        error_msg("Reading radial rab from pseudo potential file failed");
    }

    /* parse local potential */
    ezxml_t loc_pot = ezxml_get(upfFP, "PP_LOCAL", -1);
    if (loc_pot == NULL)
    {
        error_msg("Reading local part from pseudo potential file failed");
    }

    loc_pseudo->ngrid = ngrid;
    loc_pseudo->Zval = Zval;

    loc_pseudo->Vloc_atomic = malloc(sizeof(ELPH_float) * ngrid);
    CHECK_ALLOC(loc_pseudo->Vloc_atomic);

    loc_pseudo->r_grid = malloc(sizeof(ELPH_float) * ngrid);
    CHECK_ALLOC(loc_pseudo->r_grid);

    loc_pseudo->rab_grid = malloc(sizeof(ELPH_float) * ngrid);
    CHECK_ALLOC(loc_pseudo->rab_grid);

    if (parser_doubles_from_string(loc_pot->txt, loc_pseudo->Vloc_atomic) != ngrid)
    {
        error_msg("Parsing local potential from upf-1 failed");
    }
    if (parser_doubles_from_string(mesh->txt, loc_pseudo->r_grid) != ngrid)
    {
        error_msg("Parsing radial mesh from upf-1 failed");
    }
    if (parser_doubles_from_string(dmesh->txt, loc_pseudo->rab_grid) != ngrid)
    {
        error_msg("Parsing Rab from upf-1 failed");
    }

    ezxml_free(upfFP);
    free(xml_buf);
}

//----

static void get_upf2_element(FILE* fp, char* atomic_sym)
{
    /*
    atomic_sym must be atleast 3 bytes long
    */
    if (fseek(fp, 0, SEEK_SET) != 0)
    {
        error_msg("error setting the start seek for upf file");
    }

    atomic_sym[2] = '\0';

    ezxml_t upfFP = ezxml_parse_fp(fp);
    if (upfFP == NULL)
    {
        error_msg("Reading from pseudo potential file failed");
    }

    /* Read header details */

    ezxml_t header = ezxml_get(upfFP, "PP_HEADER", -1);
    if (header == NULL)
    {
        error_msg("Reading header from pseudo potential file failed");
    }

    const char* temp_upf_ptr = ezxml_attr(header, "element");

    // get the element
    memcpy(atomic_sym, temp_upf_ptr, sizeof(char) * 2);
    ezxml_free(upfFP);
}

//----

static void get_upf1_element(FILE* fp, char* atomic_sym)
{
    /*
    atomic_sym must be atleast 3 bytes long
    */
    if (fseek(fp, 0, SEEK_SET) != 0)
    {
        error_msg("error setting the start seek for upf file");
    }

    atomic_sym[2] = '\0';
    char* xml_buf = malloc(1000);
    CHECK_ALLOC(xml_buf);

    bool read_header = false;
    while (fgets(xml_buf, 1000, fp))
    {
        if (strstr(xml_buf, "<PP_HEADER>"))
        {
            read_header = true;
            break;
        }
    }
    if (!read_header)
    {
        error_msg("Unable to parse header info from upf file");
    }

    // Now read line by line
    fgets(xml_buf, 1000, fp); // version number
    fgets(xml_buf, 1000, fp); // element label
    char tmp_read[30];
    char* token_tmp = strtok(xml_buf, " ");
    sscanf(token_tmp, "%s", tmp_read);
    size_t atm_lab_size = strlen(tmp_read);
    if (atm_lab_size > 2)
    {
        error_msg("Buffer over when reading element from upf file");
    }
    strcpy(atomic_sym, tmp_read);
    if (atm_lab_size == 1)
    {
        atomic_sym[1] = ' ';
    }

    fgets(xml_buf, 1000, fp); // pseudo type
    if (!strstr(xml_buf, "NC"))
    {
        error_msg("Pseudo potential is not norm conserving");
    }

    free(xml_buf);
    xml_buf = NULL;
}

//----
