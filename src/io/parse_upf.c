/*
This file contains functions that parses upf data from pseudo pots . 
Only Local part, grids , valance electron info are read rest are available in yambo
*/

#include "ezxml/ezxml.h"
#include "io.h"


static void parse_floats_from_string(char * str,  ND_int nparse, ELPH_float * out);

void parse_upf2(const char * filename, ELPH_float * Zval, ND_array(Nd_floatS)* Vloc, \
                ND_array(Nd_floatS)* r_grid, ND_array(Nd_floatS)* rab_grid)
{
    /*
    This function initiates and allocates memory for Vloc, r_grid, rab_grid
    */
    ezxml_t upfFP = ezxml_parse_file(filename);
    if (upfFP == NULL) error_msg("Reading from pseudo potential file failed");
    
    /* Read header details */
    ezxml_t header = ezxml_get(upfFP,"PP_HEADER", -1);
    if (header == NULL) error_msg("Reading header from pseudo potential file failed");

    *Zval = (ELPH_float)atof(ezxml_attr(header, "z_valence")) ;
    ND_int ngrid = atoi(ezxml_attr(header, "mesh_size")) ;
    
    if (strcmp(ezxml_attr(header, "pseudo_type"), "NC")) \
            error_msg("Pseudo potential is not norm conserving");


    /* Parse mesh grid*/
    ezxml_t mesh    = ezxml_get(upfFP,"PP_MESH",0, "PP_R", -1) ;
    if (mesh == NULL) error_msg("Reading radial mesh from pseudo potential file failed");

    /* parse rab */
    ezxml_t dmesh   = ezxml_get(upfFP,"PP_MESH",0, "PP_RAB", -1) ;
    if (dmesh == NULL) error_msg("Reading radial rab from pseudo potential file failed");

    /* parse local potential */
    ezxml_t loc_pot = ezxml_get(upfFP,"PP_LOCAL", -1) ;
    if (loc_pot == NULL) error_msg("Reading local part from pseudo potential file failed");

    ND_function(init,Nd_floatS) (Vloc,    1, nd_idx{ngrid}); 
    ND_function(init,Nd_floatS) (r_grid,  1, nd_idx{ngrid}); 
    ND_function(init,Nd_floatS) (rab_grid,1, nd_idx{ngrid}); 

    ND_function(malloc,Nd_floatS) (Vloc);
    ND_function(malloc,Nd_floatS) (r_grid);
    ND_function(malloc,Nd_floatS) (rab_grid);

    parse_floats_from_string(loc_pot->txt,  ngrid, Vloc->data);
    parse_floats_from_string(mesh->txt,     ngrid, r_grid->data);
    parse_floats_from_string(dmesh->txt,    ngrid, rab_grid->data);
    
    ezxml_free(upfFP);

}


//----

void get_upf_element(const char * filename, char* atomic_sym)
{
    /*
    This function initiates and allocates memory for Vloc, r_grid, rab_grid
    atomic_sym must be atleast 3 bytes long
    */
    atomic_sym[2] = '\0';

    ezxml_t upfFP = ezxml_parse_file(filename);
    if (upfFP == NULL) error_msg("Reading from pseudo potential file failed");

    
    /* Read header details */

    ezxml_t header = ezxml_get(upfFP,"PP_HEADER", -1);
    if (header == NULL)  error_msg("Reading header from pseudo potential file failed");

    const char * temp_upf_ptr = ezxml_attr(header, "element");

    // get the element
    memcpy(atomic_sym,temp_upf_ptr, sizeof(char)*2);
    ezxml_free(upfFP);
}


//----

static void parse_floats_from_string(char * str,  ND_int nparse, ELPH_float * out)
{   
    char * start_ptr = str;
    char * end_ptr;
    ND_int  iparsed = 0;
    double temp = 0.0 ;
    while(true)
    {   
        if (*start_ptr == '\0') break ;
        temp = strtod(start_ptr,&end_ptr) ;
        if (start_ptr == end_ptr) break ;
        start_ptr = end_ptr;
        out[ iparsed] = temp ;
        ++ iparsed;
    }

    if ( iparsed != nparse) error_msg("str to float Conversion failed");

}

