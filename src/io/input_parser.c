#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inih/ini.h"
#include "io.h"


#define MAX_STR_INPUT 10000
#define MAX_PSEUDO_NUM 200

static void Bcast_input_data(struct usr_input *input, int root, MPI_Comm comm);
static int handler(void* user, const char* section, const char* name,
            const char* value);

// function to alloc, initiate usr_input
void init_usr_input(struct usr_input ** input)
{
    // this function also sets defaults for the user input file
    *input = malloc(sizeof(struct usr_input));
    struct usr_input * inp = *input;

    inp->save_dir = calloc(MAX_STR_INPUT,sizeof(char));
    inp->pseudos  = malloc(sizeof(char *)*MAX_PSEUDO_NUM);

    for (int i = 0; i<MAX_PSEUDO_NUM; ++i) inp->pseudos[i] = NULL;

    // defaults
    inp->nkpool = 1;
    inp->nqpool = 1;
    inp->start_bnd = 0;
    inp->end_bnd   = 0;
    strcpy(inp->save_dir,"./SAVE");
    inp->pseudo_dir = inp->save_dir + 600 ;
    inp->dvscf_file = inp->save_dir + 1200;
    strcpy(inp->pseudo_dir,"");
    strcpy(inp->dvscf_file,"");
    inp->dimension = '3';

}

// function to free usr_input struct data
void free_usr_input(struct usr_input *input)
{   
    free(input->pseudos);
    free(input->save_dir);
    free(input);
    
}

static void Bcast_input_data(struct usr_input *input, int root, MPI_Comm comm)
{
    int mpi_error;

    // all char * will be bcasted in one single call
    mpi_error = MPI_Bcast(input->save_dir, MAX_STR_INPUT, MPI_CHAR, root, comm);
    mpi_error = MPI_Bcast(&input->dimension, 1, MPI_CHAR, root, comm);
    mpi_error = MPI_Bcast(&input->nkpool, 1, MPI_INT, root, comm);
    mpi_error = MPI_Bcast(&input->nqpool, 1, MPI_INT, root, comm);
    mpi_error = MPI_Bcast(&input->start_bnd, 1, MPI_INT, root, comm);
    mpi_error = MPI_Bcast(&input->end_bnd, 1, MPI_INT, root, comm);

}

static int handler(void* user, const char* section, const char* name,
            const char* value)
{
    struct usr_input * inp = user;

    size_t nospace_len =0 ;
    size_t val_len = strlen(value);
    for (size_t i =0; i < val_len; ++i)
    {
        if (value[i] != ' ') ++nospace_len;
    }

    if (nospace_len == 0)
    {
        printf("Invalid input for %s ",name);
        error_msg("Invalid input");
    }

    if (strcmp(section, "input") == 0)
    {
        if (strcmp(name, "nkpool") == 0)
        {
            inp->nkpool =  atoi(value);
        }
        else if (strcmp(name, "nqpool") == 0)
        {
            inp->nqpool =  atoi(value);
        }
        else if (strcmp(name, "start_bnd") == 0)
        {
            inp->start_bnd = atoi(value);
        }
        else if (strcmp(name, "end_bnd") == 0)
        {
            inp->end_bnd = atoi(value);
        }
        else if (strcmp(name, "save_dir") == 0)
        {
            strcpy(inp->save_dir,value);
        }
        else if (strcmp(name, "pseudo_dir") == 0)
        {
            strcpy(inp->pseudo_dir,value);
        }
        else if (strcmp(name, "dvscf_file") == 0)
        {
            strcpy(inp->dvscf_file,value);
        }
        else if (strcmp(name, "pseudos") == 0)
        {   
            // for now we copy to memory starting from inp->save_dir + 1800 
            char * tmp_pseudo = inp->save_dir + 1800;
            // we will later extract from inp->save_dir + 1800
            strcpy(tmp_pseudo,value);
        }
        else if (strcmp(name, "dimension") == 0)
        {   
            int temp_a2i_val = atoi(value);
            if (temp_a2i_val < 2 || temp_a2i_val >3) 
            {
                error_msg("Invalid value for dimension. Only 2 or 3 is supported");
            }
            else 
            {
                inp->dimension = temp_a2i_val + '0';
            }
        }
        else
        {   
            error_msg("Invalid variable in input file.");
        }
    }    
    else
    {   
        //error_msg("Invalid variable in input file");
        error_msg("Invalid input file. input file must contain [input] at the top."); 
    }
    return 1;
}


void read_input_file(const char * input_file, struct usr_input ** input_data)
{   
    // input_data must be free outside

    init_usr_input(input_data);

    int my_rank, mpi_error;
    mpi_error = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0)
    {
        if (ini_parse(input_file, handler, *input_data) < 0) {
            error_msg("Cannot open input file input file.");
        }
        if (strlen((*input_data)->pseudo_dir) == 0 ) error_msg("Pseudo Directory not found input");
        if (strlen((*input_data)->dvscf_file) == 0 ) error_msg("dVscf file not found in input");
    }
    // broad cast;
    Bcast_input_data(*input_data, 0, MPI_COMM_WORLD);


    // parse the pseudo
    char * pseudo_parse = strtok ((*input_data)->save_dir + 1800,",");
    
    int ipot = 0;
    while (pseudo_parse != NULL) {
        while(*pseudo_parse == ' ') ++pseudo_parse;
        if (strlen(pseudo_parse) < 1) continue;
        (*input_data)->pseudos[ipot] = pseudo_parse;
        pseudo_parse = strtok (NULL, ",");
        ++ipot;
    }

    

}







