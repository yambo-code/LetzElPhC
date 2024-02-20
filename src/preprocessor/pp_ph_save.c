/*
This file contains functions which are system dependent
*/
#include "../elphC.h"
#include "../common/string_func.h"
#include "../io/ezxml/ezxml.h"
#include "../common/error.h"

#define PH_X_INP_READ_BUF_SIZE 512
#define ELPH_MAX_ENV_SIZE 64
#define PH_SAVE_DIR_NAME_DEFAULT "ph_save"


#if defined(_WIN32)
    #include <direct.h>
    #define mkdir(path, mode) _mkdir(path)
    #define ELPH_COPY_CMD_DEFAULT "copy"
#else
    #include <sys/types.h>
    #include <sys/stat.h>
    #define ELPH_COPY_CMD_DEFAULT "cp"
#endif

void create_ph_save_dir_pp_qe(const char * inp_file)
{   

    // first check if system can find a command processor
    if (!system(NULL))
    {
        fprintf(stderr, "System cannot find a command processor.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    //  parse the input file
    // open the  qe ph.x input file

    char ELPH_COPY_CMD[ELPH_MAX_ENV_SIZE];
    char PH_SAVE_DIR_NAME[ELPH_MAX_ENV_SIZE];

    // check of env exits for copy cmd and ph_save
    char* env_var_tmp = getenv("ELPH_COPY_CMD");
    if (env_var_tmp && strlen(env_var_tmp) > 0 )
    {   
        if (strlen(env_var_tmp) > (ELPH_MAX_ENV_SIZE-1))
        {
            fprintf(stderr, "Warning : length of ELPH_COPY_CMD environment "
                    "variable must be strictly less than %d\n", (int)ELPH_MAX_ENV_SIZE);
        }
        snprintf(ELPH_COPY_CMD, ELPH_MAX_ENV_SIZE, "%s", env_var_tmp);
    }
    else
    {
        strcpy(ELPH_COPY_CMD,ELPH_COPY_CMD_DEFAULT);
    }
    //
    //
    env_var_tmp = getenv("ELPH_SAVE_DIR");
    if (env_var_tmp && strlen(env_var_tmp) > 0 )
    {   
        if (strlen(env_var_tmp) > (ELPH_MAX_ENV_SIZE-1))
        {
            fprintf(stderr, "Warning : length of ELPH_SAVE_DIR environment "
                    "variable must be strictly less than %d\n",(int)ELPH_MAX_ENV_SIZE);
            
        }
        snprintf(PH_SAVE_DIR_NAME, ELPH_MAX_ENV_SIZE, "%s", env_var_tmp);
    }
    else
    {
        strcpy(PH_SAVE_DIR_NAME, PH_SAVE_DIR_NAME_DEFAULT);
    }
    //
    FILE* fp = fopen(inp_file, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Unable to open given ph.x input file : %s \n", inp_file);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    char * read_buf = malloc(4*PH_X_INP_READ_BUF_SIZE);
    char * key_str  = read_buf + PH_X_INP_READ_BUF_SIZE;
    char * val_str  = read_buf + 2*PH_X_INP_READ_BUF_SIZE;
    char * tmp_buf  = read_buf + 3*PH_X_INP_READ_BUF_SIZE;

    char * inputs_vals = malloc(5*PH_X_INP_READ_BUF_SIZE);

    char * out_dir      = inputs_vals;
    char * dyn_prefix   = inputs_vals + PH_X_INP_READ_BUF_SIZE;
    char * dvscf_prefix = inputs_vals + 2*PH_X_INP_READ_BUF_SIZE;
    char * drho_prefix  = inputs_vals + 3*PH_X_INP_READ_BUF_SIZE;
    char * scf_prefix   = inputs_vals + 4*PH_X_INP_READ_BUF_SIZE;
    
    // set defaults
    bool ldisp = false;
    bool elph_yambo = false; 
    // this is true if electron_phonon is set to either 'yambo' or 'dvscf' or 'Wannier'
    
    env_var_tmp = getenv("ESPRESSO_TMPDIR");
    if (env_var_tmp)
    {
        strcpy(out_dir,env_var_tmp);
    }
    else
    {
        strcpy(out_dir,"./");
    }

    strcpy(dyn_prefix,"matdyn");
    strcpy(dvscf_prefix,"");
    strcpy(drho_prefix,"");
    strcpy(scf_prefix,"pwscf");

    // ===== read stuff =====
    while (fgets(read_buf, PH_X_INP_READ_BUF_SIZE, fp))
    {   
        // remove comments
        str_replace_chars(read_buf, ",'\"!", "   \0");
        // lower case all the chars
        for (char* p = read_buf; *p; ++p)
        {
            *p = tolower(*p);
        }
        if (strlen(read_buf) == 0)
        {
            continue;
        }
        // now read key
        char* token = strtok(read_buf, "=");
        strcpy(key_str,token);

        //  read value
        token = strtok(NULL, "=");
        if (token)
        {
            strcpy(val_str,token);
        }
        else
        {   
            // line does not contain key value
            continue;
        }

        // remove spaces
        sscanf(key_str,"%s",tmp_buf);
        strcpy(key_str,tmp_buf);

        sscanf(val_str,"%s",tmp_buf);
        strcpy(val_str,tmp_buf);

        if (!strcmp(key_str,"ldisp"))
        {
            if (!strcmp(val_str,".true.") || !strcmp(val_str,"1") || !strcmp(val_str,"t"))
            {
                ldisp = true;
            }
        }
        //
        else if (!strcmp(key_str,"outdir"))
        {
            strcpy(out_dir,val_str);
        }
        //
        else if (!strcmp(key_str,"fildyn"))
        {
            strcpy(dyn_prefix,val_str);
        }
        //
        else if (!strcmp(key_str,"fildvscf"))
        {
            strcpy(dvscf_prefix,val_str);
        }
        //
        else if (!strcmp(key_str,"fildrho"))
        {
            strcpy(drho_prefix,val_str);
        }
        //
        else if (!strcmp(key_str,"prefix"))
        {
            strcpy(scf_prefix,val_str);
        }
        //
        else if (!strcmp(key_str,"electron_phonon"))
        {
            if (!strcmp(val_str,"yambo") || !strcmp(val_str,"dvscf") || !strcmp(val_str,"wannier"))
            {
                elph_yambo = true;
            }
        }
    }
    fclose(fp);
    
    if (strlen(dvscf_prefix) == 0)
    {
        fprintf(stderr,"Error : fildvscf is not set in ph.x input. "
        "Change in potential is not stored. "
        "You must rerun the ph.x calculation again with "
        "fildvscf flag set in ph.x input\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // if (strlen(drho_prefix) == 0)
    // {
    //     fprintf(stderr,"Error : fildrho is not set in ph.x input. "
    //     "Change in potential is not stored. "
    //     "You must rerun the ph.x calculation again with "
    //     "fildrho flag set in ph.x input\n");
    //     MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    // }
    
    if (strstr(dyn_prefix,".xml"))
    {
        fprintf(stderr,"Error : xml format for dyn files not yet supported.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!ldisp)
    {
        fprintf(stderr,"Warning : ldisp is set .false. in the ph.x calculation. "
        "Make sure dyn0 file is present with k-point compatible q-point grid.\n");
    }

    char * src_file_tmp  = read_buf;
    char * dest_file_tmp = read_buf + PH_X_INP_READ_BUF_SIZE;
    char * exe_cmd_tmp   = read_buf + 2*PH_X_INP_READ_BUF_SIZE;
    // from now we use read buffer as file_name_buf

    // now create the ph_save directory
    if (mkdir(PH_SAVE_DIR_NAME, 0777))
    {
        if (errno != EEXIST)
        {
            error_msg("Failed to create ph_save directory.");
        }
    }
    
    char prefix_dir[100];
    // 1) copy data-file-schema.xml
    snprintf(prefix_dir,100,"%s.save",scf_prefix);

    cwk_path_join_multiple( (const char * []){out_dir,prefix_dir, 
                            "data-file-schema.xml",NULL}, 
                            src_file_tmp, PH_X_INP_READ_BUF_SIZE);

    cwk_path_join_multiple( (const char * []){PH_SAVE_DIR_NAME, 
                            "data-file-schema.xml",NULL}, 
                            dest_file_tmp, PH_X_INP_READ_BUF_SIZE);
    
    snprintf(exe_cmd_tmp, 2*PH_X_INP_READ_BUF_SIZE, "%s %s %s",ELPH_COPY_CMD,src_file_tmp,dest_file_tmp);
    system(exe_cmd_tmp);

    // 2) copy pseudo pots files
    // since src_file_tmp store file path of data-file-schema.xml we 
    // can read the pseudo pots used
    FILE* fp_xml = fopen(src_file_tmp, "r");
    if (!fp_xml)
    {
        error_msg("Error opening data-file-schema.xml file");
    }

    ezxml_t qexml = ezxml_parse_fp(fp);
    if (qexml == NULL)
    {
        error_msg("Error parsing data-file-schema.xml file");
    }

    ezxml_t atom_specs = ezxml_get(qexml, "input", 0, "atomic_species", -1);
    if (atom_specs == NULL)
    {
        error_msg("error reading atomic spices from data-file-schema.xml file");
    }

    ND_int ntype = atoll(ezxml_attr(atom_specs, "ntyp"));

    for (ND_int itype = 0; itype < ntype; ++itype)
    {
        char * tmp_str = ezxml_get(atom_specs, "species", itype, "pseudo_file", -1)->txt;

        cwk_path_join_multiple( (const char * []){out_dir,prefix_dir, 
                            tmp_str,NULL}, src_file_tmp, 
                            PH_X_INP_READ_BUF_SIZE);

        cwk_path_join_multiple( (const char * []){PH_SAVE_DIR_NAME, 
                            tmp_str,NULL}, dest_file_tmp,
                            PH_X_INP_READ_BUF_SIZE);
        snprintf(exe_cmd_tmp, 2*PH_X_INP_READ_BUF_SIZE, "%s %s %s",ELPH_COPY_CMD,src_file_tmp,dest_file_tmp);

        system(exe_cmd_tmp);
    }
    ezxml_free(qexml);
    fclose(fp_xml);

    // 3) copy dyn files

    snprintf(src_file_tmp, PH_X_INP_READ_BUF_SIZE, "%s0",dyn_prefix);
    // first copy the dyn0 file
    cwk_path_join_multiple( (const char * []){PH_SAVE_DIR_NAME, 
                            "dyn0",NULL}, dest_file_tmp,
                            PH_X_INP_READ_BUF_SIZE);
    
    snprintf(exe_cmd_tmp, 2*PH_X_INP_READ_BUF_SIZE, "%s %s %s",ELPH_COPY_CMD,src_file_tmp,dest_file_tmp);

    system(exe_cmd_tmp);

    // open dyn0 file to get number of qpoints
    FILE* fp_dyn0 = fopen(src_file_tmp, "r");

    if (fp_dyn0 == NULL)
    {
        fprintf(stderr, "Opening file %s failed. "
        "This could be due to ldisp not set to .true. \n", src_file_tmp);

        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // dummy
    fgets(src_file_tmp, PH_X_INP_READ_BUF_SIZE, fp);

    int ndyn;
    // read number of dyn files
    fgets(src_file_tmp, PH_X_INP_READ_BUF_SIZE, fp);

    if (sscanf(src_file_tmp, "%d", &ndyn) != 1)
    {
        error_msg("Error reading 2nd line from dyn0 file");
    }
    fclose(fp_dyn0);

    for (int idyn = 0; idyn < ndyn; ++idyn)
    {
        snprintf(src_file_tmp, PH_X_INP_READ_BUF_SIZE, "%s%d", dyn_prefix, idyn+1);

        char dyn_tmp_str[32];
        snprintf(dyn_tmp_str, 32, "dyn%d", idyn+1);

        cwk_path_join_multiple( (const char * []){PH_SAVE_DIR_NAME, 
                            dyn_tmp_str,NULL}, dest_file_tmp,
                            PH_X_INP_READ_BUF_SIZE);
        
        snprintf(exe_cmd_tmp, 2*PH_X_INP_READ_BUF_SIZE, "%s %s %s",ELPH_COPY_CMD,src_file_tmp,dest_file_tmp);
        system(exe_cmd_tmp);
    }

    // 3) copy dvscf files
    for (int idyn = 0; idyn < ndyn; ++idyn)
    {   
        char dvscf_tmp_str[32];
        // set dyn suffix number
        // if elph_yambo == true then iq_ is added
        if (elph_yambo)
        {
            snprintf(dvscf_tmp_str, 32, "%d_1", idyn+1);
        }
        else
        {
            strcpy(dvscf_tmp_str,"1");
        }
        // set the source file
        if (idyn)
        {   
            // use dest_file_tmp and exe_cmd_tmp as tmp buffer string
            snprintf(dest_file_tmp, PH_X_INP_READ_BUF_SIZE, "%s.q_%d", scf_prefix, idyn+1);
            snprintf(exe_cmd_tmp, PH_X_INP_READ_BUF_SIZE, "%s.%s%s", 
                    scf_prefix, dvscf_prefix, dvscf_tmp_str);

            cwk_path_join_multiple( (const char * []){out_dir,"_ph0", 
                            dest_file_tmp, exe_cmd_tmp, NULL}, src_file_tmp, 
                            PH_X_INP_READ_BUF_SIZE);
        }
        else
        {
            snprintf(exe_cmd_tmp, PH_X_INP_READ_BUF_SIZE, "%s.%s%s", 
                    scf_prefix, dvscf_prefix, dvscf_tmp_str);
            
            cwk_path_join_multiple( (const char * []){out_dir,"_ph0", 
                            exe_cmd_tmp, NULL}, src_file_tmp, 
                            PH_X_INP_READ_BUF_SIZE);
        }

        snprintf(dvscf_tmp_str, 32, "dvscf%d", idyn+1);

        // set the destination
        cwk_path_join_multiple( (const char * []){PH_SAVE_DIR_NAME, 
                            dvscf_tmp_str,NULL}, dest_file_tmp,
                            PH_X_INP_READ_BUF_SIZE);
        
        snprintf(exe_cmd_tmp, 2*PH_X_INP_READ_BUF_SIZE, "%s %s %s",ELPH_COPY_CMD,src_file_tmp,dest_file_tmp);
        system(exe_cmd_tmp);
    }

    // 4) copy drho
    // currently not needed

    // 5) copy pattern files
    for (int idyn = 0; idyn < ndyn; ++idyn)
    {   
        char dyn_tmp_str[32];

        snprintf(dest_file_tmp, PH_X_INP_READ_BUF_SIZE, "%s.phsave", scf_prefix);
        snprintf(dyn_tmp_str, 32, "patterns.%d.xml", idyn+1);

        cwk_path_join_multiple( (const char * []){out_dir,"_ph0", 
                            dest_file_tmp, dyn_tmp_str, NULL}, src_file_tmp, 
                            PH_X_INP_READ_BUF_SIZE);

        cwk_path_join_multiple( (const char * []){PH_SAVE_DIR_NAME, 
                            dyn_tmp_str,NULL}, dest_file_tmp,
                            PH_X_INP_READ_BUF_SIZE);
        
        snprintf(exe_cmd_tmp, 2*PH_X_INP_READ_BUF_SIZE, "%s %s %s",ELPH_COPY_CMD,src_file_tmp,dest_file_tmp);
        system(exe_cmd_tmp);
    }
    
    free(read_buf);
    free(inputs_vals);

}


