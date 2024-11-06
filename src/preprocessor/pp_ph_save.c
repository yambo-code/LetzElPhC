/*
This file contains functions which are os dependent.
*/
#include <errno.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "../common/cwalk/cwalk.h"
#include "../common/error.h"
#include "../common/string_func.h"
#include "../elphC.h"
#include "../io/ezxml/ezxml.h"
#include "ELPH_copy.h"
#include "preprocessor.h"

#if defined(_WIN32)
#include <direct.h>
#define mkdir(path, mode) _mkdir(path)
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

static ND_int find_nqpools(const char* out_dir, char* buffer_tmp,
                           ND_int buffer_size);

void create_ph_save_dir_pp_qe(const char* inp_file)
{
    //  parse the input file
    // open the  qe ph.x input file

    char PH_SAVE_DIR_NAME[ELPH_MAX_ENV_SIZE];

    // check if env exits for ph_save
    char* env_var_tmp = getenv("ELPH_PH_SAVE_DIR");
    if (env_var_tmp && strlen(env_var_tmp) > 0)
    {
        if (strlen(env_var_tmp) > (ELPH_MAX_ENV_SIZE - 1))
        {
            fprintf(stderr,
                    "Warning : length of ELPH_PH_SAVE_DIR environment "
                    "variable must be strictly less than %d\n",
                    (int)ELPH_MAX_ENV_SIZE);
        }
        snprintf(PH_SAVE_DIR_NAME, ELPH_MAX_ENV_SIZE, "%s", env_var_tmp);
    }
    else
    {
        strncpy_custom(PH_SAVE_DIR_NAME, PH_SAVE_DIR_NAME_DEFAULT,
                       ELPH_MAX_ENV_SIZE);
    }
    //
    FILE* fp = fopen(inp_file, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Unable to open given ph.x input file : %s \n",
                inp_file);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    char* read_buf = malloc(4 * PH_X_INP_READ_BUF_SIZE);
    CHECK_ALLOC(read_buf);

    char* key_str = read_buf + PH_X_INP_READ_BUF_SIZE;
    char* val_str = read_buf + 2 * PH_X_INP_READ_BUF_SIZE;
    char* tmp_buf = read_buf + 3 * PH_X_INP_READ_BUF_SIZE;

    char* inputs_vals = malloc(5 * PH_X_INP_READ_BUF_SIZE);
    CHECK_ALLOC(inputs_vals);

    char* out_dir = inputs_vals;
    char* dyn_prefix = inputs_vals + PH_X_INP_READ_BUF_SIZE;
    char* dvscf_prefix = inputs_vals + 2 * PH_X_INP_READ_BUF_SIZE;
    char* drho_prefix = inputs_vals + 3 * PH_X_INP_READ_BUF_SIZE;
    char* scf_prefix = inputs_vals + 4 * PH_X_INP_READ_BUF_SIZE;

    // set defaults
    bool ldisp = false;
    bool elph_yambo = false;
    // this is true if electron_phonon is set to either 'yambo' or 'dvscf' or
    // 'Wannier'

    env_var_tmp = getenv("ESPRESSO_TMPDIR");
    if (env_var_tmp)
    {
        strncpy_custom(out_dir, env_var_tmp, PH_X_INP_READ_BUF_SIZE);
    }
    else
    {
        strcpy(out_dir, "./");
    }

    strcpy(dyn_prefix, "matdyn");
    strcpy(dvscf_prefix, "");
    strcpy(drho_prefix, "");
    strcpy(scf_prefix, "pwscf");

    // ===== read stuff =====
    while (fgets(read_buf, PH_X_INP_READ_BUF_SIZE, fp))
    {
        // remove comments
        str_replace_chars(read_buf, ",'\"!", "   \0");

        if (strlen(read_buf) == 0)
        {
            continue;
        }
        // now read key
        char* token = strtok(read_buf, "=");
        strncpy_custom(key_str, token, PH_X_INP_READ_BUF_SIZE);
        // lower case the key
        lowercase_str(key_str);
        //  read value
        token = strtok(NULL, "=");
        if (token)
        {
            strncpy_custom(val_str, token, PH_X_INP_READ_BUF_SIZE);
        }
        else
        {
            // line does not contain key value
            continue;
        }

        // remove spaces
        sscanf(key_str, "%s", tmp_buf);
        strncpy_custom(key_str, tmp_buf, PH_X_INP_READ_BUF_SIZE);

        sscanf(val_str, "%s", tmp_buf);
        strncpy_custom(val_str, tmp_buf, PH_X_INP_READ_BUF_SIZE);

        if (!strcmp(key_str, "ldisp"))
        {
            // lowercase val_str
            lowercase_str(val_str);

            if (!strcmp(val_str, ".true.") || !strcmp(val_str, "1") ||
                !strcmp(val_str, "t"))
            {
                ldisp = true;
            }
        }
        //
        else if (!strcmp(key_str, "outdir"))
        {
            strncpy_custom(out_dir, val_str, PH_X_INP_READ_BUF_SIZE);
        }
        //
        else if (!strcmp(key_str, "fildyn"))
        {
            strncpy_custom(dyn_prefix, val_str, PH_X_INP_READ_BUF_SIZE);
        }
        //
        else if (!strcmp(key_str, "fildvscf"))
        {
            strncpy_custom(dvscf_prefix, val_str, PH_X_INP_READ_BUF_SIZE);
        }
        //
        else if (!strcmp(key_str, "fildrho"))
        {
            strncpy_custom(drho_prefix, val_str, PH_X_INP_READ_BUF_SIZE);
        }
        //
        else if (!strcmp(key_str, "prefix"))
        {
            strncpy_custom(scf_prefix, val_str, PH_X_INP_READ_BUF_SIZE);
        }
        //
        else if (!strcmp(key_str, "electron_phonon"))
        {
            // lowercase val_str
            lowercase_str(val_str);

            if (!strcmp(val_str, "yambo") || !strcmp(val_str, "dvscf") ||
                !strcmp(val_str, "wannier"))
            {
                elph_yambo = true;
            }
        }
    }
    fclose(fp);

    if (strlen(dvscf_prefix) == 0)
    {
        fprintf(stderr,
                "Error : fildvscf is not set in ph.x input. "
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

    if (strstr(dyn_prefix, ".xml"))
    {
        fprintf(stderr,
                "Error : xml format for dyn files not yet supported.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!ldisp)
    {
        fprintf(stderr,
                "Warning : ldisp is set .false. in the ph.x calculation. "
                "Make sure dyn0 file is present with k-point compatible "
                "q-point grid.\n");
    }

    char* src_file_tmp = read_buf;
    char* dest_file_tmp = read_buf + PH_X_INP_READ_BUF_SIZE;
    char* tmp_file_buf = read_buf + 2 * PH_X_INP_READ_BUF_SIZE;
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
    snprintf(prefix_dir, 100, "%s.save", scf_prefix);

    cwk_path_join_multiple(
        (const char*[]){out_dir, prefix_dir, "data-file-schema.xml", NULL},
        src_file_tmp, PH_X_INP_READ_BUF_SIZE);

    cwk_path_join_multiple(
        (const char*[]){PH_SAVE_DIR_NAME, "data-file-schema.xml", NULL},
        dest_file_tmp, PH_X_INP_READ_BUF_SIZE);

    if (0 != copy_files(src_file_tmp, dest_file_tmp))
    {
        printf("Warning : Error copying file %s.\n", src_file_tmp);
    }

    // 2) copy pseudo pots files
    // since src_file_tmp store file path of data-file-schema.xml we
    // can read the pseudo pots used
    FILE* fp_xml = fopen(src_file_tmp, "r");
    if (!fp_xml)
    {
        error_msg("Error opening data-file-schema.xml file");
    }

    ezxml_t qexml = ezxml_parse_fp(fp_xml);
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
        char* tmp_str =
            ezxml_get(atom_specs, "species", itype, "pseudo_file", -1)->txt;

        cwk_path_join_multiple(
            (const char*[]){out_dir, prefix_dir, tmp_str, NULL}, src_file_tmp,
            PH_X_INP_READ_BUF_SIZE);

        cwk_path_join_multiple((const char*[]){PH_SAVE_DIR_NAME, tmp_str, NULL},
                               dest_file_tmp, PH_X_INP_READ_BUF_SIZE);

        if (0 != copy_files(src_file_tmp, dest_file_tmp))
        {
            printf("Warning : Error copying file %s.\n", src_file_tmp);
        }
    }
    ezxml_free(qexml);
    fclose(fp_xml);

    // 3) copy dyn files

    snprintf(src_file_tmp, PH_X_INP_READ_BUF_SIZE, "%s0", dyn_prefix);
    // first copy the dyn0 file
    cwk_path_join_multiple((const char*[]){PH_SAVE_DIR_NAME, "dyn0", NULL},
                           dest_file_tmp, PH_X_INP_READ_BUF_SIZE);

    if (0 != copy_files(src_file_tmp, dest_file_tmp))
    {
        printf("Warning : Error copying file %s.\n", src_file_tmp);
    }

    // open dyn0 file to get number of qpoints
    FILE* fp_dyn0 = fopen(src_file_tmp, "r");

    if (fp_dyn0 == NULL)
    {
        fprintf(stderr,
                "Opening file %s failed. "
                "This could be due to ldisp not set to .true. \n",
                src_file_tmp);

        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // dummy
    fgets(src_file_tmp, PH_X_INP_READ_BUF_SIZE, fp_dyn0);

    int ndyn;
    // read number of dyn files
    fgets(src_file_tmp, PH_X_INP_READ_BUF_SIZE, fp_dyn0);

    if (sscanf(src_file_tmp, "%d", &ndyn) != 1)
    {
        error_msg("Error reading 2nd line from dyn0 file");
    }
    fclose(fp_dyn0);

    for (int idyn = 0; idyn < ndyn; ++idyn)
    {
        snprintf(src_file_tmp, PH_X_INP_READ_BUF_SIZE, "%s%d", dyn_prefix,
                 idyn + 1);

        char dyn_tmp_str[32];
        snprintf(dyn_tmp_str, 32, "dyn%d", idyn + 1);

        cwk_path_join_multiple(
            (const char*[]){PH_SAVE_DIR_NAME, dyn_tmp_str, NULL}, dest_file_tmp,
            PH_X_INP_READ_BUF_SIZE);

        if (0 != copy_files(src_file_tmp, dest_file_tmp))
        {
            printf("Warning : Error copying file %s.\n", src_file_tmp);
        }
    }

    // Find the number of qpools used in ph.x calculation
    ND_int ph_x_qpools =
        find_nqpools(out_dir, tmp_file_buf, PH_X_INP_READ_BUF_SIZE);

    // 3) copy dvscf files
    for (int idyn = 0; idyn < ndyn; ++idyn)
    {
        char dvscf_tmp_str[32];
        // set dyn suffix number
        // if elph_yambo == true then iq_ is added
        if (elph_yambo)
        {
            snprintf(dvscf_tmp_str, 32, "%d_1", idyn + 1);
        }
        else
        {
            strcpy(dvscf_tmp_str, "1");
        }
        // set the source file
        for (ND_int iphx_pool = 0; iphx_pool < ph_x_qpools; ++iphx_pool)
        {
            char ph0_tmp[16];
            snprintf(ph0_tmp, sizeof(ph0_tmp) - 1, "_ph%d", (int)iphx_pool);
            if (idyn)
            {
                // use dest_file_tmp and tmp_file_buf as tmp buffer string
                snprintf(dest_file_tmp, PH_X_INP_READ_BUF_SIZE, "%s.q_%d",
                         scf_prefix, idyn + 1);
                snprintf(tmp_file_buf, PH_X_INP_READ_BUF_SIZE, "%s.%s%s",
                         scf_prefix, dvscf_prefix, dvscf_tmp_str);

                cwk_path_join_multiple(
                    (const char*[]){out_dir, ph0_tmp, dest_file_tmp,
                                    tmp_file_buf, NULL},
                    src_file_tmp, PH_X_INP_READ_BUF_SIZE);
            }
            else
            {
                snprintf(tmp_file_buf, PH_X_INP_READ_BUF_SIZE, "%s.%s%s",
                         scf_prefix, dvscf_prefix, dvscf_tmp_str);

                cwk_path_join_multiple(
                    (const char*[]){out_dir, ph0_tmp, tmp_file_buf, NULL},
                    src_file_tmp, PH_X_INP_READ_BUF_SIZE);
            }

            snprintf(tmp_file_buf, PH_X_INP_READ_BUF_SIZE - 1, "dvscf%d",
                     idyn + 1);

            // set the destination
            cwk_path_join_multiple(
                (const char*[]){PH_SAVE_DIR_NAME, tmp_file_buf, NULL},
                dest_file_tmp, PH_X_INP_READ_BUF_SIZE);

            if (0 == copy_files(src_file_tmp, dest_file_tmp))
            {
                break;
            }
        }
    }

    // 4) copy drho
    // currently not needed

    // 5) copy pattern files
    for (int idyn = 0; idyn < ndyn; ++idyn)
    {
        char dyn_tmp_str[32];

        snprintf(dest_file_tmp, PH_X_INP_READ_BUF_SIZE, "%s.phsave",
                 scf_prefix);
        snprintf(dyn_tmp_str, 32, "patterns.%d.xml", idyn + 1);

        cwk_path_join_multiple(
            (const char*[]){out_dir, "_ph0", dest_file_tmp, dyn_tmp_str, NULL},
            src_file_tmp, PH_X_INP_READ_BUF_SIZE);

        cwk_path_join_multiple(
            (const char*[]){PH_SAVE_DIR_NAME, dyn_tmp_str, NULL}, dest_file_tmp,
            PH_X_INP_READ_BUF_SIZE);

        if (0 != copy_files(src_file_tmp, dest_file_tmp))
        {
            printf("Warning : Error copying file %s.\n", src_file_tmp);
        }
    }

    free(read_buf);
    free(inputs_vals);
}

static ND_int find_nqpools(const char* out_dir, char* buffer_tmp,
                           ND_int buffer_size)
{
    // find number of qpools used in ph.x
    ND_int counter = 0;

    char ph0_tmp[16];

    while (true)
    {
        snprintf(ph0_tmp, sizeof(ph0_tmp) - 1, "_ph%d", (int)counter);
        cwk_path_join_multiple((const char*[]){out_dir, ph0_tmp, NULL},
                               buffer_tmp, buffer_size);
        if (0 != check_dir_exists(buffer_tmp))
        {
            // issue with dir (either cannot access or does not exist)
            break;
        }
        ++counter;
    }
    return counter;
}