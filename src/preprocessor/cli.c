#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../elphC.h"
#include "ELPH_getopt.h"
#include "preprocessor.h"

void ELPH_cli_parser(int argc, char* argv[], struct calc_details* calc_info)
{
    // set default
    calc_info->calc = CALC_ELPH;
    calc_info->code = DFT_CODE_QE;
    strcpy(calc_info->input_file, "");
    /*
     * Here are the options
     * --help (help)
     * --version (version)
     * -pp (run preprocessor)
     * --code=qe
     * -F file name
     */
    bool help_cmd = false;
    bool ver_cmd = false;
    bool pp_cmd = false;

    const struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'v'},
        {"pp", no_argument, NULL, 'p'},
        {"code", required_argument, NULL, 'c'},
        {"F", required_argument, NULL, 'f'}};

    int ch;
    int longindex;
    while ((ch = ELPH_getopt_long_only(argc, argv, "", long_options,
                                       &longindex)) != -1)
    {
        switch (ch)
        {
            case 'h':
                help_cmd = true;
                break;
            case 'v':
                ver_cmd = true;
                break;
            case 'p':
                pp_cmd = true;
                break;
            case 'c':
                // get the code
                if (strstr(optarg, "qe"))
                {
                    calc_info->code = DFT_CODE_QE;
                }
                else
                {
                    fprintf(stderr, "only --code=qe is supported\n");
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }
                break;
            case 'f':
                strncpy(calc_info->input_file, optarg, 512 - 1);
                break;
            case '?':
                fprintf(stderr, "Unsupported argument given.\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    if (help_cmd)
    {
        calc_info->calc = CALC_HELP;
        return;
    }
    else if (ver_cmd)
    {
        calc_info->calc = CALC_VERSION;
        return;
    }
    else if (pp_cmd)
    {
        calc_info->calc = CALC_PH_SAVE_CREATE;
        return;
    }
    else
    {
        calc_info->calc = CALC_ELPH;
        return;
    }
}
