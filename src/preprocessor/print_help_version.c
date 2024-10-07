/*
 * This file contains functions to print the helper information and
 * version number
 *
 */
#include "../elphC.h"
#include "preprocessor.h"
#include <stdio.h>

void ELPH_print_version(void)
{
    // we print to stdout
    fprintf(stdout, "LetzElPhC, ");
    fprintf(stdout, "Version : Beta\n");
    fprintf(stdout, "Build details :\n");
// print information about the precession
#if defined(COMPILE_ELPH_DOUBLE)
    fprintf(stdout, "PRECESSION : Double\n");
#else
    fprintf(stdout, "PRECESSION : Single\n");
#endif
// print information about the openmp
#if defined(ELPH_OMP_PARALLEL_BUILD)
    fprintf(stdout, "OPENMP : True\n");
#else
    fprintf(stdout, "OPENMP : False\n");
#endif
}

void ELPH_print_help(void)
{
    fprintf(stdout, "LetzElPhC help\n");
    fprintf(stdout, "--help     : to print the help\n");
    fprintf(stdout, "--version  : to print the version\n");
    fprintf(stdout, "-F         : Input file.\n");
    fprintf(stdout, "--code     : DFT code. currently only \"qe\" is supported."
                    "Default is \"qe\"\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "-pp        : activates preprocessor. used for ph_save creation.\n"
                    "If -pp is absent, the code assumes that user requested elph calculation.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Example: \n");
    fprintf(stdout, "To create ph_save folder, run the folowing command:\n");
    fprintf(stdout, "lelphc -pp -F DFPT_input_file\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "To perform the elph calculation, an example command looks like:\n");
    fprintf(stdout, "mpirun -n 4 lelphc -F ELPH_input_file\n");

    fprintf(stdout, "\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Set the following environment variables to"
                    " change the default behaviour of preprocessor\n");

    fprintf(stdout, "\n");
    fprintf(stdout, "ELPH_COPY_CMD : Default (\"cp\") In case if you want to \n"
                    "create a softlink to dfpt data, use can set this varibale\n");
    fprintf(stdout, "To use symlinks (instead of copy), use export ELPH_COPY_CMD=\"ln -sr\". Symlinks work only on linux \n");
    fprintf(stdout, "\n");
    fprintf(stdout, "ELPH_SAVE_DIR : Default (ph_save). In case if you want to"
                    "change the name of the ph_save folder \n");
    fprintf(stdout, " Use export ELPH_SAVE_DIR=ph_save_name \n");

    fprintf(stdout, "\n");
}
