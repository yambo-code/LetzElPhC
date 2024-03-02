#pragma once

#include "../elphC.h"
#include "../common/dtypes.h"

#define PH_X_INP_READ_BUF_SIZE 512
#define ELPH_MAX_ENV_SIZE 64
#define PH_SAVE_DIR_NAME_DEFAULT "ph_save"

enum calc_type
{
    CALC_ELPH, // initiate elph calculation
    CALC_PH_SAVE_CREATE, // preprocess (creating ph_save_dir)
    CALC_HELP, // help
    CALC_VERSION // print version
};

struct calc_details
{
    enum calc_type calc;
    enum ELPH_dft_code code;
    char input_file[512];
    // name of the input file.
    // can be elph input or DFT-Phonon input
    // base on calc_type
};

void ELPH_cli_parser(int argc, char* argv[], struct calc_details* calc_info);

void create_ph_save_dir_pp_qe(const char* inp_file);

void ELPH_print_version();
void ELPH_print_help();
