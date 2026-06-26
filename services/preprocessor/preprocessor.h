#pragma once

#include "common/dtypes.h"
#include "elphC.h"

#define PH_X_INP_READ_BUF_SIZE 512
#define ELPH_MAX_ENV_SIZE 64
#define PH_SAVE_DIR_NAME_DEFAULT "ph_save"

void ELPH_cli_parser(int argc, char* argv[], struct calc_details* calc_info);

#if defined _Y6_LETZ
void create_ph_save_dir_pp_qe(const char* inp_file,char* ph_path,char* scf_path);
#else
void create_ph_save_dir_pp_qe(const char* inp_file,const char* ph_path, const char* scf_path);
#endif

void ELPH_print_version(void);
void ELPH_print_help(void);
