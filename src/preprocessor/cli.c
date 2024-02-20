#include "../elphC.h"

enum calc_type
{
    CALC_ELPH, // initiate elph calculation
    CALC_PH_SAVE_CREATE, // preprocess (creating ph_save_dir)
    CALC_HELP, // help
    CALC_VERSION // print version
};

enum calc_type ELPH_calc_type(int narg, const char ** argv)
{
    // first check if the use requested help
    enum calc_type ctype = CALC_ELPH;

    for (int i = 0; i<narg; ++i)
    {   
        const char * cli_str = argv[i];
        if (strstr(cli_str,"--help"))
        {   
            ctype = CALC_HELP;
            break;
        }
        else if (strstr(cli_str,"-pp"))
        {
            ctype = CALC_PH_SAVE_CREATE;
            break;
        }
        else if (strstr(cli_str,"--version"))
        {
            ctype = CALC_VERSION;
            break;
        }
    }
    return ctype;
}







