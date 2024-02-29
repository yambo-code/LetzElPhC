/*
 * This is the starting point of the program
*/

#include "elphC.h"
#include "preprocessor/preprocessor.h"
#include "elph/elph.h"

int main(int argc, char* argv[])
{
#if defined(ELPH_OMP_PARALLEL_BUILD)
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
    MPI_Init(&argc, &argv);
#endif
    
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    struct calc_details * calc_info = malloc(sizeof(struct calc_details)); 

    ND_int dft_code_tmp;
    // this is used to broad cast the enum
    bool run_elph_calc = false;
    if (my_rank == 0)
    {
        ELPH_cli_parser(argc, argv, calc_info);
        
        dft_code_tmp = calc_info->code;

        if (calc_info->calc == CALC_ELPH)
        {
            run_elph_calc = true;
        }
    }
    
    MPI_Bcast(&dft_code_tmp,  1, ELPH_MPI_ND_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&run_elph_calc, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(calc_info->input_file, sizeof(calc_info->input_file), MPI_CHAR, 0, MPI_COMM_WORLD);
    
    calc_info->code = dft_code_tmp;
    
    if (run_elph_calc)
    {
        elph_driver(calc_info->input_file, calc_info->code, MPI_COMM_WORLD);
    }
    else 
    {
        // run preprocessor
        if (!my_rank)
        {
            if (calc_info->calc == CALC_HELP)
            {
                printf("help\n");
            }
            else if (calc_info->calc == CALC_VERSION)
            {
                printf("version\n");
            }
            else if (calc_info->calc == CALC_PH_SAVE_CREATE)
            {
                if (calc_info->code == DFT_CODE_QE)
                {
                    create_ph_save_dir_pp_qe(calc_info->input_file);
                }
            }
        }
    }

    free(calc_info);

    MPI_Finalize();

    return 0;
}





