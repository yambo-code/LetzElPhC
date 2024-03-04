#include "elph.h"


void init_kernel(struct kernel_info * kernel)
{
    strcpy(kernel->name_str,"dfpt");
    kernel->bare_loc = true;
    kernel->non_loc = true;
    kernel->screening = ELPH_DFPT_SCREENING; 
}

void set_kernel(const char * kernel_str, struct kernel_info * kernel)
{
    size_t name_str_size =  sizeof(kernel->name_str)-1;
    kernel->name_str[name_str_size] = '\0';
    strncpy(kernel->name_str, kernel_str, name_str_size);
    
    if (!strcmp(kernel_str,"dfpt"))
    {
        kernel->bare_loc = true;
        kernel->non_loc = true;
        kernel->screening = ELPH_DFPT_SCREENING; 
    }
    else if (!strcmp(kernel_str,"dfpt_local"))
    {
        kernel->bare_loc = true;
        kernel->non_loc = false;
        kernel->screening = ELPH_DFPT_SCREENING;
    }
    else if (!strcmp(kernel_str,"bare"))
    {
        kernel->bare_loc = true;
        kernel->non_loc = true;
        kernel->screening = ELPH_NO_SCREENING;
    }
    else if (!strcmp(kernel_str,"bare_local"))
    {
        kernel->bare_loc = true;
        kernel->non_loc = false;
        kernel->screening = ELPH_NO_SCREENING;
    }
    else 
    {
        error_msg("invalid value for kernel in the input file");
    }
    
}

