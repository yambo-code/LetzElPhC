#pragma once
#include <stdio.h>
#include <stdlib.h>

/* Macro to get the function name */
#ifndef Current_Function_NAME
    #ifdef WIN32
        #define Current_Function_NAME __FUNCTION__
    #else
        #define Current_Function_NAME __func__
    #endif
#endif

/* Function to print error message to the the file */
#define error_msg(print_str) nd_error_msg(print_str, Current_Function_NAME)
void nd_error_msg(const char *error_msg, const char *func_name);
