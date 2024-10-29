#include "string_func.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
/*
This file contains some useful string functions
*/

void lowercase_str(char* str)
{
    // lower case all the chars in a string
    for (char* p = str; *p; ++p)
    {
        *p = tolower(*p);
    }
}

ND_int parser_doubles_from_string(char* str, ELPH_float* out)
{
    /*
    Extract all float values from given string

    if out == NULL, it return number of float it parsed
    */
    char* p = str;
    char* q;
    double temp_val;
    ND_int count = 0;

    while (*p)
    {
        if (isdigit(*p) || ((*p == '-' || *p == '+') && isdigit(*(p + 1))))
        {
            temp_val = strtod(p, &q);

            if (p == q)
            {
                break;
            }
            else
            {
                p = q;
            }

            if (out != NULL)
            {
                out[count] = temp_val;
            }
            ++count;
        }
        else
        {
            ++p;
        }
    }
    return count;
}

bool string_start_with(char* str, char* compare_str, bool trim)
{
    /*
    Check if given string starts with a substring
    */
    char* a;
    char* b;
    a = str;
    b = compare_str;
    if (trim)
    {
        while (isspace(*a))
        {
            ++a;
        }
        while (isspace(*b))
        {
            ++b;
        }
    }
    if (b[0] != a[0])
    {
        return false;
    }
    return !strncmp(a, b, strlen(b));
}

bool string_end_with(char* str, char* compare_str, bool trim)
{
    /*
    Check if given string ends with a substring
    */
    char* temp_str =
        malloc(sizeof(char) * (strlen(str) + strlen(compare_str) + 2));
    CHECK_ALLOC(temp_str);

    char* a = temp_str;
    char* b = temp_str + strlen(str) + 1;

    strcpy(a, str);
    strcpy(b, compare_str);

    str_reverse_in_place(a);
    str_reverse_in_place(b);
    if (trim)
    {
        while (isspace(*a))
        {
            ++a;
        }
        while (isspace(*b))
        {
            ++b;
        }
    }

    bool ret_value = !strncmp(a, b, strlen(b));
    if (b[0] != a[0])
    {
        ret_value = false;
    }
    free(temp_str);
    return ret_value;
}

char* str_reverse_in_place(char* str)
{
    /*
    Do a inplace reversing of string
    */
    ND_int len = strlen(str);

    if (len == 0)
    {
        return str;
    }

    char* p1 = str;

    char* p2 = str + len - 1;

    while (p1 < p2)
    {
        char tmp = *p1;
        *p1++ = *p2;
        *p2-- = tmp;
    }
    return str;
}

void str_replace_chars(char* str_in, const char* delimters,
                       const char* replace_chars)
{
    ND_int ndelimters = strlen(delimters);
    // if  ndelimters != strlen(replace_chars) buffer overflow

    ND_int str_in_len = strlen(str_in);

    for (ND_int i = 0; i < str_in_len; ++i)
    {
        for (ND_int j = 0; j < ndelimters; ++j)
        {
            if (str_in[i] == delimters[j])
            {
                str_in[i] = replace_chars[j];
                break;
            }
        }
    }
}
