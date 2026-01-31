#include "string_func.h"

#include <ctype.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
/*
This file contains some useful string functions
*/

char* strncpy_custom(char* dest, const char* src, size_t count)
{
    // this does strncpy(dest,src,n-1) and sets nth element as \0
    if (1 == count)
    {
        dest[0] = '\0';
    }
    else if (count > 1)
    {
        strncpy(dest, src, count - 1);
        dest[count - 1] = '\0';
    }
    return dest;
}

void lowercase_str(char* str)
{
    // lower case all the chars in a string
    for (char* p = str; *p; ++p)
    {
        *p = tolower((unsigned char)(*p));
    }
}

ND_int parse_floats_from_string(const char* str, ELPH_float* out,
                                ND_int out_size)
{
    /*
    Extract atmost out_size float values from given string

    out_size is size of the out buffer. In case more floats
    found than the buffer size, only the first out_size floats
    are written.

    Returns: Total number of floats found in the string (may exceed out_size).
    Compare return value with out_size to detect truncation.
    if out == NULL, it returns number of floats found in the string.
    */
    const char* p = str;
    char* q;
    double temp_val;
    ND_int count = 0;

    while (*p)
    {
        if (isdigit((unsigned char)(*p)) ||
            ((*p == '-' || *p == '+') && isdigit((unsigned char)(*(p + 1)))))
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

            if (out != NULL && count < out_size)
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
    if (str == NULL || compare_str == NULL)
    {
        return false;
    }
    //
    char* a;
    char* b;
    a = str;
    b = compare_str;
    if (trim)
    {
        while (isspace((unsigned char)(*a)))
        {
            ++a;
        }
        while (isspace((unsigned char)(*b)))
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
    if (str == NULL || compare_str == NULL)
    {
        return false;
    }
    //
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
        while (isspace((unsigned char)(*a)))
        {
            ++a;
        }
        while (isspace((unsigned char)(*b)))
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
    if (str == NULL)
    {
        return NULL;
    }
    //
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

static int my_strcasecmp(const char* a, const char* b)
{
    if (!a || !b)
    {
        return (a == b) ? 0 : -1;
    }

    while (*a && *b)
    {
        int diff = tolower((unsigned char)*a) - tolower((unsigned char)*b);
        if (diff != 0)
        {
            return diff;
        }
        a++;
        b++;
    }
    return tolower((unsigned char)*a) - tolower((unsigned char)*b);
}

bool parse_bool_input(const char* str)
{
    if (!str)
    {
        return false;
    }

    const char* s = str;

    // Skip leading whitespace
    while (*s && isspace((unsigned char)*s))
    {
        s++;
    }

    // Empty string after trimming
    if (*s == '\0')
    {
        return false;
    }

    // Check for Fortran-style logicals: .true. or .false. (case-insensitive)
    if (*s == '.')
    {
        s++;  // s now points to the first character after the dot

        // Check for .true. (5 characters: t, r, u, e, .)
        if (tolower((unsigned char)s[0]) == 't' &&
            tolower((unsigned char)s[1]) == 'r' &&
            tolower((unsigned char)s[2]) == 'u' &&
            tolower((unsigned char)s[3]) == 'e' && s[4] == '.' &&
            s[5] == '\0')  // Ensure it is the end of the string
        {
            return true;
        }
        // Check for .false. (6 characters: f, a, l, s, e, .)
        else if (tolower((unsigned char)s[0]) == 'f' &&
                 tolower((unsigned char)s[1]) == 'a' &&
                 tolower((unsigned char)s[2]) == 'l' &&
                 tolower((unsigned char)s[3]) == 's' &&
                 tolower((unsigned char)s[4]) == 'e' && s[5] == '.' &&
                 s[6] == '\0')  // Ensure it is the end of the string
        {
            return false;
        }

        error_msg(
            ".true./true/yes/on/1 or .false./false/no/off/0 (case insensitive) "
            "are only accepted.");
        return false;
    }

    // Check for common true values (case-insensitive)
    // *Note: The '1' check is also fixed to ensure no trailing characters.
    if (my_strcasecmp(s, "true") == 0 || my_strcasecmp(s, "yes") == 0 ||
        my_strcasecmp(s, "on") == 0 || (*s == '1' && *(s + 1) == '\0'))
    {
        return true;
    }

    // Check for common false values (case-insensitive)
    // *Note: The '0' check is also fixed to ensure no trailing characters.
    else if (my_strcasecmp(s, "false") == 0 || my_strcasecmp(s, "no") == 0 ||
             my_strcasecmp(s, "off") == 0 || (*s == '0' && *(s + 1) == '\0'))
    {
        return false;
    }

    error_msg(
        ".true./true/yes/on/1 or .false./false/no/off/0 (case insensitive) are "
        "only accepted.");
    // Default to false for unrecognized input
    return false;
}

void strip_quotes_in_string(char* s)
{
    char* start = s;
    char* end;
    char quote;

    if (!s || !*s)
    {
        return;
    }

    /* Must start and end with same quote */
    quote = s[0];
    end = s + strlen(s) - 1;

    if ((quote == '"' || quote == '\'') && *end == quote)
    {
        /* Move inside the quotes */
        start++;
        end--;

        /* Trim leading whitespace inside quotes */
        while (start <= end && isspace((unsigned char)*start))
        {
            start++;
        }

        /* Trim trailing whitespace inside quotes */
        while (end >= start && isspace((unsigned char)*end))
        {
            end--;
        }

        /* Copy cleaned string back */
        memmove(s, start, (size_t)(end - start + 1));
        s[end - start + 1] = '\0';
    }
}
