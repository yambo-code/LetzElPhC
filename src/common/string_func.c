#include "string_func.h"

#include <ctype.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
/*
This file contains some useful string functions
*/

char* strlcpy_custom(char* restrict dest, const char* restrict src,
                     size_t count)
{
    // copies MAX(count,strlen(src)) characters into dest pointer
    if (count == 0)
    {
        return dest;
    }

    char* d = dest;
    const char* s = src;
    size_t n = count;

    while (--n > 0 && *s)
    {
        *d++ = *s++;
        // Copy until count runs out or src ends
    }

    *d = '\0';
    // Always ensure null termination

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
    char* end;
    ND_int count = 0;

    while (*p)
    {
        double val = strtod(p, &end);
        // Attempt to parse a double

        if (p == end)
        {
            p++;
            // No conversion occurred, advance manually to skip non-numeric char
        }
        else
        {
            if (out != NULL && count < out_size)
            {
                out[count] = val;
            }
            count++;
            p = end;
            // Advance p to where strtod stopped parsing
        }
    }
    return count;
}

bool string_start_with(char* str, char* compare_str, bool trim)
{
    if (str == NULL || compare_str == NULL)
    {
        return false;
    }

    char* a = str;
    char* b = compare_str;

    if (trim)
    {
        // Trim leading whitespace
        while (isspace((unsigned char)*a))
        {
            ++a;
        }
        while (isspace((unsigned char)*b))
        {
            ++b;
        }
    }

    // Determine length of compare_str without trailing whitespace
    size_t blen = strlen(b);
    if (trim)
    {
        while (blen > 0 && isspace((unsigned char)b[blen - 1]))
        {
            --blen;
        }
    }

    // If compare string is empty after trimming
    if (blen == 0)
    {
        return false;
    }

    // Compare only up to trimmed length
    return strncmp(a, b, blen) == 0;
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

    size_t s_len = strlen(str);
    size_t c_len = strlen(compare_str);

    if (c_len == 0)
    {
        return true;
    }
    if (s_len == 0)
    {
        return false;
    }

    const char* s_end = str + s_len - 1;
    // Point to the last character of str

    const char* c_start = compare_str;

    const char* c_end = compare_str + c_len - 1;
    // Point to the last character of compare_str

    if (trim)
    {
        while (s_end >= str && isspace((unsigned char)*s_end))
        {
            s_end--;
            // Backtrack past trailing whitespace in str
        }
        while (c_end >= compare_str && isspace((unsigned char)*c_end))
        {
            c_end--;
            // Backtrack past trailing whitespace in compare_str
        }
        while (c_start <= c_end && isspace((unsigned char)*c_start))
        {
            c_start++;
            // trim leading whitespace in compare_str
        }
    }

    while (c_end >= c_start)
    {
        if (s_end < str || *s_end != *c_end)
        {
            return false;
        }

        s_end--;
        c_end--;
    }

    return true;
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
    // It is responsiblity of the user to give same number of
    // delimters and replace_chars. Note that '\0' is not allowed
    // in delimters but is allowed in replace_chars.
    if (!str_in || !delimters || !replace_chars)
    {
        return;
    }

    unsigned char map[256];
    // Create a generic map where each char maps to itself initially

    for (int i = 0; i < 256; ++i)
    {
        map[i] = (unsigned char)i;
    }

    const char* d = delimters;
    const char* r = replace_chars;

    while (*d)
    {
        map[(unsigned char)*d] = (unsigned char)*r;
        // Update the map: delimiter character points to replacement character
        d++;
        r++;
    }

    char* p = str_in;
    while (*p)
    {
        *p = (char)map[(unsigned char)*p];
        // Replace in O(1) using the lookup table
        p++;
    }
}

int my_strcasecmp(const char* a, const char* b)
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

    size_t str_len = strlen(s);
    // return incase of single quotes
    if (str_len == 1)
    {
        return;
    }
    /* Must start and end with same quote */
    quote = s[0];
    end = s + str_len - 1;

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
