#pragma once
#include <stdbool.h>
#include <stddef.h>

#include "elphC.h"

// safer strncpy version
char* strlcpy_custom(char* restrict dest, const char* restrict src,
                     size_t count);

// lower case a string
void lowercase_str(char* str);

// Extract all float values from given string
// if out == NULL, it return number of float it parsed
ND_int parse_floats_from_string(const char* str, ELPH_float* out,
                                ND_int out_size);

// Check if given string starts with a substring
bool string_start_with(char* str, char* compare_str, bool trim);

// Check if given string ends with a substring
bool string_end_with(char* str, char* compare_str, bool trim);

// Do a inplace reversing of string
char* str_reverse_in_place(char* str);

void str_replace_chars(char* str_in, const char* delimters,
                       const char* replace_chars);

// given a string, returns bool
bool parse_bool_input(const char* str);

void strip_quotes_in_string(char* s);
