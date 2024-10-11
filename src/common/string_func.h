#pragma once
#include <stdbool.h>

#include "../elphC.h"

// lower case a string
void lowercase_str(char* str);

// Extract all float values from given string
// if out == NULL, it return number of float it parsed
ND_int parser_doubles_from_string(char* str, ELPH_float* out);

// Check if given string starts with a substring
bool string_start_with(char* str, char* compare_str, bool trim);

// Check if given string ends with a substring
bool string_end_with(char* str, char* compare_str, bool trim);

// Do a inplace reversing of string
char* str_reverse_in_place(char* str);

void str_replace_chars(char* str_in, const char* delimters,
                       const char* replace_chars);
