/*
 * unit_test_string_func.c
 *
 * Compile example:
 *   gcc -Wall -Wextra -std=c11 string_func.c unit_test_string_func.c -o test
 */

#include "common/string_func.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------
   Minimal test framework
   ------------------------------------------------------- */

static int tests_run = 0;
static int tests_failed = 0;

#define ASSERT_TRUE(cond)                                           \
    do                                                              \
    {                                                               \
        tests_run++;                                                \
        if (!(cond))                                                \
        {                                                           \
            tests_failed++;                                         \
            printf("FAIL: %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        }                                                           \
    } while (0)

#define ASSERT_STR_EQ(a, b) ASSERT_TRUE(strcmp((a), (b)) == 0)

#define ASSERT_INT_EQ(a, b) ASSERT_TRUE((a) == (b))

#define ASSERT_FLOAT_EQ(a, b, eps) ASSERT_TRUE(fabs((a) - (b)) < (eps))

/* -------------------------------------------------------
   Tests
   ------------------------------------------------------- */

void test_strlcpy_custom(void)
{
    char buf[16];

    strlcpy_custom(buf, "hello", sizeof(buf));
    ASSERT_STR_EQ(buf, "hello");

    strlcpy_custom(buf, "toolongstring", 5);
    ASSERT_STR_EQ(buf, "tool");

    strlcpy_custom(buf, "x", 1);
    ASSERT_STR_EQ(buf, "");
}

void test_lowercase_str(void)
{
    char s1[] = "HELLO World!";
    lowercase_str(s1);
    ASSERT_STR_EQ(s1, "hello world!");

    char s2[] = "";
    lowercase_str(s2);
    ASSERT_STR_EQ(s2, "");
}

void test_parse_floats_from_string(void)
{
    const char* str = "1.5 2 -.25 abc 0.4e1";
    ELPH_float out[4];

    ND_int count = parse_floats_from_string(str, out, 4);

    ASSERT_INT_EQ(count, 4);
    ASSERT_FLOAT_EQ(out[0], 1.5, 1e-6);
    ASSERT_FLOAT_EQ(out[1], 2.0, 1e-6);
    ASSERT_FLOAT_EQ(out[2], -0.25, 1e-6);
    ASSERT_FLOAT_EQ(out[3], 4.0, 1e-6);

    /* Count-only mode */
    ASSERT_INT_EQ(parse_floats_from_string(str, NULL, 0), 4);
}

void test_string_start_with(void)
{
    ASSERT_TRUE(string_start_with("hello world", "hello", false));
    ASSERT_TRUE(string_start_with("   hello", "hello", true));
    ASSERT_TRUE(!string_start_with("world", "hello", false));
}

void test_string_end_with(void)
{
    ASSERT_TRUE(string_end_with("hello.txt", ".txt", false));
    ASSERT_TRUE(string_end_with("hello   ", "hello", true));
    ASSERT_TRUE(!string_end_with("hello", "world", false));
}

void test_str_reverse_in_place(void)
{
    char s1[] = "abcd";
    str_reverse_in_place(s1);
    ASSERT_STR_EQ(s1, "dcba");

    char s2[] = "a";
    str_reverse_in_place(s2);
    ASSERT_STR_EQ(s2, "a");

    char s3[] = "";
    str_reverse_in_place(s3);
    ASSERT_STR_EQ(s3, "");

    char s4[] = " abc";
    str_reverse_in_place(s4);
    ASSERT_STR_EQ(s4, "cba ");
}

void test_str_replace_chars(void)
{
    char s[] = "a,b;c#";
    str_replace_chars(s, ",;#", "..\0");
    puts(s);
    ASSERT_STR_EQ(s, "a.b.c");
}

void test_my_strcasecmp(void)
{
    ASSERT_INT_EQ(my_strcasecmp("Hello", "helLo"), 0);
    ASSERT_TRUE(my_strcasecmp("abc", "abd") < 0);
    ASSERT_TRUE(my_strcasecmp("abd", "abc") > 0);
}

void test_parse_bool_input(void)
{
    ASSERT_TRUE(parse_bool_input("true"));
    ASSERT_TRUE(parse_bool_input("YES"));
    ASSERT_TRUE(parse_bool_input("1"));

    ASSERT_TRUE(!parse_bool_input("false"));
    ASSERT_TRUE(!parse_bool_input("off"));
    ASSERT_TRUE(!parse_bool_input("0"));

    ASSERT_TRUE(parse_bool_input(".true."));
    ASSERT_TRUE(!parse_bool_input(".false."));
}

void test_strip_quotes_in_string(void)
{
    char s1[] = "\" hello world \"";
    strip_quotes_in_string(s1);
    ASSERT_STR_EQ(s1, "hello world");

    char s2[] = "'test'";
    strip_quotes_in_string(s2);
    ASSERT_STR_EQ(s2, "test");

    char s3[] = "noquotes";
    strip_quotes_in_string(s3);
    ASSERT_STR_EQ(s3, "noquotes");

    char s4[] = "'    noquotes'";
    strip_quotes_in_string(s4);
    ASSERT_STR_EQ(s4, "noquotes");

    char s5[] = "'noquotes    '";
    strip_quotes_in_string(s5);
    ASSERT_STR_EQ(s5, "noquotes");

    char s6[] = "\"   noquotes  \"";
    strip_quotes_in_string(s6);
    ASSERT_STR_EQ(s6, "noquotes");
}

/* -------------------------------------------------------
   Test runner
   ------------------------------------------------------- */

int main(void)
{
    test_strlcpy_custom();
    test_lowercase_str();
    test_parse_floats_from_string();
    test_string_start_with();
    test_string_end_with();
    test_str_reverse_in_place();
    test_str_replace_chars();
    test_my_strcasecmp();
    test_parse_bool_input();
    test_strip_quotes_in_string();

    printf("\nTests run: %d\n", tests_run);
    printf("Failures : %d\n", tests_failed);

    if (tests_failed == 0)
    {
        printf("All tests passed!\n");
    }

    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
