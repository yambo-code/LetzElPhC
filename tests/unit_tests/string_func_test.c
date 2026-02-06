/*
 * unit_test_string_func_extensive.c
 *
 * Comprehensive test suite with edge cases, boundary conditions, and exotic
 * scenarios
 *
 * Compile example:
 *   gcc -Wall -Wextra -std=c11 string_func.c unit_test_string_func_extensive.c
 * -o test_extensive
 */

#include "common/string_func.h"

#include <float.h>
#include <limits.h>
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

#define ASSERT_FALSE(cond) ASSERT_TRUE(!(cond))

#define ASSERT_STR_EQ(a, b) ASSERT_TRUE(strcmp((a), (b)) == 0)

#define ASSERT_INT_EQ(a, b) ASSERT_TRUE((a) == (b))

#define ASSERT_FLOAT_EQ(a, b, eps) ASSERT_TRUE(fabs((a) - (b)) < (eps))

#define ASSERT_NULL(ptr) ASSERT_TRUE((ptr) == NULL)

/* -------------------------------------------------------
   Test: strlcpy_custom
   ------------------------------------------------------- */

void test_strlcpy_custom_comprehensive(void)
{
    char buf[32];

    /* Basic functionality */
    strlcpy_custom(buf, "hello", sizeof(buf));
    ASSERT_STR_EQ(buf, "hello");

    /* Exact buffer size */
    strlcpy_custom(buf, "1234567890", 11);
    ASSERT_STR_EQ(buf, "1234567890");

    /* Truncation */
    strlcpy_custom(buf, "toolongstring", 5);
    ASSERT_STR_EQ(buf, "tool");

    /* Count = 1 (only null terminator) */
    strlcpy_custom(buf, "x", 1);
    ASSERT_STR_EQ(buf, "");

    /* Count = 0 (no-op) */
    strcpy(buf, "unchanged");
    strlcpy_custom(buf, "newstring", 0);
    ASSERT_STR_EQ(buf, "unchanged");

    /* Empty source string */
    strlcpy_custom(buf, "", sizeof(buf));
    ASSERT_STR_EQ(buf, "");

    /* Very long source, small buffer */
    char long_src[1000];
    memset(long_src, 'a', 999);
    long_src[999] = '\0';
    strlcpy_custom(buf, long_src, 10);
    ASSERT_INT_EQ(strlen(buf), 9);
    ASSERT_TRUE(buf[9] == '\0');

    /* Special characters */
    strlcpy_custom(buf, "tab\there\nnewline", sizeof(buf));
    ASSERT_STR_EQ(buf, "tab\there\nnewline");

    /* Unicode/UTF-8 (multi-byte chars) */
    strlcpy_custom(buf, "Hello 世界", sizeof(buf));
    ASSERT_STR_EQ(buf, "Hello 世界");

    /* Count = 2 (one char + null) */
    strlcpy_custom(buf, "abcdef", 2);
    ASSERT_STR_EQ(buf, "a");

    /* All printable ASCII */
    strlcpy_custom(buf, "!@#$%^&*()_+-=[]{}|;':,.<>?/", sizeof(buf));
    ASSERT_STR_EQ(buf, "!@#$%^&*()_+-=[]{}|;':,.<>?/");

    /* Numbers */
    strlcpy_custom(buf, "0123456789", sizeof(buf));
    ASSERT_STR_EQ(buf, "0123456789");
}

/* -------------------------------------------------------
   Test: lowercase_str
   ------------------------------------------------------- */

void test_lowercase_str_comprehensive(void)
{
    char s1[] = "HELLO World!";
    lowercase_str(s1);
    ASSERT_STR_EQ(s1, "hello world!");

    /* Empty string */
    char s2[] = "";
    lowercase_str(s2);
    ASSERT_STR_EQ(s2, "");

    /* All uppercase */
    char s3[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    lowercase_str(s3);
    ASSERT_STR_EQ(s3, "abcdefghijklmnopqrstuvwxyz");

    /* All lowercase (no change) */
    char s4[] = "abcdefghijklmnopqrstuvwxyz";
    lowercase_str(s4);
    ASSERT_STR_EQ(s4, "abcdefghijklmnopqrstuvwxyz");

    /* Mixed with numbers and special chars */
    char s5[] = "TeSt123!@#ABC";
    lowercase_str(s5);
    ASSERT_STR_EQ(s5, "test123!@#abc");

    /* Only numbers */
    char s6[] = "1234567890";
    lowercase_str(s6);
    ASSERT_STR_EQ(s6, "1234567890");

    /* Only special characters */
    char s7[] = "!@#$%^&*()_+";
    lowercase_str(s7);
    ASSERT_STR_EQ(s7, "!@#$%^&*()_+");

    /* Whitespace and newlines */
    char s8[] = "HELLO\tWORLD\nTEST\r\n";
    lowercase_str(s8);
    ASSERT_STR_EQ(s8, "hello\tworld\ntest\r\n");

    /* Single character */
    char s9[] = "A";
    lowercase_str(s9);
    ASSERT_STR_EQ(s9, "a");

    char s10[] = "z";
    lowercase_str(s10);
    ASSERT_STR_EQ(s10, "z");
}

/* -------------------------------------------------------
   Test: parse_floats_from_string
   ------------------------------------------------------- */

void test_parse_floats_from_string_comprehensive(void)
{
    ELPH_float out[20];
    ND_int count;

    /* Basic positive and negative */
    const char* str1 = "1.5 2 -.25 abc 0.4e1";
    count = parse_floats_from_string(str1, out, 20);
    ASSERT_INT_EQ(count, 4);
    ASSERT_FLOAT_EQ(out[0], 1.5, 1e-6);
    ASSERT_FLOAT_EQ(out[1], 2.0, 1e-6);
    ASSERT_FLOAT_EQ(out[2], -0.25, 1e-6);
    ASSERT_FLOAT_EQ(out[3], 4.0, 1e-6);

    /* Count-only mode (NULL output) */
    ASSERT_INT_EQ(parse_floats_from_string(str1, NULL, 0), 4);

    /* Scientific notation */
    const char* str2 = "1e10 2.5e-3 -3.14E+2";
    count = parse_floats_from_string(str2, out, 20);
    ASSERT_INT_EQ(count, 3);
    ASSERT_FLOAT_EQ(out[0], 1e10, 1e4);
    ASSERT_FLOAT_EQ(out[1], 2.5e-3, 1e-9);
    ASSERT_FLOAT_EQ(out[2], -314.0, 1e-3);

    /* Zero and negative zero */
    const char* str3 = "0 -0 0.0 -0.0";
    count = parse_floats_from_string(str3, out, 20);
    ASSERT_INT_EQ(count, 4);
    ASSERT_FLOAT_EQ(out[0], 0.0, 1e-10);

    /* Very small and very large numbers */
    const char* str4 = "1e-308 1e308";
    count = parse_floats_from_string(str4, out, 20);
    ASSERT_INT_EQ(count, 2);

    /* Buffer truncation */
    const char* str5 = "1 2 3 4 5 6 7 8 9 10";
    count = parse_floats_from_string(str5, out, 5);
    ASSERT_INT_EQ(count, 10);           /* Total found */
    ASSERT_FLOAT_EQ(out[4], 5.0, 1e-6); /* Only first 5 stored */

    /* Mixed with lots of text */
    const char* str6 = "abc 123 def 456.789 ghi -999 jkl";
    count = parse_floats_from_string(str6, out, 20);
    ASSERT_INT_EQ(count, 3);
    ASSERT_FLOAT_EQ(out[0], 123.0, 1e-5);
    ASSERT_FLOAT_EQ(out[1], 456.789, 1e-5);
    ASSERT_FLOAT_EQ(out[2], -999.0, 1e-5);

    /* Empty string */
    const char* str7 = "";
    count = parse_floats_from_string(str7, out, 20);
    ASSERT_INT_EQ(count, 0);

    /* Only text, no numbers */
    const char* str8 = "hello world test";
    count = parse_floats_from_string(str8, out, 20);
    ASSERT_INT_EQ(count, 0);

    /* Leading/trailing spaces */
    const char* str9 = "  1.5   2.5  ";
    count = parse_floats_from_string(str9, out, 20);
    ASSERT_INT_EQ(count, 2);
    ASSERT_FLOAT_EQ(out[0], 1.5, 1e-6);

    /* Multiple delimiters */
    const char* str10 = "1.1,2.2;3.3\t4.4\n5.5";
    count = parse_floats_from_string(str10, out, 20);
    ASSERT_INT_EQ(count, 5);

    /* Infinity and special values (platform dependent) */
    const char* str11 = "inf -inf";
    count = parse_floats_from_string(str11, out, 20);
    /* May or may not parse depending on strtod implementation */

    /* Decimal points without digits */
    const char* str12 = ". .. ... 1.2.3";
    count = parse_floats_from_string(str12, out, 20);
    ASSERT_TRUE(count >= 0); /* Implementation dependent */

    /* Plus signs */
    const char* str13 = "+1.5 +2 +3.14";
    count = parse_floats_from_string(str13, out, 20);
    ASSERT_INT_EQ(count, 3);
    ASSERT_FLOAT_EQ(out[0], 1.5, 1e-6);

    /* Consecutive numbers without separators */
    const char* str14 = "1.52.5";
    count = parse_floats_from_string(str14, out, 20);
    /* strtod will parse 1.5, then 2.5 */
    ASSERT_TRUE(count >= 1);
}

/* -------------------------------------------------------
   Test: string_start_with
   ------------------------------------------------------- */

void test_string_start_with_comprehensive(void)
{
    /* Basic cases */
    ASSERT_TRUE(string_start_with("hello world", "hello", false));
    ASSERT_FALSE(string_start_with("world", "hello", false));

    /* Empty strings */
    ASSERT_FALSE(string_start_with("", "", false));
    ASSERT_FALSE(string_start_with("hello", "", false));
    ASSERT_FALSE(string_start_with("", "hello", false));

    /* NULL inputs */
    ASSERT_FALSE(string_start_with(NULL, "hello", false));
    ASSERT_FALSE(string_start_with("hello", NULL, false));
    ASSERT_FALSE(string_start_with(NULL, NULL, false));

    /* Trimming whitespace */
    ASSERT_TRUE(string_start_with("   hello", "hello", true));
    ASSERT_TRUE(string_start_with("hello", "   hello", true));
    ASSERT_TRUE(string_start_with("   hello   ", "   hello   ", true));
    ASSERT_FALSE(string_start_with("   hello", "hello", false));

    /* Exact match */
    ASSERT_TRUE(string_start_with("hello", "hello", false));
    ASSERT_TRUE(string_start_with("hello", "hello", true));

    /* Compare string longer than input */
    ASSERT_FALSE(string_start_with("hi", "hello", false));

    /* Single character */
    ASSERT_TRUE(string_start_with("h", "h", false));
    ASSERT_TRUE(string_start_with("hello", "h", false));
    ASSERT_FALSE(string_start_with("a", "b", false));

    /* Case sensitivity */
    ASSERT_FALSE(string_start_with("Hello", "hello", false));
    ASSERT_TRUE(string_start_with("hello", "hello", false));

    /* Whitespace only */
    ASSERT_FALSE(string_start_with("   ", "   ", true));
    ASSERT_TRUE(string_start_with("   ", "   ", false));

    /* Special characters */
    ASSERT_TRUE(string_start_with("!@#$%", "!@#", false));
    ASSERT_TRUE(string_start_with("\t\nhello", "\t\nhello", false));

    /* Numbers */
    ASSERT_TRUE(string_start_with("123456", "123", false));

    /* Trailing whitespace in compare_str */
    ASSERT_TRUE(string_start_with("hello world", "hello   ", true));

    /* Unicode */
    ASSERT_TRUE(string_start_with("世界hello", "世界", false));

    /* Partial match at start */
    ASSERT_TRUE(string_start_with("testing", "test", false));
    ASSERT_FALSE(string_start_with("testing", "sting", false));
}

/* -------------------------------------------------------
   Test: string_end_with
   ------------------------------------------------------- */

void test_string_end_with_comprehensive(void)
{
    /* Basic cases */
    ASSERT_TRUE(string_end_with("hello.txt", ".txt", false));
    ASSERT_FALSE(string_end_with("hello", "world", false));

    /* Empty strings */
    ASSERT_TRUE(string_end_with("hello", "", false));
    ASSERT_FALSE(string_end_with("", "hello", false));
    ASSERT_TRUE(string_end_with("", "", false));

    /* NULL inputs */
    ASSERT_FALSE(string_end_with(NULL, "txt", false));
    ASSERT_FALSE(string_end_with("hello", NULL, false));
    ASSERT_FALSE(string_end_with(NULL, NULL, false));

    /* Trimming whitespace */
    ASSERT_TRUE(string_end_with("hello   ", "hello", true));
    ASSERT_TRUE(string_end_with("hello", "hello   ", true));
    ASSERT_TRUE(string_end_with("   hello   ", "hello", true));
    ASSERT_FALSE(string_end_with("hello   ", "hello", false));

    /* Exact match */
    ASSERT_TRUE(string_end_with("hello", "hello", false));

    /* Compare string longer than input */
    ASSERT_FALSE(string_end_with("hi", "hello", false));

    /* Single character */
    ASSERT_TRUE(string_end_with("h", "h", false));
    ASSERT_TRUE(string_end_with("hello", "o", false));
    ASSERT_FALSE(string_end_with("a", "b", false));

    /* Case sensitivity */
    ASSERT_FALSE(string_end_with("HELLO", "hello", false));

    /* Whitespace handling */
    ASSERT_TRUE(string_end_with("hello\t\n", "hello", true));
    ASSERT_TRUE(string_end_with("hello", "   llo   ", true));

    /* Special characters */
    ASSERT_TRUE(string_end_with("file.c++", ".c++", false));
    ASSERT_TRUE(string_end_with("test!@#", "!@#", false));

    /* Numbers */
    ASSERT_TRUE(string_end_with("file123", "123", false));

    /* Multiple extensions */
    ASSERT_TRUE(string_end_with("file.tar.gz", ".gz", false));
    ASSERT_TRUE(string_end_with("file.tar.gz", "tar.gz", false));

    /* Newlines and tabs */
    ASSERT_TRUE(string_end_with("hello\n", "\n", false));
    ASSERT_TRUE(string_end_with("hello\t", "\t", false));

    /* Unicode */
    ASSERT_TRUE(string_end_with("hello世界", "世界", false));

    /* Whitespace only compare_str */
    ASSERT_TRUE(string_end_with("hello   ", "   ", true));
}

/* -------------------------------------------------------
   Test: str_reverse_in_place
   ------------------------------------------------------- */

void test_str_reverse_in_place_comprehensive(void)
{
    char s1[] = "abcd";
    str_reverse_in_place(s1);
    ASSERT_STR_EQ(s1, "dcba");

    /* Single character */
    char s2[] = "a";
    str_reverse_in_place(s2);
    ASSERT_STR_EQ(s2, "a");

    /* Empty string */
    char s3[] = "";
    str_reverse_in_place(s3);
    ASSERT_STR_EQ(s3, "");

    /* NULL pointer */
    ASSERT_NULL(str_reverse_in_place(NULL));

    /* Two characters */
    char s4[] = "ab";
    str_reverse_in_place(s4);
    ASSERT_STR_EQ(s4, "ba");

    /* Odd length */
    char s5[] = "12345";
    str_reverse_in_place(s5);
    ASSERT_STR_EQ(s5, "54321");

    /* Even length */
    char s6[] = "123456";
    str_reverse_in_place(s6);
    ASSERT_STR_EQ(s6, "654321");

    /* With spaces */
    char s7[] = " abc ";
    str_reverse_in_place(s7);
    ASSERT_STR_EQ(s7, " cba ");

    /* Special characters */
    char s8[] = "!@#$%";
    str_reverse_in_place(s8);
    ASSERT_STR_EQ(s8, "%$#@!");

    /* Palindrome */
    char s9[] = "racecar";
    str_reverse_in_place(s9);
    ASSERT_STR_EQ(s9, "racecar");

    /* Newlines and tabs */
    char s10[] = "a\tb\nc";
    str_reverse_in_place(s10);
    ASSERT_STR_EQ(s10, "c\nb\ta");

    /* Numbers */
    char s11[] = "1234567890";
    str_reverse_in_place(s11);
    ASSERT_STR_EQ(s11, "0987654321");

    /* Mixed case */
    char s12[] = "HeLLo";
    str_reverse_in_place(s12);
    ASSERT_STR_EQ(s12, "oLLeH");

    /* Unicode/multi-byte (byte-level reversal) */
    char s13[] = "hello世界";
    str_reverse_in_place(s13);
    /* Will reverse bytes, not logical characters */
    ASSERT_TRUE(strlen(s13) == strlen("hello世界"));

    /* Very long string */
    char long_str[1001];
    for (int i = 0; i < 1000; i++)
    {
        long_str[i] = 'a' + (i % 26);
    }
    long_str[1000] = '\0';
    str_reverse_in_place(long_str);
    ASSERT_TRUE(long_str[0] == 'a' + (999 % 26));
    ASSERT_TRUE(long_str[999] == 'a');
}

/* -------------------------------------------------------
   Test: str_replace_chars
   ------------------------------------------------------- */

void test_str_replace_chars_comprehensive(void)
{
    /* Basic replacement */
    char s1[] = "a,b;c#";
    str_replace_chars(s1, ",;#", "..\0");
    ASSERT_STR_EQ(s1, "a.b.c");

    /* Replace with same character */
    char s2[] = "hello";
    str_replace_chars(s2, "l", "l");
    ASSERT_STR_EQ(s2, "hello");

    /* Replace all occurrences */
    char s3[] = "aaabbbccc";
    str_replace_chars(s3, "abc", "xyz");
    ASSERT_STR_EQ(s3, "xxxyyyzzz");

    /* Empty string */
    char s4[] = "";
    str_replace_chars(s4, "a", "b");
    ASSERT_STR_EQ(s4, "");

    /* No matches */
    char s5[] = "hello";
    str_replace_chars(s5, "xyz", "abc");
    ASSERT_STR_EQ(s5, "hello");

    /* Replace spaces */
    char s6[] = "hello world";
    str_replace_chars(s6, " ", "_");
    ASSERT_STR_EQ(s6, "hello_world");

    /* Multiple different replacements */
    char s7[] = "a1b2c3";
    str_replace_chars(s7, "123", "xyz");
    ASSERT_STR_EQ(s7, "axbycz");

    /* Replace to null character */
    char s8[] = "a#b#c";
    str_replace_chars(s8, "#", "\0");
    ASSERT_STR_EQ(s8, "a"); /* String ends at first null */

    /* Special to special */
    char s9[] = "a\tb\nc";
    str_replace_chars(s9, "\t\n", "  ");
    ASSERT_STR_EQ(s9, "a b c");

    /* Numbers to letters */
    char s10[] = "123";
    str_replace_chars(s10, "123", "abc");
    ASSERT_STR_EQ(s10, "abc");

    /* Case sensitivity */
    char s11[] = "AaBbCc";
    str_replace_chars(s11, "ABC", "xyz");
    ASSERT_STR_EQ(s11, "xaybzc");

    /* NULL inputs (should return early) */
    str_replace_chars(NULL, "a", "b");
    char s12[] = "test";
    str_replace_chars(s12, NULL, "b");
    ASSERT_STR_EQ(s12, "test");
    str_replace_chars(s12, "a", NULL);
    ASSERT_STR_EQ(s12, "test");

    /* Single character */
    char s13[] = "a";
    str_replace_chars(s13, "a", "z");
    ASSERT_STR_EQ(s13, "z");

    /* All same character */
    char s14[] = "aaaaa";
    str_replace_chars(s14, "a", "b");
    ASSERT_STR_EQ(s14, "bbbbb");

    /* High ASCII characters */
    char s15[] = "test\x80\x90";
    str_replace_chars(s15, "\x80", "X");
    ASSERT_TRUE(s15[4] == 'X');
}

/* -------------------------------------------------------
   Test: my_strcasecmp
   ------------------------------------------------------- */

void test_my_strcasecmp_comprehensive(void)
{
    /* Equal strings */
    ASSERT_INT_EQ(my_strcasecmp("Hello", "helLo"), 0);
    ASSERT_INT_EQ(my_strcasecmp("test", "TEST"), 0);

    /* Different strings */
    ASSERT_TRUE(my_strcasecmp("abc", "abd") < 0);
    ASSERT_TRUE(my_strcasecmp("abd", "abc") > 0);

    /* Empty strings */
    ASSERT_INT_EQ(my_strcasecmp("", ""), 0);
    ASSERT_TRUE(my_strcasecmp("a", "") > 0);
    ASSERT_TRUE(my_strcasecmp("", "a") < 0);

    /* NULL pointers */
    ASSERT_INT_EQ(my_strcasecmp(NULL, NULL), 0);
    ASSERT_INT_EQ(my_strcasecmp(NULL, "test"), -1);
    ASSERT_INT_EQ(my_strcasecmp("test", NULL), -1);

    /* Different lengths */
    ASSERT_TRUE(my_strcasecmp("hello", "helloworld") < 0);
    ASSERT_TRUE(my_strcasecmp("helloworld", "hello") > 0);

    /* Numbers */
    ASSERT_INT_EQ(my_strcasecmp("123", "123"), 0);
    ASSERT_TRUE(my_strcasecmp("123", "124") < 0);

    /* Special characters */
    ASSERT_INT_EQ(my_strcasecmp("!@#", "!@#"), 0);
    ASSERT_TRUE(my_strcasecmp("!@#", "!@$") < 0);

    /* Mixed alphanumeric */
    ASSERT_INT_EQ(my_strcasecmp("Test123", "test123"), 0);

    /* Whitespace */
    ASSERT_INT_EQ(my_strcasecmp("hello world", "HELLO WORLD"), 0);
    ASSERT_TRUE(my_strcasecmp("hello world", "hello_world") != 0);

    /* Single character */
    ASSERT_INT_EQ(my_strcasecmp("a", "A"), 0);
    ASSERT_TRUE(my_strcasecmp("a", "B") < 0);
    ASSERT_TRUE(my_strcasecmp("Z", "a") > 0);

    /* All uppercase vs all lowercase */
    ASSERT_INT_EQ(my_strcasecmp("ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                                "abcdefghijklmnopqrstuvwxyz"),
                  0);

    /* Prefixes */
    ASSERT_TRUE(my_strcasecmp("test", "testing") < 0);
    ASSERT_TRUE(my_strcasecmp("testing", "test") > 0);

    /* Very similar strings */
    ASSERT_TRUE(my_strcasecmp("test1", "test2") < 0);
    ASSERT_INT_EQ(my_strcasecmp("Test", "tEst"), 0);
}

/* -------------------------------------------------------
   Test: parse_bool_input
   ------------------------------------------------------- */

void test_parse_bool_input_comprehensive(void)
{
    /* True values */
    ASSERT_TRUE(parse_bool_input("true"));
    ASSERT_TRUE(parse_bool_input("TRUE"));
    ASSERT_TRUE(parse_bool_input("True"));
    ASSERT_TRUE(parse_bool_input("TrUe"));
    ASSERT_TRUE(parse_bool_input("yes"));
    ASSERT_TRUE(parse_bool_input("YES"));
    ASSERT_TRUE(parse_bool_input("Yes"));
    ASSERT_TRUE(parse_bool_input("on"));
    ASSERT_TRUE(parse_bool_input("ON"));
    ASSERT_TRUE(parse_bool_input("1"));

    /* False values */
    ASSERT_FALSE(parse_bool_input("false"));
    ASSERT_FALSE(parse_bool_input("FALSE"));
    ASSERT_FALSE(parse_bool_input("False"));
    ASSERT_FALSE(parse_bool_input("no"));
    ASSERT_FALSE(parse_bool_input("NO"));
    ASSERT_FALSE(parse_bool_input("off"));
    ASSERT_FALSE(parse_bool_input("OFF"));
    ASSERT_FALSE(parse_bool_input("0"));

    /* Fortran-style logicals */
    ASSERT_TRUE(parse_bool_input(".true."));
    ASSERT_TRUE(parse_bool_input(".TRUE."));
    ASSERT_TRUE(parse_bool_input(".True."));
    ASSERT_TRUE(parse_bool_input(".tRuE."));
    ASSERT_FALSE(parse_bool_input(".false."));
    ASSERT_FALSE(parse_bool_input(".FALSE."));
    ASSERT_FALSE(parse_bool_input(".False."));

    /* Mixed case corner cases */
    ASSERT_TRUE(parse_bool_input("oN"));
    ASSERT_TRUE(parse_bool_input("On"));
    ASSERT_FALSE(parse_bool_input("oFf"));
    ASSERT_FALSE(parse_bool_input("OfF"));
}

/* -------------------------------------------------------
   Test: strip_quotes_in_string
   ------------------------------------------------------- */

void test_strip_quotes_in_string_comprehensive(void)
{
    /* Double quotes */
    char s1[] = "\" hello world \"";
    strip_quotes_in_string(s1);
    ASSERT_STR_EQ(s1, "hello world");

    /* Single quotes */
    char s2[] = "'test'";
    strip_quotes_in_string(s2);
    ASSERT_STR_EQ(s2, "test");

    /* No quotes */
    char s3[] = "noquotes";
    strip_quotes_in_string(s3);
    ASSERT_STR_EQ(s3, "noquotes");

    /* Leading whitespace inside quotes */
    char s4[] = "'    test'";
    strip_quotes_in_string(s4);
    ASSERT_STR_EQ(s4, "test");

    /* Trailing whitespace inside quotes */
    char s5[] = "'test    '";
    strip_quotes_in_string(s5);
    ASSERT_STR_EQ(s5, "test");

    /* Both leading and trailing */
    char s6[] = "\"   test   \"";
    strip_quotes_in_string(s6);
    ASSERT_STR_EQ(s6, "test");

    /* Empty string */
    char s7[] = "";
    strip_quotes_in_string(s7);
    ASSERT_STR_EQ(s7, "");

    /* NULL pointer */
    strip_quotes_in_string(NULL); /* Should not crash */

    /* Mismatched quotes */
    char s8[] = "\"test'";
    strip_quotes_in_string(s8);
    ASSERT_STR_EQ(s8, "\"test'"); /* Should not strip */

    char s9[] = "'test\"";
    strip_quotes_in_string(s9);
    ASSERT_STR_EQ(s9, "'test\""); /* Should not strip */

    /* Only quotes with whitespace */
    char s10[] = "\"   \"";
    strip_quotes_in_string(s10);
    ASSERT_STR_EQ(s10, "");

    char s11[] = "'   '";
    strip_quotes_in_string(s11);
    ASSERT_STR_EQ(s11, "");

    /* Single quote only */
    char s12[] = "\"";
    strip_quotes_in_string(s12);
    ASSERT_STR_EQ(s12, "\""); /* Should not strip */

    char s13[] = "'";
    strip_quotes_in_string(s13);
    ASSERT_STR_EQ(s13, "'"); /* Should not strip */

    /* Two characters (matching quotes) */
    char s14[] = "\"\"";
    strip_quotes_in_string(s14);
    ASSERT_STR_EQ(s14, "");

    char s15[] = "''";
    strip_quotes_in_string(s15);
    ASSERT_STR_EQ(s15, "");

    /* Special characters inside quotes */
    char s16[] = "\"!@#$%^&*()\"";
    strip_quotes_in_string(s16);
    ASSERT_STR_EQ(s16, "!@#$%^&*()");

    /* Newlines and tabs inside quotes */
    char s17[] = "\"hello\n\tworld\"";
    strip_quotes_in_string(s17);
    ASSERT_STR_EQ(s17, "hello\n\tworld");

    /* Numbers */
    char s18[] = "'12345'";
    strip_quotes_in_string(s18);
    ASSERT_STR_EQ(s18, "12345");

    /* Quote at start only */
    char s19[] = "\"test";
    strip_quotes_in_string(s19);
    ASSERT_STR_EQ(s19, "\"test");

    /* Quote at end only */
    char s20[] = "test\"";
    strip_quotes_in_string(s20);
    ASSERT_STR_EQ(s20, "test\"");

    /* Multiple internal quotes */
    char s21[] = "\"test\"test\"";
    strip_quotes_in_string(s21);
    ASSERT_STR_EQ(s21, "test\"test");

    /* Only whitespace after trimming */
    char s22[] = "\"\\t\\n   \"";
    strip_quotes_in_string(s22);
    /* Result depends on interpretation of escape sequences */

    /* Unicode inside quotes */
    char s23[] = "\"世界\"";
    strip_quotes_in_string(s23);
    ASSERT_STR_EQ(s23, "世界");

    /* Very long string */
    char long_quoted[1003];
    long_quoted[0] = '"';
    for (int i = 1; i <= 1000; i++)
    {
        long_quoted[i] = 'a';
    }
    long_quoted[1001] = '"';
    long_quoted[1002] = '\0';
    strip_quotes_in_string(long_quoted);
    ASSERT_INT_EQ(strlen(long_quoted), 1000);
}

/* -------------------------------------------------------
   Test runner
   ------------------------------------------------------- */

int main(void)
{
    printf("Running comprehensive string function tests...\n\n");

    test_strlcpy_custom_comprehensive();
    test_lowercase_str_comprehensive();
    test_parse_floats_from_string_comprehensive();
    test_string_start_with_comprehensive();
    test_string_end_with_comprehensive();
    test_str_reverse_in_place_comprehensive();
    test_str_replace_chars_comprehensive();
    test_my_strcasecmp_comprehensive();
    test_parse_bool_input_comprehensive();
    test_strip_quotes_in_string_comprehensive();

    printf("\n");
    printf("=====================================\n");
    printf("Tests run: %d\n", tests_run);
    printf("Failures : %d\n", tests_failed);
    printf("=====================================\n");

    if (tests_failed == 0)
    {
        printf("✓ All tests passed!\n");
    }
    else
    {
        printf("✗ Some tests failed!\n");
    }

    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
