#pragma once 

typedef long long int ND_int; /* type used for all indices internally. */
#define nd_idx (ND_int []) /* macro for giving indices to ND_functions*/
#define ND_ALL ((void *)0)
#define ND_END 0
#define PRI_NdInt "%lld"

// Note : if we change ND_function and ND_array macros, change also in nd_array_src.h
#define ND_function(FUN_NAME, TYPE_SMALL)         ND_function_HIDDEN(FUN_NAME, TYPE_SMALL)
#define ND_function_HIDDEN(FUN_NAME, TYPE_SMALL)  nd_ ## FUN_NAME ## _ ## TYPE_SMALL

#define ND_array(TYPE_SMALL)         ND_array_HIDDEN(TYPE_SMALL)
#define ND_array_HIDDEN(TYPE_SMALL)  nd_arr ## _ ## TYPE_SMALL

