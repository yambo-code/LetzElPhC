// we define some constants here
//
//
#pragma once
#include "../elphC.h"

#if defined(COMPILE_ELPH_DOUBLE)

#define ELPH_EPS 1e-6
#define ELPH_PI 3.1415927
#define ELPH_SQRT2 1.4142136
#define ELPH_e2 2.0  // in Ry units

#else

#define ELPH_EPS 1e-6
#define ELPH_PI 3.1415927f
#define ELPH_SQRT2 1.4142136f
#define ELPH_e2 2.0f  // in Ry units
#endif