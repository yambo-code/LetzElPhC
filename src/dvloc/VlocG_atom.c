#include <math.h>

#include "../common/constants.h"
#include "../common/numerical_func.h"
#include "../elphC.h"
#include "dvloc.h"

/*Function to compute local short range bare potential of Each atom placed at
 * origin */
ELPH_float Vloc_Gspace(ELPH_float* work_arr, const char cutoff,
                       const ELPH_float Gnorm, const ND_int ngrid,
                       const ELPH_float* Vloc_atomic, const ELPH_float* r_grid,
                       const ELPH_float* rab_grid, const ELPH_float Zval,
                       const ELPH_float eta, const ELPH_float volume)
{
    /* work_arr is the work array. dims = (ngrid)
    We do not call malloc/free inside this function as this function is called
    with a very high frequency and can lead to fragmentation of memory. Gnorm ==
    0 only whene q = 0 and G = 0

    */
    ELPH_float VlocGq;
    /* PseudoPot is spilt into long range and short range */
    /*
    Lt x-> 0 erf(r)/r = 2/sqrt(pi)
    */
    ELPH_float sqrt_eta = sqrt(eta);
    /* Compute the integrate i.e V_{SR}(r)*/
    for (ND_int imesh = 0; imesh < ngrid; ++imesh)
    {
        ELPH_float fac, radius;
        radius = r_grid[imesh];

        if (Gnorm < ELPH_EPS)
        {
            fac = radius;
        }
        else
        {
            fac = sin(radius * Gnorm) / Gnorm;
        }

        if (cutoff == '3' && Gnorm < ELPH_EPS)
        {
            fac *= (radius * Vloc_atomic[imesh] + ELPH_e2 * Zval);
        }
        else
        {
            fac *= (radius * Vloc_atomic[imesh] +
                    ELPH_e2 * Zval * erf(radius * sqrt_eta));
        }
        work_arr[imesh] = fac;
    }
    VlocGq = simpson(work_arr, rab_grid, ngrid);

    return 4 * ELPH_PI * VlocGq / volume;
}
