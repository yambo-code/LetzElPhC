#include "dvloc.h"

/*Function to compute local bare potential of Each atom placed at origin */
ELPH_float Vloc_Gspace(ELPH_float * work_arr, const char cutoff, const ELPH_float Gnorm, \
            const ND_array(Nd_floatS)* Vloc_atomic, const ND_array(Nd_floatS)* r_grid, \
            const ND_array(Nd_floatS)* rab_grid, const ELPH_float Zval,const ELPH_float eta, \
            const ELPH_float cutoff_fac)
{   
    /* work_arr is the work array. dims = (ngrid)
    We do not call malloc/free inside this function as this function is called with a very high
    frequency and can lead to fragmentation of memory.
    Gnorm == 0 only whene q = 0 and G = 0

    Note the output must be divided by volume of the cell !
    */
    ELPH_float VlocGq;
    ND_int ngrid = r_grid->dims[0];

    /* PseudoPot is spilt into long range and short range */

    /*
    Lt x-> 0 erf(r)/r = 2/sqrt(pi)
    */
    ELPH_float sqrt_eta = sqrt(eta);
    /* Compute the integrate i.e V_{SR}(r)*/
    for (ND_int imesh = 0; imesh<ngrid; ++imesh)
    {   
        ELPH_float fac, radius;
        radius = r_grid->data[imesh];
        
        if (Gnorm < ELPH_EPS) fac = radius;
        else fac = sin(radius*Gnorm)/Gnorm ;

        if (cutoff == '3' && Gnorm < ELPH_EPS)
        {
            fac *= (radius*Vloc_atomic->data[imesh] + ELPH_e2*Zval);
        }
        else fac *= (radius*Vloc_atomic->data[imesh] + ELPH_e2*Zval*erf(radius*sqrt_eta));
        
        work_arr[imesh] = fac;
    }
    VlocGq = simpson(work_arr, rab_grid->data, ngrid);

    /* Add the long range part again */
    /* We set V_LR(G=0) = 0  as it is cancelled by the hartree term */
    if (Gnorm >= ELPH_EPS)
    {   
        VlocGq -= ELPH_e2*Zval*exp(-Gnorm*Gnorm*0.25/eta)*cutoff_fac/(Gnorm*Gnorm) ;
    }
    
    /* Note the output must be multiplied with volume */
    return 4*ELPH_PI*VlocGq ;

}