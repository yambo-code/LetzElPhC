#include <netcdf.h>
#include <netcdf_par.h>
#include <stdio.h>
#include "elphC.h"
#include "common/error.h"
#include "nc_utils.h"

void quick_read(const int ncid, char* var_name, void* data_out)
{
    int varid, retval;
    if ((retval = nc_inq_varid(ncid, var_name, &varid)))
        ERR(retval);
    if ((retval = nc_get_var(ncid, varid, data_out)))
        ERR(retval);
}

void quick_read_float(const int ncid, char* var_name, ELPH_float* data_out)
{
    int varid, retval;
    if ((retval = nc_inq_varid(ncid, var_name, &varid)))
        ERR(retval);

#if defined(COMPILE_ELPH_DOUBLE)
    if ((retval = nc_get_var_double(ncid, varid, data_out)))
        ERR(retval);
#else
    if ((retval = nc_get_var_float(ncid, varid, data_out)))
        ERR(retval);
#endif
}

void quick_read_sub(const int ncid, char* var_name, const size_t* startp,
                    const size_t* countp, void* data_out)
{
    int varid, retval;
    if ((retval = nc_inq_varid(ncid, var_name, &varid)))
        ERR(retval);

    if ((retval = nc_var_par_access(ncid, varid, NC_COLLECTIVE)))
        ERR(retval);

    if ((retval = nc_get_vara(ncid, varid, startp, countp, data_out)))
        ERR(retval);
}
