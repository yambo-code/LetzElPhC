/*
 * NetCDF utility functions for reading variables
 */
#pragma once

#include <netcdf.h>
#include "elphC.h"

/* Read entire variable from NetCDF file */
void quick_read(const int ncid, char* var_name, void* data_out);

/* Read entire ELPH_float variable from NetCDF file */
void quick_read_float(const int ncid, char* var_name, ELPH_float* data_out);

/* Read slice of variable from NetCDF file (collective I/O) */
void quick_read_sub(const int ncid, char* var_name, const size_t* startp,
                    const size_t* countp, void* data_out);
