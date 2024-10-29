#pragma once

/* Function to print error message to the the file */
#define MPI_error_msg(err_code)                                         \
    {                                                                   \
        if (err_code != MPI_SUCCESS)                                    \
        {                                                               \
            ELPH_MPI_error_msg(err_code, __FILE__, __LINE__, __func__); \
        }                                                               \
    }

#define error_msg(print_str) \
    elph_error_msg(print_str, __FILE__, __LINE__, __func__)

#define CHECK_ALLOC(ptr)                                     \
    {                                                        \
        if (!(ptr))                                          \
        {                                                    \
            error_msg("Failed to allocate " #ptr " buffer"); \
        }                                                    \
    }

// Netcdf error macro
#define ERR(e)                                          \
    {                                                   \
        fprintf(stderr, "Error: %s\n", nc_strerror(e)); \
        error_msg("netcdf_error");                      \
    }

void elph_error_msg(const char* error_msg, const char* file,
                    const long long int line, const char* func_name);

void ELPH_MPI_error_msg(int err_code, const char* file,
                        const long long int line, const char* func_name);

// MPI_error_msg(mpi_error);
// if (buf == NULL) {error_msg("Failed to allocate buf array");}
