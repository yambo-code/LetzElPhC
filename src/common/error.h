#pragma once

// some error macros

#define ERR_DIR_DOES_NOT_EXIST 1
#define ERR_NOT_A_DIRECTORY -1

#define ERR_FILE_OPEN_READ 1   // cannot open file to read
#define ERR_FILE_OPEN_WRITE 2  // cannot open file to write
#define ERR_FILE_COPY_FAIL 3   // copying files failed

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
