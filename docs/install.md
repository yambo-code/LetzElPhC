# Installing the Code

## Mandatory Requirements

*   **GNU Make**
*   **C99 compiler** with complex number support (e.g., GCC, Clang, ICC, AMD C-Compiler, MinGW for Windows, PGI, or Arm C compilers).
*   **MPI implementation** supporting at least MPI-standard 2.1 (e.g., Open-MPI, MPICH, Intel MPI, Microsoft MPI).
*   **FFTW-3** or **Intel-MKL**.
*   **HDF5** and **NetCDF-C** libraries with **Parallel I/O support** (compiled with MPI).
*   A **BLAS library** (e.g., OpenBLAS, BLIS, Intel-MKL, or Atlas).

## Installation Process

LetzElPhC employs a standard make build system.

1.  **Prepare Makefile**:
    *   Navigate to the `sample_config` directory.
    *   Copy a sample makefile to the `src` directory and rename it to `make.inc`.
    *   Example: `cp sample_config/make_gcc_omp.inc src/make.inc`

2.  **Edit Configuration**:
    *   Navigate to the `src` directory.
    *   Edit `make.inc` to set your compiler, flags, and library paths.

3.  **Compile**:
    In the `src` directory, execute:
    ```bash
    make
    # To compile in parallel (e.g., 4 processes):
    make -j 4
    ```

Upon successful compilation, the **`lelphc`** executable will be created in the `src` directory.

### `make.inc` Variables Explanation

Below is an annotated example of a `make.inc` file. Use this to configure your build.

```make
CC                  :=  mpicc
# The MPI C compiler wrapper (e.g., mpicc, mpiicc).

CFLAGS              := -O3
# Compiler flags. -O3 activates high-level optimization.

LD_FLAGS            :=
# Linker flags. Add any necessary library paths or link options here.

# **** OPENMP BUILD (Optional) ***
# To enable multi-threading via OpenMP (recommended for performance):
# 1. Uncomment OPENMP_FLAGS line to define the macro.
# 2. Add the OpenMP compiler flag (e.g., -fopenmp for GCC, -qopenmp for Intel) to both CFLAGS and LD_FLAGS.
# OPENMP_FLAGS   	:= -DELPH_OMP_PARALLEL_BUILD
# CFLAGS            += -fopenmp 
# LD_FLAGS          += -fopenmp 

# **** DOUBLE PRECISION BUILD (Optional) ***
# By default, LetzElPhC runs in single precision to save memory.
# To run in double precision (higher accuracy, 2x memory usage):
# Uncomment the following line:
# CFLAGS              += -DCOMPILE_ELPH_DOUBLE

# FFTW3 include and libs
FFTW_INC 	        :=  -I/opt/homebrew/include
FFTW3_LIB           :=  -L/opt/homebrew/lib -lfftw3_threads -lfftw3f -lfftw3f_omp -lfftw3_omp -lfftw3
# Note: If compiling in single precision (default), link against the single precision FFTW libraries (often suffixed with 'f', e.g., -lfftw3f).
# If compiling in double precision, link against the standard double precision FFTW libraries (e.g., -lfftw3).

# BLAS and LAPACK
BLAS_LIB 	        :=  -L/opt/homebrew/opt/openblas/lib -lopenblas
# Link your BLAS/LAPACK provider. If separate, add both: -lblas -llapack

# NetCDF (Parallel IO required)
NETCDF_INC          :=  -I/path/to/netcdf/include
NETCDF_LIB 	        :=  -L/path/to/netcdf/lib -lnetcdf
# Ensure your NetCDF library was compiled with parallel I/O support (usually requires HDF5).

# HDF5
HDF5_LIB            :=  -L/path/to/hdf5/lib -lhdf5

# Extra
INC_DIRS            :=
LIBS                :=

```

!!! tip "Tip"
    If you have difficulty locating the required libraries, check your **YAMBO** installation directory. Open the `report` file located in the `config` directory; it lists all necessary libraries and include paths used by Yambo, which are generally compatible with LetzElPhC.
