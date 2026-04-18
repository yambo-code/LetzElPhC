# Library paths and compiler are taken directly from yambo's configure output
# (config/setup, two levels above this plugin directory).
# Run `bash yambo_build_rr.sh` (or your configure script) before `make` here.
YAMBO_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../../..)
include $(YAMBO_ROOT)/config/setup

CC           := $(cc)

### if you are using yambo <= 5.1.2, add -DYAMBO_LT_5_1 to CFLAGS
### add -DCOMPILE_ELPH_DOUBLE to compile in double precision
CFLAGS       := -O2 -Wall -Wextra
LD_FLAGS     :=
OPENMP_FLAGS :=
#-DELPH_OMP_PARALLEL_BUILD  ## uncomment for openmp build

# FFTW: yambo's lfft/ifft cover fftw3; add fftw3f (single-precision) needed by LetzElPhC
FFTW_INC     := $(ifft)
FFTW3_LIB    := $(lfft) -lfftw3f

BLAS_LIB     := $(lblas)

# NetCDF: yambo's lnetcdf/inetcdf point to the parallel-enabled build (nc_open_par)
NETCDF_INC   := $(inetcdf)
NETCDF_LIB   := $(lnetcdf)

# HDF5: yambo links -lhdf5_fortran; LetzElPhC needs the C high-level API (-lhdf5_hl).
# Extract the -L flags from yambo's hdf5 config and substitute C-appropriate libraries.
HDF5_LIB     := $(filter -L%,$(lhdf5)) -lhdf5_hl -lhdf5 -lz

INC_DIRS     :=
LIBS         :=


#### Notes / Extra CFLAGS
### for openmp: add -DELPH_OMP_PARALLEL_BUILD to CFLAGS and -fopenmp to LD_FLAGS and CFLAGS
