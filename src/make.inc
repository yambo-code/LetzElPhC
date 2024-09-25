CC                  :=  mpicc

### add -DCOMPILE_ELPH_DOUBLE if you want to compile the code in double precession to cflags
#CFLAGS             := -O3 -march=native  -fopenmp
#LD_FLAGS           := -fopenmp 

##Debug flags, incase of testing always turn on the address sanitizer ! (clang is recommanded for this)
CFLAGS       	    := -g -fsanitize=address,undefined -fno-omit-frame-pointer # for mkl add 
LD_FLAGS            := -g -fsanitize=address,undefined -Wl,-no_compact_unwind

OPENMP_FLAGS   	    :=  #-DELPH_OMP_PARALLEL_BUILD ## uncomment for openmp build

FFTW_INC 	        :=  -I/usr/include 
FFTW3_LIB           :=  -L/usr/lib -lfftw3_threads -lfftw3f -lfftw3f_omp -lfftw3_omp -lfftw3

BLAS_LIB 	        :=  -L/usr/lib -llapack -lblas

NETCDF_INC          :=  -I/usr/include 
NETCDF_LIB 	        :=  -L/usr/lib -lnetcdf
HDF5_LIB            :=  -L/usr/lib  -lhdf5

INC_DIRS            := 
LIBS                := 


#### Notes Extra CFLAGS
### add -DCOMPILE_ELPH_DOUBLE if you want to compile the code in double precession
### for openmp use -DELPH_OMP_PARALLEL_BUILD in CFLAGS and set -fopenmp in LD_FLAGS and CFLAGS
