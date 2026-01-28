CC                  :=  mpicc

### if you are using yambo <= 5.1.2, you need to add "-DYAMBO_LT_5_1" to cflags
### add -DCOMPILE_ELPH_DOUBLE if you want to compile the code in double precession to cflags
#CFLAGS             := -O3 -Wall -Wextra -march=native  -fopenmp
#LD_FLAGS           := -fopenmp 

##Debug flags, incase of testing always turn on the address sanitizer ! (clang is recommanded for this)
CFLAGS       	    := -g -Wall -Wextra # -fsanitize=address,undefined -fno-omit-frame-pointer # for mkl add 
LD_FLAGS            := -g  # -fsanitize=address,undefined  -Wl,-no_compact_unwind

OPENMP_FLAGS   	    :=  
#-DELPH_OMP_PARALLEL_BUILD ## uncomment for openmp build

FFTW_INC 	    :=  -I${YAMBO_EXT_LIBS}/gfortran/mpifort/include/
FFTW3_LIB           :=  -L${YAMBO_EXT_LIBS}/gfortran/mpifort/lib -lfftw3f -lfftw3
#FFTW_INC 	    :=  -I/opt/homebrew/include 
#FFTW3_LIB           :=  -L/opt/homebrew/lib -lfftw3_threads -lfftw3f -lfftw3f_omp -lfftw3_omp -lfftw3

#BLAS_LIB 	        :=  -L/opt/homebrew/opt/openblas/lib -lopenblas 
BLAS_LIB 	        :=  -L${YAMBO_EXT_LIBS}/gfortran/mpifort/lib -lopenblas 

NETCDF_INC          :=  -I${YAMBO_EXT_LIBS}/gfortran/mpifort/v4/parallel/include
NETCDF_LIB 	    :=  -L${YAMBO_EXT_LIBS}/gfortran/mpifort/v4/parallel/lib -lnetcdf
HDF5_LIB            :=  -L${YAMBO_EXT_LIBS}/gfortran/mpifort/v4/parallel/lib -lhdf5_hl -lhdf5 -lz
#NETCDF_INC          :=  -I/Users/murali/softwares/core/include 
#NETCDF_LIB 	    :=  -L/Users/murali/softwares/core/lib -lnetcdf
#HDF5_LIB            :=  -L/opt/homebrew/lib  -lhdf5

INC_DIRS            := 
LIBS                := 


#### Notes Extra CFLAGS
### add -DCOMPILE_ELPH_DOUBLE if you want to compile the code in double precession
### if you are using yambo <= 5.1.2, you need to add "-DYAMBO_LT_5_1" to cflags
### for openmp use -DELPH_OMP_PARALLEL_BUILD in CFLAGS and set -fopenmp in LD_FLAGS and CFLAGS
