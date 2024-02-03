CC               :=  mpicc
CPP              :=  gcc -E -P

### if you are using yambo <= 5.1.2, you need to add "-DYAMBO_LT_5_1" to cflags
#CFLAGS          :=  -O3 -march=native # -DND_MKL_BLAS -fopenmp

##Debug
CFLAGS       	 := -O1 -g -fsanitize=address -fno-omit-frame-pointer # for mkl add 

OPENMP_FLAGS   	 :=  #-DELPH_OMP_PARALLEL_BUILD ## uncomment for openmp build

FFTW_INC 	 :=  -I/opt/homebrew/include 
FFTW3_LIB        :=  -L/opt/homebrew/lib -lfftw3_threads -lfftw3f -lfftw3f_omp -lfftw3_omp -lfftw3

### If you donot have cblas.h header leave BLAS_INC empty and add "-DND_DEFAULT_BLAS_HEADER" to cflags
BLAS_INC 	 :=  -I/opt/homebrew/Cellar/openblas/0.3.24/include   
BLAS_LIB 	 :=  -L/opt/homebrew/Cellar/openblas/0.3.24/lib -lopenblas 

NETCDF_INC       :=  -I/Users/murali/softwares/core/include 
NETCDF_LIB 	 :=  -L/Users/murali/softwares/core/lib -lnetcdf
HDF5_LIB         :=  -L/opt/homebrew/lib  -lhdf5

INC_DIRS         := # -I/Users/murali/softwares/tblis/include
LIBS             :=  #-L/Users/murali/softwares/tblis/include -ltblis

#LD_FLAGS         := #-fopenmp 

#Debug
LD_FLAGS         := -g -fsanitize=address -Wl,-no_compact_unwind


#### Notes Extra CFLAGS
### if you are using yambo <= 5.1.2, you need to add "-DYAMBO_LT_5_1" to cflags
## if using intel mkl use -DND_MKL_BLAS in clags
### for openmp use -DELPH_OMP_PARALLEL_BUILD in CFLAGS and set -fopenmp in LD_FLAGS and CFLAGS
