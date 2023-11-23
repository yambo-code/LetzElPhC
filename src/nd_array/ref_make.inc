CC               :=  gcc-13
CPP              :=  $(CC) -E -P
CFLAGS           :=  -O3 -march=native #-O1 -g  -fsanitize=address -fno-omit-frame-pointer 
# for mkl add -DND_MKL_BLAS

OPENMP_FLAGS     :=  # 

FFTW_INC 	 :=  -I/opt/homebrew/include 
FFTW3_LIB        :=  -L/opt/homebrew/lib -lfftw3_threads -lfftw3f -lfftw3f_omp -lfftw3_omp -lfftw3

BLAS_INC 	 :=  -I/opt/homebrew/Cellar/openblas/0.3.24/include 
BLAS_LIB 	 :=  -L/opt/homebrew/Cellar/openblas/0.3.24/lib -lopenblas 

NETCDF_INC       :=  -I/Users/murali/softwares/core/include 
NETCDF_LIB 	 :=  -L/Users/murali/softwares/core/lib -lnetcdf
HDF5_LIB         :=  -L/opt/homebrew/lib  -lhdf5

INC_DIRS         := # -I/Users/murali/softwares/tblis/include
LIBS             :=  #-L/Users/murali/softwares/tblis/include -ltblis

LD_FLAGS         :=  #-g -fsanitize=address -Wl,-no_compact_unwind



#### Extra CFLAGS
## -DCOMPILE_ND_TBLIS for tblis and einsum support
## if using intel mkl use -DND_MKL_BLAS in clags
