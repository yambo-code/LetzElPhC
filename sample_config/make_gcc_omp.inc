CC                  :=  mpicc
CFLAGS              := -O3 

### if you are using yambo <= 5.1.2, Uncomment the below 
#CFLAGS             += -DYAMBO_LT_5_1

### If you need to compile with double precession support, uncomment the below
#CFLAGS             += -DCOMPILE_ELPH_DOUBLE


### Comment the below 3 lines if you "donot" want openmp build
CFLAGS              += -fopenmp
LD_FLAGS            += -fopenmp 
OPENMP_FLAGS   	    := -DELPH_OMP_PARALLEL_BUILD 

## set the following env variables FFTW3_DIR, BLAS_LIB_DIR, NETCDF_DIR, LAPACK_LIB_DIR, HDF5_DIR 


FFTW_INC 	        :=  -I${FFTW3_DIR}/include
FFTW3_LIB           :=  -L${FFTW3_DIR}/lib -lfftw3_threads -lfftw3f -lfftw3f_omp -lfftw3_omp -lfftw3
## you can remove -lfftw3_threads -lfftw3f_omp -lfftw3_omp if you donot 
## want to build with openmp support
## refer https://www.fftw.org/fftw3_doc/Precision.html

BLAS_LIB 	        :=  -L${BLAS_LIB_DIR} -lopenblas 
# BLAS_LIB 	        :=  -L${BLAS_LIB_DIR} -lblas -L${LAPACK_LIB_DIR} -llapack
NETCDF_INC          :=  -I${NETCDF_DIR}/include 
NETCDF_LIB 	        :=  -L${NETCDF_DIR}/lib -lnetcdf
HDF5_LIB            :=  -L${HDF5_DIR}/lib -lhdf5

