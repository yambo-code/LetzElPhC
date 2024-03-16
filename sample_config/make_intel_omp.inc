CC                  :=  mpiicc
CFLAGS              := -O3 

### if you are using yambo <= 5.1.2, Uncomment the below 
#CFLAGS             += -DYAMBO_LT_5_1

### If you need to compile with double precession support, uncomment the below
#CFLAGS             += -DCOMPILE_ELPH_DOUBLE


### Comment the below 3 lines if you "donot" want openmp build
CFLAGS              += -qopenmp
LD_FLAGS            += -qopenmp 
OPENMP_FLAGS   	    := -DELPH_OMP_PARALLEL_BUILD 


## set the following env variables FFTW3_DIR, MKLROOT, NETCDF_DIR, HDF5_DIR 

FFTW_INC 	        :=  -I${MKLROOT}/include/fftw
FFTW3_LIB           :=  

BLAS_LIB 	        :=  -L${MKLROOT}/lib -lmkl_rt -lpthread -lm -ldl
NETCDF_INC          :=  -I${NETCDF_DIR}/include 
NETCDF_LIB 	        :=  -L${NETCDF_DIR}/lib -lnetcdf
HDF5_LIB            :=  -L${HDF5_DIR}/lib -lhdf5

