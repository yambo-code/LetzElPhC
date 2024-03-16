## brew install openmpi hdf5-mpi fftw 
# install netcdf manually. See at the end
#
CC                  :=  mpicc
CFLAGS              := -O3 

### if you are using yambo <= 5.1.2, Uncomment the below 
#CFLAGS             += -DYAMBO_LT_5_1

### If you need to compile with double precession support, uncomment the below
#CFLAGS             += -DCOMPILE_ELPH_DOUBLE

## set the following env variables FFTW3_DIR, BLAS_LIB_DIR, NETCDF_DIR, LAPACK_LIB_DIR, HDF5_DIR 
FFTW_INC 	        :=  -I/opt/homebrew/include 
FFTW3_LIB           :=  -L/opt/homebrew/lib -lfftw3_threads -lfftw3f -lfftw3f_omp -lfftw3_omp -lfftw3

BLAS_LIB 	        :=   -lblas -llapack
NETCDF_INC          :=  -I${NETCDF_DIR}/include 
NETCDF_LIB 	        :=  -L${NETCDF_DIR}/lib -lnetcdf
HDF5_LIB            :=  -L/opt/homebrew/lib -lhdf5



## to install netcdf manually
# Download the tar ball from https://github.com/Unidata/netcdf-c/releases
# untar and go to the folder
# mkdir install_dir
# ./configure CC=mpicc --prefix=$(pwd)/install_dir LIBS="-L/opt/homebrew/lib -lhdf5 -lz -lm -ldl -lcurl" CFLAGS="-I/opt/homebrew/include"; make install
#
# you NETCDF_DIR is $(pwd)/install_dir 
