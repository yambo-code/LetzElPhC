TARGET 		 :=  elphc

CC               :=  mpicc
CFLAGS           :=  -O1 -g  -fsanitize=address -fno-omit-frame-pointer #-O3 #-march=native
OPENMP_FLAGS     :=  #-DELPH_OMP_PARALLEL_BUILD

FFTW_INC 	 :=  -I/opt/homebrew/include 
FFTW3_LIB        :=  -L/opt/homebrew/lib -lfftw3_threads -lfftw3f -lfftw3f_omp -lfftw3_omp -lfftw3

BLAS_INC 	 :=  -I/opt/homebrew/Cellar/openblas/0.3.24/include 
BLAS_LIB 	 :=  -L/opt/homebrew/Cellar/openblas/0.3.24/lib -lopenblas 

NETCDF_INC       :=  -I/Users/murali/softwares/core/include 
NETCDF_LIB 	 :=  -L/Users/murali/softwares/core/lib -lnetcdf
HDF5_LIB         :=  -L/opt/homebrew/lib  -lhdf5

LD_FLAGS         :=  -g -fsanitize=address -Wl,-no_compact_unwind



#### make start here

SUBDIR = nonloc io/ezxml io common dvloc elph wfc symmetries fft
ND_ARR_LIB = -L./nd_array/src -lnd_array


INCS = $(wildcard *.h $(foreach fd, $(SUBDIR), $(fd)/*.h))
SRCS = $(wildcard *.c $(foreach fd, $(SUBDIR), $(fd)/*.c))
OBJS = $(addprefix ./, $(SRCS:c=o)) 

INC_DIRS = $(FFTW_INC) $(BLAS_INC) $(NETCDF_INC)
LIBS     = $(ND_ARR_LIB) $(FFTW3_LIB) $(BLAS_LIB) $(NETCDF_LIB) $(HDF5_LIB) -lm


PHONY := $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o ./$@ $(OBJS) $(LIBS) $(LD_FLAGS)

%.o: %.c $(INCS)
	$(CC) $(CFLAGS) -c $< $(INC_DIRS) $(OPENMP_FLAGS) -o $@

PHONY += clean
clean:
	rm ./$(TARGET) $(OBJS)

PHONY += echoes
echoes:
	@echo "INC files: $(INCS)"
	@echo "SRC files: $(SRCS)"
	@echo "OBJ files: $(OBJS)"
	@echo "LIB files: $(LIBS)"
	@echo "INC DIR: $(INC_DIRS)"


.PHONY = $(PHONY)
