TARGET 		 :=  elphc

CC               :=  mpicc
CFLAGS           :=  -O3 #-march=native
OPENMP_FLAGS     :=  #-DELPH_OMP_PARALLEL_BUILD

FFTW_INC 	 :=  -I/opt/homebrew/include 
FFTW3_LIB        :=  -L/opt/homebrew/lib -lfftw3_threads -lfftw3f -lfftw3f_omp -lfftw3_omp -lfftw3

BLAS_INC 	 :=  -I/opt/homebrew/Cellar/openblas/0.3.23/include 
BLAS_LIB 	 :=  -L/opt/homebrew/Cellar/openblas/0.3.23/lib -lopenblas 

NETCDF_INC       :=  -I/Users/murali/softwares/core/include 
NETCDF_LIB 	 :=  -L/Users/murali/softwares/core/lib -lnetcdf
HDF5_LIB         :=  -L/opt/homebrew/lib  -lhdf5


SUBDIR = nonloc io common dvloc elph wfc

INCS = $(wildcard *.h $(foreach fd, $(SUBDIR), $(fd)/*.h))
SRCS = $(wildcard *.c $(foreach fd, $(SUBDIR), $(fd)/*.c))
OBJS = $(addprefix ./, $(SRCS:c=o)) 

INC_DIRS = $(FFTW_INC) $(BLAS_INC) $(NETCDF_INC)
LIBS     = $(FFTW3_LIB) $(BLAS_LIB) $(NETCDF_LIB) $(HDF5_LIB) -lm


PHONY := $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $(OUTDIR)/$@ $(OBJS) $(LIBS)

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
