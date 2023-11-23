include make.inc

TARGET 		 :=  elphc

#### make start here

SUBDIR = nonloc io/ezxml io common dvloc elph wfc symmetries fft
ND_ARR_LIB = -L./nd_array/src -lnd_array


INCS = $(wildcard *.h $(foreach fd, $(SUBDIR), $(fd)/*.h))
SRCS = $(wildcard *.c $(foreach fd, $(SUBDIR), $(fd)/*.c))
OBJS = $(addprefix ./, $(SRCS:c=o)) 

INC_DIRS = $(FFTW_INC) $(BLAS_INC) $(NETCDF_INC)
LIBS     = $(ND_ARR_LIB) $(FFTW3_LIB) $(BLAS_LIB) $(NETCDF_LIB) $(HDF5_LIB) -lm


libnd_array = ./nd_array/src/libnd_array.a


PHONY := $(TARGET)
$(TARGET): $(libnd_array) $(OBJS) 
	$(CC) -o ./$@ $(OBJS) $(LIBS) $(LD_FLAGS)

%.o: %.c $(INCS) 
	$(CC) $(CFLAGS) -c $< $(INC_DIRS) $(OPENMP_FLAGS) -o $@

$(libnd_array):
	$(MAKE) -C nd_array

PHONY += clean
clean:
	rm ./$(TARGET) $(OBJS)
	$(MAKE) -C nd_array clean

PHONY += echoes
echoes:
	@echo "INC files: $(INCS)"
	@echo "SRC files: $(SRCS)"
	@echo "OBJ files: $(OBJS)"
	@echo "LIB files: $(LIBS)"
	@echo "INC DIR: $(INC_DIRS)"


.PHONY = $(PHONY)
