include make.inc

TARGET 		 :=  lelphc

#### make start here

SUBDIR = nonloc io/ezxml io/inih io io/qe common dvloc elph wfc symmetries fft
ND_ARR_LIB = -L./nd_array -lnd_array


INCS = $(wildcard *.h $(foreach fd, $(SUBDIR), $(fd)/*.h))
SRCS = $(wildcard *.c $(foreach fd, $(SUBDIR), $(fd)/*.c))
OBJS = $(addprefix ./, $(SRCS:c=o)) 

INC_DIRS = $(FFTW_INC) $(BLAS_INC) $(NETCDF_INC)
LIBS     = $(ND_ARR_LIB) $(FFTW3_LIB) $(BLAS_LIB) $(NETCDF_LIB) $(HDF5_LIB) -lm


libnd_array = ./nd_array/libnd_array.a
libnd_make = ./nd_array/make.inc

PHONY := $(TARGET)
$(TARGET): $(libnd_make)  $(libnd_array) $(OBJS) 
	$(CC) -o ./$@ $(OBJS) $(LIBS) $(LD_FLAGS)

%.o: %.c $(INCS) make.inc
	$(CC) $(CFLAGS) -c $< $(INC_DIRS) $(OPENMP_FLAGS) -o $@

$(libnd_make): make.inc
	cp make.inc nd_array

$(libnd_array): $(libnd_make)
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
