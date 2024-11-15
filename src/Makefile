include make.inc

TARGET 		 :=  lelphc

#### make start here

SUBDIR = nonloc io/ezxml io/inih io \
		 io/qe common dvloc elph wfc \
		 symmetries fft common/cwalk \
		 preprocessor interpolation

.SUFFIXES :

INCS = $(wildcard *.h $(foreach fd, $(SUBDIR), $(fd)/*.h))
SRCS = $(wildcard *.c $(foreach fd, $(SUBDIR), $(fd)/*.c))
OBJS = $(addprefix ./, $(SRCS:c=o)) 

INC_DIRS += $(FFTW_INC) $(NETCDF_INC)
LIBS     += $(FFTW3_LIB) $(BLAS_LIB) $(NETCDF_LIB) $(HDF5_LIB) -lm


PHONY := $(TARGET)
$(TARGET): $(OBJS) 
	$(CC) -o ./$@ $(OBJS) $(LIBS) $(LD_FLAGS)

%.o: %.c $(INCS) make.inc
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
