include ../make.inc

TARGET 		 :=  libnd_array.a

#### make start here

SUBDIR = src

INC_DIRS += $(FFTW_INC) $(BLAS_INC) $(NETCDF_INC)

### Any includes must be included here
INCS = $(SUBDIR)/nd_error.h $(SUBDIR)/nd_array.h $(SUBDIR)/common_def.h

## sources will be automatically detected
SRCS = $(wildcard *.c $(foreach fd, $(SUBDIR), $(fd)/*.c))
HEADER_GEN_SRC = $(SUBDIR)/nd_array_src.h

strip_src := $(foreach F,$(SRCS),$(lastword $(subst /, ,$F)))

OBJ_I =  $(addprefix ND_INT_, $(strip_src:c=o)) 
OBJ_F = $(addprefix  ND_FLOAT_, $(strip_src:c=o))
OBJ_D = $(addprefix  ND_DOUBLE_, $(strip_src:c=o)) 
OBJ_FC = $(addprefix ND_SINGLE_COMPLEX_, $(strip_src:c=o)) 
OBJ_DC = $(addprefix ND_DOUBLE_COMPLEX_, $(strip_src:c=o)) 

OBJS = $(OBJ_I) $(OBJ_F) $(OBJ_D) $(OBJ_FC) $(OBJ_DC)

ND_INT_%.o: $(SUBDIR)/%.c $(SRCS) $(HEADER_GEN_SRC)
	$(CC) $(CFLAGS) -c ./$< $(INC_DIRS) -DCOMPILE_ND_INT $(OPENMP_FLAGS) -o ./$(SUBDIR)/$@

ND_FLOAT_%.o: $(SUBDIR)/%.c $(SRCS) $(HEADER_GEN_SRC)
	$(CC) $(CFLAGS) -c ./$< $(INC_DIRS) -DCOMPILE_ND_FLOAT $(OPENMP_FLAGS) -o ./$(SUBDIR)/$@
	
ND_DOUBLE_%.o: $(SUBDIR)/%.c $(SRCS) $(HEADER_GEN_SRC)
	$(CC) $(CFLAGS) -c ./$< $(INC_DIRS) -DCOMPILE_ND_DOUBLE -DCOMPILE_ONCE $(OPENMP_FLAGS) -o ./$(SUBDIR)/$@
	
ND_SINGLE_COMPLEX_%.o: $(SUBDIR)/%.c $(SRCS) $(HEADER_GEN_SRC)
	$(CC) $(CFLAGS) -c ./$< $(INC_DIRS) -DCOMPILE_ND_SINGLE_COMPLEX $(OPENMP_FLAGS) -o ./$(SUBDIR)/$@
	
ND_DOUBLE_COMPLEX_%.o: $(SUBDIR)/%.c $(SRCS) $(HEADER_GEN_SRC)
	$(CC) $(CFLAGS) -c ./$< $(INC_DIRS) -DCOMPILE_ND_DOUBLE_COMPLEX $(OPENMP_FLAGS) -o ./$(SUBDIR)/$@

PHONY := $(TARGET)

ND_HEADERS = nd_INT.h nd_FLOAT.h nd_DOUBLE.h nd_SINGLE_COMPLEX.h nd_DOUBLE_COMPLEX.h

$(TARGET): $(OBJS) $(INCS)
	$(CPP) $(HEADER_GEN_SRC) -DCOMPILE_ND_INT -DGENERATE_HEADERS > ./$(SUBDIR)/nd_INT.h
	$(CPP) $(HEADER_GEN_SRC) -DCOMPILE_ND_FLOAT -DGENERATE_HEADERS > ./$(SUBDIR)/nd_FLOAT.h
	$(CPP) $(HEADER_GEN_SRC) -DCOMPILE_ND_DOUBLE -DGENERATE_HEADERS -DCOMPILE_ONCE > ./$(SUBDIR)/nd_DOUBLE.h
	$(CPP) $(HEADER_GEN_SRC) -DCOMPILE_ND_SINGLE_COMPLEX -DGENERATE_HEADERS > ./$(SUBDIR)/nd_SINGLE_COMPLEX.h
	$(CPP) $(HEADER_GEN_SRC) -DCOMPILE_ND_DOUBLE_COMPLEX -DGENERATE_HEADERS > ./$(SUBDIR)/nd_DOUBLE_COMPLEX.h
	ar rcs ./$(SUBDIR)/$@ ./$(SUBDIR)/*.o

clean_file = $(TARGET) $(OBJS) $(ND_HEADERS)
clean_file := $(addprefix ./$(SUBDIR)/,$(clean_file)) 

PHONY += clean
clean:
	rm $(clean_file)

.PHONY = $(PHONY)
