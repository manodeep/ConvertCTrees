### Set the default compiler -- options are icc/gcc/clang. 
ifneq (USE_MPI, $(findstring USE_MPI,$(OPT)))
CC:=gcc
else
CC:=mpicc
endif

#### Add any compiler specific flags you want
CFLAGS:=-ggdb 

#### Add any compiler specific link flags you want
CLINK:=

### Add include paths here
INCLUDE:=

### The POSIX_SOURCE flag is required to get the definition of strtok_r
CFLAGS += -Wsign-compare -Wall -Wextra -Wshadow -Wpadded -Wunused -std=c99 -g -m64 -fPIC  -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_BSD_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE -D_DARWIN_C_SOURCE -O3 -Ofast
GSL_CFLAGS := $(shell gsl-config --cflags) 
GSL_LIBDIR := $(shell gsl-config --prefix)/lib
GSL_LINK   := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR)

ifneq (DOUBLE_PREC,$(findstring DOUBLE_PREC,$(OPT)))
	VECTOR_TYPE:=float
else
	VECTOR_TYPE:=double
endif


ifeq (icc,$(findstring icc,$(CC)))
  CFLAGS += -xhost -opt-prefetch -opt-prefetch-distance=16 #-vec-report6  
  ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
		CFLAGS += -openmp
		CLINK  += -openmp 
  endif
else

  ### compiler specific flags for gcc
  ifeq (gcc,$(findstring gcc,$(CC)))
		CFLAGS += -ftree-vectorize -flto -funroll-loops #-ftree-vectorizer-verbose=6 -fopt-info-vec-missed #-fprofile-use -fprofile-correction
        CLINK += -flto
    ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
			CFLAGS += -fopenmp
			CLINK  += -fopenmp
    endif
  endif

  ### compiler specific flags for clang
  ifeq (clang,$(findstring clang,$(CC)))
		CFLAGS += -funroll-loops
		ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
      CLANG_VERSION:=$(shell $(CC) -dumpversion 2>&1)
      ifeq ($(CLANG_OMP_AVAIL),1)
			  CFLAGS += -fopenmp
			  CLINK  += -fopenmp=libomp
      else
        $(warning clang does not support OpenMP - please use gcc/icc for compiling with openmp. Removing USE_OMP from compile options)
        OPT:=$(filter-out -DUSE_OMP,$(OPT))
      endif
    endif
  endif

  ifeq (USE_AVX,$(findstring USE_AVX,$(OPT)))
    CFLAGS  +=  -mavx -mpopcnt
  endif

  #### common options for gcc and clang
  CFLAGS  += -march=native -fno-strict-aliasing
  CFLAGS  += -Wformat=2  -Wpacked  -Wnested-externs -Wpointer-arith  -Wredundant-decls  -Wfloat-equal -Wcast-qual  
  CFLAGS  +=  -Wcast-align #-Wmissing-declarations #-Wmissing-prototypes
  CFLAGS  += -Wnested-externs -Wstrict-prototypes  #-D_POSIX_C_SOURCE=2 -Wpadded -Wconversion
  CLINK += -lm
endif

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
 ## use the clang assembler instead of GNU assembler
 ## http://stackoverflow.com/questions/10327939/erroring-on-no-such-instruction-while-assembling-project-on-mac-os-x-lion
 ifeq (gcc,$(findstring gcc,$(CC)))
	 CFLAGS += -Wa,-q
 endif
endif

