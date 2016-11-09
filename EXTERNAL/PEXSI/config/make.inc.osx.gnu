#!/usr/bin/bash
COMPILE_MODE     = release
USE_PROFILE      = 0
PAR_ND_LIBRARY   = ptscotch
SEQ_ND_LIBRARY   = scotch
USE_SYMPACK      = 0

# Different compiling and linking options.
SUFFIX       = osx_v0.10.0

#compiler, options and tools
################################################################
CC           = mpicc 
CXX          = mpic++
FC           = mpif90
LOADER       = mpic++


AR           = ar 
ARFLAGS      = rvcu
# For System V based machine without ranlib, like Cray and SGI,
# use touch instead.
#RANLIB      = touch
RANLIB       = ranlib

CP           = cp
RM           = rm
RMFLAGS      = -f
################################################################


#pexsi directory
PEXSI_DIR     = $(HOME)/Projects/pexsi
PEXSI_BUILD_DIR = $(PEXSI_DIR)/build

# Required libraries directories
DSUPERLU_DIR  = $(HOME)/Documents/Software/SuperLU_DIST_5.1.2
METIS_DIR     = $(HOME)/Software/metis-5.1.2/build_release
PARMETIS_DIR  = $(HOME)/Software/parmetis-4.0.2/build/Darwin-x86_64
PTSCOTCH_DIR  = $(HOME)/Software/scotch_6.0.0

#useful on some systems (OSX using homebrew compilers)
GFORTRAN_DIR  = /usr/local/Cellar/gfortran/4.8.2/gfortran/lib

# Graph partitioning libraries
SEQ_ND_DIR  = $(HOME)/software/release/scotch_6.0.0
SEQ_ND_LIB     = -L${SEQ_ND_DIR}/lib -lscotchmetis -lscotch -lscotcherr
PAR_ND_DIR  = $(HOME)/software/release/scotch_6.0.0
PAR_ND_LIB     = -L${PAR_ND_DIR}/lib -lptscotchparmetis -lptscotch -lptscotcherr -lscotch

# Includes
PEXSI_INCLUDE    = -I${PEXSI_DIR}/include 
DSUPERLU_INCLUDE = -I${DSUPERLU_DIR}/SRC
INCLUDES         = ${PEXSI_INCLUDE} ${DSUPERLU_INCLUDE} 

# Libraries
CPP_LIB          = -lstdc++ -lmpi -lmpi_cxx
GFORTRAN_LIB     = /usr/local/lib/libgfortran.dylib 
LAPACK_LIB       = -llapack
BLAS_LIB         = -lblas
DSUPERLU_LIB     = ${DSUPERLU_DIR}/build_release/libsuperlu_dist_5.1.2.a
PEXSI_LIB        = ${PEXSI_DIR}/src/libpexsi_${SUFFIX}.a


################ End of configuration #####################


# Different compiling and linking options.
ifeq (${COMPILE_MODE}, release)
  COMPILE_DEF    = -DRELEASE
  COMPILE_FLAG   = $(RELEASE_COMPILE_FLAG)
endif
ifeq (${COMPILE_MODE}, debug)
  COMPILE_DEF    = -DDEBUG=1
  COMPILE_FLAG   = -O0 -w -g
endif

LIBS  = ${PEXSI_LIB} ${DSUPERLU_LIB} ${PAR_ND_LIB} ${SEQ_ND_LIB} ${LAPACK_LIB} ${BLAS_LIB} ${GFORTRAN_LIB}

COMPILE_DEF += -DAdd_ -std=c++11

LIBS  = ${PEXSI_LIB} ${DSUPERLU_LIB} ${PAR_ND_LIB} ${SEQ_ND_LIB} ${LAPACK_LIB} ${BLAS_LIB} ${GFORTRAN_LIB}
COMPILE_DEF  += -DAdd_
CPPFLAG = -std=c++11


ifeq (${USE_SYMPACK}, 1)
  #symPACK related definitions
  SYMPACK_DIR = ${HOME}/sympack_install
  include ${SYMPACK_DIR}/include/sympack.mak
  INCLUDES += ${SYMPACK_INCLUDE} 
  LIBS+= ${SYMPACK_LIB} ${LAPACK_LIB} ${BLAS_LIB} ${GFORTRAN_LIB}
  COMPILE_DEF  += -DWITH_SYMPACK
  CPPFLAG += -std=c++11
endif



CFLAGS       = ${COMPILE_FLAG} ${PROFILE_FLAG} ${INCLUDES} -std=c99
FFLAGS       = ${COMPILE_FLAG} ${PROFILE_FLAG} ${INCLUDES}
CXXFLAGS     = ${COMPILE_FLAG} ${CPPFLAG} ${PROFILE_FLAG} ${INCLUDES} 
CCDEFS       = ${COMPILE_DEF} 
CPPDEFS      = ${COMPILE_DEF} 
LOADOPTS     = ${PROFILE_FLAG} ${LIBS} -Wl,--allow-multiple-definition
#FLOADOPTS    = ${LIBS} -L/usr/local/lib  -lstdc++ -Wl,--allow-multiple-definition 


# Generate auto-dependencies 
%.d: %.c
	@set -e; rm -f $@; \
	$(CC) -M $(CCDEFS) $(CFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CXX) -M $(CPPDEFS) $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$
