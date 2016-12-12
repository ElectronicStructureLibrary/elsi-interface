#!/usr/bin/bash
COMPILE_MODE = release
USE_PROFILE  = 0

# Compiler and tools
CC      = $(MPICC)
CXX     = $(MPICXX)
FC      = $(MPIFC)
LOADER  = $(MPICXX)

AR      = $(ARCH)
ARFLAGS = $(ARCHFLAGS)
RANLIB  = ranlib

CP      = cp
RM      = rm
RMFLAGS = -f

# Includes
PEXSI_INCLUDE = -I${PEXSI_DIR}/include
INCLUDES      = ${PEXSI_INCLUDE} ${SUPERLU_INCLUDE}

# Libraries
CPP_LIB       = $(CXX_LDFLAGS)
SCALAPACK_LIB = $(SCALAPACK_LDFLAGS)
PEXSI_LIB     = ${PEXSI_DIR}/lib/libpexsi.a

# Different compiling and linking options
ifeq (${COMPILE_MODE}, release)
  COMPILE_DEF  = -DRELEASE
  COMPILE_FLAG = -O3 -w -g
endif
ifeq (${COMPILE_MODE}, debug)
  COMPILE_DEF  = -DDEBUG=1
  COMPILE_FLAG = -O2 -w -g
endif

ifeq (${PAR_ND_LIBRARY}, ptscotch)
  PAR_ND_LIB = ${PTSCOTCH_LIB}
else
  PAR_ND_LIB = ${PARMETIS_LIB}
endif 

ifeq (${SEQ_ND_LIBRARY}, scotch)
  SEQ_ND_LIB = ${SCOTCH_LIB}
else
  SEQ_ND_LIB = ${METIS_LIB}
endif 

ifeq (${USE_PROFILE}, 1)
  PROFILE_FLAG = -DPROFILE
endif

LIBS = ${PEXSI_LIB} ${SUPERLU_LIB} ${PARMETIS_LIB} ${METIS_LIB} ${SCALAPACK_LIB}

COMPILE_DEF += -DAdd_

CPPFLAG = -std=c++11

CFLAGS    = ${CFLAGS_P} ${PROFILE_FLAG} ${INCLUDES} -std=c99
FFLAGS    = ${FFLAGS_P} ${PROFILE_FLAG} ${INCLUDES}
CXXFLAGS  = ${CXXFLAGS_P} ${CPPFLAG} ${PROFILE_FLAG} ${INCLUDES}
CCDEFS    = ${COMPILE_DEF}
CPPDEFS   = ${COMPILE_DEF}
LOADOPTS  = ${PROFILE_FLAG} ${LIBS} -Wl,--allow-multiple-definition
FLOADOPTS = ${PROFILE_FLAG} ${LIBS} ${CPP_LIB} -Wl,--allow-multiple-definition

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
