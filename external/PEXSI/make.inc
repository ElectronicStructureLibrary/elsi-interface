#!/usr/bin/bash
COMPILE_MODE     = release
USE_PROFILE      = 0
#PAR_ND_LIBRARY   = ptscotch
#SEQ_ND_LIBRARY   = scotch

# Compiler and tools
################################################################
CC           = $(MPICC)
CXX          = $(MPICXX)
FC           = $(MPIFC)

AR           = $(ARCH)
ARFLAGS      = $(ARCHFLAGS)

CP           = cp
RM           = rm
RMFLAGS      = -f
################################################################

# PEXSI directory
PEXSI_DIR     = $(THIS_DIR)/external/PEXSI

# Required libraries directories
DSUPERLU_DIR  = $(THIS_DIR)/external/PEXSI/external
PARMETIS_DIR  = $(THIS_DIR)/external/PEXSI/external

# Includes
PEXSI_INCLUDE    = -I$(PEXSI_DIR)/include 
DSUPERLU_INCLUDE = -I$(DSUPERLU_DIR)/include
INCLUDES         = $(PEXSI_INCLUDE) $(DSUPERLU_INCLUDE)

# Libraries
SCALAPACK_LIB    = $(SCALAPACK) $(LAPACK) $(BLAS)
METIS_LIB        = -L$(PARMETIS_DIR)/lib -lmetis
PARMETIS_LIB     = -L$(PARMETIS_DIR)/lib -lparmetis
#SCOTCH_LIB       = -L$(SCOTCH_DIR)/lib -lscotchmetis -lscotcherrexit \
#                   -lscotcherr -lscotch
DSUPERLU_LIB     = $(DSUPERLU_DIR)/lib/libsuperlu.a
PEXSI_LIB        = $(PEXSI_DIR)/lib/libpexsi.a

# Different compiling and linking options.
ifeq (${COMPILE_MODE}, release)
  COMPILE_DEF    = -DDEBUG=0 -DRELEASE
  COMPILE_FLAG   = -O3 -w 
endif
ifeq (${COMPILE_MODE}, debug)
  COMPILE_DEF    = -DDEBUG=1
  COMPILE_FLAG   = -O0 -g -w
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
  PROFILE_FLAG  = -DPROFILE
endif

LIBS  = ${PEXSI_LIB} ${DSUPERLU_LIB} ${PAR_ND_LIB} ${SEQ_ND_LIB} ${SCALAPACK_LIB} ${BLACS_LIB} ${FORTRAN_LIB}

COMPILE_DEF  += -DAdd_
LOADER = $(MPICXX)
CPPFLAG = -std=c++0x
FLOADOPTS    = ${LIBS} -lstdc++

CFLAGS       = ${CFLAGS_P} ${PROFILE_FLAG} ${INCLUDES}
FFLAGS       = ${FFLAGS_P} ${PROFILE_FLAG} ${INCLUDES}
CXXFLAGS     = ${CXXFLAGS_P} ${CPPFLAG} ${PROFILE_FLAG} ${INCLUDES} 
CCDEFS       = ${COMPILE_DEF} 
CPPDEFS      = ${COMPILE_DEF} 
LOADOPTS     = ${PROFILE_FLAG} ${LIBS}

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
