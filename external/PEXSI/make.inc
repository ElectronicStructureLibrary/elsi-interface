#!/usr/bin/bash
COMPILE_MODE     = release
USE_PROFILE      = 0

# Different compiling and linking options
#SUFFIX       = linux_v0.9.2

# Compiler and tools
CC           = $(MPICC)
CXX          = $(MPICXX)
FC           = $(MPIFC)
LOADER       = $(MPICXX)

AR           = $(ARCH)
ARFLAGS      = $(ARCHFLAGS)
RANLIB       = ranlib

CP           = cp
RM           = rm
RMFLAGS      = -f

# PEXSI directory
#PEXSI_DIR = $(THIS_DIR)/external/PEXSI

# Required libraries directories
#DSUPERLU_DIR  = /home/wy29/yu/SuperLU_DIST_4.3
#METIS_DIR     = $(HOME)/Software/metis-5.1.0/build_release
#PARMETIS_DIR  = /home/wy29/yu/parmetis-4.0.3
#PTSCOTCH_DIR  = $(HOME)/Software/scotch_6.0.0/build_release
#LAPACK_DIR    = $(HOME)/Software/lapack-3.5.0
#BLAS_DIR      = $(HOME)/Software/OpenBLAS/build_release

# Includes
PEXSI_INCLUDE    = -I${PEXSI_DIR}/include
#DSUPERLU_INCLUDE = -I${DSUPERLU_DIR}/SRC
INCLUDES         = ${PEXSI_INCLUDE} ${SUPERLU_INCLUDE}

# Libraries
CPP_LIB          = -lstdc++ -lmpi #-lmpi_cxx
#GFORTRAN_LIB     = /usr/lib/gcc/x86_64-linux-gnu/4.8/libgfortran.a
#LAPACK_LIB       = -L${MKLROOT}/lib/intel64 -lmkl_lapack95_lp64
#BLAS_LIB         = -L${MKLROOT}/lib/intel64 -lmkl_blas95_lp64 \
                  # -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
                  # -lpthread -lm
#DSUPERLU_LIB     = ${DSUPERLU_DIR}/lib/libsuperlu_dist_4.3.a
SCALAPACK_LIB    = $(SCALAPACK) $(LAPACK) $(BLAS)
PEXSI_LIB        = ${PEXSI_DIR}/lib/libpexsi.a

# Graph partitioning libraries
#METIS_LIB        = -L${PARMETIS_DIR}/build/Linux-x86_64/libmetis -lmetis
#PARMETIS_LIB     = -L${PARMETIS_DIR}/build/Linux-x86_64/libparmetis -lparmetis
#SCOTCH_LIB       = -L${PTSCOTCH_DIR}/lib -lscotchmetis -lscotch -lscotcherr
#PTSCOTCH_LIB     = -L${PTSCOTCH_DIR}/lib -lptscotchparmetis -lptscotch -lptscotcherr -lscotch

# Different compiling and linking options
ifeq (${COMPILE_MODE}, release)
  COMPILE_DEF    = -DRELEASE
  COMPILE_FLAG   = -O3 -w -g
endif
ifeq (${COMPILE_MODE}, debug)
  COMPILE_DEF    = -DDEBUG=1
  COMPILE_FLAG   = -O2 -w -g
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

#LIBS  = ${PEXSI_LIB} ${DSUPERLU_LIB} ${PAR_ND_LIB} ${SEQ_ND_LIB} ${LAPACK_LIB} ${BLAS_LIB} ${GFORTRAN_LIB}
LIBS  = ${PEXSI_LIB} ${SUPERLU_LIB} ${PARMETIS_LIB} ${METIS_LIB} ${SCALAPACK_LIB}

COMPILE_DEF  += -DAdd_

CPPFLAG = -std=c++11

CFLAGS       = ${CFLAGS_P} ${PROFILE_FLAG} ${INCLUDES} -std=c99
FFLAGS       = ${FFLAGS_P} ${PROFILE_FLAG} ${INCLUDES}
CXXFLAGS     = ${CXXFLAGS_P} ${CPPFLAG} ${PROFILE_FLAG} ${INCLUDES}
CCDEFS       = ${COMPILE_DEF}
CPPDEFS      = ${COMPILE_DEF}
LOADOPTS     = ${PROFILE_FLAG} ${LIBS} -Wl,--allow-multiple-definition
FLOADOPTS    = ${PROFILE_FLAG} ${LIBS} ${CPP_LIB} -Wl,--allow-multiple-definition \
#               -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 \
#               -lmkl_lapack95_lp64 -lmkl_blas95_lp64 \
#               -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

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
