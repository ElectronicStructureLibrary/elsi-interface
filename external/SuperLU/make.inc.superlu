# The name of the libraries to be created/linked to
DSuperLUroot = $(THIS_DIR)/external/SuperLU/Super/LU
#DSUPERLULIB  = $(DSuperLUroot)/lib/libsuperlu.a
DSUPERLULIB  = $(THIS_DIR)/external/lib/libsuperlu.a
BLASDEF      = -DUSE_VENDOR_BLAS
BLASLIB      = -L${MKLROOT}/lib/intel64 \
               -I${MKLROOT}/include/intel64/lp64 -lmkl_intel_lp64 \
               -lmkl_sequential -lmkl_core -lpthread
# parmetis 4.0.2, 32-bit integer
METISLIB     = -L$(THIS_DIR)/external/lib -lmetis
PARMETISLIB  = -L$(THIS_DIR)/external/lib -lparmetis

# Define all the libraries
LIBS         = $(DSUPERLULIB) $(BLASLIB) $(PARMETISLIB) $(METISLIB) $(FLIBS)

# C compiler setup
CC           = $(MPICC)
# CFLAGS should be set to be the C flags that include optimization
CFLAGS       = -O2 -DDEBUGlevel=0 -DPRNTlevel=1
# NOOPTS should be set to be the C flags that turn off any optimization
NOOPTS       = -O0

# FORTRAN compiler setup
FORTRAN      = $(MPIFC)
F90FLAGS     = -O2

LOADER       = $(MPICC)

# C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
# Need follow the convention of how C calls a Fortran routine.
CDEFS        = -DAdd_

ARCH         = $(AR)
ARCHFLAGS    = cr
