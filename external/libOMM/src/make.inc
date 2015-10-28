SHELL     = /bin/sh
FORTRAN   = $(MPIFC)
OPTS      = $(FFLAGS_O)
MSLIBPATH = $(libOMM_DIR)/src/MatrixSwitch-0.1.2/src/
LIBNAME   = libOMM.a

ifeq ($(strip $(XLF)),yes)
  XLFPP = -WF,
endif

ifeq ($(strip $(MPI)),yes)
  FPPFLAGS += ${XLFPP}-DMPI
endif
ifeq ($(strip $(LAP)),yes)
  FPPFLAGS += ${XLFPP}-DLAP
endif
ifeq ($(strip $(SLAP)),yes)
  FPPFLAGS += ${XLFPP}-DSLAP
endif
