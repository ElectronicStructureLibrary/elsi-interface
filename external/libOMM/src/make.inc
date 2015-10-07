SHELL     = /bin/sh
FORTRAN   = $(MPIFC)
OPTS      = $(FFLAGS_O)
MSLIBPATH = $(libOMM_DIR)/src/MatrixSwitch-0.1.2/src/
LIBNAME   = libOMM.a

ifeq ($(strip $(MPI)),yes)
  FPPFLAGS += -DMPI
endif
ifeq ($(strip $(LAP)),yes)
  FPPFLAGS += -DLAP
endif
ifeq ($(strip $(SLAP)),yes)
  FPPFLAGS += -DSLAP
endif
