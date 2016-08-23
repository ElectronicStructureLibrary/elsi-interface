SHELL = /bin/sh
OPTS  = $(FFLAGS_O)

ifeq ($(strip $(XLF)),yes)
  XLFPP = -WF,
endif

FPPFLAGS  = ${XLFPP}-DCONV
FPPFLAGS += ${XLFPP}-DPSP

ifeq ($(strip $(MPI)),yes)
  FPPFLAGS += ${XLFPP}-DMPI
endif
ifeq ($(strip $(LAP)),yes)
  FPPFLAGS += ${XLFPP}-DLAP
endif
ifeq ($(strip $(SLAP)),yes)
  FPPFLAGS += ${XLFPP}-DSLAP
endif

PSPLIBNAME = libpspblas.a
PSPBUILDPATH = $(ROOTDIR)/build_psp

MSLIBNAME = libMatrixSwitch.a
MSBUILDPATH = $(ROOTDIR)/build_MatrixSwitch

OMMLIBNAME = libOMM.a
OMMBUILDPATH = $(ROOTDIR)/build_libOMM

BLASLIB = $(SCALAPACK_LDFLAGS)
LAPLIB = $(SCALAPACK_LDFLAGS)
SLAPLIB = $(SCALAPACK_LDFLAGS)
