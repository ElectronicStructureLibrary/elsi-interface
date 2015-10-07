#-include ../src/make.inc
#-include make.inc

ifeq ($(strip $(MPI)),yes)
  FPPFLAGS += -DMPI
endif
ifeq ($(strip $(LAP)),yes)
  FPPFLAGS += -DLAP
endif
ifeq ($(strip $(SLAP)),yes)
  FPPFLAGS += -DSLAP
endif

.SUFFIXES:
.SUFFIXES: .x .f90 .F90

all : example1.x example2.x example3.x example4.x example5.x example6.x example7.x

example1.x : example1.F90

example2.x : example2.F90

example3.x : example3.F90

example4.x : example4.F90

example5.x : example5.F90

example6.x : example6.F90

example7.x : example7.F90

clean : 
	rm -f *.x

LINK_MACRO = $< -I../src -I../src/MatrixSwitch-0.1.2/src -L../src -lOMM -L../src/MatrixSwitch-0.1.2/src -lmatrixswitch $(SCALAPACK) -lblas -llapack -o $@

.F90.x : 
	$(MPIFC) $(FCFLAGS_O) $(FPPFLAGS) $(LINK_MACRO)

.f90.x : 
	$(MPIFC) $(FCFLAGS_O) $(LINK_MACRO)
