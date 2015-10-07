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
.SUFFIXES: .o .a .f90 .F90

all : libmatrixswitch.a

ALLOBJ = MatrixSwitch.o

libmatrixswitch.a : $(ALLOBJ)
	ar cr $@ $(ALLOBJ)
	ranlib $@

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.mod

.F90.o : 
	$(MPIFC) $(FCFLAGS_O) $(FPPFLAGS) -DCONV -c $< -o $@

.f90.o : 
	$(MPIFC) $(FCFLAGS_O) -c $< -o $@
