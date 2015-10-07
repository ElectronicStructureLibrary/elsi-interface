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

all : MatrixSwitch libOMM

MatrixSwitch :
	cd MatrixSwitch-0.1.2/src; $(MAKE)

libOMM : libOMM.a

ALLOBJ = omm_params.o omm_rand.o omm_quartic.o omm_ops.o omm.o omm_wrappers.o

libOMM.a : omm_params.o $(ALLOBJ)
	ar cr $@ $(ALLOBJ)
	ranlib $@

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.mod
	cd MatrixSwitch-0.1.2/src; rm -f *.a; rm -f *.o; rm -f *.mod

.F90.o : 
	$(MPIFC) $(FCFLAGS_O) $(FPPFLAGS) -IMatrixSwitch-0.1.2/src/ -c $< -o $@

.f90.o : 
	$(MPIFC) $(FCFLAGS_O) -IMatrixSwitch-0.1.2/src/ -c $< -o $@
