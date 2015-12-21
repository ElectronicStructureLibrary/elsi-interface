include make.inc

.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90

all : MatrixSwitch libOMM

MatrixSwitch :
	cd $(MSLIBPATH); $(MAKE)

libOMM : $(LIBNAME)

ALLOBJ = omm_params.o omm_rand.o omm_quartic.o omm_ops.o omm.o omm_wrappers.o

$(LIBNAME) : omm_params.o $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.mod
	cd $(MSLIBPATH); rm -f *.a; rm -f *.o; rm -f *.mod

.F90.o : 
	$(MPIFC) $(FFLAGS_O) $(FPPFLAGS) -I$(MSLIBPATH) -c $< -o $@

.f90.o : 
	$(MPIFC) $(FFLAGS_O) -I$(MSLIBPATH) -c $< -o $@
