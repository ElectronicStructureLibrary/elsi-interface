include ../../make.inc

.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90

all : libmatrixswitch.a

ALLOBJ = MatrixSwitch.o

libmatrixswitch.a : $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.mod

.F90.o : 
	$(MPIFC) $(FCFLAGS_O) $(FPPFLAGS) -DCONV -c $< -o $@

.f90.o : 
	$(MPIFC) $(FCFLAGS_O) -c $< -o $@
