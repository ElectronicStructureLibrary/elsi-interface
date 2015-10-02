include make.inc

.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90

all : $(LIBNAME)

ALLOBJ = MatrixSwitch.o

$(LIBNAME) : $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.mod

.F90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $< -o $@

.f90.o : 
	$(FORTRAN) $(OPTS) -c $< -o $@
