include ../src/make.inc
include make.inc

.SUFFIXES:
.SUFFIXES: .x .f90 .F90

all : example1.x example2.x example3.x example4.x

example1.x : example1.F90

example2.x : example2.F90

example3.x : example3.F90

example4.x : example4.F90

clean : 
	rm -f *.x

.F90.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -I$(LIBPATH) $< $(LIBPATH)/$(LIBNAME) $(SLAPLIB) $(LAPLIB) $(BLASLIB) -o $@

.f90.x : 
	$(FORTRAN) $(OPTS) -I$(LIBPATH) $< $(LIBPATH)/$(LIBNAME) $(SLAPLIB) $(LAPLIB) $(BLASLIB) -o $@
