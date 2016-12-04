include ../src/make.inc
include make.inc

.SUFFIXES:
.SUFFIXES: .x .f90 .F90

all : example1.x example2.x example3.x example4.x example5.x example6.x example7.x example8.x example9.x

example1.x : example1.F90

example2.x : example2.F90

example3.x : example3.F90

example4.x : example4.F90

example5.x : example5.F90

example6.x : example6.F90

example7.x : example7.F90

example8.x : example8.F90

example9.x : example9.F90

example10.x : example10.F90

example11.x : example11.F90

example12.x : example12.F90

example13.x : example13.F90

clean : 
	rm -f *.x

LINK_MACRO = $< -I$(BUILDPATH)/include -I$(MSINCLUDEPATH) -I$(PSPINCLUDEPATH) $(BUILDPATH)/lib/$(LIBNAME) $(MSLIB) $(PSPLIB) $(SLAPLIB) $(LAPLIB) $(BLASLIB) -o $@

.F90.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) $(LINK_MACRO)

.f90.x : 
	$(FORTRAN) $(OPTS) $(LINK_MACRO)
