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

LINK_MACRO = $< -I$(BUILDPATH)/include -I$(PSPINCLUDEPATH) $(BUILDPATH)/lib/$(LIBNAME) $(PSPLIB) $(SLAPLIB) $(LAPLIB) $(BLASLIB) -o $@

.F90.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) $(LINK_MACRO)

.f90.x : 
	$(FORTRAN) $(OPTS) $(LINK_MACRO)
