PSPDIR = ..
include $(PSPDIR)/make.inc
include make.inc

.SUFFIXES:
.SUFFIXES: .x .f90 .F90

all : pdgemmScaling.x pzgemmScaling.x

pdgemmScaling.x : pdgemmScaling.F90

pzgemmScaling.x : pzgemmScaling.F90

clean :
	rm -f *.x

LINK_MACRO = $< -I$(BUILDPATH)/include $(BUILDPATH)/lib/$(LIBNAME) $(SLAPLIB) $(LAPLIB) $(BLASLIB) -o $@

.F90.x :
	$(FORTRAN) $(OPTS) $(FPPFLAGS) $(LINK_MACRO)

.f90.x :
	$(FORTRAN) $(OPTS) $(LINK_MACRO)
