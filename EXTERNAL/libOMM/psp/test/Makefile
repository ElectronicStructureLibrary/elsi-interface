PSPDIR = ..
include $(PSPDIR)/make.inc
include make.inc

.SUFFIXES:
.SUFFIXES: .x .f90 .F90

all : pdgemmScaling.x pzgemmScaling.x testList.x matSum.x testMPI.x testMPI2.x testMPI3.x pdgemm.x spBLAS.x 

pdgemmScaling.x : pdgemmScaling.F90

pzgemmScaling.x : pzgemmScaling.F90

testList.x : testList.F90

testMPI.x : testMPI.F90

testMPI2.x : testMPI2.F90

testMPI3.x : testMPI3.F90

matSum.x : matSum.F90

pdgemm.x : pdgemm.F90

spBLAS.x : spBLAS.F90

clean :
	rm -f *.x

LINK_MACRO = $< -I$(BUILDPATH)/include $(BUILDPATH)/lib/$(LIBNAME) $(SLAPLIB) $(LAPLIB) $(BLASLIB) -o $@

.F90.x :
	$(FORTRAN) $(OPTS) $(FPPFLAGS) $(LINK_MACRO)

.f90.x :
	$(FORTRAN) $(OPTS) $(LINK_MACRO)
