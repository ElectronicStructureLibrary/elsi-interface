-include make.sys

all: readwrite

readwrite: readwrite.o elsi.o
	$(MPIFC) -o $@ readwrite.o elsi.o $(ELPALIB) $(SCALAPACK)

readwrite.o: readwrite.f90 elsi.o
	$(MPIFC) -c $< elsi.o

elsi.o: elsi.f90
	$(MPIFC) -c $< $(ELPALIB)

clean:
	rm -rf *.o *.mod readwrite
