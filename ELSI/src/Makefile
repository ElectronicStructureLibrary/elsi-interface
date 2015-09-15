-include make.sys

all: readwrite

readwrite: readwrite.o elsi.o
	$(MPIFC) -o $@ readwrite.o elsi.o mpi_tools.o $(ELPA) $(SCALAPACK) $(HDF5)

readwrite.o: readwrite.f90 elsi.o
	$(MPIFC) -c $< elsi.o $(ELPA) $(SCALAPACK) $(HDF5)

elsi.o: elsi.f90 mpi_tools.o
	$(MPIFC) -c $< mpi_tools.o $(ELPA) $(SCALAPACK) $(HDF5)

mpi_tools.o: mpi_tools.f90
	$(MPIFC) -c $<

clean:
	rm -rf *.o *.mod readwrite
