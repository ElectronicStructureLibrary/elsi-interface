-include make.sys

all: readwrite

readwrite: readwrite.o elsi.o hdf5_tools.o
	$(MPIFC) $(FFLAGS) -o $@ readwrite.o elsi.o mpi_tools.o hdf5_tools.o $(ELPA) $(SCALAPACK) $(HDF5)

readwrite.o: readwrite.f90 elsi.o
	$(MPIFC) $(FFLAGS)-c $< elsi.o $(ELPA) $(SCALAPACK) $(HDF5)

elsi.o: elsi.f90 mpi_tools.o hdf5_tools.o
	$(MPIFC) $(FFLAGS) -c $< mpi_tools.o hdf5_tools.o $(ELPA) $(SCALAPACK) $(HDF5)

mpi_tools.o: mpi_tools.f90
	$(MPIFC) $(FFLAGS) -c $< $(SCALAPACK)

hdf5_tools.o: hdf5_tools.f90
	$(MPIFC) $(FFLAGS) -c $< $(HDF5)

clean:
	rm -rf *.o *.mod readwrite
