-include make.sys

all: write read_and_write read_and_solve

read_and_solve: read_and_solve.o elsi.o hdf5_tools.o
	$(MPIFC) $(FFLAGS) -o $@ read_and_solve.o elsi.o mpi_tools.o hdf5_tools.o $(ELPA) $(SCALAPACK) $(HDF5)

read_and_write: read_and_write.o elsi.o hdf5_tools.o
	$(MPIFC) $(FFLAGS) -o $@ read_and_write.o elsi.o mpi_tools.o hdf5_tools.o $(ELPA) $(SCALAPACK) $(HDF5)

write: write.o elsi.o hdf5_tools.o
	$(MPIFC) $(FFLAGS) -o $@ write.o elsi.o mpi_tools.o hdf5_tools.o $(ELPA) $(SCALAPACK) $(HDF5)

read_and_solve.o: read_and_solve.f90 elsi.o
	$(MPIFC) $(FFLAGS)-c $< elsi.o $(ELPA) $(SCALAPACK) $(HDF5)

read_and_write.o: read_and_write.f90 elsi.o
	$(MPIFC) $(FFLAGS)-c $< elsi.o $(ELPA) $(SCALAPACK) $(HDF5)

write.o: write.f90 elsi.o
	$(MPIFC) $(FFLAGS)-c $< elsi.o $(ELPA) $(SCALAPACK) $(HDF5)

elsi.o: elsi.f90 mpi_tools.o hdf5_tools.o
	$(MPIFC) $(FFLAGS) -c $< mpi_tools.o hdf5_tools.o $(ELPA) $(SCALAPACK) $(HDF5)

mpi_tools.o: mpi_tools.f90
	$(MPIFC) $(FFLAGS) -c $< $(SCALAPACK)

hdf5_tools.o: hdf5_tools.f90
	$(MPIFC) $(FFLAGS) -c $< $(HDF5)

clean:
	rm -rf *.o *.mod 
