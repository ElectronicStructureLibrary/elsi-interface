-include make.sys

all: construct_and_write read_and_write read_and_solve

read_and_solve: read_and_solve.o elsi.o mpi_tools.o hdf5_tools.o dimensions.o
	$(MPIFC) $(FFLAGS) -o $@ $^ $(ELPA) $(SCALAPACK) $(HDF5)

read_and_write: read_and_write.o elsi.o mpi_tools.o hdf5_tools.o dimensions.o
	$(MPIFC) $(FFLAGS) -o $@ $^ $(ELPA) $(SCALAPACK) $(HDF5)

construct_and_write: construct_and_write.o elsi.o mpi_tools.o hdf5_tools.o dimensions.o
	$(MPIFC) $(FFLAGS) -o $@ $^ $(ELPA) $(SCALAPACK) $(HDF5)

read_and_solve.o: read_and_solve.f90 elsi.o
	$(MPIFC) $(FFLAGS)-c $^ elsi.o $(ELPA) $(SCALAPACK) $(HDF5)

read_and_write.o: read_and_write.f90 elsi.o
	$(MPIFC) $(FFLAGS)-c $^ $(ELPA) $(SCALAPACK) $(HDF5)

construct_and_write.o: construct_and_write.f90 elsi.o
	$(MPIFC) $(FFLAGS)-c $^ $(ELPA) $(SCALAPACK) $(HDF5)

elsi.o: elsi.f90 mpi_tools.o hdf5_tools.o dimensions.o
	$(MPIFC) $(FFLAGS) -c $^ $(ELPA) $(SCALAPACK) $(HDF5)

mpi_tools.o: mpi_tools.f90 dimensions.o
	$(MPIFC) $(FFLAGS) -c $^ $(SCALAPACK)

hdf5_tools.o: hdf5_tools.f90 dimensions.o
	$(MPIFC) $(FFLAGS) -c $^ $(HDF5)

dimensions.o: dimensions.f90
	$(MPIFC) $(FFLAGS) -c $<

clean:
	rm -rf *.o *.mod construct_and_write read_and_write read_and_solve 
