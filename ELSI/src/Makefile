-include make.sys

all: construct_and_write read_and_write read_and_solve

read_and_solve: read_and_solve.o elsi_interface.o elsi_mpi_tools.o \
   elsi_hdf5_tools.o elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -o $@ $^ $(ELPA) $(SCALAPACK) $(HDF5)

read_and_write: read_and_write.o elsi_interface.o elsi_mpi_tools.o\
   elsi_hdf5_tools.o elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -o $@ $^ $(ELPA) $(SCALAPACK) $(HDF5)

construct_and_write: construct_and_write.o elsi_interface.o \
   elsi_mpi_tools.o elsi_hdf5_tools.o elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -o $@ $^ $(ELPA) $(SCALAPACK) $(HDF5)

read_and_solve.o: read_and_solve.f90 elsi_interface.o
	$(MPIFC) $(FFLAGS)-c $^ $(ELPA) $(SCALAPACK) $(HDF5)

read_and_write.o: read_and_write.f90 elsi_interface.o
	$(MPIFC) $(FFLAGS)-c $^ $(ELPA) $(SCALAPACK) $(HDF5)

construct_and_write.o: construct_and_write.f90 elsi_interface.o
	$(MPIFC) $(FFLAGS)-c $^ $(ELPA) $(SCALAPACK) $(HDF5)

elsi_interface.o: elsi_interface.f90 elsi_mpi_tools.o \
   elsi_hdf5_tools.o elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -c $^ $(ELPA) $(SCALAPACK) $(HDF5)

elsi_mpi_tools.o: elsi_mpi_tools.f90 elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -c $^ $(SCALAPACK)

elsi_hdf5_tools.o: elsi_hdf5_tools.f90 elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -c $^ $(HDF5)

elsi_dimensions.o: elsi_dimensions.f90
	$(MPIFC) $(FFLAGS) -c $<

clean:
	rm -rf *.o *.mod construct_and_write read_and_write read_and_solve 
