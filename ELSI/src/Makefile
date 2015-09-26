-include make.sys

all: construct_and_write.x read_and_write.x read_and_solve.x

read_and_solve.x: read_and_solve.o elsi_interface.o elsi_mpi_tools.o \
   elsi_hdf5_tools.o elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -o $@ $^ $(ELPA) $(SCALAPACK) $(HDF5)

read_and_write.x: read_and_write.o elsi_interface.o elsi_mpi_tools.o\
   elsi_hdf5_tools.o elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -o $@ $^ $(ELPA) $(SCALAPACK) $(HDF5)

construct_and_write.x: construct_and_write.o elsi_interface.o \
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
	rm -rf *.o *.mod *.x ../lib ../include ../bin

install:
	ar r libelsi.a elsi*.o
	mkdir -p ../lib
	mv libelsi.a ../lib/.
	mkdir -p ../include
	mv *.mod ../include/.
	mkdir -p ../bin
	mv construct_and_write.x ../bin/construct_and_write
	mv read_and_write.x ../bin/read_and_write
	mv read_and_solve.x ../bin/read_and_solve
