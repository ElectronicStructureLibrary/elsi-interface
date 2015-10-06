-include make.sys

all: example/construct_and_write_elpa.x \
      example/set_and_write_elpa.x \
      example/read_and_write_elpa.x \
      example/read_and_solve_elpa.x

example/%.x: example/%.o src/elsi_interface.o \
   src/elsi_mpi_tools.o src/elsi_hdf5_tools.o src/elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -o $@ $^ -Isrc/ $(ELPA) $(OMM) $(SCALAPACK) $(HDF5)

example/%.o: example/%.f90 src/elsi_interface.o
	$(MPIFC) $(FFLAGS) -o $@ -c $^ -Isrc/ $(ELPA) $(OMM) $(SCALAPACK) $(HDF5)

src/elsi_interface.o: src/elsi_interface.f90 src/elsi_mpi_tools.o \
   src/elsi_hdf5_tools.o src/elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -module src/ -o $@ -c $^ $(ELPA) $(OMM) $(SCALAPACK) \
      $(HDF5)

src/elsi_mpi_tools.o: src/elsi_mpi_tools.f90 src/elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -module src/ -o $@ -c $^ $(OMM) $(SCALAPACK)

src/elsi_hdf5_tools.o: src/elsi_hdf5_tools.f90 src/elsi_dimensions.o
	$(MPIFC) $(FFLAGS) -module src/ -o $@ -c $^ $(HDF5)

src/elsi_dimensions.o: src/elsi_dimensions.f90
	$(MPIFC) $(FFLAGS) -module src/ -o $@ -c $<

clean:
	rm -rf `find . -name "*.o"`
	rm -rf `find . -name "*.mod"`
	rm -rf `find . -name "*.x"`
	rm -rf lib 
	rm -rf include 
	rm -rf bin

install:
	ar r src/libelsi.a src/elsi*.o
	mkdir -p lib
	mv src/libelsi.a lib/libelsi.a
	mkdir -p include
	mv src/*.mod include/.
	mkdir -p bin
	mv example/construct_and_write_elpa.x bin/construct_and_write_elpa
	mv example/set_and_write_elpa.x bin/set_and_write_elpa
	mv example/read_and_write_elpa.x bin/read_and_write_elpa
	mv example/read_and_solve_elpa.x bin/read_and_solve_elpa

check:
	cd regression && ./regression.sh
