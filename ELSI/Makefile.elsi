#-include make.sys

all: example/construct_and_write_elpa.x \
     example/set_and_write_elpa.x \
     example/read_and_write_elpa.x \
     example/read_and_solve_elpa.x \
     example/construct_and_write_omm.x \
     example/set_and_write_omm.x \
     example/read_and_write_omm.x \
     example/read_and_write_pexsi.x \
     example/read_and_solve_omm.x \
     example/read_and_solve_pexsi.x \
     example/set_and_convert_pexsi.x

example/%.x: example/%.o src/elsi_interface.o src/elsi_mpi_tools.o src/elsi_hdf5_tools.o src/elsi_matrix_conversion.o src/elsi_timings.o src/elsi_dimensions.o
	$(MPIFC) $(FFLAGS_I) -o $@ $^ -Isrc/ $(ELPA) $(OMM) $(PEXSI) $(SCALAPACK) $(LAPACK) $(BLAS) $(HDF5)

example/%.o: example/%.f90 src/elsi_interface.o
	$(MPIFC) $(FFLAGS_I) -o $@ -c $^ -Isrc/ $(ELPA) $(OMM) $(PEXSI) $(SCALAPACK) $(HDF5)

src/elsi_interface.o: src/elsi_interface.f90 src/elsi_mpi_tools.o src/elsi_hdf5_tools.o src/elsi_matrix_conversion.o src/elsi_timings.o src/elsi_dimensions.o
	$(MPIFC) $(FFLAGS_I) -o $@ -c $^ $(ELPA) $(OMM) $(PEXSI) $(SCALAPACK) $(HDF5)

src/elsi_mpi_tools.o: src/elsi_mpi_tools.f90 src/elsi_dimensions.o
	$(MPIFC) $(FFLAGS_I) -o $@ -c $^ $(OMM) $(PEXSI) $(SCALAPACK)

src/elsi_hdf5_tools.o: src/elsi_hdf5_tools.f90 src/elsi_dimensions.o
	$(MPIFC) $(FFLAGS_I) -o $@ -c $^ $(HDF5) $(PEXSI)

src/elsi_matrix_conversion.o: src/elsi_matrix_conversion.f90 
	$(MPIFC) $(FFLAGS_I) -o $@ -c $< $(HDF5) $(PEXSI)

src/elsi_timings.o: src/elsi_timings.f90 
	$(MPIFC) $(FFLAGS_I) -o $@ -c $< $(PEXSI)

src/elsi_dimensions.o: src/elsi_dimensions.f90 
	$(MPIFC) $(FFLAGS_I) -o $@ -c $< $(PEXSI)

clean:
	rm -rf `find . -name "*.o"`
	rm -rf `find . -name "*.mod"`
	rm -rf `find . -name "*.x"`
	rm -rf lib 
	rm -rf include 
	rm -rf bin
	rm -f src/libelsi.a

install:
	$(ARCH) $(ARCHFLAGS)  src/libelsi.a src/elsi*.o
	$(RANLIB) src/libelsi.a
	mkdir -p lib
	cp src/libelsi.a lib/libelsi.a
	mkdir -p include
	cp *.mod include/.
	mkdir -p bin
	cp example/construct_and_write_elpa.x bin/construct_and_write_elpa
	cp example/set_and_write_elpa.x bin/set_and_write_elpa
	cp example/read_and_write_elpa.x bin/read_and_write_elpa
	cp example/read_and_solve_elpa.x bin/read_and_solve_elpa
	cp example/construct_and_write_omm.x bin/construct_and_write_omm
	cp example/set_and_write_omm.x bin/set_and_write_omm
	cp example/read_and_write_omm.x bin/read_and_write_omm
	cp example/read_and_write_pexsi.x bin/read_and_write_pexsi
	cp example/read_and_solve_omm.x bin/read_and_solve_omm
	cp example/read_and_solve_pexsi.x bin/read_and_solve_pexsi
	cp example/set_and_convert_pexsi.x bin/set_and_convert_pexsi

check:
	cd regression && ./regression.sh
