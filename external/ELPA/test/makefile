# This is an example of Makefile.
-include make.sys
#
PREPROCESSOR = -DINSTALLER -DHAVE_ENVIRONMENT_CHECKING
ifeq ($(strip $(ARCHITECTURE)),Generic)
  PREPROCESSOR += -DWITH_COMPLEX_GENERIC_KERNEL -DWITH_REAL_GENERIC_KERNEL
endif
ifeq ($(strip $(ARCHITECTURE)),Simple)
  PREPROCESSOR += -DWITH_COMPLEX_GENERIC_SIMPLE_KERNEL -DWITH_REAL_GENERIC_SIMPLE_KERNEL
endif
ifeq ($(strip $(ARCHITECTURE)),BlueGeneQ)
  PREPROCESSOR += -DWITH_COMPLEX_GENERIC_KERNEL -DWITH_REAL_GENERIC_KERNEL
endif
ifeq ($(strip $(ARCHITECTURE)),SSE)
  PREPROCESSOR += -DWITH_COMPLEX_SSE_KERNEL -DWITH_REAL_SSE_KERNEL
endif
ifeq ($(strip $(ARCHITECTURE)),AVX)
  ifeq ($(strip $(BLOCK_REAL)),2)
    PREPROCESSOR += -DWITH_REAL_AVX_BLOCK2_KERNEL
  endif
  ifeq ($(strip $(BLOCK_REAL)),4)
    PREPROCESSOR += -DWITH_REAL_AVX_BLOCK4_KERNEL
  endif
  ifeq ($(strip $(BLOCK_REAL)),6)
    PREPROCESSOR += -DWITH_REAL_AVX_BLOCK6_KERNEL
  endif
  ifeq ($(strip $(BLOCK_COMPLEX)),1)
    PREPROCESSOR += -DWITH_COMPLEX_AVX_BLOCK1_KERNEL
  endif
  ifeq ($(strip $(BLOCK_COMPLEX)),2)
    PREPROCESSOR += -DWITH_COMPLEX_AVX_BLOCK2_KERNEL
  endif
endif
ifeq ($(strip $(OPENMP)),yes)
  PREPROCESSOR += -DWITH_OPENMP
endif
#
ifeq ($(strip $(ARCHITECTURE)),)
  $(info Please set the variable ARCHITECTURE.)
  $(info How to choose appropriate kernels is documented in ../src/elpa2_kernels/README_elpa2_kernels.)
  $(info Often, generic kernels and mpi-only may be a good first step,)
  $(info but considering your specific computer platform may make a big difference.)
  $(error ARCHITECTURE not set.)
endif
ifeq ($(strip $(ARCHITECTURE)),Generic)
  ELPA2_KERNELS_O ?= elpa2_kernels_real.o elpa2_kernels_complex.o 
endif
ifeq ($(strip $(ARCHITECTURE)),Simple)
  ELPA2_KERNELS_O ?= elpa2_kernels_real_simple.o elpa2_kernels_complex_simple.o
endif
ifeq ($(strip $(ARCHITECTURE)),SSE)
  ELPA2_KERNELS_O ?= elpa2_kernels_asm_x86_64.o
endif
ifeq ($(strip $(ARCHITECTURE)),BlueGeneQ)
  ELPA2_KERNELS_O ?= elpa2_kernels_real_bgq.o elpa2_kernels_complex.o
endif
ifeq ($(strip $(ARCHITECTURE)),AVX)
  FCFLAGS += -lstdc++
  ELPA2_KERNELS_O ?= elpa2_kernels_real_sse-avx_2hv.o elpa2_kernels_complex_sse-avx_1hv.o
  ifeq ($(strip $(BLOCK_REAL)),2)
    ELPA2_KERNELS_O +=
  endif
  ifeq ($(strip $(BLOCK_REAL)),4)
    ELPA2_KERNELS_O += elpa2_kernels_real_sse-avx_4hv.o
  endif
  ifeq ($(strip $(BLOCK_REAL)),6)
    ELPA2_KERNELS_O += elpa2_kernels_real_sse-avx_4hv.o elpa2_kernels_real_sse-avx_6hv.o
  endif
  ifeq ($(strip $(BLOCK_COMPLEX)),1)
    ELPA2_KERNELS_O +=
  endif
  ifeq ($(strip $(BLOCK_COMPLEX)),2)
    ELPA2_KERNELS_O += elpa2_kernels_complex_sse-avx_2hv.o
  endif
endif
ifeq ($(strip $(ELPA2_KERNELS_O)),)
  $(info Please check your choice for ARCHITECTURE: $(strip $(ARCHITECTURE)).)
  $(info How to choose appropriate kernels is documented in ../src/elpa2_kernels/README_elpa2_kernels.)
  $(info Often, generic kernels and mpi-only may be a good first step,)
  $(info but considering your specific computer platform may make a big difference.)
  $(error Invalid choice for ARCHITECTURE.)
endif
#
ifeq ($(strip $(OPENMP)),yes)
  FCFLAGS += -openmp -lstdc++
  FTIMINGS_O ?= ftimings.o ftimings_type.o ftimings_value.o time.o highwater_mark.o virtual_memory.o resident_set_size.o
else
  ifneq ($(strip $(OPENMP)),no)
    $(info Please set OPENMP to yes/no.)
    $(error Invalid choice for OPENMP.)
  endif
endif
#
SHARED_SOURCES_O = util.o blacs_infrastructure.o check_correctnes.o prepare_matrix.o read_input_parameters.o setup_mpi.o
ELPA_QR_O = elpa_pdgeqrf.o elpa_pdlarfb.o elpa_qrkernels.o qr_utils.o
$(info Building ELPA with $(strip $(ARCHITECTURE)) architecture.)
$(info OpenMP -- $(strip $(OPENMP)).)
#
all: test_real test_complex test_real2 test_complex2

test_real: test_real.o elpa1.o elpa_utilities.o $(SHARED_SOURCES_O) $(FTIMINGS_O)
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -o $@ $^ $(SCALAPACK)

test_complex: test_complex.o elpa1.o elpa_utilities.o $(SHARED_SOURCES_O) $(FTIMINGS_O)
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -o $@ $^ $(SCALAPACK)

test_real2: test_real2.o elpa1.o elpa2.o elpa_utilities.o elpa2_utilities.o $(ELPA_QR_O) $(SHARED_SOURCES_O) $(FTIMINGS) $(ELPA2_KERNELS_O)
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -o $@ $^ $(SCALAPACK)

test_complex2: test_complex2.o elpa1.o elpa2.o elpa_utilities.o elpa2_utilities.o $(ELPA_QR_O) $(SHARED_SOURCES_O) $(FTIMINGS) $(ELPA2_KERNELS_O)
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -o $@ $^ $(SCALAPACK)

test_real.o: test_real.f90 elpa1.o
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -c $<

test_complex.o: test_complex.f90 elpa1.o
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -c $<

test_real2.o: test_real2.F90 elpa2.o $(ELPA2_KERNELS_O)
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -c $<

test_complex2.o: test_complex2.F90 elpa2.o $(ELPA2_KERNELS_O)
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -c $<

elpa1.o: ../src/elpa1.f90 elpa_utilities.o $(SHARED_SOURCES_O)
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -c $<

elpa2.o: ../src/elpa2.F90 elpa1.o elpa2_utilities.o elpa_pdgeqrf.o $(ELPA2_KERNELS_O)
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -c $<

ifeq ($(OPENMP),yes)
ftimings.o : ../src/ftimings/ftimings.F90 ftimings_type.o ftimings_value.o
	$(MPIFC) $(FCFLAGS) -c $<

ftimings_type.o : ../src/ftimings/ftimings_type.F90
	$(MPIFC) $(FCFLAGS) -c $<

ftimings_value.o : ../src/ftimings/ftimings_value.F90
	$(MPIFC) $(FCFLAGS) -c $<

time.o : ../src/ftimings/time.c
	$(MPICC) $(CCFLAGS) -c $<

highwater_mark.o : ../src/ftimings/highwater_mark.c
	$(MPICC) $(CCFLAGS) -c $<

virtual_memory.o : ../src/ftimings/virtual_memory.c
	$(MPICC) $(CCFLAGS) -c $<

resident_set_size.o : ../src/ftimings/resident_set_size.c
	$(MPICC) $(CCFLAGS) -c $<
endif

ifeq ($(ARCHITECTURE),Generic)
elpa2_kernels_real.o : ../src/elpa2_kernels/elpa2_kernels_real.f90
	$(MPIFC) $(FCFLAGS) -c $<

elpa2_kernels_complex.o : ../src/elpa2_kernels/elpa2_kernels_complex.f90
	$(MPIFC) $(FCFLAGS) -c $<
endif

ifeq ($(ARCHITECTURE),Simple)
elpa2_kernels_real_simple.o : ../src/elpa2_kernels/elpa2_kernels_real_simple.f90
	$(MPIFC) $(FCFLAGS) -c $<

elpa2_kernels_complex_simple.o : ../src/elpa2_kernels/elpa2_kernels_complex_simple.f90
	$(MPIFC) $(FCFLAGS) -c $<
endif

ifeq ($(ARCHITECTURE),AVX)
elpa2_kernels_complex_sse-avx_1hv.o : ../src/elpa2_kernels/elpa2_kernels_complex_sse-avx_1hv.cpp
	$(MPICC) $(CCFLAGS) -c $<

elpa2_kernels_real_sse-avx_2hv.o : ../src/elpa2_kernels/elpa2_kernels_real_sse-avx_2hv.c
	$(MPICC) $(CCFLAGS) -c $<

ifeq ($(BLOCK_COMPLEX),2)
elpa2_kernels_complex_sse-avx_2hv.o : ../src/elpa2_kernels/elpa2_kernels_complex_sse-avx_2hv.cpp
	$(MPICC) $(CCFLAGS) -c $<
endif

ifeq ($(BLOCK_REAL),4)
elpa2_kernels_real_sse-avx_4hv.o : ../src/elpa2_kernels/elpa2_kernels_real_sse-avx_4hv.c
	$(MPICC) $(CCFLAGS) -c $<
endif

ifeq ($(BLOCK_REAL),6)
elpa2_kernels_real_sse-avx_4hv.o : ../src/elpa2_kernels/elpa2_kernels_real_sse-avx_4hv.c
	$(MPICC) $(CCFLAGS) -c $<

elpa2_kernels_real_sse-avx_6hv.o : ../src/elpa2_kernels/elpa2_kernels_real_sse-avx_6hv.c
	$(MPICC) $(CCFLAGS) -c $<
endif
endif

ifeq ($(ARCHITECTURE),SSE)
elpa2_kernels_asm_x86_64.o : ../src/elpa2_kernels/elpa2_kernels_asm_x86_64.s
	$(MPIFC) $(FCFLAGS) -xSSE4.2 -c $<
endif

ifeq ($(ARCHITECTURE),BlueGeneQ)
elpa2_kernels_complex.o : ../src/elpa2_kernels/elpa2_kernels_complex.f90
	$(MPIFC) $(FCFLAGS) -c $<

elpa2_kernels_real_bgq.o : ../src/elpa2_kernels/elpa2_kernels_real_bgq.f90
	$(MPIFC) $(FCFLAGS) -c $<
endif

elpa_utilities.o: ../src/elpa_utilities.F90
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -c $<

elpa2_utilities.o: ../src/elpa2_utilities.F90
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -c $<

util.o: ./shared_sources/util.F90
	$(MPIFC) $(FCFLAGS) -c $<

blacs_infrastructure.o: ./shared_sources/blacs_infrastructure.F90
	$(MPIFC) $(FCFLAGS) -c $<

check_correctnes.o: ./shared_sources/check_correctnes.F90
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -c $<

read_input_parameters.o: ./shared_sources/read_input_parameters.F90
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -c $<

setup_mpi.o: ./shared_sources/setup_mpi.F90
	$(MPIFC) $(FCFLAGS) $(PREPROCESSOR) -c $<

prepare_matrix.o: ./shared_sources/prepare_matrix.F90
	$(MPIFC) $(FCFLAGS) -c $<

elpa_pdgeqrf.o: ../src/elpa_qr/elpa_pdgeqrf.f90 elpa_pdlarfb.o qr_utils.o elpa_qrkernels.o
	$(MPIFC) $(FCFLAGS) -c $<

elpa_pdlarfb.o: ../src/elpa_qr/elpa_pdlarfb.f90 qr_utils.o
	$(MPIFC) $(FCFLAGS) -c $<

qr_utils.o: ../src/elpa_qr/qr_utils.f90
	$(MPIFC) $(FCFLAGS) -c $<

elpa_qrkernels.o: ../src/elpa_qr/elpa_qrkernels.f90
	$(MPIFC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod test_real test_complex test_real2 test_complex2

install:
	ar r libelpa.a *.o
	mkdir -p ~/lib
	mv libelpa.a ~/lib/.
	mkdir -p ~/include/elpa/modules
	mv *.mod ~/include/elpa/modules/.
#	mkdir -p ~/bin
#	mv *.x ~/bin/.
