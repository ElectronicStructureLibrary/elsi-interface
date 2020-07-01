### Draco ###

SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -ip -xcore-avx2 -fp-model precise" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -ip -xcore-avx2 -fp-model precise -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -ip -xcore-avx2 -fp-model precise -std=c++11" CACHE STRING "C++ flags")

SET(USE_EXTERNAL_ELPA ON CACHE BOOL "Use external ELPA")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")

SET(INC_PATHS "$ENV{ELPA_HOME}/include/elpa-2020.05.001/modules" CACHE STRING "External library include paths")
SET(LIB_PATHS "$ENV{ELPA_HOME}/lib $ENV{MKLROOT}/lib/intel64" CACHE STRING "External library paths")
SET(LIBS "elpa mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "External libraries")
