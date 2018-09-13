### Cobra ###

SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -xcore-avx512" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -xcore-avx512 -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -xcore-avx512 -std=c++11" CACHE STRING "C++ flags")

SET(ELPA2_KERNEL "AVX512" CACHE STRING "Use ELPA AVX512 kernel")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")

SET(LIB_PATHS "$ENV{MKLROOT}/lib/intel64" CACHE STRING "External library paths")
SET(LIBS "mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "External libraries")
