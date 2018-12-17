### Cobra ###

SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -xcore-avx512" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -xcore-avx512 -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -xcore-avx512 -std=c++11" CACHE STRING "C++ flags")

SET(USE_EXTERNAL_ELPA ON CACHE BOOL "Use external ELPA")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")

SET(INC_PATHS "/mpcdf/soft/SLE_12_SP3/packages/skylake/elpa/intel_18_0_5-impi_2018_4/2018.05.001-standard/include/elpa-2018.05.001/modules" CACHE STRING "External library include paths")
SET(LIB_PATHS "/mpcdf/soft/SLE_12_SP3/packages/skylake/elpa/intel_18_0_5-impi_2018_4/2018.05.001-standard/lib $ENV{MKLROOT}/lib/intel64" CACHE STRING "External library paths")
SET(LIBS "elpa mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "External libraries")
