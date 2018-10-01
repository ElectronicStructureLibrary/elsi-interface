### Generic Intel ###

SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_Fortran_FLAGS "-O3 -fp-model precise" CACHE STRING "Fortran flags")

SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")

SET(LIB_PATHS "$ENV{MKLROOT}/lib/intel64" CACHE STRING "External library paths")
SET(LIBS "mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "External libraries")
