### Generic Intel ###

SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -xAVX -fp-model precise" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -xAVX -fp-model precise" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -xAVX -fp-model precise" CACHE STRING "C++ flags")

SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ELPA2_KERNEL "AVX" CACHE STRING "Use ELPA AVX kernel")
SET(USE_EXTERNAL_ELPA ON CACHE BOOL "Use external ELPA")

SET(INC_PATHS "$ENV{ELPAROOT}/include/elpa-2018.05.001/modules" CACHE STRING "Include paths")
SET(LIB_PATHS "$ENV{MKLROOT}/lib/intel64 $ENV{ELPAROOT}/lib" CACHE STRING "External library paths")
SET(LIBS "mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core elpa" CACHE STRING "External libraries")
