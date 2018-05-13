### Generic Intel ###

SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -xAVX -fp-model precise" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -xAVX -fp-model precise -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -xAVX -fp-model precise -std=c++11" CACHE STRING "C++ flags")

#SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
#SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ELPA2_KERNEL "AVX" CACHE STRING "Use ELPA AVX kernel")
