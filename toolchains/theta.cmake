### Theta ###

SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "cc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "CC" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-fast -no-ipo" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-fast -no-ipo -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-fast -no-ipo -std=c++11" CACHE STRING "C++ flags")

SET(ELPA2_KERNEL "AVX512" CACHE STRING "Use ELPA AVX512 kernel")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(SCOTCH_LAST_RESORT "aprun -n 1 -N 1" CACHE STRING "Command to run PT-SCOTCH header generation")
