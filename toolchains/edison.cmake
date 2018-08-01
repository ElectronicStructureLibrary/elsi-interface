### Edison ###

SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "cc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "CC" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-fast -no-ipo" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-fast -no-ipo" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-fast -no-ipo" CACHE STRING "C++ flags")

SET(ELPA2_KERNEL "AVX" CACHE STRING "Use ELPA AVX kernel")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")
