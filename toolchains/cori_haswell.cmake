### Cori-Haswell ###

SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "cc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "CC" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -ip -xcore-avx2 -fp-model precise" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -ip -xcore-avx2 -fp-model precise -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -ip -xcore-avx2 -fp-model precise -std=c++11" CACHE STRING "C++ flags")

SET(ELPA2_KERNEL "AVX2" CACHE STRING "Use ELPA AVX2 kernel")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")
