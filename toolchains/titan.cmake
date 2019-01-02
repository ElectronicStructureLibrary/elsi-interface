### Titan ###

# module load cmake3/3.11.3
# module load gcc/8.2.0
# module load pgi/18.10.0

SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "cc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "CC" CACHE STRING "MPI C++ compiler")

SET(CMAKE_C_FLAGS "-c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "--c++11" CACHE STRING "C++ flags")

SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")
