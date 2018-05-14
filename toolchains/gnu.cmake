### Generic GNU ###

SET(CMAKE_Fortran_COMPILER "mpif90" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -mavx" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -mavx -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -mavx -std=c++11" CACHE STRING "C++ flags")

SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ELPA2_KERNEL "AVX" CACHE STRING "Use ELPA AVX kernel")

SET(LIB_PATHS "/home/wy29/opt/scalapack-2.0.2;/home/wy29/opt/lapack-3.8.0" CACHE STRING "External library paths")
SET(LIBS "libscalapack.a;liblapack.a;librefblas.a" CACHE STRING "External libraries")
