### Summit ###

# module load cmake/3.15.2
# module load xl/16.1.1-3
# module load spectrum-mpi/10.3.0.0-20190419
# module load essl/6.2.0-20190419
# module load netlib-lapack/3.8.0
# module load netlib-scalapack/2.0.2

SET(CMAKE_Fortran_COMPILER "mpixlf" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpixlc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpixlC" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-Ofast -qarch=pwr9 -qstrict" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-Ofast -qarch=pwr9 -qstrict -qlanglvl=stdc99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-Ofast -qarch=pwr9 -qstrict -qlanglvl=extended0x" CACHE STRING "C++ flags")

SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")
SET(ADD_UNDERSCORE OFF CACHE BOOL "Do not suffix C functions with an underscore")

SET(LIB_PATHS "$ENV{OLCF_ESSL_ROOT}/lib64" CACHE STRING "External library paths")
SET(LIBS "essl;scalapack;lapack;blas" CACHE STRING "External libraries")
