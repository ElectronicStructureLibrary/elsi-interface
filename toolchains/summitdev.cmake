### Summitdev ###

# module load essl/5.5.0-20161110
# module load netlib-lapack
# module load netlib-scalapack/2.0.2
# module load xl
# Need CMake 3.11 or newer (see http://gitlab.kitware.com/cmake/cmake/issues/17784)

SET(CMAKE_Fortran_COMPILER "mpixlf" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpixlc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpixlC" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-Ofast -qarch=pwr8 -qstrict -qxlf2003=polymorphic" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-Ofast -qarch=pwr8 -qstrict -qlanglvl=stdc99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-Ofast -qarch=pwr8 -qstrict -qlanglvl=extended0x" CACHE STRING "C++ flags")

SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")
SET(ADD_UNDERSCORE OFF CACHE BOOL "Do not suffix C functions with an underscore")

SET(LIB_PATHS "/sw/summitdev/essl/5.5.0-20161110/lib64" CACHE STRING "External library paths")
SET(LIBS "essl;scalapack;lapack;blas" CACHE STRING "External libraries")
