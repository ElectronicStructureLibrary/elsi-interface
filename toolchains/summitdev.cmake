### Summitdev ###
### Need CMake 3.11.0 (see http://gitlab.kitware.com/cmake/cmake/issues/17784)

SET(CMAKE_Fortran_COMPILER "mpixlf" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpixlc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpixlC" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-Ofast -qarch=pwr8 -qstrict -qxlf2003=polymorphic" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-Ofast -qarch=pwr8 -qstrict -qlanglvl=stdc99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-Ofast -qarch=pwr8 -qstrict -qlanglvl=extended0x" CACHE STRING "C++ flags")

SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
