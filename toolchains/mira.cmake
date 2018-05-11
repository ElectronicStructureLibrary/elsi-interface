### Mira ###

SET(CMAKE_Fortran_COMPILER "mpixlf2003_r" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpixlc_r" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpixlcxx_r" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -qstrict -qmaxmem=-1 -qsmp=noomp -qxlf90=autodealloc -qfree=f90 -qarch=qp -qtune=qp -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes -qessl" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -qstrict -qlanglvl=stdc99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -qstrict -qlanglvl=extended0x" CACHE STRING "C++ flags")

SET(ELPA2_KERNEL "BGQ" CACHE STRING "Use ELPA BlueGene Q kernel")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
