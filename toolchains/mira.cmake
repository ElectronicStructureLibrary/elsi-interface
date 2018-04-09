### Mira ###

SET(CMAKE_Fortran_COMPILER "mpixlf2003_r" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpixlc_r" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpixlcxx_r" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -qstrict -qmaxmem=-1 -qsmp=noomp -qxlf90=autodealloc -qfree=f90 -qarch=qp -qtune=qp -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes -qessl" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -qstrict -qlanglvl=stdc99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -qstrict -qlanglvl=extended0x" CACHE STRING "C++ flags")

SET(ELPA2_KERNEL "BGQ" CACHE STRING "Use ELPA BlueGene Q kernel")
SET(ENABLE_PEXSI "ON" CACHE BOOL "Enable PEXSI")

SET(PTSCOTCH_DIR "/home/wzvyu/mira/scotch_6.0.5a" CACHE PATH "PtScotch directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PtScotch libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PtScotch include directory")
