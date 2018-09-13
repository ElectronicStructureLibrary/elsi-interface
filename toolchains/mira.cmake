### Mira ###

SET(CMAKE_Fortran_COMPILER "mpixlf2003_r" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpixlc_r" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpixlcxx_r" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -qstrict -qmaxmem=-1 -qsmp=noomp -qxlf90=autodealloc -qfree=f90 -qarch=qp -qtune=qp -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes -qessl" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -qstrict -qlanglvl=stdc99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -qstrict -qlanglvl=extended0x" CACHE STRING "C++ flags")
SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--allow-multiple-definition" CACHE STRING "Linker flags")

SET(ELPA2_KERNEL "BGQ" CACHE STRING "Use ELPA BlueGene Q kernel")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")

SET(IBM_MAIN_DIR "/soft/compilers/ibmcmp-oct2017" CACHE PATH "IBM main directory")
SET(LIB_PATHS "/soft/libraries/alcf/current/xl/SCALAPACK/lib;/soft/libraries/alcf/current/xl/LAPACK/lib;/soft/libraries/essl/current/essl/5.1/lib64;${IBM_MAIN_DIR}/xlf/bg/14.1/bglib64;${IBM_MAIN_DIR}/xlsmp/bg/3.1/bglib64;${IBM_MAIN_DIR}/xlmass/bg/7.3/bglib64" CACHE STRING "External library paths")
SET(LIBS "libscalapack.a;liblapack.a;libesslsmpbg.a;libxlf90_r.a;libxlopt.a;libxlfmath.a;libxl.a;libxlsmp.a;libmass.a" CACHE STRING "External libraries")
