### Generic Intel (debug) ###

SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -fp-model precise -g -traceback -check bounds -check uninit -check pointers -fpe0" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -fp-model precise -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -fp-model precise -std=c++11" CACHE STRING "C++ flags")

SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(USE_EXTERNAL_ELPA ON CACHE BOOL "Use external ELPA")
SET(USE_EXTERNAL_OMM ON CACHE BOOL "Use external libOMM")
SET(USE_EXTERNAL_SUPERLU ON CACHE BOOL "Use external SuperLU_DIST")
SET(USE_EXTERNAL_NTPOLY ON CACHE BOOL "Use external NTPoly")

SET(LIB_PATHS "$ENV{MKLROOT}/lib/intel64 /home/wy29/opt/elpa-2018.11.001/lib /home/wy29/opt/omm/lib /home/wy29/opt/SuperLU_DIST_6.1.0/lib /home/wy29/opt/scotch_6.0.0/lib /home/wy29/opt/NTPoly/lib" CACHE STRING "External library paths")
SET(INC_PATHS "/home/wy29/opt/elpa-2018.11.001/include /home/wy29/opt/omm/include /home/wy29/opt/SuperLU_DIST_6.1.0/SRC /home/wy29/opt/scotch_6.0.0/include /home/wy29/opt/NTPoly/include" CACHE STRING "External library include paths")
SET(LIBS "OMM MatrixSwitch elpa superlu_dist ptscotchparmetis ptscotch ptscotcherr scotchmetis scotch scotcherr NTPoly mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "External libraries")
