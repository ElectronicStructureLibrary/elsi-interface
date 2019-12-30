### Generic Intel (debug) ###

SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -ip -fp-model precise -g -traceback -check bounds -check uninit -check pointers -fpe0 -warn all -init=snan,arrays" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -ip -fp-model precise -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -ip -fp-model precise -std=c++11" CACHE STRING "C++ flags")
SET(CMAKE_EXE_LINKER_FLAGS "-lstdc++" CACHE STRING "Linker flags")

SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_EIGENEXA ON CACHE BOOL "Enable EigenExa")
SET(USE_EXTERNAL_ELPA ON CACHE BOOL "Use external ELPA")
SET(USE_EXTERNAL_OMM ON CACHE BOOL "Use external libOMM")
SET(USE_EXTERNAL_PEXSI ON CACHE BOOL "Use external PEXSI")
SET(USE_EXTERNAL_NTPOLY ON CACHE BOOL "Use external NTPoly")

SET(LIB_PATHS "$ENV{MKLROOT}/lib/intel64 /home/wy29/opt/elpa/lib /home/wy29/opt/omm/lib /home/wy29/opt/pexsi/src /home/wy29/opt/SuperLU_DIST_6.1.1/lib /home/wy29/opt/scotch_6.0.0/lib /home/wy29/opt/NTPoly/lib /home/wy29/opt/EigenExa-2.4b/lib" CACHE STRING "External library paths")
SET(INC_PATHS "/home/wy29/opt/elpa/include /home/wy29/opt/omm/include /home/wy29/opt/pexsi/src /home/wy29/opt/SuperLU_DIST_6.1.1/SRC /home/wy29/opt/scotch_6.0.0/include /home/wy29/opt/NTPoly/include /home/wy29/opt/EigenExa-2.4b/include" CACHE STRING "External library include paths")
SET(LIBS "OMM MatrixSwitch elpa pexsi superlu_dist ptscotchparmetis ptscotch ptscotcherr scotchmetis scotch scotcherr NTPoly EigenExa mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "External libraries")
