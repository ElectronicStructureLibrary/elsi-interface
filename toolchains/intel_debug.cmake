### Generic Intel (debug) ###

SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -xAVX -fp-model precise -g -traceback -check bounds -check uninit -check pointers -fpe0" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -xAVX -fp-model precise -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -xAVX -fp-model precise -std=c++11" CACHE STRING "C++ flags")

SET(ENABLE_TESTS "ON" CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS "ON" CACHE BOOL "Enable C tests")
SET(ENABLE_PEXSI "ON" CACHE BOOL "Enable PEXSI")
SET(ENABLE_MKL "ON" CACHE BOOL "Pass -mkl=cluster to linker")
# ENABLE_MLK is valid only if linking with Intel compiler and MKL libraries
# exist; use "MATH_LIB" otherwise

SET(ELPA_LIB "/home/wy29/opt/elpa-2017.11.001/lib/libelpa.a" CACHE STRING "ELPA libraries")
SET(ELPA_INC "/home/wy29/opt/elpa-2017.11.001/include" CACHE PATH "ELPA include directory")
SET(OMM_LIB "/home/wy29/opt/omm/build/lib/libOMM.a;/home/wy29/opt/omm/build/lib/libMatrixSwitch.a" CACHE STRING "libOMM libraries")
SET(OMM_INC "/home/wy29/opt/omm/build/include" CACHE PATH "libOMM include directory")
SET(SUPERLU_LIB "/home/wy29/opt/SuperLU_DIST_5.3.0/lib/libsuperlu_dist.a" CACHE STRING "SuperLU_DIST libraries")
SET(SUPERLU_INC "/home/wy29/opt/SuperLU_DIST_5.3.0/SRC" CACHE PATH "SuperLU_DIST include directory")
SET(PTSCOTCH_DIR "/home/wy29/opt/scotch_6.0.5a" CACHE PATH "PtScotch directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PtScotch libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PtScotch include directory")
