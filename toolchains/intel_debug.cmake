### Intel compilers + MPI + MKL (debug) ###

SET(CMAKE_Fortran_COMPILER "ifort")
SET(CMAKE_C_COMPILER "icc")
SET(CMAKE_CXX_COMPILER "icpc")

SET(MPI_Fortran_COMPILER "mpiifort")
SET(MPI_C_COMPILER "mpiicc")
SET(MPI_CXX_COMPILER "mpiicpc")

SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -xAVX -fp-model precise -g -traceback -check bounds -check uninit -check pointers -fpe0" CACHE STRING "Fortran flags" FORCE)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -xAVX -fp-model precise -g -Wall -Wno-unused-parameter" CACHE STRING "C flags" FORCE)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -xAVX -fp-model precise -g -Wall -Wno-unused-parameter" CACHE STRING "C++ flags" FORCE)

SET(ENABLE_TESTS "ON" CACHE BOOL "Enable Fortran tests")
SET(ENABLE_PEXSI "ON" CACHE BOOL "Enable PEXSI")
SET(ENABLE_MKL "ON" CACHE BOOL "Use Intel MKL libraries")

SET(ELPA_LIB "/home/wy29/opt/elpa-2017.11.001.rc1/lib/libelpa.a" CACHE STRING "ELPA libraries")
SET(ELPA_INC "/home/wy29/opt/elpa-2017.11.001.rc1/include" CACHE PATH "ELPA include directory")
SET(OMM_LIB "/home/wy29/opt/omm/build/lib/libOMM.a;/home/wy29/opt/omm/build/lib/libMatrixSwitch.a" CACHE STRING "libOMM libraries")
SET(OMM_INC "/home/wy29/opt/omm/build/include" CACHE PATH "libOMM include directory")
SET(SUPERLU_LIB "/home/wy29/opt/SuperLU_DIST_5.1.3/lib/libsuperlu_dist.a" CACHE STRING "SuperLU_DIST libraries")
SET(SUPERLU_INC "/home/wy29/opt/SuperLU_DIST_5.1.3/SRC" CACHE PATH "SuperLU_DIST include directory")
SET(PTSCOTCH_DIR "/home/wy29/opt/scotch_6.0.4" CACHE PATH "PtScotch directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PtScotch Libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PtScotch include directory")
