### Intel compilers + MPI + MKL ###

SET(CMAKE_Fortran_COMPILER "mpiifort")
SET(CMAKE_C_COMPILER "mpiicc")
SET(CMAKE_CXX_COMPILER "mpiicpc")

SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -xAVX -fp-model precise" CACHE STRING "Fortran flags" FORCE)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -xAVX -fp-model precise -std=c99" CACHE STRING "C flags" FORCE)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -xAVX -fp-model precise -std=c++11" CACHE STRING "C++ flags" FORCE)

SET(ENABLE_TESTS "ON" CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS "ON" CACHE BOOL "Enable C tests")
SET(ENABLE_ELPA_AVX "ON" CACHE BOOL "Use ELPA AVX kernel")
SET(ENABLE_PEXSI "ON" CACHE BOOL "Enable PEXSI")
SET(ENABLE_MKL "ON" CACHE BOOL "Use Intel MKL libraries")

SET(PTSCOTCH_DIR "/home/wy29/opt/scotch_6.0.5a" CACHE PATH "PtScotch directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PtScotch Libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PtScotch include directory")
