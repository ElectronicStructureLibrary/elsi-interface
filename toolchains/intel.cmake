### Generic Intel ###

SET(CMAKE_Fortran_COMPILER "mpiifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpiicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -xAVX -fp-model precise" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -xAVX -fp-model precise -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -xAVX -fp-model precise -std=c++11" CACHE STRING "C++ flags")

SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ELPA2_KERNEL "AVX" CACHE STRING "Use ELPA AVX kernel")
SET(ENABLE_MKL ON CACHE BOOL "Pass -mkl=cluster to linker")
# ENABLE_MLK is valid only if linking with Intel compiler and MKL libraries
# exist; use "MATH_LIB" otherwise

SET(PTSCOTCH_DIR "/home/wy29/opt/scotch_6.0.5a" CACHE PATH "PT-SCOTCH directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PT-SCOTCH libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PT-SCOTCH include directory")
