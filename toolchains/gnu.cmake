### Generic GNU ###

SET(CMAKE_Fortran_COMPILER "mpif90" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O3 -mavx" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O3 -mavx -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-O3 -mavx -std=c++11" CACHE STRING "C++ flags")

SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ADD_UNDERSCORE ON CACHE BOOL "Add underscore")
SET(ELPA2_KERNEL "AVX" CACHE STRING "Use ELPA AVX kernel")

SET(PTSCOTCH_DIR "/home/wy29/opt/scotch_6.0.5a_gcc" CACHE PATH "PT-SCOTCH directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PT-SCOTCH libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PT-SCOTCH include directory")
SET(MATH_LIB "/home/wy29/opt/scalapack-2.0.2/libscalapack.a;/home/wy29/opt/lapack-3.8.0/liblapack.a;/home/wy29/opt/lapack-3.8.0/librefblas.a" CACHE STRING "Linear algebra libraries")