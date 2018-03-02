### Intel compilers + OpenMPI + Intel MKL + Netlib ScaLAPACK ###

SET(CMAKE_Fortran_COMPILER "mpif90")
SET(CMAKE_C_COMPILER "mpicc")
SET(CMAKE_CXX_COMPILER "mpicxx")

SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fast -no-ipo" CACHE STRING "Fortran flags" FORCE)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fast -no-ipo" CACHE STRING "C flags" FORCE)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -no-ipo" CACHE STRING "C++ flags" FORCE)

SET(ENABLE_TESTS "ON" CACHE BOOL "Enable Fortran tests")
SET(ENABLE_ELPA_AVX "ON" CACHE BOOL "Use ELPA AVX kernel")
SET(ENABLE_PEXSI "ON" CACHE BOOL "Enable PEXSI")

SET(PTSCOTCH_DIR "/Users/vyu/Soft/scotch_6.0.5a" CACHE PATH "PtScotch directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PtScotch Libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PtScotch include directory")
SET(MATH_LIB "/Users/vyu/Soft/scalapack-2.0.2/libscalapack.a;/opt/intel/mkl/lib/libmkl_intel_lp64.a;/opt/intel/mkl/lib/libmkl_sequential.a;/opt/intel/mkl/lib/libmkl_core.a" CACHE STRING "Linear algebra libraries")
