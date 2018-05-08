### Cori-Haswell ###

SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "cc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "CC" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-fast -no-ipo" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-fast -no-ipo -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-fast -no-ipo -std=c++11" CACHE STRING "C++ flags")

SET(ELPA2_KERNEL "AVX2" CACHE STRING "Use ELPA AVX2 kernel")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")

SET(PTSCOTCH_DIR "/global/homes/v/vwzyu/cori_haswell/soft/scotch_6.0.5a" CACHE PATH "PT-SCOTCH directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PT-SCOTCH libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PT-SCOTCH include directory")
