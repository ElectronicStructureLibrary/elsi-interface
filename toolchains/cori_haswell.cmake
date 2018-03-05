### Intel compilers + MPI + MKL (Cori-Haswell) ###

SET(CMAKE_Fortran_COMPILER "ftn")
SET(CMAKE_C_COMPILER "cc")
SET(CMAKE_CXX_COMPILER "CC")

SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fast -no-ipo" CACHE STRING "Fortran flags" FORCE)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fast -no-ipo -std=c99" CACHE STRING "C flags" FORCE)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -no-ipo -std=c++11" CACHE STRING "C++ flags" FORCE)

SET(ENABLE_ELPA_AVX2 "ON" CACHE BOOL "Use ELPA AVX2 kernel")
SET(ENABLE_PEXSI "ON" CACHE BOOL "Enable PEXSI")

SET(PTSCOTCH_DIR "/global/homes/v/vwzyu/cori_haswell/soft/scotch_6.0.5a" CACHE PATH "PtScotch directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PtScotch Libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PtScotch include directory")
