### Titan ###

SET(CMAKE_Fortran_COMPILER "ftn")
SET(CMAKE_C_COMPILER "cc")
SET(CMAKE_CXX_COMPILER "CC")

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99" CACHE STRING "C flags" FORCE)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --c++11" CACHE STRING "C++ flags" FORCE)

SET(ENABLE_PEXSI "ON" CACHE BOOL "Enable PEXSI")
SET(ADD_UNDERSCORE "ON" CACHE BOOL "Add underscore")

SET(PTSCOTCH_DIR "/ccs/home/vwzyu/titan/scotch_6.0.5a" CACHE PATH "PtScotch directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PtScotch Libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PtScotch include directory")
