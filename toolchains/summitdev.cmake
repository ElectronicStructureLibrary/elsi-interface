### Summitdev ###

SET(CMAKE_Fortran_COMPILER "mpixlf")
SET(CMAKE_C_COMPILER "mpixlc")
SET(CMAKE_CXX_COMPILER "mpixlC")

SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -qarch=pwr8 -qstrict -qxlf2003=polymorphic" CACHE STRING "Fortran flags" FORCE)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -qarch=pwr8 -qstrict -qlanglvl=stdc99" CACHE STRING "C flags" FORCE)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -qarch=pwr8 -qstrict -qlanglvl=extended0x" CACHE STRING "C++ flags" FORCE)

SET(ENABLE_PEXSI "ON" CACHE BOOL "Enable PEXSI")

SET(PTSCOTCH_DIR "/ccs/home/vwzyu/summitdev/scotch_6.0.5a" CACHE PATH "PtScotch directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PtScotch Libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PtScotch include directory")
