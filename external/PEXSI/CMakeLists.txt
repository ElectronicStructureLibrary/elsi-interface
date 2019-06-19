### Source files ###
LIST(APPEND pexsi_src
  src/TreeBcast.cpp
  src/blas.cpp
  src/f_interface.f90
  src/getPole.cpp
  src/get_perm_c_parmetis.c
  src/global.cpp
  src/interface.cpp
  src/lapack.cpp
  src/mpi_interf.cpp
  src/pdsymbfact.c
  src/pole.cpp
  src/ppexsi.cpp
  src/pzsymbfact.c
  src/superlu_dist_internal_complex.cpp
  src/superlu_dist_internal_real.cpp
  src/timer.cpp
  src/utility.cpp)

ADD_LIBRARY(pexsi ${pexsi_src})

TARGET_COMPILE_DEFINITIONS(pexsi PRIVATE RELEASE)
IF(ADD_UNDERSCORE)
  TARGET_COMPILE_DEFINITIONS(pexsi PRIVATE Add_)
ENDIF()

ADD_DEPENDENCIES(pexsi superlu_dist)

TARGET_LINK_LIBRARIES(pexsi PRIVATE
  ${LIBS}
  superlu_dist
  ptscotchparmetis
  ptscotch
  ptscotcherr
  scotchmetis
  scotch
  scotcherr)

TARGET_INCLUDE_DIRECTORIES(pexsi PRIVATE
  include
  ${INC_PATHS}
  ${PROJECT_SOURCE_DIR}/external/SuperLU_DIST/src
  ${PROJECT_SOURCE_DIR}/external/SCOTCH/src/libscotch
  ${PROJECT_SOURCE_DIR}/external/SCOTCH/src/libscotchmetis
  ${PROJECT_BINARY_DIR}/generated)

IF(NOT CMAKE_VERSION VERSION_LESS 3.8.2)
  TARGET_COMPILE_FEATURES(pexsi PRIVATE c_std_99 cxx_std_11)
ENDIF()

INSTALL(TARGETS pexsi
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
