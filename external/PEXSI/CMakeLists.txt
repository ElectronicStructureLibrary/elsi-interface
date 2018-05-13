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
  src/utility.cpp
)

ADD_LIBRARY(pexsi ${pexsi_src})

TARGET_COMPILE_DEFINITIONS(pexsi PRIVATE RELEASE)
IF(ADD_UNDERSCORE)
  TARGET_COMPILE_DEFINITIONS(pexsi PRIVATE Add_)
ENDIF()

IF(NOT USE_EXTERNAL_SUPERLU)
  ADD_DEPENDENCIES(pexsi superlu_dist)
ENDIF()

IF(NOT USE_EXTERNAL_SUPERLU)
  TARGET_LINK_LIBRARIES(pexsi PRIVATE
    superlu_dist
    ptscotchparmetis
    ptscotch
    ptscotcherr
    scotchmetis
    scotch
    scotcherr
  )
ELSE()
  TARGET_LINK_LIBRARIES(pexsi PRIVATE ${SUPERLU_LIB})
ENDIF()
TARGET_LINK_LIBRARIES(pexsi PRIVATE ${MATH_LIB})

TARGET_INCLUDE_DIRECTORIES(pexsi PRIVATE include)
TARGET_INCLUDE_DIRECTORIES(pexsi PRIVATE ${INC_PATH})
IF(NOT USE_EXTERNAL_SUPERLU)
  TARGET_INCLUDE_DIRECTORIES(pexsi PRIVATE ${PROJECT_SOURCE_DIR}/external/SCOTCH/src/libscotch)
  TARGET_INCLUDE_DIRECTORIES(pexsi PRIVATE ${PROJECT_SOURCE_DIR}/external/SCOTCH/src/libscotchmetis)
  TARGET_INCLUDE_DIRECTORIES(pexsi PRIVATE ${PROJECT_BINARY_DIR}/generated)
  TARGET_INCLUDE_DIRECTORIES(pexsi PRIVATE ${PROJECT_SOURCE_DIR}/external/SuperLU_DIST/src)
ENDIF()

IF(NOT "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" STREQUAL "${PROJECT_BINARY_DIR}/lib")
  INSTALL(TARGETS pexsi
    DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
ENDIF()
