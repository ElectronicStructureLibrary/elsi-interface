INCLUDE_DIRECTORIES(include)

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

IF(NOT SUPERLU_FOUND)
  ADD_DEPENDENCIES(pexsi superlu_dist)
ENDIF()

IF(BUILD_SHARED_LIBS)
  TARGET_LINK_LIBRARIES(pexsi ${PTSCOTCH_LIB} ${MATH_LIB})
  IF(ENABLE_MKL)
    SET_TARGET_PROPERTIES(pexsi
      PROPERTIES LINK_FLAGS "-mkl=cluster")
  ENDIF()
  IF(NOT SUPERLU_FOUND)
    TARGET_LINK_LIBRARIES(pexsi superlu_dist)
  ELSE()
    TARGET_LINK_LIBRARIES(pexsi ${SUPERLU_LIB})
  ENDIF()
ENDIF()

INSTALL(TARGETS pexsi
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
