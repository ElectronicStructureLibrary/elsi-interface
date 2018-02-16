### Link line ###
SET(LINK_LIB "${OMM_LIB} ${ELPA_LIB} ${PEXSI_LIB} ${SUPERLU_LIB} ${PTSCOTCH_LIB} ${MATH_LIB}")

### Core Fortran tests ###
LIST(APPEND ftest_src
  Fortran/elsi_test.f90
  Fortran/test_dm_complex.f90
  Fortran/test_dm_complex_sparse.f90
  Fortran/test_dm_real.f90
  Fortran/test_dm_real_sparse.f90
  Fortran/test_ev_complex.f90
  Fortran/test_ev_complex_sparse.f90
  Fortran/test_ev_real.f90
  Fortran/test_ev_real_sparse.f90
)

ADD_EXECUTABLE(elsi_test ${ftest_src})
ADD_DEPENDENCIES(elsi_test elsi)
TARGET_LINK_LIBRARIES(elsi_test elsi ${LINK_LIB})
INSTALL(TARGETS elsi_test
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

### Matrix reading/writing tests ###
LIST(APPEND ftestrw_src
  Fortran/elsi_test_rw.f90
  Fortran/test_rw_complex.f90
  Fortran/test_rw_real.f90
)

ADD_EXECUTABLE(elsi_test_rw ${ftestrw_src})
ADD_DEPENDENCIES(elsi_test_rw elsi)
TARGET_LINK_LIBRARIES(elsi_test_rw elsi ${LINK_LIB})
INSTALL(TARGETS elsi_test_rw
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

### Core C tests ###
LIST(APPEND ctest_src
  C/elsi_test_c.c
  C/test_dm_complex_c.c
  C/test_dm_real_c.c
  C/test_ev_complex_c.c
  C/test_ev_real_c.c
)

ADD_EXECUTABLE(elsi_test_c ${ctest_src})
ADD_DEPENDENCIES(elsi_test_c elsi)
TARGET_LINK_LIBRARIES(elsi_test_c elsi ${LINK_LIB})
INSTALL(TARGETS elsi_test_c
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

### Test matrices ###
SET(h_real "${CMAKE_SOURCE_DIR}/test/matrices/H_real.csc")
SET(s_real "${CMAKE_SOURCE_DIR}/test/matrices/S_real.csc")
SET(h_cmplx "${CMAKE_SOURCE_DIR}/test/matrices/H_complex.csc")
SET(s_cmplx "${CMAKE_SOURCE_DIR}/test/matrices/S_complex.csc")
