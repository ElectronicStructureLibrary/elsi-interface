### Core Fortran tests ###
list(APPEND ftest_src
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

add_executable(elsi_test ${ftest_src})
add_dependencies(elsi_test elsi)
target_link_libraries(elsi_test ${ALL_LIBS})
install(TARGETS elsi_test
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

### Matrix reading/writing tests ###
list(APPEND ftestrw_src
  Fortran/elsi_test_rw.f90
  Fortran/test_rw_complex.f90
  Fortran/test_rw_real.f90
)

add_executable(elsi_test_rw ${ftestrw_src})
add_dependencies(elsi_test_rw elsi)
target_link_libraries(elsi_test_rw ${ALL_LIBS})
install(TARGETS elsi_test_rw
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

### Core C tests ###
list(APPEND ctest_src
  C/elsi_test_c.c
  C/test_dm_complex_c.c
  C/test_dm_real_c.c
  C/test_ev_complex_c.c
  C/test_ev_real_c.c
)

add_executable(elsi_test_c ${ctest_src})
add_dependencies(elsi_test_c elsi)
target_link_libraries(elsi_test_c ${ALL_LIBS})
install(TARGETS elsi_test_c
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

### Test matrices ###
set(h_real "${CMAKE_SOURCE_DIR}/test/matrices/H_real.csc")
set(s_real "${CMAKE_SOURCE_DIR}/test/matrices/S_real.csc")
set(h_cmplx "${CMAKE_SOURCE_DIR}/test/matrices/H_complex.csc")
set(s_cmplx "${CMAKE_SOURCE_DIR}/test/matrices/S_complex.csc")
