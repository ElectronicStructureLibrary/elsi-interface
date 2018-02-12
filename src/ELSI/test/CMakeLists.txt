list(APPEND ftest_src
  ${elsi_dir}/test/Fortran/elsi_test.f90
  ${elsi_dir}/test/Fortran/test_dm_complex.f90
  ${elsi_dir}/test/Fortran/test_dm_complex_sparse.f90
  ${elsi_dir}/test/Fortran/test_dm_real.f90
  ${elsi_dir}/test/Fortran/test_dm_real_sparse.f90
  ${elsi_dir}/test/Fortran/test_ev_complex.f90
  ${elsi_dir}/test/Fortran/test_ev_complex_sparse.f90
  ${elsi_dir}/test/Fortran/test_ev_real.f90
  ${elsi_dir}/test/Fortran/test_ev_real_sparse.f90
  )

add_executable(elsi_test ${ftest_src})
add_dependencies(elsi_test elsi)
target_link_libraries(elsi_test ${LIBS})

list(APPEND ftestrw_src
  ${elsi_dir}/test/Fortran/elsi_test_rw.f90
  ${elsi_dir}/test/Fortran/test_rw_complex.f90
  ${elsi_dir}/test/Fortran/test_rw_real.f90
  )

add_executable(elsi_test_rw ${ftestrw_src})
add_dependencies(elsi_test_rw elsi)
target_link_libraries(elsi_test_rw ${LIBS})

if(ENABLE_C_INTERFACE)
  list(APPEND ctest_src
    ${elsi_dir}/test/C/elsi_test_c.c
    ${elsi_dir}/test/C/test_dm_complex_c.c
    ${elsi_dir}/test/C/test_dm_real_c.c
    ${elsi_dir}/test/C/test_ev_complex_c.c
    ${elsi_dir}/test/C/test_ev_real_c.c
    )
  add_executable(elsi_test_c ${ctest_src})
  add_dependencies(elsi_test_c elsi)
  target_link_libraries(elsi_test_c ${LIBS})
endif()
