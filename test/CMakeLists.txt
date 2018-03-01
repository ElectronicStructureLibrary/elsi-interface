### Link line ###
SET(LINK_LIB
  ${OMM_LIB}
  ${ELPA_LIB}
  ${PEXSI_LIB}
  ${SIPS_LIB}
  ${SUPERLU_LIB}
  ${PTSCOTCH_LIB}
  ${MATH_LIB}
  ${MPI_Fortran_LIBRARIES}
  ${MPI_C_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
)

### Fortran tests ###
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
IF(ENABLE_MKL)
  SET_TARGET_PROPERTIES(elsi_test
    PROPERTIES LINK_FLAGS "-mkl=cluster")
ENDIF()
INSTALL(TARGETS elsi_test
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

### Matrix IO tests ###
LIST(APPEND ftestrw_src
  Fortran/elsi_test_rw.f90
  Fortran/test_rw_complex.f90
  Fortran/test_rw_real.f90
)

ADD_EXECUTABLE(elsi_test_rw ${ftestrw_src})
ADD_DEPENDENCIES(elsi_test_rw elsi)
TARGET_LINK_LIBRARIES(elsi_test_rw elsi ${LINK_LIB})
IF(ENABLE_MKL)
  SET_TARGET_PROPERTIES(elsi_test_rw
    PROPERTIES LINK_FLAGS "-mkl=cluster")
ENDIF()
INSTALL(TARGETS elsi_test_rw
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

### C tests ###
IF(ENABLE_C_TESTS)
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
  IF(ENABLE_MKL)
    SET_TARGET_PROPERTIES(elsi_test_c
      PROPERTIES LINK_FLAGS "-mkl=cluster")
  ENDIF()
  INSTALL(TARGETS elsi_test_c
    DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
ENDIF()

### Test matrices ###
SET(h_real "${CMAKE_SOURCE_DIR}/test/matrices/H_real.csc")
SET(s_real "${CMAKE_SOURCE_DIR}/test/matrices/S_real.csc")
SET(h_cmplx "${CMAKE_SOURCE_DIR}/test/matrices/H_complex.csc")
SET(s_cmplx "${CMAKE_SOURCE_DIR}/test/matrices/S_complex.csc")

### CTest ###
INCLUDE(CTest)

MACRO(test_fortran arg1 arg2 arg3 arg4 arg5 arg6 arg7)
  ADD_TEST(test_fortran_${arg1}
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
    ${EXECUTABLE_OUTPUT_PATH}/elsi_test
    ${arg2} ${arg3} ${arg4} ${arg5} ${arg6} ${arg7}
  )
  SET_TESTS_PROPERTIES(test_fortran_${arg1}
    PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
  )
ENDMACRO()

### Parallel Fortran tests: ev + ELPA ###
test_fortran(01 e d r 1 ${h_real} ${s_real})
test_fortran(02 e d c 1 ${h_cmplx} ${s_cmplx})
test_fortran(03 e s r 1 ${h_real} ${s_real})
test_fortran(04 e s c 1 ${h_cmplx} ${s_cmplx})

### Parallel Fortran tests: dm + ELPA ###
test_fortran(05 d d r 1 ${h_real} ${s_real})
test_fortran(06 d d c 1 ${h_cmplx} ${s_cmplx})
test_fortran(07 d s r 1 ${h_real} ${s_real})
test_fortran(08 d s c 1 ${h_cmplx} ${s_cmplx})

### Parallel Fortran tests: dm + libOMM ###
test_fortran(09 d d r 2 ${h_real} ${s_real})
test_fortran(10 d d c 2 ${h_cmplx} ${s_cmplx})
test_fortran(11 d s r 2 ${h_real} ${s_real})
test_fortran(12 d s c 2 ${h_cmplx} ${s_cmplx})

### Parallel Fortran tests: dm + PEXSI ###
IF(ENABLE_PEXSI)
  test_fortran(13 d d r 3 ${h_real} ${s_real})
  test_fortran(14 d d c 3 ${h_cmplx} ${s_cmplx})
  test_fortran(15 d s r 3 ${h_real} ${s_real})
  test_fortran(16 d s c 3 ${h_cmplx} ${s_cmplx})
ENDIF()

### Parallel Fortran tests: ev + SIPs ###
IF(ENABLE_SIPS)
  test_fortran(17 e d r 5 ${h_real} ${s_real})
ENDIF()

### Serial Fortran tests: ev + ELPA (LAPACK) ###
ADD_TEST(test_serial_01
  ${EXECUTABLE_OUTPUT_PATH}/elsi_test
  e d r 1 ${h_real} ${s_real}
)
SET_TESTS_PROPERTIES(test_serial_01
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

ADD_TEST(test_serial_02
  ${EXECUTABLE_OUTPUT_PATH}/elsi_test
  e d c 1 ${h_cmplx} ${s_cmplx}
)
SET_TESTS_PROPERTIES(test_serial_02
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

### Fortran tests: matrix IO ###
ADD_TEST(test_matio_01
  ${EXECUTABLE_OUTPUT_PATH}/elsi_test_rw
  r ${h_real} ${s_real}
)
SET_TESTS_PROPERTIES(test_matio_01
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

ADD_TEST(test_matio_02
  ${EXECUTABLE_OUTPUT_PATH}/elsi_test_rw
  c ${h_cmplx} ${s_cmplx}
)
SET_TESTS_PROPERTIES(test_matio_02
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

ADD_TEST(test_matio_03
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ${EXECUTABLE_OUTPUT_PATH}/elsi_test_rw
  r ${h_real} ${s_real}
)
SET_TESTS_PROPERTIES(test_matio_03
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

ADD_TEST(test_matio_04
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
  ${EXECUTABLE_OUTPUT_PATH}/elsi_test_rw
  c ${h_cmplx} ${s_cmplx}
)
SET_TESTS_PROPERTIES(test_matio_04
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

IF(ENABLE_C_TESTS)
  MACRO(test_c arg1 arg2 arg3 arg4 arg5 arg6 arg7)
    ADD_TEST(test_c_${arg1}
      ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
      ${EXECUTABLE_OUTPUT_PATH}/elsi_test_c
      ${arg2} ${arg3} ${arg4} ${arg5} ${arg6} ${arg7}
    )
    SET_TESTS_PROPERTIES(test_c_${arg1}
      PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
    )
  ENDMACRO()

  ### Parallel C tests: ev + ELPA ###
  test_c(01 e d r 1 ${h_real} ${s_real})
  test_c(02 e d c 1 ${h_cmplx} ${s_cmplx})

  ### Parallel C tests: dm + ELPA ###
  test_c(03 d d r 1 ${h_real} ${s_real})
  test_c(04 d d c 1 ${h_cmplx} ${s_cmplx})

  ### Parallel C tests: dm + libOMM ###
  test_c(05 d d r 2 ${h_real} ${s_real})
  test_c(06 d d c 2 ${h_cmplx} ${s_cmplx})

  ### Parallel C tests: dm + PEXSI ###
  IF(ENABLE_PEXSI)
    test_c(07 d d r 3 ${h_real} ${s_real})
    test_c(08 d d c 3 ${h_cmplx} ${s_cmplx})
  ENDIF()
ENDIF()
