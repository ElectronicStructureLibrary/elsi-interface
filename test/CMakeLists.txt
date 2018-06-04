### Fortran tests ###
LIST(APPEND ftest_src
  Fortran/elsi_test.f90
  Fortran/test_dm_cmplx_den.f90
  Fortran/test_dm_cmplx_csc1.f90
  Fortran/test_dm_cmplx_csc2.f90
  Fortran/test_dm_real_den.f90
  Fortran/test_dm_real_csc1.f90
  Fortran/test_dm_real_csc2.f90
  Fortran/test_ev_cmplx_den.f90
  Fortran/test_ev_cmplx_csc1.f90
  Fortran/test_ev_cmplx_csc2.f90
  Fortran/test_ev_real_den.f90
  Fortran/test_ev_real_csc1.f90
  Fortran/test_ev_real_csc2.f90
)

MACRO(add_test_target arg1)
  ADD_DEPENDENCIES(${arg1} elsi)
  TARGET_LINK_LIBRARIES(${arg1} PRIVATE elsi)
  TARGET_INCLUDE_DIRECTORIES(${arg1} PRIVATE ${PROJECT_SOURCE_DIR}/src)
  IF(NOT USE_EXTERNAL_OMM)
    TARGET_LINK_LIBRARIES(${arg1} PRIVATE OMM MatrixSwitch)
  ENDIF()
  IF(NOT USE_EXTERNAL_ELPA)
    TARGET_LINK_LIBRARIES(${arg1} PRIVATE elpa)
  ENDIF()
  TARGET_LINK_LIBRARIES(${arg1} PRIVATE fortjson)
  IF(ENABLE_SIPS)
    TARGET_LINK_LIBRARIES(${arg1} PRIVATE sips)
  ENDIF()
  IF(ENABLE_PEXSI)
    TARGET_LINK_LIBRARIES(${arg1} PRIVATE pexsi)
    IF(NOT USE_EXTERNAL_SUPERLU)
      TARGET_LINK_LIBRARIES(${arg1} PRIVATE
        superlu_dist
        ptscotchparmetis
        ptscotch
        ptscotcherr
        scotchmetis
        scotch
        scotcherr
      )
    ENDIF()
  ENDIF()
  TARGET_LINK_LIBRARIES(${arg1} PRIVATE ${LIBS})
  IF(NOT "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}" STREQUAL "${PROJECT_BINARY_DIR}/bin")
    INSTALL(TARGETS ${arg1}
      DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
  ENDIF()
ENDMACRO()

ADD_EXECUTABLE(elsi_test ${ftest_src})
add_test_target(elsi_test)
SET_TARGET_PROPERTIES(elsi_test
  PROPERTIES LINKER_LANGUAGE Fortran)

### Matrix IO tests ###
LIST(APPEND ftestrw_src
  Fortran/elsi_test_rw.f90
  Fortran/test_rw_cmplx.f90
  Fortran/test_rw_real.f90
)

ADD_EXECUTABLE(elsi_test_rw ${ftestrw_src})
add_test_target(elsi_test_rw)
SET_TARGET_PROPERTIES(elsi_test_rw
  PROPERTIES LINKER_LANGUAGE Fortran)

### C tests ###
IF(ENABLE_C_TESTS)
  LIST(APPEND ctest_src
    C/elsi_test_c.c
    C/test_dm_cmplx_den_c.c
    C/test_dm_real_den_c.c
    C/test_ev_cmplx_den_c.c
    C/test_ev_real_den_c.c
  )

  ADD_EXECUTABLE(elsi_test_c ${ctest_src})
  add_test_target(elsi_test_c)
  SET_TARGET_PROPERTIES(elsi_test_c
    PROPERTIES LINKER_LANGUAGE C)
ENDIF()

### Test matrices ###
SET(h_real "${PROJECT_SOURCE_DIR}/test/matrices/H_real.csc")
SET(s_real "${PROJECT_SOURCE_DIR}/test/matrices/S_real.csc")
SET(h_cmplx "${PROJECT_SOURCE_DIR}/test/matrices/H_complex.csc")
SET(s_cmplx "${PROJECT_SOURCE_DIR}/test/matrices/S_complex.csc")

### CTest ###
INCLUDE(CTest)

SEPARATE_ARGUMENTS(MPIEXEC_NP UNIX_COMMAND "${MPIEXEC_NP}")
SEPARATE_ARGUMENTS(MPIEXEC_1P UNIX_COMMAND "${MPIEXEC_1P}")

MACRO(test_fortran arg1 arg2 arg3 arg4 arg5 arg6 arg7)
  ADD_TEST(test_fortran_${arg1}
    ${MPIEXEC_NP} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test
    ${arg2} ${arg3} ${arg4} ${arg5} ${arg6} ${arg7}
  )
  SET_TESTS_PROPERTIES(test_fortran_${arg1}
    PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
  )
ENDMACRO()

### Parallel Fortran tests: ev + ELPA ###
test_fortran(01 e 0 r 1 ${h_real} ${s_real})
test_fortran(02 e 0 c 1 ${h_cmplx} ${s_cmplx})
test_fortran(03 e 1 r 1 ${h_real} ${s_real})
test_fortran(04 e 1 c 1 ${h_cmplx} ${s_cmplx})
test_fortran(05 e 2 r 1 ${h_real} ${s_real})
test_fortran(06 e 2 c 1 ${h_cmplx} ${s_cmplx})

### Parallel Fortran tests: dm + ELPA ###
test_fortran(07 d 0 r 1 ${h_real} ${s_real})
test_fortran(08 d 0 c 1 ${h_cmplx} ${s_cmplx})
test_fortran(09 d 1 r 1 ${h_real} ${s_real})
test_fortran(10 d 1 c 1 ${h_cmplx} ${s_cmplx})
test_fortran(11 d 2 r 1 ${h_real} ${s_real})
test_fortran(12 d 2 c 1 ${h_cmplx} ${s_cmplx})

### Parallel Fortran tests: dm + libOMM ###
test_fortran(13 d 0 r 2 ${h_real} ${s_real})
test_fortran(14 d 0 c 2 ${h_cmplx} ${s_cmplx})
test_fortran(15 d 1 r 2 ${h_real} ${s_real})
test_fortran(16 d 1 c 2 ${h_cmplx} ${s_cmplx})
test_fortran(17 d 2 r 2 ${h_real} ${s_real})
test_fortran(18 d 2 c 2 ${h_cmplx} ${s_cmplx})

### Parallel Fortran tests: dm + PEXSI ###
IF(ENABLE_PEXSI)
  test_fortran(19 d 0 r 3 ${h_real} ${s_real})
  test_fortran(20 d 0 c 3 ${h_cmplx} ${s_cmplx})
  test_fortran(21 d 1 r 3 ${h_real} ${s_real})
  test_fortran(22 d 1 c 3 ${h_cmplx} ${s_cmplx})
  test_fortran(23 d 2 r 3 ${h_real} ${s_real})
  test_fortran(24 d 2 c 3 ${h_cmplx} ${s_cmplx})
ENDIF()

### Parallel Fortran tests: ev + SIPS ###
IF(ENABLE_SIPS)
  test_fortran(25 e 0 r 5 ${h_real} ${s_real})
  test_fortran(26 e 1 r 5 ${h_real} ${s_real})
  test_fortran(27 e 2 r 5 ${h_real} ${s_real})
ENDIF()

### Parallel Fortran tests: dm + SIPS ###
IF(ENABLE_SIPS)
  test_fortran(28 d 0 r 5 ${h_real} ${s_real})
  test_fortran(29 d 1 r 5 ${h_real} ${s_real})
  test_fortran(30 d 2 r 5 ${h_real} ${s_real})
ENDIF()

### Serial Fortran tests: ev + ELPA (LAPACK) ###
ADD_TEST(test_serial_01
  ${MPIEXEC_1P} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test
  e 0 r 1 ${h_real} ${s_real}
)
SET_TESTS_PROPERTIES(test_serial_01
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

ADD_TEST(test_serial_02
  ${MPIEXEC_1P} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test
  e 0 c 1 ${h_cmplx} ${s_cmplx}
)
SET_TESTS_PROPERTIES(test_serial_02
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

### Fortran tests: matrix IO ###
ADD_TEST(test_matio_01
  ${MPIEXEC_1P} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test_rw
  r ${h_real} ${s_real}
)
SET_TESTS_PROPERTIES(test_matio_01
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

ADD_TEST(test_matio_02
  ${MPIEXEC_1P} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test_rw
  c ${h_cmplx} ${s_cmplx}
)
SET_TESTS_PROPERTIES(test_matio_02
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

ADD_TEST(test_matio_03
  ${MPIEXEC_NP} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test_rw
  r ${h_real} ${s_real}
)
SET_TESTS_PROPERTIES(test_matio_03
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

ADD_TEST(test_matio_04
  ${MPIEXEC_NP} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test_rw
  c ${h_cmplx} ${s_cmplx}
)
SET_TESTS_PROPERTIES(test_matio_04
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
)

IF(ENABLE_C_TESTS)
  MACRO(test_c arg1 arg2 arg3 arg4 arg5 arg6 arg7)
    ADD_TEST(test_c_${arg1}
      ${MPIEXEC_NP} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test_c
      ${arg2} ${arg3} ${arg4} ${arg5} ${arg6} ${arg7}
    )
    SET_TESTS_PROPERTIES(test_c_${arg1}
      PROPERTIES PASS_REGULAR_EXPRESSION "Passed"
    )
  ENDMACRO()

  ### Parallel C tests: ev + ELPA ###
  test_c(01 e 0 r 1 ${h_real} ${s_real})
  test_c(02 e 0 c 1 ${h_cmplx} ${s_cmplx})

  ### Parallel C tests: dm + ELPA ###
  test_c(03 d 0 r 1 ${h_real} ${s_real})
  test_c(04 d 0 c 1 ${h_cmplx} ${s_cmplx})

  ### Parallel C tests: dm + libOMM ###
  test_c(05 d 0 r 2 ${h_real} ${s_real})
  test_c(06 d 0 c 2 ${h_cmplx} ${s_cmplx})

  ### Parallel C tests: dm + PEXSI ###
  IF(ENABLE_PEXSI)
    test_c(07 d 0 r 3 ${h_real} ${s_real})
    test_c(08 d 0 c 3 ${h_cmplx} ${s_cmplx})
  ENDIF()
ENDIF()
