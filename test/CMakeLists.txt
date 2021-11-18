### Fortran tests ###
LIST(APPEND ftest_src
  Fortran/elsi_test.f90
  Fortran/test_bse_cmplx_den.f90
  Fortran/test_bse_real_den.f90
  Fortran/test_dm_cmplx_den.f90
  Fortran/test_dm_cmplx_coo.f90
  Fortran/test_dm_cmplx_csc1.f90
  Fortran/test_dm_cmplx_csc2.f90
  Fortran/test_dm_real_den.f90
  Fortran/test_dm_real_coo.f90
  Fortran/test_dm_real_csc1.f90
  Fortran/test_dm_real_csc2.f90
  Fortran/test_ev_cmplx_den.f90
  Fortran/test_ev_cmplx_coo.f90
  Fortran/test_ev_cmplx_csc1.f90
  Fortran/test_ev_cmplx_csc2.f90
  Fortran/test_ev_real_den.f90
  Fortran/test_ev_real_coo.f90
  Fortran/test_ev_real_csc1.f90
  Fortran/test_ev_real_csc2.f90
  Fortran/test_occ_normal.f90
  Fortran/test_occ_non_aufbau.f90)

MACRO(add_test_target arg1)
  ADD_DEPENDENCIES(${arg1} elsi)

  TARGET_LINK_LIBRARIES(${arg1} PRIVATE elsi)

  TARGET_INCLUDE_DIRECTORIES(${arg1} PRIVATE
    ${INC_PATHS}
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    ${PROJECT_SOURCE_DIR}/src)
ENDMACRO()

ADD_EXECUTABLE(elsi_test ${ftest_src})
add_test_target(elsi_test)
SET_TARGET_PROPERTIES(elsi_test
  PROPERTIES LINKER_LANGUAGE Fortran)

### Matrix IO tests ###
LIST(APPEND ftestrw_src
  Fortran/elsi_test_rw.f90
  Fortran/test_rw_cmplx.f90
  Fortran/test_rw_real.f90)

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
    C/test_ev_real_den_c.c)

  ADD_EXECUTABLE(elsi_test_c ${ctest_src})
  add_test_target(elsi_test_c)
  SET_TARGET_PROPERTIES(elsi_test_c
    PROPERTIES LINKER_LANGUAGE C)
ENDIF()

### Test matrices ###
SET(h_real "${PROJECT_SOURCE_DIR}/test/matrices/H_real.csc")
SET(s_real "${PROJECT_SOURCE_DIR}/test/matrices/S_real.csc")
SET(h_cmplx "${PROJECT_SOURCE_DIR}/test/matrices/H_cmplx.csc")
SET(s_cmplx "${PROJECT_SOURCE_DIR}/test/matrices/S_cmplx.csc")
SET(a_real "${PROJECT_SOURCE_DIR}/test/matrices/A_real.csc")
SET(b_real "${PROJECT_SOURCE_DIR}/test/matrices/B_real.csc")
SET(a_cmplx "${PROJECT_SOURCE_DIR}/test/matrices/A_cmplx.csc")
SET(b_cmplx "${PROJECT_SOURCE_DIR}/test/matrices/B_cmplx.csc")

### CTest ###
INCLUDE(CTest)

SEPARATE_ARGUMENTS(MPIEXEC_NP UNIX_COMMAND "${MPIEXEC_NP}")
SEPARATE_ARGUMENTS(MPIEXEC_1P UNIX_COMMAND "${MPIEXEC_1P}")

MACRO(test_fortran arg1 arg2 arg3 arg4 arg5 arg6 arg7)
  ADD_TEST(test_fortran_${arg1}
    ${MPIEXEC_NP} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test
    ${arg2} ${arg3} ${arg4} ${arg5} ${arg6} ${arg7})
  SET_TESTS_PROPERTIES(test_fortran_${arg1}
    PROPERTIES PASS_REGULAR_EXPRESSION "Passed")
ENDMACRO()

### Parallel Fortran tests: ev + ELPA ###
test_fortran(01_elpa e 0 r 1 ${h_real} ${s_real})
test_fortran(02_elpa e 0 c 1 ${h_cmplx} ${s_cmplx})
test_fortran(03_elpa e 1 r 1 ${h_real} ${s_real})
test_fortran(04_elpa e 1 c 1 ${h_cmplx} ${s_cmplx})
test_fortran(05_elpa e 2 r 1 ${h_real} ${s_real})
test_fortran(06_elpa e 2 c 1 ${h_cmplx} ${s_cmplx})
test_fortran(07_elpa e 3 r 1 ${h_real} ${s_real})
test_fortran(08_elpa e 3 c 1 ${h_cmplx} ${s_cmplx})

### Parallel Fortran tests: dm + ELPA ###
test_fortran(09_elpa d 0 r 1 ${h_real} ${s_real})
test_fortran(10_elpa d 0 c 1 ${h_cmplx} ${s_cmplx})
test_fortran(11_elpa d 1 r 1 ${h_real} ${s_real})
test_fortran(12_elpa d 1 c 1 ${h_cmplx} ${s_cmplx})
test_fortran(13_elpa d 2 r 1 ${h_real} ${s_real})
test_fortran(14_elpa d 2 c 1 ${h_cmplx} ${s_cmplx})
test_fortran(15_elpa d 3 r 1 ${h_real} ${s_real})
test_fortran(16_elpa d 3 c 1 ${h_cmplx} ${s_cmplx})

### Parallel Fortran tests: dm + libOMM ###
test_fortran(01_omm d 0 r 2 ${h_real} ${s_real})
test_fortran(02_omm d 0 c 2 ${h_cmplx} ${s_cmplx})
test_fortran(03_omm d 1 r 2 ${h_real} ${s_real})
test_fortran(04_omm d 1 c 2 ${h_cmplx} ${s_cmplx})
test_fortran(05_omm d 2 r 2 ${h_real} ${s_real})
test_fortran(06_omm d 2 c 2 ${h_cmplx} ${s_cmplx})
test_fortran(07_omm d 3 r 2 ${h_real} ${s_real})
test_fortran(08_omm d 3 c 2 ${h_cmplx} ${s_cmplx})

### Parallel Fortran tests: dm + PEXSI ###
IF(ENABLE_PEXSI)
  test_fortran(01_pexsi d 0 r 3 ${h_real} ${s_real})
  test_fortran(02_pexsi d 0 c 3 ${h_cmplx} ${s_cmplx})
  test_fortran(03_pexsi d 1 r 3 ${h_real} ${s_real})
  test_fortran(04_pexsi d 1 c 3 ${h_cmplx} ${s_cmplx})
  test_fortran(05_pexsi d 2 r 3 ${h_real} ${s_real})
  test_fortran(06_pexsi d 2 c 3 ${h_cmplx} ${s_cmplx})
  test_fortran(07_pexsi d 3 r 3 ${h_real} ${s_real})
  test_fortran(08_pexsi d 3 c 3 ${h_cmplx} ${s_cmplx})
ENDIF()

### Parallel Fortran tests: ev + EigenExa ###
IF(ENABLE_EIGENEXA)
  test_fortran(01_eigenexa e 0 r 4 ${h_real} ${s_real})
  test_fortran(02_eigenexa e 1 r 4 ${h_real} ${s_real})
  test_fortran(03_eigenexa e 2 r 4 ${h_real} ${s_real})
  test_fortran(04_eigenexa e 3 r 4 ${h_real} ${s_real})
ENDIF()

### Parallel Fortran tests: dm + EigenExa ###
IF(ENABLE_EIGENEXA)
  test_fortran(05_eigenexa d 0 r 4 ${h_real} ${s_real})
  test_fortran(06_eigenexa d 1 r 4 ${h_real} ${s_real})
  test_fortran(07_eigenexa d 2 r 4 ${h_real} ${s_real})
  test_fortran(08_eigenexa d 3 r 4 ${h_real} ${s_real})
ENDIF()

### Parallel Fortran tests: ev + SLEPc-SIPs ###
IF(ENABLE_SIPS)
  test_fortran(01_sips e 0 r 5 ${h_real} ${s_real})
  test_fortran(02_sips e 1 r 5 ${h_real} ${s_real})
  test_fortran(03_sips e 2 r 5 ${h_real} ${s_real})
  test_fortran(04_sips e 3 r 5 ${h_real} ${s_real})
ENDIF()

### Parallel Fortran tests: dm + SLEPc-SIPs ###
IF(ENABLE_SIPS)
  test_fortran(05_sips d 0 r 5 ${h_real} ${s_real})
  test_fortran(06_sips d 1 r 5 ${h_real} ${s_real})
  test_fortran(07_sips d 2 r 5 ${h_real} ${s_real})
  test_fortran(08_sips d 3 r 5 ${h_real} ${s_real})
ENDIF()

### Parallel Fortran tests: dm + NTPoly ###
test_fortran(01_ntpoly d 0 r 6 ${h_real} ${s_real})
test_fortran(02_ntpoly d 0 c 6 ${h_cmplx} ${s_cmplx})
test_fortran(03_ntpoly d 1 r 6 ${h_real} ${s_real})
test_fortran(04_ntpoly d 1 c 6 ${h_cmplx} ${s_cmplx})
test_fortran(05_ntpoly d 2 r 6 ${h_real} ${s_real})
test_fortran(06_ntpoly d 2 c 6 ${h_cmplx} ${s_cmplx})
test_fortran(07_ntpoly d 3 r 6 ${h_real} ${s_real})
test_fortran(08_ntpoly d 3 c 6 ${h_cmplx} ${s_cmplx})

### Parallel Fortran tests: bse + BSEPACK ###
IF(ENABLE_BSEPACK)
  test_fortran(01_bsepack b 0 r 8 ${a_real} ${b_real})
  test_fortran(02_bsepack b 0 c 8 ${a_cmplx} ${b_cmplx})
ENDIF()

### Serial Fortran tests: ev + LAPACK ###
ADD_TEST(test_serial_01_lapack
  ${MPIEXEC_1P} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test e 0 r 1 ${h_real} ${s_real})
SET_TESTS_PROPERTIES(test_serial_01_lapack
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed")

ADD_TEST(test_serial_02_lapack
  ${MPIEXEC_1P} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test e 0 c 1 ${h_cmplx} ${s_cmplx})
SET_TESTS_PROPERTIES(test_serial_02_lapack
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed")

### Serial Fortran tests: ev + MAGMA ###
IF(ENABLE_MAGMA)
  ADD_TEST(test_serial_01_magma
    ${MPIEXEC_1P} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test e 0 r 7 ${h_real} ${s_real})
  SET_TESTS_PROPERTIES(test_serial_01_magma
    PROPERTIES PASS_REGULAR_EXPRESSION "Passed")

  ADD_TEST(test_serial_02_magma
    ${MPIEXEC_1P} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test e 0 c 7 ${h_cmplx} ${s_cmplx})
  SET_TESTS_PROPERTIES(test_serial_02_magma
    PROPERTIES PASS_REGULAR_EXPRESSION "Passed")
ENDIF()

### Fortran tests: matrix IO ###
ADD_TEST(test_matio_01
  ${MPIEXEC_1P} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test_rw r ${h_real} ${s_real})
SET_TESTS_PROPERTIES(test_matio_01
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed")

ADD_TEST(test_matio_02
  ${MPIEXEC_1P} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test_rw c ${h_cmplx} ${s_cmplx})
SET_TESTS_PROPERTIES(test_matio_02
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed")

ADD_TEST(test_matio_03
  ${MPIEXEC_NP} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test_rw r ${h_real} ${s_real})
SET_TESTS_PROPERTIES(test_matio_03
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed")

ADD_TEST(test_matio_04
  ${MPIEXEC_NP} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test_rw c ${h_cmplx} ${s_cmplx})
SET_TESTS_PROPERTIES(test_matio_04
  PROPERTIES PASS_REGULAR_EXPRESSION "Passed")

IF(ENABLE_C_TESTS)
  MACRO(test_c arg1 arg2 arg3 arg4 arg5 arg6 arg7)
    ADD_TEST(test_c_${arg1}
      ${MPIEXEC_NP} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/elsi_test_c
      ${arg2} ${arg3} ${arg4} ${arg5} ${arg6} ${arg7})
    SET_TESTS_PROPERTIES(test_c_${arg1}
      PROPERTIES PASS_REGULAR_EXPRESSION "Passed")
  ENDMACRO()

  ### Parallel C tests: ev + ELPA ###
  test_c(01_elpa e 0 r 1 ${h_real} ${s_real})
  test_c(02_elpa e 0 c 1 ${h_cmplx} ${s_cmplx})

  ### Parallel C tests: dm + ELPA ###
  test_c(03_elpa d 0 r 1 ${h_real} ${s_real})
  test_c(04_elpa d 0 c 1 ${h_cmplx} ${s_cmplx})

  ### Parallel C tests: dm + libOMM ###
  test_c(01_omm d 0 r 2 ${h_real} ${s_real})
  test_c(02_omm d 0 c 2 ${h_cmplx} ${s_cmplx})

  ### Parallel C tests: dm + PEXSI ###
  IF(ENABLE_PEXSI)
    test_c(01_pexsi d 0 r 3 ${h_real} ${s_real})
    test_c(02_pexsi d 0 c 3 ${h_cmplx} ${s_cmplx})
  ENDIF()

  ### Parallel C tests: dm + NTPoly ###
  test_c(01_ntpoly d 0 r 6 ${h_real} ${s_real})
  test_c(02_ntpoly d 0 c 6 ${h_cmplx} ${s_cmplx})
ENDIF()
