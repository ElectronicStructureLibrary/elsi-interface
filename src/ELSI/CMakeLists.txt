list(APPEND elsi_src
  ${elsi_dir}/src/elsi_version.f90
  ${elsi_dir}/src/elsi_mu.f90
  ${elsi_dir}/src/elsi_solver.f90
  ${elsi_dir}/src/elsi_utils.f90
  ${elsi_dir}/src/elsi_interface.f90
  ${elsi_dir}/src/elsi_datatype.f90
  ${elsi_dir}/src/elsi_mutator.f90
  ${elsi_dir}/src/elsi_matconv.f90
  ${elsi_dir}/src/elsi_setup.f90
  ${elsi_dir}/src/elsi_io.f90
  ${elsi_dir}/src/elsi_constants.f90
  ${elsi_dir}/src/elsi_matio.f90
  ${elsi_dir}/src/elsi_chess.f90
  ${elsi_dir}/src/elsi_dmp.f90
  ${elsi_dir}/src/elsi_lapack.f90
  ${elsi_dir}/src/elsi_malloc.f90
  ${elsi_dir}/src/elsi_mpi.f90
  ${elsi_dir}/src/elsi_omm.f90
  ${elsi_dir}/src/elsi_pexsi.f90
  ${elsi_dir}/src/elsi_precision.f90
  ${elsi_dir}/src/elsi_sips.f90
  ${elsi_dir}/src/elsi_sort.f90
  ${elsi_dir}/src/elsi_timings.f90
  ${elsi_dir}/src/stub_chess.f90
  )

if(NOT ENABLE_PEXSI)
  list(APPEND elsi_src
    ${elsi_dir}/src/stub_pexsi.f90
    )
endif()
if(NOT ENABLE_SIPS)
  list(APPEND elsi_src
    ${elsi_dir}/src/stub_sips.f90
    )
endif()

if(NOT EXTERNAL_ELPA)
  list(APPEND elsi_src
    ${elsi_dir}/src/elsi_elpa.f90
    )
else()
  list(APPEND elsi_src
    ${elsi_dir}/src/elsi_aeo.f90
    )
endif()

if(ENABLE_C_INTERFACE)
  list(APPEND elsi_src
    ${elsi_dir}/c_binding/elsi_c_interface.f90
    ${elsi_dir}/c_binding/elsi_c2f.f90
    )
endif()

add_library(elsi ${elsi_src})

if(INTERNAL_ELPA)
  add_dependencies(elsi elpa)
endif()

if(INTERNAL_OMM)
  add_dependencies(elsi omm)
endif()

install(TARGETS elsi
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(DIRECTORY ${CMAKE_BINARY_DIR}/include
  DESTINATION ${CMAKE_INSTALL_PREFIX})
if(ENABLE_C_INTERFACE)
  install(FILES ${elsi_dir}/c_binding/elsi.h
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
endif()

# For tests
#add_executable(test_dm_complex ${elsi_dir}/test/test_dm_complex.f90)
#target_link_libraries(test_dm_complex elsi elpa omm MatrixSwitch ${LIBS})
#install(TARGETS test_dm_complex DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
#include(CTest)
#add_test(TutorialRuns Tutorial 25)
#add_test(TutorialComp25 Tutorial 25)
#set_tests_properties(TutorialComp25 PROPERTIES PASS_REGULAR_EXPRESSION "25 is 5")

#define a macro to simplify adding tests, then use it
#macro (do_test arg result)
#  add_test(TutorialComp${arg} Tutorial ${arg})
#  set_tests_properties(TutorialComp${arg}
#    PROPERTIES PASS_REGULAR_EXPRESSION ${result})
#endmacro (do_test)
#do_test(25 "25 is 5")
#do_test(-25 "-25 is 0")
