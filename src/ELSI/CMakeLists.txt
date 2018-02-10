enable_language(Fortran)

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

if(NOT enable_pexsi)
  list(APPEND elsi_src ${elsi_dir}/src/stub_pexsi.f90)
endif()
if(NOT enable_sips)
  list(APPEND elsi_src ${elsi_dir}/src/stub_sips.f90)
endif()

if(NOT external_elpa)
  list(APPEND elsi_src ${elsi_dir}/src/elsi_elpa.f90)
else()
  list(APPEND elsi_src ${elsi_dir}/src/elsi_aeo.f90)
endif()

add_library(elsi ${elsi_src})

target_include_directories(elsi PUBLIC)
target_compile_options(elsi PUBLIC)

install(TARGETS elsi DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

add_dependencies(elsi elpa omm)
