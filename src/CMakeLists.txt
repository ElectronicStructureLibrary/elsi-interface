### Source files ###
list(APPEND elsi_src
  elsi_c_interface.f90
  elsi_c2f.f90
  elsi_chess.f90
  elsi_constants.f90
  elsi_datatype.f90
  elsi_dmp.f90
  elsi_interface.f90
  elsi_io.f90
  elsi_lapack.f90
  elsi_malloc.f90
  elsi_matconv.f90
  elsi_matio.f90
  elsi_mpi.f90
  elsi_mu.f90
  elsi_mutator.f90
  elsi_omm.f90
  elsi_pexsi.f90
  elsi_precision.f90
  elsi_setup.f90
  elsi_sips.f90
  elsi_solver.f90
  elsi_sort.f90
  elsi_timings.f90
  elsi_version.f90
  elsi_utils.f90
  stub_chess.f90
  )

if(NOT ENABLE_PEXSI)
  list(APPEND elsi_src stub_pexsi.f90)
endif()

if(NOT ENABLE_SIPS)
  list(APPEND elsi_src stub_sips.f90)
endif()

if(NOT EXTERNAL_ELPA)
  list(APPEND elsi_src elsi_elpa.f90)
else()
  list(APPEND elsi_src elsi_aeo.f90)
endif()

add_library(elsi ${elsi_src})

if(NOT EXTERNAL_ELPA)
  add_dependencies(elsi elpa)
endif()

if(NOT EXTERNAL_OMM)
  add_dependencies(elsi OMM)
endif()

### Installation ###
install(TARGETS elsi
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

install(FILES elsi.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
  PATTERN "*.mod")
