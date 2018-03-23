### Generate version info ###
MESSAGE(STATUS "Getting ELSI version")
EXECUTE_PROCESS(COMMAND rm -f ${PROJECT_SOURCE_DIR}/src/elsi_version.f90)
EXECUTE_PROCESS(COMMAND ${PROJECT_SOURCE_DIR}/src/elsi_version_writer.pl ${elsi_DATESTAMP}
  OUTPUT_FILE ${PROJECT_SOURCE_DIR}/src/elsi_version.f90
  OUTPUT_QUIET
)

### Source files ###
LIST(APPEND elsi_src
  elsi_c_interface.f90
  elsi_c2f.f90
  elsi_constants.f90
  elsi_datatype.f90
  elsi_dmp.f90
  elsi_interface.f90
  elsi_io.f90
  elsi_lapack.f90
  elsi_malloc.f90
  elsi_mat_io.f90
  elsi_mat_redist.f90
  elsi_mpi.f90
  elsi_mutator.f90
  elsi_occ.f90
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
)

IF(NOT ENABLE_PEXSI)
  LIST(APPEND elsi_src stub_pexsi.f90)
ENDIF()

IF(NOT ENABLE_SIPS)
  LIST(APPEND elsi_src stub_sips.f90)
ENDIF()

IF(NOT ELPA_FOUND)
  LIST(APPEND elsi_src elsi_elpa.f90)
ELSE()
  LIST(APPEND elsi_src elsi_aeo.f90)
ENDIF()

ADD_LIBRARY(elsi ${elsi_src})

IF(NOT ELPA_FOUND)
  ADD_DEPENDENCIES(elsi elpa)
ENDIF()

IF(NOT OMM_FOUND)
  ADD_DEPENDENCIES(elsi OMM)
ENDIF()

IF(ENABLE_PEXSI)
  ADD_DEPENDENCIES(elsi pexsi)
ENDIF()

IF(ENABLE_SIPS)
  ADD_DEPENDENCIES(elsi sips)
ENDIF()

IF(BUILD_SHARED_LIBS)
  IF(NOT OMM_FOUND)
    TARGET_LINK_LIBRARIES(elsi OMM MatrixSwitch)
  ELSE()
    TARGET_LINK_LIBRARIES(elsi ${OMM_LIB})
  ENDIF()
  IF(NOT ELPA_FOUND)
    TARGET_LINK_LIBRARIES(elsi elpa)
  ELSE()
    TARGET_LINK_LIBRARIES(elsi ${ELPA_LIB})
  ENDIF()
  IF(ENABLE_SIPS)
    TARGET_LINK_LIBRARIES(elsi sips ${SLEPC_LIB})
  ENDIF()
  IF(ENABLE_PEXSI)
    IF(NOT SUPERLU_FOUND)
      TARGET_LINK_LIBRARIES(elsi pexsi superlu_dist ${PTSCOTCH_LIB})
    ELSE()
      TARGET_LINK_LIBRARIES(elsi pexsi ${SUPERLU_LIB} ${PTSCOTCH_LIB})
    ENDIF()
  ENDIF()
  TARGET_LINK_LIBRARIES(elsi ${MATH_LIB})
  IF(ENABLE_MKL)
    SET_TARGET_PROPERTIES(elsi
      PROPERTIES LINK_FLAGS "-mkl=cluster")
  ENDIF()
ENDIF()

### Installation ###
INSTALL(TARGETS elsi
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

INSTALL(FILES elsi.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})

INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
  PATTERN "*.mod")

INSTALL(DIRECTORY ${ARCHIVE_OUTPUT_PATH}/pkgconfig
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
