### Generate version info ###
ADD_CUSTOM_COMMAND(OUTPUT ${PROJECT_BINARY_DIR}/generated/elsi_version.f90
  COMMAND perl ${PROJECT_SOURCE_DIR}/src/elsi_version_writer.pl ${elsi_DATESTAMP} > ${PROJECT_BINARY_DIR}/generated/elsi_version.f90
  VERBATIM
)

### Source files ###
LIST(APPEND elsi_src
  ${PROJECT_BINARY_DIR}/generated/elsi_version.f90
  elsi_c_interface.f90
  elsi_c2f.f90
  elsi_constants.f90
  elsi_datatype.f90
  elsi_dmp.f90
  elsi_interface.f90
  elsi_io.f90
  elsi_json.f90
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
  elsi_utils.f90
)

SET_SOURCE_FILES_PROPERTIES(${PROJECT_BINARY_DIR}/generated/elsi_version.f90
  PROPERTIES GENERATED TRUE)

IF(NOT ENABLE_PEXSI)
  LIST(APPEND elsi_src stub_pexsi.f90)
ENDIF()

IF(NOT ENABLE_SIPS)
  LIST(APPEND elsi_src stub_sips.f90)
ENDIF()

IF(NOT USE_EXTERNAL_ELPA)
  LIST(APPEND elsi_src
    elsi_elpa.f90
    stub_aeo.f90
  )
ELSE()
  LIST(APPEND elsi_src elsi_aeo.f90)
ENDIF()

ADD_LIBRARY(elsi ${elsi_src})

IF(NOT USE_EXTERNAL_ELPA)
  ADD_DEPENDENCIES(elsi elpa)
ENDIF()
IF(NOT USE_EXTERNAL_OMM)
  ADD_DEPENDENCIES(elsi OMM)
ENDIF()
ADD_DEPENDENCIES(elsi fortjson)
IF(ENABLE_PEXSI)
  ADD_DEPENDENCIES(elsi pexsi)
ENDIF()
IF(ENABLE_SIPS)
  ADD_DEPENDENCIES(elsi sips)
ENDIF()

IF(NOT USE_EXTERNAL_OMM)
  TARGET_LINK_LIBRARIES(elsi PRIVATE OMM MatrixSwitch)
ENDIF()
IF(NOT USE_EXTERNAL_ELPA)
  TARGET_LINK_LIBRARIES(elsi PRIVATE elpa)
ENDIF()
TARGET_LINK_LIBRARIES(elsi PRIVATE fortjson)
IF(ENABLE_SIPS)
  TARGET_LINK_LIBRARIES(elsi PRIVATE sips)
ENDIF()
IF(ENABLE_PEXSI)
  TARGET_LINK_LIBRARIES(elsi PRIVATE pexsi)
  IF(NOT USE_EXTERNAL_SUPERLU)
    TARGET_LINK_LIBRARIES(elsi PRIVATE
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
TARGET_LINK_LIBRARIES(elsi PRIVATE ${LIBS})

TARGET_INCLUDE_DIRECTORIES(elsi PRIVATE ${LIB_PATH})

### Installation ###
IF(NOT "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" STREQUAL "${PROJECT_BINARY_DIR}/lib")
  INSTALL(TARGETS elsi
    DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
ENDIF()

IF(NOT "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}" STREQUAL "${PROJECT_BINARY_DIR}/include")
  INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    PATTERN "*.mod")
ENDIF()

INSTALL(FILES elsi.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
