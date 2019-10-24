### Generate version info ###
INCLUDE(GetGitRevision)
get_git_head_revision(GIT_COMMIT)

STRING(TIMESTAMP DATETIME UTC)

CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/src/elsi_version.f90.in
  ${PROJECT_BINARY_DIR}/generated/elsi_version.f90
  @ONLY)

### Source files ###
LIST(APPEND elsi_src
  ${PROJECT_BINARY_DIR}/generated/elsi_version.f90
  elsi_c_interface.f90
  elsi_constant.f90
  elsi_datatype.f90
  elsi_decision.f90
  elsi_geo.f90
  elsi_get.f90
  elsi_input.f90
  elsi_interface.f90
  elsi_lapack.f90
  elsi_malloc.f90
  elsi_mpi.f90
  elsi_ntpoly.f90
  elsi_occ.f90
  elsi_omm.f90
  elsi_output.f90
  elsi_pexsi.f90
  elsi_precision.f90
  elsi_redist.f90
  elsi_rw.f90
  elsi_set.f90
  elsi_setup.f90
  elsi_solver.f90
  elsi_sort.f90
  elsi_util.f90)

SET_SOURCE_FILES_PROPERTIES(${PROJECT_BINARY_DIR}/generated/elsi_version.f90
  PROPERTIES GENERATED TRUE)

IF(NOT ENABLE_PEXSI)
  LIST(APPEND elsi_src stub_pexsi.f90)
ENDIF()
IF(NOT ENABLE_EIGENEXA)
  LIST(APPEND elsi_src stub_eigenexa.f90)
ELSE()
  LIST(APPEND elsi_src elsi_eigenexa.f90)
ENDIF()
IF(NOT ENABLE_SIPS)
  LIST(APPEND elsi_src stub_sips.f90)
ELSE()
  LIST(APPEND elsi_src elsi_sips.f90)
ENDIF()
IF(NOT ENABLE_MAGMA)
  LIST(APPEND elsi_src stub_magma.f90)
ELSE()
  LIST(APPEND elsi_src elsi_magma.f90)
ENDIF()
IF(NOT USE_EXTERNAL_ELPA)
  LIST(APPEND elsi_src
    elsi_elpa.f90
    stub_aeo.f90)
ELSE()
  LIST(APPEND elsi_src elsi_aeo.f90)
ENDIF()
IF(USE_EXTERNAL_NTPOLY)
  LIST(APPEND elsi_src ${PROJECT_SOURCE_DIR}/external/NTPoly/src/Interface.f90)
ENDIF()

ADD_LIBRARY(elsi ${elsi_src})

IF(NOT USE_EXTERNAL_ELPA)
  ADD_DEPENDENCIES(elsi elpa)
ENDIF()
IF(NOT USE_EXTERNAL_OMM)
  ADD_DEPENDENCIES(elsi OMM)
ENDIF()
IF(ENABLE_PEXSI)
  ADD_DEPENDENCIES(elsi pexsi)
ENDIF()
IF(ENABLE_SIPS)
  ADD_DEPENDENCIES(elsi sips)
ENDIF()
IF(NOT USE_EXTERNAL_NTPOLY)
  ADD_DEPENDENCIES(elsi NTPoly)
ENDIF()
ADD_DEPENDENCIES(elsi fortjson)

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
  IF(NOT USE_EXTERNAL_PEXSI)
    TARGET_LINK_LIBRARIES(elsi PRIVATE
      pexsi
      superlu_dist
      ptscotchparmetis
      ptscotch
      ptscotcherr
      scotchmetis
      scotch
      scotcherr)
  ENDIF()
ENDIF()
IF(NOT USE_EXTERNAL_NTPOLY)
  TARGET_LINK_LIBRARIES(elsi PRIVATE NTPoly)
ENDIF()

TARGET_LINK_LIBRARIES(elsi PRIVATE ${LIBS})

TARGET_INCLUDE_DIRECTORIES(elsi PRIVATE
  ${INC_PATHS}
  ${CMAKE_Fortran_MODULE_DIRECTORY})

### Generate pkg-config file ###
GET_TARGET_PROPERTY(PKG_LIBS elsi LINK_LIBRARIES)
SET(PKG_LIBS ";elsi;${PKG_LIBS}")
FOREACH(_lib ${PKG_LIBS})
  GET_FILENAME_COMPONENT(_lib_we ${_lib} NAME_WE)
  STRING(REPLACE "${_lib}" "${_lib_we}" PKG_LIBS "${PKG_LIBS}")
ENDFOREACH()
STRING(REPLACE ";lib" ";" PKG_LIBS "${PKG_LIBS}")
STRING(REPLACE ";" " -l" PKG_LIBS "${PKG_LIBS}")
IF(LIB_PATHS)
  SET(PKG_LIB_PATHS ";${LIB_PATHS}")
  STRING(REPLACE ";" " -L" PKG_LIB_PATHS "${PKG_LIB_PATHS}")
ENDIF()
IF(INC_PATHS)
  SET(PKG_INC_PATHS ";${INC_PATHS}")
  STRING(REPLACE ";" " -I" PKG_INC_PATHS "${PKG_INC_PATHS}")
ENDIF()

CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/.elsi.pc.in
  ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/pkgconfig/elsi.pc
  @ONLY)

### Installation ###
# Libraries
INSTALL(TARGETS elsi
  EXPORT elsiConfig
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
# Header files
INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PATTERN "*.mod")
INSTALL(FILES elsi.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# pkg-config file
INSTALL(FILES ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/pkgconfig/elsi.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
# elsiConfig.cmake
INSTALL(EXPORT elsiConfig
  NAMESPACE elsi::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/elsi)
# elsiConfigVersion.cmake
INCLUDE(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(
  ${PROJECT_BINARY_DIR}/elsiConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
INSTALL(FILES ${PROJECT_BINARY_DIR}/elsiConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/elsi)
# Export ELSI
TARGET_INCLUDE_DIRECTORIES(elsi INTERFACE
  $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
EXPORT(EXPORT elsiConfig NAMESPACE elsi::
  FILE ${PROJECT_BINARY_DIR}/elsiConfig.cmake)
