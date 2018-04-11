### Source file ###
LIST(APPEND fortjson_src
  src/fortjson.f90
  src/fjson_constants.f90
  src/fjson_datatype.f90
  src/fjson_precision.f90
  src/fjson_rw.f90
  src/fjson_setup.f90
  src/fjson_utilities.f90
)

ADD_LIBRARY(fortjson ${fortjson_src})

IF(NOT "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" STREQUAL "${PROJECT_BINARY_DIR}/lib")
  INSTALL(TARGETS fortjson
    DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
ENDIF()
