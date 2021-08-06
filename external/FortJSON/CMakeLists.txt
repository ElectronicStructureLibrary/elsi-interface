### Source files ###
LIST(APPEND fortjson_src
  src/fortjson.f90
  src/fjson_constants.f90
  src/fjson_datatype.f90
  src/fjson_precision.f90
  src/fjson_rw.f90
  src/fjson_setup.f90
  src/fjson_utilities.f90)

ADD_LIBRARY(fortjson ${fortjson_src})

IF(ELSI_USE_FIND_PACKAGE_MPI)
  find_package(MPI)
  TARGET_LINK_LIBRARIES(fortjson PUBLIC MPI::MPI_C MPI::MPI_Fortran)
ENDIF()

TARGET_INCLUDE_DIRECTORIES(fortjson PRIVATE
  ${INC_PATHS}
  ${CMAKE_Fortran_MODULE_DIRECTORY})

INSTALL(TARGETS fortjson
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
