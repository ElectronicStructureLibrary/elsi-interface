SET(ELPA_LIB $ENV{ELPA_LIB})
SET(ELPA_INC $ENV{ELPA_INC})

IF(ELPA_INC AND ELPA_LIB)
  MESSAGE(STATUS "Using external ELPA")
  MESSAGE(STATUS "  ELPA_LIB: ${ELPA_LIB}")
  MESSAGE(STATUS "  ELPA_INC: ${ELPA_INC}")

  FOREACH(usr_lib ${ELPA_LIB})
    IF(NOT EXISTS ${usr_lib})
      MESSAGE(FATAL_ERROR "User provided ELPA library not found: ${usr_lib}")
    ENDIF()
  ENDFOREACH()

  FOREACH(usr_dir ${ELPA_INC})
    IF(NOT EXISTS ${usr_dir})
      MESSAGE(FATAL_ERROR "User provided ELPA include path not found: ${usr_dir}")
    ENDIF()

    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${usr_dir}")
  ENDFOREACH()

  SET(ELPA_FOUND TRUE)
ELSE()
  MESSAGE(STATUS "Enabling internal ELPA")

  SET(ELPA_FOUND FALSE)
ENDIF()
