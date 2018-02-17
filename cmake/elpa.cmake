SET(ELPA_LIB $ENV{ELPA_LIB})
SET(ELPA_INC $ENV{ELPA_INC})

IF(ELPA_INC AND ELPA_LIB)
  MESSAGE(STATUS "Using external ELPA")
  MESSAGE(STATUS "  ${GREEN}ELPA_LIB${COLORRESET}: ${ELPA_LIB}")
  MESSAGE(STATUS "  ${GREEN}ELPA_INC${COLORRESET}: ${ELPA_INC}")

  FOREACH(usr_lib ${ELPA_LIB})
    IF(NOT EXISTS ${usr_lib})
      MESSAGE(FATAL_ERROR "${MAGENTA}User provided ELPA library not found: ${usr_lib}${COLORRESET}")
    ENDIF()
  ENDFOREACH()

  FOREACH(usr_dir ${ELPA_INC})
    IF(NOT EXISTS ${usr_dir})
      MESSAGE(FATAL_ERROR "${MAGENTA}User provided ELPA include path not found: ${usr_dir}${COLORRESET}")
    ENDIF()

    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${usr_dir}")
  ENDFOREACH()

  SET(ELPA_FOUND TRUE)
ELSE()
  MESSAGE(STATUS "Enabling internal ELPA")

  SET(ELPA_FOUND FALSE)
ENDIF()
