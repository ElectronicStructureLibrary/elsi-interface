SET(PTSCOTCH_LIB $ENV{PTSCOTCH_LIB})
SET(PTSCOTCH_INC $ENV{PTSCOTCH_INC})

IF(PTSCOTCH_INC AND PTSCOTCH_LIB)
  MESSAGE(STATUS "Using external PtScotch")
  MESSAGE(STATUS "  PTSCOTCH_LIB: ${PTSCOTCH_LIB}")
  MESSAGE(STATUS "  PTSCOTCH_INC: ${PTSCOTCH_INC}")

  FOREACH(usr_lib ${PTSCOTCH_LIB})
    IF(NOT EXISTS ${usr_lib})
      MESSAGE(FATAL_ERROR "User provided PtScotch library not found: ${usr_lib}")
    ENDIF()
  ENDFOREACH()

  FOREACH(usr_dir ${PTSCOTCH_INC})
    IF(NOT EXISTS ${usr_dir})
      MESSAGE(FATAL_ERROR "User provided PtScotch include path not found: ${usr_dir}")
    ENDIF()

    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I${usr_dir}")
  ENDFOREACH()

  SET(PTSCOTCH_FOUND TRUE)
ELSE()
  MESSAGE(FATAL_ERROR "PtScotch not provided")

  SET(PTSCOTCH_FOUND FALSE)
ENDIF()
