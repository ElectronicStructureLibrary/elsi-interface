SET(OMM_LIB $ENV{OMM_LIB})
SET(OMM_INC $ENV{OMM_INC})

IF(OMM_INC AND OMM_LIB)
  MESSAGE(STATUS "Using external libOMM")
  MESSAGE(STATUS "${GREEN}OMM_LIB${COLORRESET}: ${OMM_LIB}")
  MESSAGE(STATUS "${GREEN}OMM_INC${COLORRESET}: ${OMM_INC}")

  FOREACH(usr_lib ${OMM_LIB})
    IF(NOT EXISTS ${usr_lib})
      MESSAGE(FATAL_ERROR "${MAGENTA}User provided libOMM library not found: ${usr_lib}${COLORRESET}")
    ENDIF()
  ENDFOREACH()

  FOREACH(usr_dir ${OMM_INC})
    IF(NOT EXISTS ${usr_dir})
      MESSAGE(FATAL_ERROR "${MAGENTA}User provided libOMM include path not found: ${usr_dir}${COLORRESET}")
    ENDIF()

    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${usr_dir}")
  ENDFOREACH()

  SET(OMM_FOUND TRUE)
ELSE()
  MESSAGE(STATUS "Enabling internal libOMM")

  SET(OMM_FOUND FALSE)
ENDIF()
