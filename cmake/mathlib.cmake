SET(MATH_LIB $ENV{MATH_LIB})
SET(MATH_INC $ENV{MATH_INC})

IF(MATH_LIB)
  MESSAGE(STATUS "Using linear algebra libraries")
  MESSAGE(STATUS "  ${GREEN}MATH_LIB${COLORRESET}: ${MATH_LIB}")

  FOREACH(usr_lib ${MATH_LIB})
    IF(NOT EXISTS ${usr_lib})
      MESSAGE(FATAL_ERROR "${MAGENTA}User provided linear algebra library not found: ${usr_lib}${COLORRESET}")
    ENDIF()
  ENDFOREACH()

  IF(MATH_INC)
    MESSAGE(STATUS "  ${GREEN}MATH_INC${COLORRESET}: ${MATH_INC}")
    FOREACH(usr_dir ${MATH_INC})
      IF(NOT EXISTS ${usr_dir})
        MESSAGE(FATAL_ERROR "${MAGENTA}User provided linear algebra include path not found: ${usr_dir}${COLORRESET}")
      ENDIF()

      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${usr_dir}")
    ENDFOREACH()
  ENDIF()

  SET(MATH_FOUND TRUE)
ELSE()
  MESSAGE(WARNING "${MAGENTA}Linear algebra libraries not provided${COLORRESET}")

  SET(MATH_FOUND FALSE)
ENDIF()
