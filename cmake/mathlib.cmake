IF(MATH_LIB)
  MESSAGE(STATUS "Using linear algebra libraries")
  MESSAGE(STATUS "${GREEN}MATH_LIB${COLORRESET}: ${MATH_LIB}")

  FOREACH(usr_lib ${MATH_LIB})
    IF(NOT EXISTS ${usr_lib})
      MESSAGE(FATAL_ERROR "${MAGENTA}User provided linear algebra library not found: ${usr_lib}${COLORRESET}")
    ENDIF()
  ENDFOREACH()

  SET(MATH_FOUND TRUE)
ELSE()
  IF(NOT ENABLE_MKL)
    MESSAGE(WARNING "${MAGENTA}Linear algebra libraries not provided${COLORRESET}")
  ENDIF()

  SET(MATH_FOUND FALSE)
ENDIF()
