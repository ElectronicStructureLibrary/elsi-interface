SET(MATH_LIB $ENV{MATH_LIB})
SET(MATH_INC $ENV{MATH_INC})

IF(MATH_LIB)
  MESSAGE(STATUS "Using linear algebra libraries")
  MESSAGE(STATUS "  MATH_LIB: ${MATH_LIB}")
  MESSAGE(STATUS "  MATH_INC: ${MATH_INC}")

  FOREACH(usr_lib ${MATH_LIB})
    IF(NOT EXISTS ${usr_lib})
      MESSAGE(FATAL_ERROR "User provided linear algebra library not found: ${usr_lib}")
    ENDIF()
  ENDFOREACH()

  FOREACH(usr_dir ${MATH_INC})
    IF(NOT EXISTS ${usr_dir})
      MESSAGE(FATAL_ERROR "User provided linear algebra include path not found: ${usr_dir}")
    ENDIF()

    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${usr_dir}")
  ENDFOREACH()

  SET(MATH_FOUND TRUE)
ELSE()
  MESSAGE(WARNING "Linear algebra libraries not provided")

  SET(MATH_FOUND FALSE)
ENDIF()
