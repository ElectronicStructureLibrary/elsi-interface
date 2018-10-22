# Taken from http://stackoverflow.com/a/19578320
STRING(ASCII 27 Esc)

SET(COLORRESET "${Esc}[m")
SET(COLORBOLD "${Esc}[1m")
SET(RED "${Esc}[31m")
SET(GREEN "${Esc}[32m")
SET(YELLOW "${Esc}[33m")
SET(BLUE "${Esc}[34m")
SET(MAGENTA "${Esc}[35m")
SET(CYAN "${Esc}[36m")
SET(WHITE "${Esc}[37m")
SET(BOLDRED "${Esc}[1;31m")
SET(BOLDGREEN "${Esc}[1;32m")
SET(BOLDYELLOW "${Esc}[1;33m")
SET(BOLDBLUE "${Esc}[1;34m")
SET(BOLDMAGENTA "${Esc}[1;35m")
SET(BOLDCYAN "${Esc}[1;36m")
SET(BOLDWHITE "${Esc}[1;37m")

# Taken from FHI-aims with permission of copyright holders
# Convert a string into a list
# Do nothing if input is already a list
MACRO(convert_to_list var)
  LIST(LENGTH ${var} _length)
  IF(_length EQUAL 1)
    STRING(REPLACE " " ";" ${var} ${${var}})
  ENDIF()
ENDMACRO()

# Taken from FHI-aims with permission of copyright holders
# Go through directories listed in LIB_PATHS and turn entries of LIBS into actual targets
FUNCTION(generate_library_targets _PATHS _LIBRARIES)
  FOREACH(LIB ${${_LIBRARIES}})
    # If the target already exists, skip it.
    # This can occur when the same library is set multiple times in LIBS due to cyclic dependencies
    IF(NOT TARGET ${LIB})
      FIND_LIBRARY(LIB_FULLPATH ${LIB} HINTS ${${_PATHS}})
      IF(LIB_FULLPATH)
        MESSAGE(STATUS "Found ${LIB_FULLPATH}")
        ADD_LIBRARY(${LIB} UNKNOWN IMPORTED)
        SET_TARGET_PROPERTIES(${LIB} PROPERTIES IMPORTED_LOCATION ${LIB_FULLPATH})
        UNSET(LIB_FULLPATH CACHE)
      ELSE()
        MESSAGE(FATAL_ERROR "${Magenta}Could not find ${LIB}${ColorReset}")
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDFUNCTION()
