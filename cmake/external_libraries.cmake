# Taken from FHI-aims with permission of copyright holders

# Convert a string into a list.
# Do nothing if input is already a list.
MACRO(convert_to_list var)
  LIST(LENGTH ${var} _length)
  IF(_length EQUAL 1)
    STRING(REPLACE " " ";" ${var} ${${var}})
  ENDIF()
ENDMACRO()

# Go through directories listed in LIB_PATHS and turn entries of LIBS into actual targets.
FUNCTION(generate_library_targets _PATHS _LIBRARIES)
  FOREACH(LIB ${${_LIBRARIES}})
    # If the target already exists, skip it.
    # This can occur when the same library is set multiple times in LIBS due to cyclic dependencies.
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
