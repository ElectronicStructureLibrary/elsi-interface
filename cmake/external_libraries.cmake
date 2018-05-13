# Taken from FHI-aims with permission of copyright holders

# Convert a string into a list.
# Do nothing if input is already a list.
MACRO(convert_to_list var)
  LIST(LENGTH ${var} _length)
  IF(_length EQUAL 1)
    STRING(REPLACE " " ";" ${var} ${${var}})
  ENDIF()
ENDMACRO()

# Go through directories listed in LIBRARY_PATHS and turn the entries of
# LIBRARIES into actual targets
FUNCTION(generate_library_targets _PATHS _LIBRARIES)
  FOREACH(LIB ${${_LIBRARIES}})
    FIND_LIBRARY(LIB_FULLPATH ${LIB} PATHS ${${_PATHS}})

    IF(LIB_FULLPATH)
      MESSAGE(STATUS "Found ${LIB_FULLPATH}")
      ADD_LIBRARY(${LIB} UNKNOWN IMPORTED)
      SET_TARGET_PROPERTIES(${LIB} PROPERTIES IMPORTED_LOCATION ${LIB_FULLPATH})
    ELSE()
      MESSAGE(FATAL_ERROR "${Magenta}Could not find ${LIB}${ColorReset}")
    ENDIF()

    UNSET(LIB_FULLPATH CACHE)
  ENDFOREACH()
ENDFUNCTION()
