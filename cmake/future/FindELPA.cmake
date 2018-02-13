find_path(ELPA_INC_DIR
  NAMES elpa.mod elpa1.mod elpa2.mod
  HINTS ${ELPA_INC_DIR} ${ELPA_DIR}
  PATH_SUFFIXES include
  )

find_library(ELPA_LIB
  NAMES elpa
  HINTS ${ELPA_LIB_DIR} ${ELPA_DIR}
  PATH_SUFFIXES lib
  )

find_package_handle_standard_args(ELPA
  REQUIRED_VARS ELPA_INC_DIR ELPA_LIB)

if(ELPA_FOUND)
  set(ELPA_LIBRARIES ${ELPA_LIB})
  set(ELPA_INCLUDE_DIRS ${ELPA_INC_DIR})
endif()
