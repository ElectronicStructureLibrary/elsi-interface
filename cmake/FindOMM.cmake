find_path(OMM_INC_DIR
  NAMES MatrixSwitch.mod
  HINTS ${OMM_INC_DIR} ${OMM_DIR}
  PATH_SUFFIXES include
  )

find_library(OMM_LIB
  NAMES omm MatrixSwitch
  HINTS ${OMM_LIB_DIR} ${OMM_DIR}
  PATH_SUFFIXES lib
  )

find_package_handle_standard_args(OMM
  REQUIRED_VARS OMM_INC_DIR OMM_LIB)

if(OMM_FOUND)
  set(OMM_LIBRARIES ${OMM_LIB})
  set(OMM_INCLUDE_DIRS ${OMM_INC_DIR})
endif()
