### MatrixSwitch source files ###
LIST(APPEND ms_src
  src/MatrixSwitch/MatrixSwitch_ops.f90
  src/MatrixSwitch/ms_mm_multiply.f90
  src/MatrixSwitch/ms_m_add.f90
  src/MatrixSwitch/ms_m_set.f90
  src/MatrixSwitch/ms_m_copy.f90
  src/MatrixSwitch/ms_m_register.f90
  src/MatrixSwitch/MatrixSwitch.f90)

### libOMM source files ###
LIST(APPEND omm_src
  src/libOMM/omm_params.f90
  src/libOMM/omm_rand.f90
  src/libOMM/omm_quartic.f90
  src/libOMM/omm_ops.f90
  src/libOMM/omm.f90
  src/libOMM/omm_callback.f90)

ADD_LIBRARY(OMM ${omm_src})
ADD_LIBRARY(MatrixSwitch ${ms_src})

ADD_DEPENDENCIES(OMM MatrixSwitch)

IF(NOT USE_EXTERNAL_ELPA)
  ADD_DEPENDENCIES(OMM elpa)
ENDIF()

TARGET_LINK_LIBRARIES(OMM PRIVATE MatrixSwitch ${LIBS})
IF(NOT USE_EXTERNAL_ELPA)
  TARGET_LINK_LIBRARIES(OMM PRIVATE elpa)
ENDIF()
TARGET_LINK_LIBRARIES(MatrixSwitch PRIVATE ${LIBS})

TARGET_INCLUDE_DIRECTORIES(OMM PRIVATE ${INC_PATHS})

INSTALL(TARGETS MatrixSwitch
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
INSTALL(TARGETS OMM
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
