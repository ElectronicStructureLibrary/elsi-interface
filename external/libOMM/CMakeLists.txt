### MatrixSwitch source files ###
LIST(APPEND ms_src
  MatrixSwitch/MatrixSwitch_ops.f90
  MatrixSwitch/ms_mm_multiply.f90
  MatrixSwitch/ms_m_add.f90
  MatrixSwitch/ms_m_set.f90
  MatrixSwitch/ms_m_copy.f90
  MatrixSwitch/ms_m_register.f90
  MatrixSwitch/MatrixSwitch.f90
)

### libOMM source files ###
LIST(APPEND omm_src
  libOMM/omm_params.f90
  libOMM/omm_rand.f90
  libOMM/omm_quartic.f90
  libOMM/omm_ops.f90
  libOMM/omm.f90
  libOMM/omm_callback.f90
)

ADD_LIBRARY(OMM ${omm_src})
ADD_LIBRARY(MatrixSwitch ${ms_src})

ADD_DEPENDENCIES(OMM MatrixSwitch)

IF(NOT ELPA_FOUND)
  ADD_DEPENDENCIES(OMM elpa)
ENDIF()

INSTALL(TARGETS MatrixSwitch
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

INSTALL(TARGETS OMM
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
