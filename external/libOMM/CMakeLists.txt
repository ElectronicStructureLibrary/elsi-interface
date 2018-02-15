### MatrixSwitch source files ###
list(APPEND ms_src
  MatrixSwitch/MatrixSwitch_ops.f90
  MatrixSwitch/ms_mm_multiply.f90
  MatrixSwitch/ms_m_add.f90
  MatrixSwitch/ms_m_set.f90
  MatrixSwitch/ms_m_copy.f90
  MatrixSwitch/ms_m_register.f90
  MatrixSwitch/MatrixSwitch.f90
  )

### libOMM source files ###
list(APPEND omm_src
  libOMM/omm_params.f90
  libOMM/omm_rand.f90
  libOMM/omm_quartic.f90
  libOMM/omm_ops.f90
  libOMM/omm.f90
  libOMM/omm_callback.f90
  )

add_library(OMM ${omm_src})
add_library(MatrixSwitch ${ms_src})

add_dependencies(OMM MatrixSwitch)

if(NOT EXTERNAL_ELPA)
  add_dependencies(OMM elpa)
endif()

install(TARGETS MatrixSwitch
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

install(TARGETS OMM
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
