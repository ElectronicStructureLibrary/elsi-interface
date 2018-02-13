list(APPEND ms_src
  ${omm_dir}/MatrixSwitch/src/MatrixSwitch_ops.f90
  ${omm_dir}/MatrixSwitch/src/ms_mm_multiply.f90
  ${omm_dir}/MatrixSwitch/src/ms_m_add.f90
  ${omm_dir}/MatrixSwitch/src/ms_m_set.f90
  ${omm_dir}/MatrixSwitch/src/ms_m_copy.f90
  ${omm_dir}/MatrixSwitch/src/ms_m_register.f90
  ${omm_dir}/MatrixSwitch/src/MatrixSwitch.f90
  )

list(APPEND omm_src
  ${omm_dir}/libOMM/src/omm_params.f90
  ${omm_dir}/libOMM/src/omm_rand.f90
  ${omm_dir}/libOMM/src/omm_quartic.f90
  ${omm_dir}/libOMM/src/omm_ops.f90
  ${omm_dir}/libOMM/src/omm.f90
  ${omm_dir}/libOMM/src/omm_callback.f90
  )

add_library(OMM ${omm_src})
add_library(MatrixSwitch ${ms_src})

add_dependencies(OMM MatrixSwitch)

if(INTERNAL_ELPA)
  add_dependencies(OMM elpa)
endif()

install(TARGETS MatrixSwitch DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(TARGETS OMM DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
