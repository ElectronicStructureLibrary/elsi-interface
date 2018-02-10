list(APPEND ms_src
  ${omm_dir}/MatrixSwitch/src/MatrixSwitch_ops.f90
  ${omm_dir}/MatrixSwitch/src/ms_mm_multiply.f90
  ${omm_dir}/MatrixSwitch/src/ms_m_add.f90
  ${omm_dir}/MatrixSwitch/src/ms_m_set.f90
  ${omm_dir}/MatrixSwitch/src/ms_m_copy.f90
  ${omm_dir}/MatrixSwitch/src/ms_m_register.f90
  ${omm_dir}/MatrixSwitch/src/MatrixSwitch.f90
  )

add_library(MatrixSwitch ${ms_src})

target_include_directories(MatrixSwitch PUBLIC)
target_compile_options(MatrixSwitch PUBLIC)

install(TARGETS MatrixSwitch DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

list(APPEND omm_src
  ${omm_dir}/libOMM/src/omm_params.f90
  ${omm_dir}/libOMM/src/omm_rand.f90
  ${omm_dir}/libOMM/src/omm_quartic.f90
  ${omm_dir}/libOMM/src/omm_ops.f90
  ${omm_dir}/libOMM/src/omm.f90
  ${omm_dir}/libOMM/src/omm_callback.f90
  )

add_library(omm ${omm_src})

target_include_directories(omm PUBLIC)
target_compile_options(omm PUBLIC)

install(TARGETS omm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

add_dependencies(omm elpa MatrixSwitch)
