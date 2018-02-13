list(APPEND sips_src
  ${sips_dir}/src/m_qetsc.F90
  )

add_library(sips ${sips_src})

install(TARGETS sips DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
