### Source file ###
list(APPEND sips_src src/m_qetsc.F90)

add_library(sips ${sips_src})

install(TARGETS sips
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
