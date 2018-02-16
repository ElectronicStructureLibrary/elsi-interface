### Source file ###
LIST(APPEND sips_src src/m_qetsc.F90)

ADD_LIBRARY(sips ${sips_src})

INSTALL(TARGETS sips
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
