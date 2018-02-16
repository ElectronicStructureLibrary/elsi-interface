enable_language(CXX)
include_directories(include)

### Source files ###
list(APPEND pexsi_src
  src/TreeBcast.cpp
  src/blas.cpp
  src/f_interface.f90
  src/getPole.cpp
  src/get_perm_c_parmetis.c
  src/global.cpp
  src/interface.cpp
  src/lapack.cpp
  src/mpi_interf.cpp
  src/pdsymbfact.c
  src/pole.cpp
  src/ppexsi.cpp
  src/pzsymbfact.c
  src/superlu_dist_internal_complex.cpp
  src/superlu_dist_internal_real.cpp
  src/timer.cpp
  src/utility.cpp
)

add_library(pexsi ${pexsi_src})

if(NOT EXTERNAL_SUPERLU)
  add_dependencies(pexsi superlu_dist)
endif()

install(TARGETS pexsi
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
