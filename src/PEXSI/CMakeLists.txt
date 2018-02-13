enable_language(CXX)
include_directories(${pexsi_dir}/include)

list(APPEND pexsi_src
  ${pexsi_dir}/src/TreeBcast.cpp
  ${pexsi_dir}/src/blas.cpp
  ${pexsi_dir}/src/f_interface.f90
  ${pexsi_dir}/src/getPole.cpp
  ${pexsi_dir}/src/get_perm_c_parmetis.c
  ${pexsi_dir}/src/global.cpp
  ${pexsi_dir}/src/interface.cpp
  ${pexsi_dir}/src/lapack.cpp
  ${pexsi_dir}/src/mpi_interf.cpp
  ${pexsi_dir}/src/pdsymbfact.c
  ${pexsi_dir}/src/pole.cpp
  ${pexsi_dir}/src/ppexsi.cpp
  ${pexsi_dir}/src/pzsymbfact.c
  ${pexsi_dir}/src/superlu_dist_internal_complex.cpp
  ${pexsi_dir}/src/superlu_dist_internal_real.cpp
  ${pexsi_dir}/src/timer.cpp
  ${pexsi_dir}/src/utility.cpp
  )

add_library(pexsi ${pexsi_src})

install(TARGETS pexsi DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
