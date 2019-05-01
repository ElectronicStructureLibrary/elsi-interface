### Source files ###
LIST(APPEND elpa_src
  src/aligned_mem.f90
  src/check_pd.f90
  src/elpa1_auxiliary.f90
  src/elpa1_compute_private.f90
  src/elpa1.f90
  src/elpa2_compute.f90
  src/elpa2.f90
  src/elpa_utilities.f90
  src/mod_mpi.f90
  src/mod_pack_unpack_complex.f90
  src/mod_pack_unpack_real.f90
  src/mod_precision.f90
  src/mod_single_hh_trafo_real.f90)

### Kernels ###
IF(ELPA2_KERNEL STREQUAL "AVX512")
  MESSAGE(STATUS "Using ELPA AVX512 kernel")
  LIST(APPEND elpa_src
    src/mod_compute_hh_trafo_real_avx512.f90
    src/mod_compute_hh_trafo_complex_avx512.f90
    src/mod_fortran_interfaces_avx512.f90
    src/elpa2_utilities_avx512.f90
    src/kernel_real_avx512_2hv.c
    src/kernel_complex_avx512_1hv.c)
ELSEIF(ELPA2_KERNEL STREQUAL "AVX2")
  MESSAGE(STATUS "Using ELPA AVX2 kernel")
  LIST(APPEND elpa_src
    src/mod_compute_hh_trafo_real_avx2.f90
    src/mod_compute_hh_trafo_complex_avx2.f90
    src/mod_fortran_interfaces_avx2.f90
    src/elpa2_utilities_avx2.f90
    src/kernel_real_avx2_2hv.c
    src/kernel_complex_avx2_1hv.c)
ELSEIF(ELPA2_KERNEL STREQUAL "AVX")
  MESSAGE(STATUS "Using ELPA AVX kernel")
  LIST(APPEND elpa_src
    src/mod_compute_hh_trafo_real_avx.f90
    src/mod_compute_hh_trafo_complex_avx.f90
    src/mod_fortran_interfaces_avx.f90
    src/elpa2_utilities_avx.f90
    src/kernel_real_avx_2hv.c
    src/kernel_complex_avx_1hv.c)
ELSEIF(ELPA2_KERNEL STREQUAL "BGQ")
  MESSAGE(STATUS "Using ELPA Blue Gene Q kernel")
  LIST(APPEND elpa_src
    src/mod_compute_hh_trafo_real_bgq.f90
    src/mod_compute_hh_trafo_complex.f90
    src/elpa2_utilities_bgq.f90
    src/kernel_real_bgq.f90
    src/kernel_complex.f90)
ELSE()
  MESSAGE(STATUS "Using ELPA generic kernel")
  LIST(APPEND elpa_src
    src/mod_compute_hh_trafo_real.f90
    src/mod_compute_hh_trafo_complex.f90
    src/elpa2_utilities.f90
    src/kernel_real.f90
    src/kernel_complex.f90)
ENDIF()

ADD_LIBRARY(elpa ${elpa_src})

TARGET_LINK_LIBRARIES(elpa PRIVATE ${LIBS})

TARGET_INCLUDE_DIRECTORIES(elpa PRIVATE
  ${INC_PATHS}
  ${CMAKE_Fortran_MODULE_DIRECTORY})

IF(ELPA2_KERNEL STREQUAL "AVX512" OR ELPA2_KERNEL STREQUAL "AVX2" OR ELPA2_KERNEL STREQUAL "AVX")
  IF(NOT CMAKE_VERSION VERSION_LESS 3.8.2)
    TARGET_COMPILE_FEATURES(elpa PRIVATE c_std_99)
  ENDIF()
ENDIF()

INSTALL(TARGETS elpa
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
