### Source files ###
LIST(APPEND elpa_src
  src/aligned_mem.f90
  src/check_for_gpu.f90
  src/check_singularity.f90
  src/elpa1_auxiliary.f90
  src/elpa1_compute_private.f90
  src/elpa1.f90
  src/elpa1_utilities.f90
  src/elpa2_compute.f90
  src/elpa2.f90
  src/elpa_pdgeqrf.f90
  src/elpa_pdlarfb.f90
  src/elpa_qrkernels.f90
  src/elpa_utilities.f90
  src/interface_c_kernel.f90
  src/mod_cuda.f90
  src/mod_mpi.f90
  src/mod_pack_unpack_complex.f90
  src/mod_pack_unpack_real_gpu.f90
  src/mod_pack_unpack_real.f90
  src/mod_precision.f90
  src/mod_single_hh_trafo_real.f90
  src/qr_utils.f90
  src/timer_dummy.f90
)

### Kernels ###
IF(ENABLE_ELPA_AVX512)
  MESSAGE(STATUS "Enabling ELPA AVX 512 kernel")
  LIST(APPEND elpa_src
    src/mod_compute_hh_trafo_real_avx512.f90
    src/mod_compute_hh_trafo_complex_avx512.f90
    src/mod_fortran_interfaces_avx512.f90
    src/elpa2_utilities_avx512.f90
    src/elpa2_kernels_real_avx512_2hv.c
    src/elpa2_kernels_complex_avx512_1hv.c
  )
ELSEIF(ENABLE_ELPA_AVX2)
  MESSAGE(STATUS "Enabling ELPA AVX 2 kernel")
  LIST(APPEND elpa_src
    src/mod_compute_hh_trafo_real_avx2.f90
    src/mod_compute_hh_trafo_complex_avx2.f90
    src/mod_fortran_interfaces_avx2.f90
    src/elpa2_utilities_avx2.f90
    src/elpa2_kernels_real_avx2_2hv.c
    src/elpa2_kernels_complex_avx2_1hv.c
  )
ELSEIF(ENABLE_ELPA_AVX)
  MESSAGE(STATUS "Enabling ELPA AVX kernel")
  LIST(APPEND elpa_src
    src/mod_compute_hh_trafo_real_avx.f90
    src/mod_compute_hh_trafo_complex_avx.f90
    src/mod_fortran_interfaces_avx.f90
    src/elpa2_utilities_avx.f90
    src/elpa2_kernels_real_avx_2hv.c
    src/elpa2_kernels_complex_avx_1hv.c
  )
ELSEIF(ENABLE_ELPA_BGQ)
  MESSAGE(STATUS "Enabling ELPA Blue Gene Q kernel")
  LIST(APPEND elpa_src
    src/mod_compute_hh_trafo_bgq.f90
    src/mod_compute_hh_trafo_complex.f90
    src/mod_fortran_interfaces.f90
    src/elpa2_utilities_bgq.f90
    src/elpa2_kernels_real_bgq.f90
    src/elpa2_kernels_complex.f90
  )
ELSE()
  MESSAGE(STATUS "Using ELPA generic kernel")
  LIST(APPEND elpa_src
    src/mod_compute_hh_trafo_real.f90
    src/mod_compute_hh_trafo_complex.f90
    src/mod_fortran_interfaces.f90
    src/elpa2_utilities.f90
    src/elpa2_kernels_real.f90
    src/elpa2_kernels_complex.f90
  )
ENDIF()

ADD_LIBRARY(elpa ${elpa_src})

INSTALL(TARGETS elpa
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

SET(ELPA_LIB "${CMAKE_BINARY_DIR}/lib/libelpa.a")
