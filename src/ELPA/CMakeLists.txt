list(APPEND elpa_src
  ${elpa_dir}/src/aligned_mem.f90
  ${elpa_dir}/src/check_for_gpu.f90
  ${elpa_dir}/src/check_singularity.f90
  ${elpa_dir}/src/elpa1_auxiliary.f90
  ${elpa_dir}/src/elpa1_compute_private.f90
  ${elpa_dir}/src/elpa1.f90
  ${elpa_dir}/src/elpa1_utilities.f90
  ${elpa_dir}/src/elpa2_compute.f90
  ${elpa_dir}/src/elpa2.f90
  ${elpa_dir}/src/elpa_pdgeqrf.f90
  ${elpa_dir}/src/elpa_pdlarfb.f90
  ${elpa_dir}/src/elpa_qrkernels.f90
  ${elpa_dir}/src/elpa_utilities.f90
  ${elpa_dir}/src/interface_c_kernel.f90
  ${elpa_dir}/src/mod_cuda.f90
  ${elpa_dir}/src/mod_mpi.f90
  ${elpa_dir}/src/mod_pack_unpack_complex.f90
  ${elpa_dir}/src/mod_pack_unpack_real_gpu.f90
  ${elpa_dir}/src/mod_pack_unpack_real.f90
  ${elpa_dir}/src/mod_precision.f90
  ${elpa_dir}/src/mod_single_hh_trafo_real.f90
  ${elpa_dir}/src/qr_utils.f90
  ${elpa_dir}/src/timer_dummy.f90
  )

if(ENABLE_ELPA_AVX512)
  list(APPEND elpa_src
    ${elpa_dir}/src/mod_compute_hh_trafo_real_avx512.f90
    ${elpa_dir}/src/mod_compute_hh_trafo_complex_avx512.f90
    ${elpa_dir}/src/mod_fortran_interfaces_avx512.f90
    ${elpa_dir}/src/elpa2_utilities_avx512.f90
    ${elpa_dir}/src/elpa2_kernels_real_avx512_2hv.c
    ${elpa_dir}/src/elpa2_kernels_complex_avx512_1hv.c
    )
elseif(ENABLE_ELPA_AVX2)
  list(APPEND elpa_src
    ${elpa_dir}/src/mod_compute_hh_trafo_real_avx2.f90
    ${elpa_dir}/src/mod_compute_hh_trafo_complex_avx2.f90
    ${elpa_dir}/src/mod_fortran_interfaces_avx2.f90
    ${elpa_dir}/src/elpa2_utilities_avx2.f90
    ${elpa_dir}/src/elpa2_kernels_real_avx2_2hv.c
    ${elpa_dir}/src/elpa2_kernels_complex_avx2_1hv.c
    )
elseif(ENABLE_ELPA_AVX)
  list(APPEND elpa_src
    ${elpa_dir}/src/mod_compute_hh_trafo_real_avx.f90
    ${elpa_dir}/src/mod_compute_hh_trafo_complex_avx.f90
    ${elpa_dir}/src/mod_fortran_interfaces_avx.f90
    ${elpa_dir}/src/elpa2_utilities_avx.f90
    ${elpa_dir}/src/elpa2_kernels_real_avx_2hv.c
    ${elpa_dir}/src/elpa2_kernels_complex_avx_1hv.c
    )
elseif(ENABLE_ELPA_BGQ)
  list(APPEND elpa_src
    ${elpa_dir}/src/mod_compute_hh_trafo_bgq.f90
    ${elpa_dir}/src/mod_compute_hh_trafo_complex.f90
    ${elpa_dir}/src/mod_fortran_interfaces.f90
    ${elpa_dir}/src/elpa2_utilities_bgq.f90
    ${elpa_dir}/src/elpa2_kernels_real_bgq.f90
    ${elpa_dir}/src/elpa2_kernels_complex.f90
    )
else()
  list(APPEND elpa_src
    ${elpa_dir}/src/mod_compute_hh_trafo_real.f90
    ${elpa_dir}/src/mod_compute_hh_trafo_complex.f90
    ${elpa_dir}/src/mod_fortran_interfaces.f90
    ${elpa_dir}/src/elpa2_utilities.f90
    ${elpa_dir}/src/elpa2_kernels_real.f90
    ${elpa_dir}/src/elpa2_kernels_complex.f90
    )
endif()

add_library(elpa ${elpa_src})

install(TARGETS elpa DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
