### Source files ###
LIST(APPEND elpa_src
  src/aligned_mem.f90
  src/check_for_gpu.f90
  src/elpa1_auxiliary.f90
  src/elpa1_compute_private.f90
  src/elpa1.f90
  src/elpa2_compute.f90
  src/elpa2_determine_workload.f90
  src/elpa2.f90
  src/elpa_abstract_impl.f90
  src/elpa_api.f90
  src/elpa_autotune_impl.f90
  src/elpa_c_interface.c
  src/elpa_constants.f90
  src/elpa.f90
  src/elpa_impl.f90
  src/elpa_index.c
  src/elpa_utilities.f90
  src/ftimings.f90
  src/ftimings_type.f90
  src/ftimings_value.f90
  src/highwater_mark.c
  src/mod_blas_interfaces.f90
  src/mod_elpa_skewsymmetric_blas.f90
  src/mod_omp.f90
  src/mod_pack_unpack_cpu.f90
  src/mod_pack_unpack_gpu.f90
  src/mod_precision.f90
  src/mod_redist_band.f90
  src/mod_scalapack_interfaces.f90
  src/resident_set_size.c
  src/single_hh_trafo_real.f90
  src/mod_vendor_agnostic_layer.f90
  src/matrix_plot.f90
  src/merge_recursive.f90
  src/merge_systems.f90
  src/distribute_global_column.f90
  src/time.c
  src/virtual_memory.c
  src/global_product.f90                    
  src/global_gather.f90                     
  src/resort_ev.f90                         
  src/transform_columns.f90                 
  src/check_monotony.f90                    
  src/add_tmp.f90                           
  src/v_add_s.f90                           
  src/solve_secular_equation.f90
  src/elpa_cholesky.f90                     
  src/elpa_invert_trm.f90                   
  src/solve_tridi.f90
  src/cholesky_cuda_stub.f90 
  src/invert_trm_cuda_stub.f90
  src/mod_thread_affinity.f90
  src/elpa_pdgeqrf.f90                  
  src/elpa_pdlarfb.f90
  src/qr_utils.f90
  src/elpa_qrkernels.f90
  src/tests_variable_definitions.f90
  src/cannon.c 
  src/check_thread_affinity.c)

IF(USE_MPI_MODULE)
  LIST(APPEND elpa_src src/mod_mpi.f90)
ELSE()
  LIST(APPEND elpa_src src/mod_mpifh.f90)
ENDIF()

### Kernels ###
SET(HAVE_GENERIC_KERNELS "0")
SET(HAVE_AVX_KERNELS "0")
SET(HAVE_AVX2_KERNELS "0")
SET(HAVE_AVX512_KERNELS "0")
SET(HAVE_GPU_KERNELS "0")

IF(USE_GPU_CUDA)
  MESSAGE(STATUS "Enabling GPU (CUDA) support in ELPA")
  LIST(APPEND elpa_src
    src/cudaFunctions.cu
    src/cuUtils.cu
    src/elpa_generated_fortran_interfaces.f90
    src/ev_tridi_band_gpu_complex.cu
    src/ev_tridi_band_gpu_real.cu
    src/interface_c_kernel.f90
    src/kernels/kernels_complex.f90
    src/kernels/kernels_real.f90
    src/mod_compute_hh_trafo_cuda.f90
    src/mod_cuda.f90)

  SET(ELPA_REAL_DEFAULT "ELPA_2STAGE_REAL_GENERIC")
  SET(ELPA_COMPLEX_DEFAULT "ELPA_2STAGE_COMPLEX_GENERIC")
  SET(HAVE_GENERIC_KERNELS "1")
  SET(HAVE_GPU_KERNELS "1")
ELSE()
  LIST(APPEND elpa_src
    src/mod_cuda_stub.f90 
    src/mod_mkl_offload_stub.f90 
    src/mod_hip_stub.f90 
    src/interface_c_cuda_kernel_stub.f90 
    src/interface_c_gpu_kernel_stub.f90 
    src/interface_c_hip_kernel_stub.f90 
    src/test_gpu_vendor_agnostic_layer_stub.f90)

  IF(ELPA2_KERNEL STREQUAL "AVX512")
    MESSAGE(STATUS "Using ELPA AVX512 kernel")
    LIST(APPEND elpa_src
      src/elpa_generated_fortran_interfaces_avx512.f90
      src/kernels/complex_avx512_1hv_double_precision.c
      src/kernels/complex_avx512_1hv_single_precision.c
      src/kernels/real_avx512_2hv_double_precision.c
      src/kernels/real_avx512_2hv_single_precision.c
      src/mod_compute_hh_trafo_avx512.f90)

    SET(ELPA_REAL_DEFAULT "ELPA_2STAGE_REAL_AVX512_BLOCK2")
    SET(ELPA_COMPLEX_DEFAULT "ELPA_2STAGE_COMPLEX_AVX512_BLOCK1")
    SET(HAVE_AVX512_KERNELS "1")
  ELSEIF(ELPA2_KERNEL STREQUAL "AVX2")
    MESSAGE(STATUS "Using ELPA AVX2 kernel")
    LIST(APPEND elpa_src
      src/elpa_generated_fortran_interfaces_avx2.f90
      src/kernels/complex_avx2_1hv_double_precision.c
      src/kernels/complex_avx2_1hv_single_precision.c
      src/kernels/real_avx2_2hv_double_precision.c
      src/kernels/real_avx2_2hv_single_precision.c
      src/mod_compute_hh_trafo_avx2.f90)

    SET(ELPA_REAL_DEFAULT "ELPA_2STAGE_REAL_AVX2_BLOCK2")
    SET(ELPA_COMPLEX_DEFAULT "ELPA_2STAGE_COMPLEX_AVX2_BLOCK1")
    SET(HAVE_AVX2_KERNELS "1")
  ELSEIF(ELPA2_KERNEL STREQUAL "AVX")
    MESSAGE(STATUS "Using ELPA AVX kernel")
    LIST(APPEND elpa_src
      src/elpa_generated_fortran_interfaces_avx.f90
      src/kernels/complex_avx_1hv_double_precision.c
      src/kernels/complex_avx_1hv_single_precision.c
      src/kernels/real_avx_2hv_double_precision.c
      src/kernels/real_avx_2hv_single_precision.c
      src/mod_compute_hh_trafo_avx.f90)

    SET(ELPA_REAL_DEFAULT "ELPA_2STAGE_REAL_AVX_BLOCK2")
    SET(ELPA_COMPLEX_DEFAULT "ELPA_2STAGE_COMPLEX_AVX_BLOCK1")
    SET(HAVE_AVX_KERNELS "1")
  ELSE()
    MESSAGE(STATUS "Using ELPA generic kernel")
    LIST(APPEND elpa_src
      src/elpa_generated_fortran_interfaces.f90
      src/kernels/kernels_complex.f90
      src/kernels/kernels_real.f90
      src/mod_compute_hh_trafo.f90)

    SET(ELPA_REAL_DEFAULT "ELPA_2STAGE_REAL_GENERIC")
    SET(ELPA_COMPLEX_DEFAULT "ELPA_2STAGE_COMPLEX_GENERIC")
    SET(HAVE_GENERIC_KERNELS "1")
  ENDIF()
ENDIF()

CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/external/ELPA-2021.11.001/src/elpa/elpa_constants.h.in
  ${PROJECT_BINARY_DIR}/generated/elpa_constants.h
  @ONLY)

ADD_LIBRARY(elpa ${elpa_src})

IF(ELSI_USE_FIND_PACKAGE_MPI)
  find_package(MPI)
  TARGET_LINK_LIBRARIES(elpa PUBLIC MPI::MPI_C MPI::MPI_Fortran)
ENDIF()

TARGET_LINK_LIBRARIES(elpa PRIVATE ${LIBS})

TARGET_INCLUDE_DIRECTORIES(elpa PRIVATE
  ${INC_PATHS}
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${PROJECT_BINARY_DIR}/generated
  src)

IF(NOT CMAKE_VERSION VERSION_LESS 3.8.2)
  TARGET_COMPILE_FEATURES(elpa PRIVATE c_std_99)
ENDIF()

INSTALL(TARGETS elpa
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
