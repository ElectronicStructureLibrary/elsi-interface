elpa_dir=/home/yy244/elsi/elpa-2021.11.001/build_generic_2
for file in \
  aligned_mem.f90                       \
  check_for_gpu.f90                     \
  elpa1_auxiliary.f90                   \
  elpa1_compute_private.f90             \
  elpa1.f90                             \
  elpa2_compute.f90                     \
  elpa2_determine_workload.f90          \
  elpa2.f90                             \
  elpa_abstract_impl.f90                \
  elpa_api.f90                          \
  elpa_autotune_impl.f90                \
  elpa_constants.f90                    \
  elpa.f90                              \
  elpa_impl.f90                         \
  elpa_utilities.f90                    \
  ftimings.f90                          \
  ftimings_type.f90                     \
  ftimings_value.f90                    \
  mod_blas_interfaces.f90               \
  mod_elpa_skewsymmetric_blas.f90       \
  mod_omp.f90                           \
  mod_pack_unpack_cpu.f90               \
  mod_pack_unpack_gpu.f90               \
  mod_precision.f90                     \
  mod_redist_band.f90                   \
  mod_scalapack_interfaces.f90          \
  mod_vendor_agnostic_layer.f90         \
  matrix_plot.f90                       \
  merge_recursive.f90                   \
  merge_systems.f90                     \
  distribute_global_column.f90          \
  single_hh_trafo_real.f90              \
  global_product.f90                    \
  global_gather.f90                     \
  resort_ev.f90                         \
  transform_columns.f90                 \
  check_monotony.f90                    \
  add_tmp.f90                           \
  v_add_s.f90                           \
  solve_secular_equation.f90            \
  elpa_cholesky.f90                     \
  elpa_invert_trm.f90                   \
  solve_tridi.f90
do
prefix=${file%.*}
echo $prefix
origin_file=`ls $elpa_dir/*$prefix.F90* | tail -n 1`
cp $origin_file $file
done 

elpa_dir=/home/yy244/elsi/elpa-2021.11.001/build_generic_2
for file in \
  mod_cuda_stub.f90 \
  mod_mkl_offload_stub.f90 \
  mod_hip_stub.f90 \
  interface_c_cuda_kernel_stub.f90 \
  interface_c_gpu_kernel_stub.f90 \
  interface_c_hip_kernel_stub.f90 \
  test_gpu_vendor_agnostic_layer_stub.f90 
do
prefix=${file%_stub.*}
echo $prefix
origin_file=`ls $elpa_dir/*$prefix.F90* | tail -n 1`
#echo $origin_file
cp $origin_file $file
done 

elpa_dir=/home/yy244/elsi/elpa-2021.11.001/src
for file in \
  elpa_c_interface.c \
  elpa_index.c 
do
echo $file
origin_file=`ls $elpa_dir/$file | tail -n 1`
#echo $origin_file
cp $origin_file $file
done 

elpa_dir=/home/yy244/elsi/elpa-2021.11.001/src/ftimings
for file in \
  highwater_mark.c \
  resident_set_size.c \
  time.c \
  virtual_memory.c
do
echo $file
origin_file=`ls $elpa_dir/$file | tail -n 1`
#echo $origin_file
cp $origin_file $file
done 

elpa_dir=/home/yy244/elsi/elpa-2021.11.001/elpa
cp $elpa_dir/elpa_constants.h.in elpa

cp ../../ELPA/src/mod_mpi.f90 .
cp ../../ELPA/src/mod_mpifh.f90 .

elpa_dir=/home/yy244/elsi/elpa-2021.11.001/build_generic_2
elpa_dir_avx=/home/yy244/elsi/elpa-2021.11.001/build_AVX_2
elpa_dir_avx2=/home/yy244/elsi/elpa-2021.11.001/build_AVX2_2
elpa_dir_avx512=/home/yy244/elsi/elpa-2021.11.001/build_AVX512_2

origin_file=`ls $elpa_dir/*mod_compute_hh_trafo* | tail -n 1`
cp $origin_file mod_compute_hh_trafo.f90

origin_file=`ls $elpa_dir_avx/*mod_compute_hh_trafo* | tail -n 1`
cp $origin_file mod_compute_hh_trafo_avx.f90

origin_file=`ls $elpa_dir_avx2/*mod_compute_hh_trafo* | tail -n 1`
cp $origin_file mod_compute_hh_trafo_avx2.f90

origin_file=`ls $elpa_dir_avx512/*mod_compute_hh_trafo* | tail -n 1`
cp $origin_file mod_compute_hh_trafo_avx512.f90

origin_file=`ls $elpa_dir/*elpa_generated_fortran_interfaces* | tail -n 1`
cp $origin_file elpa_generated_fortran_interfaces.f90

origin_file=`ls $elpa_dir_avx/*elpa_generated_fortran_interfaces* | tail -n 1`
cp $origin_file elpa_generated_fortran_interfaces_avx.f90

origin_file=`ls $elpa_dir_avx2/*elpa_generated_fortran_interfaces* | tail -n 1`
cp $origin_file elpa_generated_fortran_interfaces_avx2.f90

origin_file=`ls $elpa_dir_avx512/*elpa_generated_fortran_interfaces* | tail -n 1`
cp $origin_file elpa_generated_fortran_interfaces_avx512.f90

elpa_dir_elpa=/home/yy244/elsi/elpa-2021.11.001/build_generic_2/elpa
cp -r $elpa_dir_elpa .

mkdir kernels

origin_file=`ls $elpa_dir/*kernels_complex* | tail -n 1`
cp $origin_file kernels/kernels_complex.f90

origin_file=`ls $elpa_dir/*kernels_real* | tail -n 1`
cp $origin_file kernels/kernels_real.f90

cp $elpa_dir_avx/src/elpa2/kernels/*.c kernels
cp $elpa_dir_avx2/src/elpa2/kernels/*.c kernels
cp $elpa_dir_avx512/src/elpa2/kernels/*.c kernels
