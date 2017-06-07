! Copyright (c) 2015-2017, the ELSI team. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!  * Neither the name of the "ELectronic Structure Infrastructure" project nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!> 
!! This module contains the C interfaces of ELSI.
!!
module ELSI_C_INTERFACE

   use, intrinsic :: iso_c_binding
   use ELSI

   implicit none

   private

   integer(kind=c_int) :: n_g_size_c
   integer(kind=c_int) :: n_l_rows_c
   integer(kind=c_int) :: n_l_cols_c
   integer(kind=c_int) :: nnz_l_pexsi_c

contains

subroutine elsi_init_c_wrapper(handle_c,solver,parallel_mode,matrix_storage_format,&
                               matrix_size,n_electrons,n_states)&
                               bind(C,name="c_elsi_init")

   implicit none

   type(c_ptr)                            :: handle_c
   integer(kind=c_int), value, intent(in) :: solver
   integer(kind=c_int), value, intent(in) :: parallel_mode
   integer(kind=c_int), value, intent(in) :: matrix_storage_format
   integer(kind=c_int), value, intent(in) :: matrix_size
   real(kind=c_double), value, intent(in) :: n_electrons
   integer(kind=c_int), value, intent(in) :: n_states

   type(elsi_handle), pointer :: handle_f

   allocate(handle_f)

   call elsi_init(handle_f,solver,parallel_mode,matrix_storage_format,&
                  matrix_size,n_electrons,n_states)

   handle_c = c_loc(handle_f)

   n_g_size_c = matrix_size

   if(parallel_mode == 0) then
      n_l_rows_c = handle_f%n_l_rows
      n_l_cols_c = handle_f%n_l_cols
   endif

end subroutine

subroutine elsi_set_mpi_c_wrapper(handle_c,mpi_comm)&
                                  bind(C,name="c_elsi_set_mpi")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: mpi_comm

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mpi(handle_f,mpi_comm)

end subroutine

subroutine elsi_set_blacs_c_wrapper(handle_c,blacs_ctxt,block_size)&
                                    bind(C,name="c_elsi_set_blacs")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: blacs_ctxt
   integer(kind=c_int), value, intent(in) :: block_size

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_blacs(handle_f,blacs_ctxt,block_size)

   n_l_rows_c = handle_f%n_l_rows
   n_l_cols_c = handle_f%n_l_cols

end subroutine

subroutine elsi_set_csc_c_wrapper(handle_c,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)&
                                  bind(C,name="c_elsi_set_csc")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: nnz_g
   integer(kind=c_int), value, intent(in) :: nnz_l
   integer(kind=c_int), value, intent(in) :: n_l_cols
   integer(kind=c_int),        intent(in) :: row_ind(nnz_l)
   integer(kind=c_int),        intent(in) :: col_ptr(n_l_cols+1)

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_csc(handle_f,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)

   nnz_l_pexsi_c = nnz_l

end subroutine

subroutine elsi_finalize_c_wrapper(handle_c) bind(C,name="c_elsi_finalize")

   implicit none

   type(c_ptr), value :: handle_c

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_finalize(handle_f)

end subroutine

subroutine elsi_customize_c_wrapper(handle_c,print_detail,overlap_is_unit,&
                                    zero_threshold,no_singularity_check,&
                                    singularity_tolerance,stop_singularity,&
                                    uplo) bind(C,name="c_elsi_customize")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: print_detail
   integer(kind=c_int), value, intent(in) :: overlap_is_unit
   real(kind=c_double), value, intent(in) :: zero_threshold
   integer(kind=c_int), value, intent(in) :: no_singularity_check
   real(kind=c_double), value, intent(in) :: singularity_tolerance
   integer(kind=c_int), value, intent(in) :: stop_singularity
   integer(kind=c_int), value, intent(in) :: uplo

   type(elsi_handle), pointer :: handle_f
   logical :: print_detail_f
   logical :: overlap_is_unit_f
   logical :: no_singularity_check_f
   logical :: stop_singularity_f

   if(print_detail == 0) then
      print_detail_f = .false.
   else
      print_detail_f = .true.
   endif

   if(overlap_is_unit == 0) then
      overlap_is_unit_f = .false.
   else
      overlap_is_unit_f = .true.
   endif

   if(no_singularity_check == 0) then
      no_singularity_check_f = .false.
   else
      no_singularity_check_f = .true.
   endif

   if(stop_singularity == 0) then
      stop_singularity_f = .false.
   else
      stop_singularity_f = .true.
   endif

   call c_f_pointer(handle_c,handle_f)

   call elsi_customize(handle_f,print_detail_f,overlap_is_unit_f,zero_threshold,&
                       no_singularity_check_f,singularity_tolerance,&
                       stop_singularity_f,uplo)

end subroutine

subroutine elsi_customize_mu_c_wrapper(handle_c,broadening_scheme,broadening_width,&
                                       occ_accuracy,mu_max_steps,spin_degeneracy)&
                                       bind(C,name="c_elsi_customize_mu")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: broadening_scheme
   real(kind=c_double), value, intent(in) :: broadening_width
   real(kind=c_double), value, intent(in) :: occ_accuracy
   integer(kind=c_int), value, intent(in) :: mu_max_steps
   real(kind=c_double), value, intent(in) :: spin_degeneracy

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_customize_mu(handle_f,broadening_scheme,broadening_width,&
                          occ_accuracy,mu_max_steps,spin_degeneracy)

end subroutine

subroutine elsi_customize_omm_c_wrapper(handle_c,n_elpa_steps,omm_flavor,eigen_shift,&
                                        omm_tolerance,use_pspblas,omm_output)&
                                        bind(C,name="c_elsi_customize_omm")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: n_elpa_steps
   integer(kind=c_int), value, intent(in) :: omm_flavor
   real(kind=c_double), value, intent(in) :: eigen_shift
   real(kind=c_double), value, intent(in) :: omm_tolerance
   integer(kind=c_int), value, intent(in) :: use_pspblas
   integer(kind=c_int), value, intent(in) :: omm_output

   type(elsi_handle), pointer :: handle_f
   logical :: use_pspblas_f
   logical :: omm_output_f

   if(use_pspblas == 0) then
      use_pspblas_f = .false.
   else
      use_pspblas_f = .true.
   endif

   if(omm_output == 0) then
      omm_output_f = .false.
   else
      omm_output_f = .true.
   endif

   call c_f_pointer(handle_c,handle_f)

   call elsi_customize_omm(handle_f,n_elpa_steps,omm_flavor,eigen_shift,&
                           omm_tolerance,use_pspblas_f,omm_output_f)

end subroutine

subroutine elsi_customize_pexsi_c_wrapper(handle_c,temperature,gap,delta_e,n_poles,&
                                          n_procs_per_pole,max_iteration,mu_min,&
                                          mu_max,mu0,mu_inertia_tolerance,&
                                          mu_inertia_expansion,mu_safeguard,&
                                          n_electron_accuracy,matrix_type,&
                                          is_symbolic_factorize,ordering,&
                                          np_symbolic_factorize,verbosity)&
                                          bind(C,name="c_elsi_customize_pexsi")

   implicit none

   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: temperature
   real(kind=c_double), value, intent(in) :: gap
   real(kind=c_double), value, intent(in) :: delta_e
   integer(kind=c_int), value, intent(in) :: n_poles
   integer(kind=c_int), value, intent(in) :: n_procs_per_pole
   integer(kind=c_int), value, intent(in) :: max_iteration
   real(kind=c_double), value, intent(in) :: mu_min
   real(kind=c_double), value, intent(in) :: mu_max
   real(kind=c_double), value, intent(in) :: mu0
   real(kind=c_double), value, intent(in) :: mu_inertia_tolerance
   real(kind=c_double), value, intent(in) :: mu_inertia_expansion
   real(kind=c_double), value, intent(in) :: mu_safeguard
   real(kind=c_double), value, intent(in) :: n_electron_accuracy
   integer(kind=c_int), value, intent(in) :: matrix_type
   integer(kind=c_int), value, intent(in) :: is_symbolic_factorize
   integer(kind=c_int), value, intent(in) :: ordering
   integer(kind=c_int), value, intent(in) :: np_symbolic_factorize
   integer(kind=c_int), value, intent(in) :: verbosity

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_customize_pexsi(handle_f,temperature,gap,delta_e,n_poles,n_procs_per_pole,&
           max_iteration,mu_min,mu_max,mu0,mu_inertia_tolerance,mu_inertia_expansion,&
           mu_safeguard,n_electron_accuracy,matrix_type,is_symbolic_factorize,&
           ordering,np_symbolic_factorize,verbosity)

end subroutine

subroutine elsi_customize_elpa_c_wrapper(handle_c,elpa_solver,elpa_output)&
                                         bind(C,name="c_elsi_customize_elpa")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: elpa_solver
   integer(kind=c_int), value, intent(in) :: elpa_output

   type(elsi_handle), pointer :: handle_f
   logical :: elpa_output_f

   if(elpa_output == 0) then
      elpa_output_f = .false.
   else
      elpa_output_f = .true.
   endif


   call c_f_pointer(handle_c,handle_f)

   call elsi_customize_elpa(handle_f,elpa_solver,elpa_output_f)

end subroutine

subroutine elsi_ev_real_c_wrapper(handle_c,H,S,e_val,e_vec) bind(C,name="c_elsi_ev_real")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: H(n_l_rows_c,n_l_cols_c)
   real(kind=c_double) :: S(n_l_rows_c,n_l_cols_c)
   real(kind=c_double) :: e_val(n_g_size_c)
   real(kind=c_double) :: e_vec(n_l_rows_c,n_l_cols_c)

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_ev_real(handle_f,H,S,e_val,e_vec)

end subroutine

subroutine elsi_ev_complex_c_wrapper(handle_c,H,S,e_val,e_vec) bind(C,name="c_elsi_ev_complex")

   implicit none

   type(c_ptr), value             :: handle_c
   complex(kind=c_double_complex) :: H(n_l_rows_c,n_l_cols_c)
   complex(kind=c_double_complex) :: S(n_l_rows_c,n_l_cols_c)
   real(kind=c_double)            :: e_val(n_g_size_c)
   complex(kind=c_double_complex) :: e_vec(n_l_rows_c,n_l_cols_c)

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_ev_complex(handle_f,H,S,e_val,e_vec)

end subroutine

subroutine elsi_ev_real_sparse_c_wrapper(handle_c,H,S,e_val,e_vec) bind(C,name="c_elsi_ev_real_sparse")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: H(nnz_l_pexsi_c)
   real(kind=c_double) :: S(nnz_l_pexsi_c)
   real(kind=c_double) :: e_val(n_g_size_c)
   real(kind=c_double) :: e_vec(n_l_rows_c,n_l_cols_c)

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_ev_real_sparse(handle_f,H,S,e_val,e_vec)

end subroutine

function elsi_dm_real_c_wrapper(handle_c,H,S,D) result(energy) bind(C,name="c_elsi_dm_real")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: H(n_l_rows_c,n_l_cols_c)
   real(kind=c_double) :: S(n_l_rows_c,n_l_cols_c)
   real(kind=c_double) :: D(n_l_rows_c,n_l_cols_c)
   real(kind=c_double) :: energy

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_dm_real(handle_f,H,S,D,energy)

end function

function elsi_dm_real_sparse_c_wrapper(handle_c,H,S,D) result(energy) bind(C,name="c_elsi_dm_real_sparse")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: H(nnz_l_pexsi_c)
   real(kind=c_double) :: S(nnz_l_pexsi_c)
   real(kind=c_double) :: D(nnz_l_pexsi_c)
   real(kind=c_double) :: energy

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_dm_real_sparse(handle_f,H,S,D,energy)

end function

subroutine elsi_collect_pexsi_c_wrapper(handle_c,mu,edm,fdm) bind(C,name="c_elsi_collect_pexsi")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: mu
   real(kind=c_double) :: edm(nnz_l_pexsi_c)
   real(kind=c_double) :: fdm(nnz_l_pexsi_c)

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_collect_pexsi(handle_f,mu,edm,fdm)

end subroutine

end module ELSI_C_INTERFACE
