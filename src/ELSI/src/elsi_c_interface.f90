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
!! This file contains the C interfaces of ELSI.
!!

subroutine elsi_init_c_wrapper(solver,parallel_mode,matrix_format,matrix_size,&
                               n_electrons_in,n_states_in)&
                               bind(C,name="c_elsi_init")

   use, intrinsic :: iso_c_binding
   use ELSI

   implicit none

   integer(kind=c_int), value, intent(in) :: solver
   integer(kind=c_int), value, intent(in) :: parallel_mode
   integer(kind=c_int), value, intent(in) :: matrix_format
   integer(kind=c_int), value, intent(in) :: matrix_size
   real(kind=c_double), value, intent(in) :: n_electrons_in
   integer(kind=c_int), value, intent(in) :: n_states_in

   call elsi_init(solver,parallel_mode,matrix_format,matrix_size,&
                  n_electrons_in,n_states_in)

end subroutine

subroutine elsi_set_mpi_c_wrapper(mpi_comm_global_in) bind(C,name="c_elsi_set_mpi")

   use, intrinsic :: iso_c_binding
   use ELSI

   implicit none

   integer(kind=c_int), value, intent(in) :: mpi_comm_global_in

   call elsi_set_mpi(mpi_comm_global_in)

end subroutine

subroutine elsi_set_blacs_c_wrapper(blacs_ctxt_in,n_b_rows_in,n_b_cols_in,&
              n_p_rows_in,n_p_cols_in) bind(C,name="c_elsi_set_blacs")

   use, intrinsic :: iso_c_binding
   use ELSI

   implicit none

   integer(kind=c_int), value, intent(in) :: blacs_ctxt_in
   integer(kind=c_int), value, intent(in) :: n_b_rows_in
   integer(kind=c_int), value, intent(in) :: n_b_cols_in
   integer(kind=c_int), value, intent(in) :: n_p_rows_in
   integer(kind=c_int), value, intent(in) :: n_p_cols_in

   call elsi_set_blacs(blacs_ctxt_in,n_b_rows_in,n_b_cols_in,&
                       n_p_rows_in,n_p_cols_in)

end subroutine

subroutine elsi_set_sparsity_c_wrapper(nnz_g_in,nnz_l_in,nnz_l_cols_in,row_ind_in,&
                                       col_ptr_in) bind(C,name="c_elsi_set_sparsity")

   use, intrinsic :: iso_c_binding
   use ELSI

   implicit none

   integer(kind=c_int), value, intent(in) :: nnz_g_in
   integer(kind=c_int), value, intent(in) :: nnz_l_in
   integer(kind=c_int), value, intent(in) :: nnz_l_cols_in
   integer(kind=c_int), intent(in)        :: row_ind_in(*)
   integer(kind=c_int), intent(in)        :: col_ptr_in(*)

   call elsi_set_sparsity(nnz_g_in,nnz_l_in,nnz_l_cols_in,row_ind_in,col_ptr_in)

end subroutine

subroutine elsi_finalize_c_wrapper() bind(C,name="c_elsi_finalize")

   use, intrinsic :: iso_c_binding
   use ELSI

   implicit none

   call elsi_finalize()

end subroutine

subroutine elsi_customize_c_wrapper(print_detail,unit_overlap,hartree_to_ev,&
                                    numerical_zero,no_check_singularity,&
                                    singularity_threshold,force_stop_singularity,&
                                    broadening_scheme,broadening_width,&
                                    mu_accuracy,mu_max_steps)&
                                    bind(C,name="c_elsi_customize")

   use, intrinsic :: iso_c_binding
   use ELSI

   implicit none

   integer(kind=c_int), value, intent(in) :: print_detail
   integer(kind=c_int), value, intent(in) :: unit_overlap
   real(kind=c_double), value, intent(in) :: hartree_to_ev
   real(kind=c_double), value, intent(in) :: numerical_zero
   integer(kind=c_int), value, intent(in) :: no_check_singularity
   real(kind=c_double), value, intent(in) :: singularity_threshold
   integer(kind=c_int), value, intent(in) :: force_stop_singularity
   integer(kind=c_int), value, intent(in) :: broadening_scheme
   real(kind=c_double), value, intent(in) :: broadening_width
   real(kind=c_double), value, intent(in) :: mu_accuracy
   integer(kind=c_int), value, intent(in) :: mu_max_steps

   logical :: print_detail_f
   logical :: unit_overlap_f
   logical :: no_check_singularity_f
   logical :: force_stop_singularity_f

   if(print_detail == 0) then
      print_detail_f = .false.
   else
      print_detail_f = .true.
   endif

   if(unit_overlap == 0) then
      unit_overlap_f = .false.
   else
      unit_overlap_f = .true.
   endif

   if(no_check_singularity == 0) then
      no_check_singularity_f = .false.
   else
      no_check_singularity_f = .true.
   endif

   if(force_stop_singularity == 0) then
      force_stop_singularity_f = .false.
   else
      force_stop_singularity_f = .true.
   endif

   call elsi_customize(print_detail_f,unit_overlap_f,hartree_to_ev,numerical_zero,&
                       no_check_singularity_f,singularity_threshold,&
                       force_stop_singularity_f,broadening_scheme,broadening_width,&
                       mu_accuracy,mu_max_steps)

end subroutine

subroutine elsi_customize_omm_c_wrapper(n_elpa_steps_omm,eigenspectrum_shift,&
                                        omm_tolerance,use_pspblas)&
                                        bind(C,name="c_elsi_customize_omm")

   use, intrinsic :: iso_c_binding
   use ELSI

   implicit none

   integer(kind=c_int), value, intent(in) :: n_elpa_steps_omm
   real(kind=c_double), value, intent(in) :: eigenspectrum_shift
   real(kind=c_double), value, intent(in) :: omm_tolerance
   integer(kind=c_int), value, intent(in) :: use_pspblas

   logical :: use_pspblas_f

   if(use_pspblas == 0) then
      use_pspblas_f = .false.
   else
      use_pspblas_f = .true.
   endif

   call elsi_customize_omm(n_elpa_steps_omm,eigenspectrum_shift,&
                           omm_tolerance,use_pspblas_f)

end subroutine

subroutine elsi_customize_pexsi_c_wrapper(temperature,gap,delta_E,n_poles,&
                                          max_iteration,mu_min,mu_max,mu0,&
                                          mu_inertia_tolerance,mu_inertia_expansion,&
                                          mu_safeguard,n_electron_accuracy,&
                                          matrix_type,is_symbolic_factorize,&
                                          ordering,np_symbolic_factorize,verbosity)&
                                          bind(C,name="c_elsi_customize_pexsi")

   use, intrinsic :: iso_c_binding
   use ELSI

   implicit none

   real(c_double), value, intent(in) :: temperature
   real(c_double), value, intent(in) :: gap
   real(c_double), value, intent(in) :: delta_E
   integer(c_int), value, intent(in) :: n_poles
   integer(c_int), value, intent(in) :: max_iteration
   real(c_double), value, intent(in) :: mu_min
   real(c_double), value, intent(in) :: mu_max
   real(c_double), value, intent(in) :: mu0
   real(c_double), value, intent(in) :: mu_inertia_tolerance
   real(c_double), value, intent(in) :: mu_inertia_expansion
   real(c_double), value, intent(in) :: mu_safeguard
   real(c_double), value, intent(in) :: n_electron_accuracy
   integer(c_int), value, intent(in) :: matrix_type
   integer(c_int), value, intent(in) :: is_symbolic_factorize
   integer(c_int), value, intent(in) :: ordering
   integer(c_int), value, intent(in) :: np_symbolic_factorize
   integer(c_int), value, intent(in) :: verbosity

   call elsi_customize_pexsi(temperature,gap,delta_E,n_poles,max_iteration,mu_min,&
                             mu_max,mu0,mu_inertia_tolerance,mu_inertia_expansion,&
                             mu_safeguard,n_electron_accuracy,matrix_type,&
                             is_symbolic_factorize,ordering,np_symbolic_factorize,&
                             verbosity)

end subroutine

subroutine elsi_customize_elpa_c_wrapper(elpa_solver)&
                                         bind(C,name="c_elsi_customize_elpa")

   use, intrinsic :: iso_c_binding
   use ELSI

   implicit none

   integer(c_int), value, intent(in) :: elpa_solver

   call elsi_customize_elpa(elpa_solver)

end subroutine

subroutine elsi_ev_real_c_wrapper(H_in,S_in,e_val_out,e_vec_out)&
                                  bind(C,name="c_elsi_ev_real")

   use, intrinsic :: iso_c_binding
   use ELSI
   use ELSI_DIMENSIONS

   implicit none

   real(kind=c_double), intent(in)  :: H_in(n_l_rows,*)
   real(kind=c_double), intent(in)  :: S_in(n_l_rows,*)
   real(kind=c_double), intent(out) :: e_val_out(n_g_size)
   real(kind=c_double), intent(out) :: e_vec_out(n_l_rows,*)

   call elsi_ev_real(H_in,S_in,e_val_out,e_vec_out)

end subroutine

subroutine elsi_ev_complex_c_wrapper(H_in,S_in,e_val_out,e_vec_out)&
                                     bind(C,name="c_elsi_ev_complex")

   use, intrinsic :: iso_c_binding
   use ELSI
   use ELSI_DIMENSIONS

   implicit none

   complex(kind=c_double_complex), intent(in)  :: H_in(n_l_rows,*)
   complex(kind=c_double_complex), intent(in)  :: S_in(n_l_rows,*)
   real(kind=c_double),            intent(out) :: e_val_out(n_g_size)
   complex(kind=c_double_complex), intent(out) :: e_vec_out(n_l_rows,*)

   call elsi_ev_complex(H_in,S_in,e_val_out,e_vec_out)

end subroutine

function elsi_dm_real_c_wrapper(H_in,S_in,D_out) result(energy_out) &
                                bind(C,name="c_elsi_dm_real")

   use, intrinsic :: iso_c_binding
   use ELSI
   use ELSI_DIMENSIONS

   implicit none

   real(kind=c_double), intent(in)  :: H_in(n_l_rows,*)
   real(kind=c_double), intent(in)  :: S_in(n_l_rows,*)
   real(kind=c_double), intent(out) :: D_out(n_l_rows,*)
   real(kind=c_double)              :: energy_out

   call elsi_dm_real(H_in,S_in,D_out,energy_out)

end function

function elsi_dm_real_sparse_c_wrapper(H_in,S_in,D_out) result(energy_out) &
                                       bind(C,name="c_elsi_dm_real_sparse")

   use, intrinsic :: iso_c_binding
   use ELSI

   implicit none

   real(kind=c_double), intent(in)  :: H_in(*)
   real(kind=c_double), intent(in)  :: S_in(*)
   real(kind=c_double), intent(out) :: D_out(*)
   real(kind=c_double)              :: energy_out

   call elsi_dm_real_sparse(H_in,S_in,D_out,energy_out)

end function
