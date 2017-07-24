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

   use, intrinsic :: ISO_C_BINDING
   use ELSI
   use ELSI_C2F

   implicit none

   private

   integer(kind=c_int) :: n_g_size_c
   integer(kind=c_int) :: n_l_rows_c
   integer(kind=c_int) :: n_l_cols_c
   integer(kind=c_int) :: nnz_l_pexsi_c

contains

subroutine elsi_init_c_wrapper(handle_c,solver,parallel_mode,&
              matrix_format,n_basis,n_electron,n_state)&
   bind(C,name="c_elsi_init")

   implicit none

   type(c_ptr)                            :: handle_c
   integer(kind=c_int), value, intent(in) :: solver
   integer(kind=c_int), value, intent(in) :: parallel_mode
   integer(kind=c_int), value, intent(in) :: matrix_format
   integer(kind=c_int), value, intent(in) :: n_basis
   real(kind=c_double), value, intent(in) :: n_electron
   integer(kind=c_int), value, intent(in) :: n_state

   type(elsi_handle), pointer :: handle_f

   allocate(handle_f)

   call elsi_init(handle_f,solver,parallel_mode,matrix_format,&
           n_basis,n_electron,n_state)

   handle_c = c_loc(handle_f)

   n_g_size_c = n_basis

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

subroutine elsi_set_mpi_global_c_wrapper(handle_c,mpi_comm_global)&
   bind(C,name="c_elsi_set_mpi_global")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: mpi_comm_global

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mpi_global(handle_f,mpi_comm_global)

end subroutine

subroutine elsi_set_spin_c_wrapper(handle_c,n_spin,i_spin)&
   bind(C,name="c_elsi_set_spin")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: n_spin
   integer(kind=c_int), value, intent(in) :: i_spin

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_spin(handle_f,n_spin,i_spin)

end subroutine

subroutine elsi_set_kpoint_c_wrapper(handle_c,n_kpt,i_kpt,weight)&
   bind(C,name="c_elsi_set_kpoint")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: n_kpt
   integer(kind=c_int), value, intent(in) :: i_kpt
   real(kind=c_double), value, intent(in) :: weight

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_kpoint(handle_f,n_kpt,i_kpt,weight)

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

subroutine elsi_set_csc_c_wrapper(handle_c,nnz_g,nnz_l,n_l_cols,&
              row_ind,col_ptr)&
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

subroutine elsi_finalize_c_wrapper(handle_c)&
   bind(C,name="c_elsi_finalize")

   implicit none

   type(c_ptr), value :: handle_c

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_finalize(handle_f)

end subroutine

subroutine elsi_customize_c_wrapper(handle_c,print_detail,&
              overlap_is_unit,zero_threshold,no_singularity_check,&
              singularity_tolerance,stop_singularity,uplo)&
   bind(C,name="c_elsi_customize")

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

   print_detail_f         = convert_c_int_to_f_logical(print_detail)
   overlap_is_unit_f      = convert_c_int_to_f_logical(overlap_is_unit)
   no_singularity_check_f = convert_c_int_to_f_logical(no_singularity_check)
   stop_singularity_f     = convert_c_int_to_f_logical(stop_singularity)

   call c_f_pointer(handle_c,handle_f)

   call elsi_customize(handle_f,print_detail_f,overlap_is_unit_f,&
           zero_threshold,no_singularity_check_f,&
           singularity_tolerance,stop_singularity_f,uplo)

end subroutine

subroutine elsi_customize_mu_c_wrapper(handle_c,broadening_scheme,&
              broadening_width,occ_accuracy,mu_max_steps,spin_degeneracy)&
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

subroutine elsi_customize_omm_c_wrapper(handle_c,n_elpa_steps,omm_flavor,&
              eigen_shift,omm_tolerance,use_pspblas,omm_output)&
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

   use_pspblas_f = convert_c_int_to_f_logical(use_pspblas)
   omm_output_f  = convert_c_int_to_f_logical(omm_output)

   call c_f_pointer(handle_c,handle_f)

   call elsi_customize_omm(handle_f,n_elpa_steps,omm_flavor,&
           eigen_shift,omm_tolerance,use_pspblas_f,omm_output_f)

end subroutine

subroutine elsi_customize_pexsi_c_wrapper(handle_c,temperature,gap,&
              delta_e,n_poles,n_procs_per_pole,max_iteration,mu_min,&
              mu_max,mu0,mu_inertia_tolerance,mu_inertia_expansion,&
              mu_safeguard,n_electron_accuracy,matrix_type,&
              is_symbolic_factorize,ordering,np_symbolic_factorize,&
              verbosity)&
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

   call elsi_customize_pexsi(handle_f,temperature,gap,delta_e,n_poles,&
           n_procs_per_pole,max_iteration,mu_min,mu_max,mu0,&
           mu_inertia_tolerance,mu_inertia_expansion,mu_safeguard,&
           n_electron_accuracy,matrix_type,is_symbolic_factorize,&
           ordering,np_symbolic_factorize,verbosity)

end subroutine

subroutine elsi_customize_elpa_c_wrapper(handle_c,elpa_solver,&
              elpa_output)&
   bind(C,name="c_elsi_customize_elpa")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: elpa_solver
   integer(kind=c_int), value, intent(in) :: elpa_output

   type(elsi_handle), pointer :: handle_f
   logical :: elpa_output_f

   elpa_output_f = convert_c_int_to_f_logical(elpa_output)

   call c_f_pointer(handle_c,handle_f)

   call elsi_customize_elpa(handle_f,elpa_solver,elpa_output_f)

end subroutine

subroutine elsi_ev_real_c_wrapper(handle_c,H,S,e_val,e_vec)&
   bind(C,name="c_elsi_ev_real")

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

subroutine elsi_ev_complex_c_wrapper(handle_c,H,S,e_val,e_vec)&
   bind(C,name="c_elsi_ev_complex")

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

subroutine elsi_ev_real_sparse_c_wrapper(handle_c,H,S,e_val,e_vec)&
   bind(C,name="c_elsi_ev_real_sparse")

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

subroutine elsi_dm_real_c_wrapper(handle_c,H,S,D,energy)&
   bind(C,name="c_elsi_dm_real")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: H(n_l_rows_c,n_l_cols_c)
   real(kind=c_double) :: S(n_l_rows_c,n_l_cols_c)
   real(kind=c_double) :: D(n_l_rows_c,n_l_cols_c)
   real(kind=c_double) :: energy

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_dm_real(handle_f,H,S,D,energy)

end subroutine

subroutine elsi_dm_complex_c_wrapper(handle_c,H,S,D,energy)&
   bind(C,name="c_elsi_dm_complex")

   implicit none

   type(c_ptr), value             :: handle_c
   complex(kind=c_double_complex) :: H(n_l_rows_c,n_l_cols_c)
   complex(kind=c_double_complex) :: S(n_l_rows_c,n_l_cols_c)
   complex(kind=c_double_complex) :: D(n_l_rows_c,n_l_cols_c)
   real(kind=c_double)            :: energy

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_dm_complex(handle_f,H,S,D,energy)

end subroutine

subroutine elsi_dm_real_sparse_c_wrapper(handle_c,H,S,D,energy)&
   bind(C,name="c_elsi_dm_real_sparse")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: H(nnz_l_pexsi_c)
   real(kind=c_double) :: S(nnz_l_pexsi_c)
   real(kind=c_double) :: D(nnz_l_pexsi_c)
   real(kind=c_double) :: energy

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_dm_real_sparse(handle_f,H,S,D,energy)

end subroutine

subroutine elsi_collect_c_wrapper(handle_c,overlap_is_singular,&
              n_singular_basis,mu)&
   bind(C,name="c_elsi_collect")

   implicit none

   type(c_ptr), value  :: handle_c
   integer(kind=c_int) :: overlap_is_singular
   integer(kind=c_int) :: n_singular_basis
   real(kind=c_double) :: mu

   type(elsi_handle), pointer :: handle_f
   logical :: overlap_is_singular_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_collect(handle_f,overlap_is_singular_f,&
           n_singular_basis,mu)

   if(overlap_is_singular_f) then
      overlap_is_singular = 1
   else
      overlap_is_singular = 0
   endif

end subroutine

subroutine elsi_collect_pexsi_c_wrapper(handle_c,mu,edm,fdm)&
   bind(C,name="c_elsi_collect_pexsi")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: mu
   real(kind=c_double) :: edm(nnz_l_pexsi_c)
   real(kind=c_double) :: fdm(nnz_l_pexsi_c)

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_collect_pexsi(handle_f,mu,edm,fdm)

end subroutine

subroutine elsi_set_output_c_wrapper(handle_c,out_level)&
   bind(C,name="c_elsi_set_output")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: out_level

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_output(handle_f,out_level)

end subroutine

subroutine elsi_set_unit_ovlp_c_wrapper(handle_c,unit_ovlp)&
   bind(C,name="c_elsi_set_unit_ovlp")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: unit_ovlp

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_unit_ovlp(handle_f,unit_ovlp)

end subroutine

subroutine elsi_set_zero_def_c_wrapper(handle_c,zero_def)&
   bind(C,name="c_elsi_set_zero_def")

   implicit none

   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: zero_def

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_zero_def(handle_f,zero_def)

end subroutine

subroutine elsi_set_sing_check_c_wrapper(handle_c,sing_check)&
   bind(C,name="c_elsi_set_sing_check")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: sing_check

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sing_check(handle_f,sing_check)

end subroutine

subroutine elsi_set_sing_tol_c_wrapper(handle_c,sing_tol)&
   bind(C,name="c_elsi_set_sing_tol")

   implicit none

   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: sing_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sing_tol(handle_f,sing_tol)

end subroutine

subroutine elsi_set_sing_stop_c_wrapper(handle_c,sing_stop)&
   bind(C,name="c_elsi_set_sing_stop")
   
   implicit none
   
   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: sing_stop
   
   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sing_stop(handle_f,sing_stop)

end subroutine

subroutine elsi_set_uplo_c_wrapper(handle_c,uplo)&
   bind(C,name="c_elsi_set_uplo")
   
   implicit none
   
   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: uplo
   
   type(elsi_handle), pointer :: handle_f
   
   call c_f_pointer(handle_c,handle_f)
   
   call elsi_set_uplo(handle_f,uplo)

end subroutine

subroutine elsi_set_elpa_solver_c_wrapper(handle_c,elpa_solver)&
   bind(C,name="c_elsi_set_elpa_solver")
   
   implicit none
   
   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: elpa_solver
   
   type(elsi_handle), pointer :: handle_f
   
   call c_f_pointer(handle_c,handle_f)
   
   call elsi_set_elpa_solver(handle_f,elpa_solver)
   
end subroutine

subroutine elsi_set_omm_flavor_c_wrapper(handle_c,omm_flavor)&
   bind(C,name="c_elsi_set_omm_flavor")
   
   implicit none
   
   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: omm_flavor
   
   type(elsi_handle), pointer :: handle_f
   
   call c_f_pointer(handle_c,handle_f)
   
   call elsi_set_omm_flavor(handle_f,omm_flavor)
   
end subroutine

subroutine elsi_set_omm_n_elpa_c_wrapper(handle_c,n_elpa)&
   bind(C,name="c_elsi_set_omm_n_elpa")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: n_elpa

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_omm_n_elpa(handle_f,n_elpa)

end subroutine

subroutine elsi_set_omm_tol_c_wrapper(handle_c,min_tol)&
   bind(C,name="c_elsi_set_omm_tol")

   implicit none

   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: min_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_omm_tol(handle_f,min_tol)

end subroutine

subroutine elsi_set_omm_psp_c_wrapper(handle_c,use_psp)&
   bind(C,name="c_elsi_set_omm_psp")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: use_psp

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_omm_psp(handle_f,use_psp)

end subroutine

subroutine elsi_set_pexsi_n_mu_c_wrapper(handle_c,n_mu)&
   bind(C,name="c_elsi_set_pexsi_n_mu")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: n_mu

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_n_mu(handle_f,n_mu)

end subroutine

subroutine elsi_set_pexsi_n_pole_c_wrapper(handle_c,n_pole)&
   bind(C,name="c_elsi_set_pexsi_n_pole")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: n_pole

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_n_pole(handle_f,n_pole)

end subroutine

subroutine elsi_set_pexsi_np_per_pole_c_wrapper(handle_c,np_per_pole)&
   bind(C,name="c_elsi_set_pexsi_np_per_pole")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: np_per_pole

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_np_per_pole(handle_f,np_per_pole)

end subroutine

subroutine elsi_set_pexsi_np_symbo_c_wrapper(handle_c,np_symbo)&
   bind(C,name="c_elsi_set_pexsi_np_symbo")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: np_symbo

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_np_symbo(handle_f,np_symbo)

end subroutine

subroutine elsi_set_pexsi_temp_c_wrapper(handle_c,temp)&
   bind(C,name="c_elsi_set_pexsi_temp")

   implicit none

   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: temp

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_temp(handle_f,temp)

end subroutine

subroutine elsi_set_pexsi_gap_c_wrapper(handle_c,gap)&
   bind(C,name="c_elsi_set_pexsi_gap")

   implicit none

   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: gap

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_gap(handle_f,gap)

end subroutine

subroutine elsi_set_pexsi_mu_min_c_wrapper(handle_c,mu_min)&
   bind(C,name="c_elsi_set_pexsi_mu_min")
   
   implicit none
   
   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: mu_min
   
   type(elsi_handle), pointer :: handle_f
   
   call c_f_pointer(handle_c,handle_f)
   
   call elsi_set_pexsi_mu_min(handle_f,mu_min)
   
end subroutine

subroutine elsi_set_pexsi_mu_max_c_wrapper(handle_c,mu_max)&
   bind(C,name="c_elsi_set_pexsi_mu_max")
   
   implicit none
   
   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: mu_max
   
   type(elsi_handle), pointer :: handle_f
   
   call c_f_pointer(handle_c,handle_f)
   
   call elsi_set_pexsi_mu_max(handle_f,mu_max)
   
end subroutine

subroutine elsi_set_pexsi_inertia_tol_c_wrapper(handle_c,inertia_tol)&
   bind(C,name="c_elsi_set_pexsi_inertia_tol")
   
   implicit none
   
   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: inertia_tol
   
   type(elsi_handle), pointer :: handle_f
   
   call c_f_pointer(handle_c,handle_f)
   
   call elsi_set_pexsi_inertia_tol(handle_f,inertia_tol)
   
end subroutine

subroutine elsi_set_sips_slice_type_c_wrapper(handle_c,slice_type)&
   bind(C,name="c_elsi_set_sips_slice_type")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: slice_type

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_slice_type(handle_f,slice_type)

end subroutine

subroutine elsi_set_sips_n_slice_c_wrapper(handle_c,n_slice)&
   bind(C,name="c_elsi_set_sips_n_slice")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: n_slice

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_n_slice(handle_f,n_slice)

end subroutine

subroutine elsi_set_sips_left_bound_c_wrapper(handle_c,left_bound)&
   bind(C,name="c_elsi_set_sips_left_bound")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: left_bound

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_left_bound(handle_f,left_bound)

end subroutine

subroutine elsi_set_sips_slice_buffer_c_wrapper(handle_c,slice_buffer)&
   bind(C,name="c_elsi_set_sips_slice_buffer")

   implicit none

   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: slice_buffer

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_slice_buffer(handle_f,slice_buffer)

end subroutine

subroutine elsi_set_mu_broaden_scheme_c_wrapper(handle_c,broaden_scheme)&
   bind(C,name="c_elsi_set_mu_broaden_scheme")

   implicit none

   type(c_ptr),         value             :: handle_c
   integer(kind=c_int), value, intent(in) :: broaden_scheme

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_broaden_scheme(handle_f,broaden_scheme)

end subroutine

subroutine elsi_set_mu_broaden_width_c_wrapper(handle_c,broaden_width)&
   bind(C,name="c_elsi_set_mu_broaden_width")
   
   implicit none
   
   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: broaden_width
   
   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_broaden_width(handle_f,broaden_width)

end subroutine

subroutine elsi_set_mu_tol_c_wrapper(handle_c,mu_tol)&
   bind(C,name="c_elsi_set_mu_tol")

   implicit none

   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: mu_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_tol(handle_f,mu_tol)

end subroutine

subroutine elsi_set_mu_spin_degen_c_wrapper(handle_c,spin_degen)&
   bind(C,name="c_elsi_set_mu_spin_degen")

   implicit none

   type(c_ptr),         value             :: handle_c
   real(kind=c_double), value, intent(in) :: spin_degen

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_spin_degen(handle_f,spin_degen)

end subroutine

subroutine elsi_get_pexsi_mu_min_c_wrapper(handle_c,mu_min)&
   bind(C,name="c_elsi_get_pexsi_mu_min")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: mu_min

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_pexsi_mu_min(handle_f,mu_min)

end subroutine

subroutine elsi_get_pexsi_mu_max_c_wrapper(handle_c,mu_max)&
   bind(C,name="c_elsi_get_pexsi_mu_max")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: mu_max

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_pexsi_mu_max(handle_f,mu_max)

end subroutine

subroutine elsi_get_ovlp_sing_c_wrapper(handle_c,ovlp_sing)&
   bind(C,name="c_elsi_get_ovlp_sing")

   implicit none

   type(c_ptr), value  :: handle_c
   integer(kind=c_int) :: ovlp_sing

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_ovlp_sing(handle_f,ovlp_sing)

end subroutine

subroutine elsi_get_mu_c_wrapper(handle_c,mu)&
   bind(C,name="c_elsi_get_mu")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: mu

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_mu(handle_f,mu)

end subroutine

subroutine elsi_get_edm_real_c_wrapper(handle_c,edm)&
   bind(C,name="c_elsi_get_edm_real")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: edm(n_l_rows_c,n_l_cols_c)

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_edm_real(handle_f,edm)

end subroutine

subroutine elsi_get_edm_complex_c_wrapper(handle_c,edm)&
   bind(C,name="c_elsi_get_edm_complex")

   implicit none

   type(c_ptr), value             :: handle_c
   complex(kind=c_double_complex) :: edm(n_l_rows_c,n_l_cols_c)

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_edm_complex(handle_f,edm)

end subroutine

end module ELSI_C_INTERFACE
