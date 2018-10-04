! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains the C interfaces of ELSI.
!!
module ELSI_C_INTERFACE

   use, intrinsic :: ISO_C_BINDING
   use ELSI

   implicit none

   private

contains

!>
!! This routine encodes the standard convention that 0 is .false., any other
!! integer is .true..
!!
function c_int_to_f_logical(int_c) result(logical_f)

   implicit none

   integer(kind=c_int), intent(in) :: int_c
   logical :: logical_f

   if(int_c == 0) then
     logical_f = .false.
   else
     logical_f = .true.
   end if

end function

!>
!! This routine converts a C string into a Fortran string. A Fortran string is
!! NOT just a character array without a NULL character. A Fortran string (i.e.
!! char(*)) is a separate data type from a character array (i.e. char,
!! dimension(*)) and they are NOT interoperable in interfaces.
!!
function c_string_to_f_string(string_c) result(string_f)

   implicit none

   character(kind=c_char,len=1), intent(in) :: string_c(*)
   character(len=:), allocatable :: string_f

   integer(kind=c_int) :: string_f_len

   string_f_len = 0

   do
     if(string_c(string_f_len+1) == C_NULL_CHAR) then
        exit
     end if

     string_f_len = string_f_len+1
   end do

   allocate(character(len=string_f_len) :: string_f)

   string_f = transfer(string_c(1:string_f_len),string_f)

end function

subroutine elsi_init_c_wrapper(handle_c,solver,parallel_mode,matrix_format,&
              n_basis,n_electron,n_state)&
   bind(C,name="c_elsi_init")

   implicit none

   type(c_ptr) :: handle_c
   integer(kind=c_int), value, intent(in) :: solver
   integer(kind=c_int), value, intent(in) :: parallel_mode
   integer(kind=c_int), value, intent(in) :: matrix_format
   integer(kind=c_int), value, intent(in) :: n_basis
   real(kind=c_double), value, intent(in) :: n_electron
   integer(kind=c_int), value, intent(in) :: n_state

   type(elsi_handle), pointer :: handle_f

   allocate(handle_f)

   call elsi_init(handle_f,solver,parallel_mode,matrix_format,n_basis,&
           n_electron,n_state)

   handle_c = c_loc(handle_f)

end subroutine

subroutine elsi_set_mpi_c_wrapper(handle_c,mpi_comm)&
   bind(C,name="c_elsi_set_mpi")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: mpi_comm

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mpi(handle_f,mpi_comm)

end subroutine

subroutine elsi_set_mpi_global_c_wrapper(handle_c,mpi_comm_global)&
   bind(C,name="c_elsi_set_mpi_global")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: mpi_comm_global

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mpi_global(handle_f,mpi_comm_global)

end subroutine

subroutine elsi_set_spin_c_wrapper(handle_c,n_spin,i_spin)&
   bind(C,name="c_elsi_set_spin")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_spin
   integer(kind=c_int), value, intent(in) :: i_spin

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_spin(handle_f,n_spin,i_spin)

end subroutine

subroutine elsi_set_kpoint_c_wrapper(handle_c,n_kpt,i_kpt,weight)&
   bind(C,name="c_elsi_set_kpoint")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
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

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: blacs_ctxt
   integer(kind=c_int), value, intent(in) :: block_size

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_blacs(handle_f,blacs_ctxt,block_size)

end subroutine

subroutine elsi_set_csc_c_wrapper(handle_c,nnz_g,nnz_l,n_lcol,row_ind,col_ptr)&
   bind(C,name="c_elsi_set_csc")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: nnz_g
   integer(kind=c_int), value, intent(in) :: nnz_l
   integer(kind=c_int), value, intent(in) :: n_lcol
   integer(kind=c_int), intent(in) :: row_ind(nnz_l)
   integer(kind=c_int), intent(in) :: col_ptr(n_lcol+1)

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_csc(handle_f,nnz_g,nnz_l,n_lcol,row_ind,col_ptr)

end subroutine

subroutine elsi_finalize_c_wrapper(handle_c)&
   bind(C,name="c_elsi_finalize")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_finalize(handle_f)

end subroutine

subroutine elsi_ev_real_c_wrapper(handle_c,ham_c,ovlp_c,eval_c,evec_c)&
   bind(C,name="c_elsi_ev_real")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: handle_f
   real(kind=c_double), pointer :: ham_f(:,:)
   real(kind=c_double), pointer :: ovlp_f(:,:)
   real(kind=c_double), pointer :: eval_f(:)
   real(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%ph%n_basis
   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_real(handle_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine elsi_ev_complex_c_wrapper(handle_c,ham_c,ovlp_c,eval_c,evec_c)&
   bind(C,name="c_elsi_ev_complex")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: handle_f
   complex(kind=c_double_complex), pointer :: ham_f(:,:)
   complex(kind=c_double_complex), pointer :: ovlp_f(:,:)
   real(kind=c_double), pointer :: eval_f(:)
   complex(kind=c_double_complex), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%ph%n_basis
   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_complex(handle_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine elsi_ev_real_sparse_c_wrapper(handle_c,ham_c,ovlp_c,eval_c,evec_c)&
   bind(C,name="c_elsi_ev_real_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: handle_f
   real(kind=c_double), pointer :: ham_f(:)
   real(kind=c_double), pointer :: ovlp_f(:)
   real(kind=c_double), pointer :: eval_f(:)
   real(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%ph%n_basis
   nnz_l = handle_f%bh%nnz_l_sp
   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_real_sparse(handle_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine elsi_ev_complex_sparse_c_wrapper(handle_c,ham_c,ovlp_c,eval_c,&
              evec_c)&
   bind(C,name="c_elsi_ev_complex_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: handle_f
   complex(kind=c_double), pointer :: ham_f(:)
   complex(kind=c_double), pointer :: ovlp_f(:)
   real(kind=c_double), pointer :: eval_f(:)
   complex(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%ph%n_basis
   nnz_l = handle_f%bh%nnz_l_sp
   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_complex_sparse(handle_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine elsi_dm_real_c_wrapper(handle_c,ham_c,ovlp_c,dm_c,energy)&
   bind(C,name="c_elsi_dm_real")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle), pointer :: handle_f
   real(kind=c_double), pointer :: ham_f(:,:)
   real(kind=c_double), pointer :: ovlp_f(:,:)
   real(kind=c_double), pointer :: dm_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%ph%n_basis
   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_dm_real(handle_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine elsi_dm_complex_c_wrapper(handle_c,ham_c,ovlp_c,dm_c,energy)&
   bind(C,name="c_elsi_dm_complex")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle), pointer :: handle_f
   complex(kind=c_double_complex), pointer :: ham_f(:,:)
   complex(kind=c_double_complex), pointer :: ovlp_f(:,:)
   complex(kind=c_double_complex), pointer :: dm_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%ph%n_basis
   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_dm_complex(handle_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine elsi_dm_real_sparse_c_wrapper(handle_c,ham_c,ovlp_c,dm_c,energy)&
   bind(C,name="c_elsi_dm_real_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle), pointer :: handle_f
   real(kind=c_double), pointer :: ham_f(:)
   real(kind=c_double), pointer :: ovlp_f(:)
   real(kind=c_double), pointer :: dm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%bh%nnz_l_sp

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(dm_c,dm_f,shape=[nnz_l])

   call elsi_dm_real_sparse(handle_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine elsi_dm_complex_sparse_c_wrapper(handle_c,ham_c,ovlp_c,dm_c,energy)&
   bind(C,name="c_elsi_dm_complex_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle), pointer :: handle_f
   complex(kind=c_double), pointer :: ham_f(:)
   complex(kind=c_double), pointer :: ovlp_f(:)
   complex(kind=c_double), pointer :: dm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%bh%nnz_l_sp

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(dm_c,dm_f,shape=[nnz_l])

   call elsi_dm_complex_sparse(handle_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine elsi_set_output_c_wrapper(handle_c,out_level)&
   bind(C,name="c_elsi_set_output")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: out_level

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_output(handle_f,out_level)

end subroutine

subroutine elsi_set_unit_ovlp_c_wrapper(handle_c,unit_ovlp)&
   bind(C,name="c_elsi_set_unit_ovlp")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: unit_ovlp

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_unit_ovlp(handle_f,unit_ovlp)

end subroutine

subroutine elsi_set_zero_def_c_wrapper(handle_c,zero_def)&
   bind(C,name="c_elsi_set_zero_def")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: zero_def

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_zero_def(handle_f,zero_def)

end subroutine

subroutine elsi_set_sing_check_c_wrapper(handle_c,sing_check)&
   bind(C,name="c_elsi_set_sing_check")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: sing_check

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sing_check(handle_f,sing_check)

end subroutine

subroutine elsi_set_sing_tol_c_wrapper(handle_c,sing_tol)&
   bind(C,name="c_elsi_set_sing_tol")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: sing_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sing_tol(handle_f,sing_tol)

end subroutine

subroutine elsi_set_sing_stop_c_wrapper(handle_c,sing_stop)&
   bind(C,name="c_elsi_set_sing_stop")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: sing_stop

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sing_stop(handle_f,sing_stop)

end subroutine

subroutine elsi_set_elpa_solver_c_wrapper(handle_c,elpa_solver)&
   bind(C,name="c_elsi_set_elpa_solver")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: elpa_solver

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_elpa_solver(handle_f,elpa_solver)

end subroutine

subroutine elsi_set_elpa_n_single_c_wrapper(handle_c,n_single)&
   bind(C,name="c_elsi_set_elpa_n_single")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_single

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_elpa_n_single(handle_f,n_single)

end subroutine

subroutine elsi_set_elpa_gpu_c_wrapper(handle_c,use_gpu)&
   bind(C,name="c_elsi_set_elpa_gpu")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: use_gpu

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_elpa_gpu(handle_f,use_gpu)

end subroutine

subroutine elsi_set_elpa_gpu_kernels_c_wrapper(handle_c,use_gpu_kernels)&
   bind(C,name="c_elsi_set_elpa_gpu_kernels")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: use_gpu_kernels

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_elpa_gpu_kernels(handle_f,use_gpu_kernels)

end subroutine

subroutine elsi_set_elpa_autotune_c_wrapper(handle_c,use_autotune)&
   bind(C,name="c_elsi_set_elpa_autotune")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: use_autotune

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_elpa_autotune(handle_f,use_autotune)

end subroutine

subroutine elsi_set_omm_flavor_c_wrapper(handle_c,omm_flavor)&
   bind(C,name="c_elsi_set_omm_flavor")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: omm_flavor

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_omm_flavor(handle_f,omm_flavor)

end subroutine

subroutine elsi_set_omm_n_elpa_c_wrapper(handle_c,n_elpa)&
   bind(C,name="c_elsi_set_omm_n_elpa")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_elpa

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_omm_n_elpa(handle_f,n_elpa)

end subroutine

subroutine elsi_set_omm_tol_c_wrapper(handle_c,omm_tol)&
   bind(C,name="c_elsi_set_omm_tol")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: omm_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_omm_tol(handle_f,omm_tol)

end subroutine

subroutine elsi_set_pexsi_n_mu_c_wrapper(handle_c,n_mu)&
   bind(C,name="c_elsi_set_pexsi_n_mu")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_mu

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_n_mu(handle_f,n_mu)

end subroutine

subroutine elsi_set_pexsi_n_pole_c_wrapper(handle_c,n_pole)&
   bind(C,name="c_elsi_set_pexsi_n_pole")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_pole

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_n_pole(handle_f,n_pole)

end subroutine

subroutine elsi_set_pexsi_np_per_pole_c_wrapper(handle_c,np_per_pole)&
   bind(C,name="c_elsi_set_pexsi_np_per_pole")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: np_per_pole

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_np_per_pole(handle_f,np_per_pole)

end subroutine

subroutine elsi_set_pexsi_np_symbo_c_wrapper(handle_c,np_symbo)&
   bind(C,name="c_elsi_set_pexsi_np_symbo")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: np_symbo

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_np_symbo(handle_f,np_symbo)

end subroutine

subroutine elsi_set_pexsi_ordering_c_wrapper(handle_c,ordering)&
   bind(C,name="c_elsi_set_pexsi_ordering")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: ordering

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_ordering(handle_f,ordering)

end subroutine

subroutine elsi_set_pexsi_temp_c_wrapper(handle_c,temp)&
   bind(C,name="c_elsi_set_pexsi_temp")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: temp

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_temp(handle_f,temp)

end subroutine

subroutine elsi_set_pexsi_gap_c_wrapper(handle_c,gap)&
   bind(C,name="c_elsi_set_pexsi_gap")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: gap

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_gap(handle_f,gap)

end subroutine

subroutine elsi_set_pexsi_delta_e_c_wrapper(handle_c,delta_e)&
   bind(C,name="c_elsi_set_pexsi_delta_e")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: delta_e

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_delta_e(handle_f,delta_e)

end subroutine

subroutine elsi_set_pexsi_mu_min_c_wrapper(handle_c,mu_min)&
   bind(C,name="c_elsi_set_pexsi_mu_min")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: mu_min

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_mu_min(handle_f,mu_min)

end subroutine

subroutine elsi_set_pexsi_mu_max_c_wrapper(handle_c,mu_max)&
   bind(C,name="c_elsi_set_pexsi_mu_max")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: mu_max

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_mu_max(handle_f,mu_max)

end subroutine

subroutine elsi_set_pexsi_inertia_tol_c_wrapper(handle_c,inertia_tol)&
   bind(C,name="c_elsi_set_pexsi_inertia_tol")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: inertia_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_inertia_tol(handle_f,inertia_tol)

end subroutine

subroutine elsi_set_sips_n_elpa_c_wrapper(handle_c,n_elpa)&
   bind(C,name="c_elsi_set_sips_n_elpa")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_elpa

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_n_elpa(handle_f,n_elpa)

end subroutine

subroutine elsi_set_sips_n_slice_c_wrapper(handle_c,n_slice)&
   bind(C,name="c_elsi_set_sips_n_slice")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_slice

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_n_slice(handle_f,n_slice)

end subroutine

subroutine elsi_set_sips_first_ev_c_wrapper(handle_c,first_ev)&
   bind(C,name="c_elsi_set_sips_first_ev")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: first_ev

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_first_ev(handle_f,first_ev)

end subroutine

subroutine elsi_set_sips_slice_type_c_wrapper(handle_c,slice_type)&
   bind(C,name="c_elsi_set_sips_slice_type")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: slice_type

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_slice_type(handle_f,slice_type)

end subroutine

subroutine elsi_set_sips_buffer_c_wrapper(handle_c,buffer)&
   bind(C,name="c_elsi_set_sips_buffer")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: buffer

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_buffer(handle_f,buffer)

end subroutine

subroutine elsi_set_sips_inertia_tol_c_wrapper(handle_c,inertia_tol)&
   bind(C,name="c_elsi_set_sips_inertia_tol")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: inertia_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_inertia_tol(handle_f,inertia_tol)

end subroutine

subroutine elsi_set_sips_interval_c_wrapper(handle_c,lower,upper)&
   bind(C,name="c_elsi_set_sips_interval")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: lower
   real(kind=c_double), value, intent(in) :: upper

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_interval(handle_f,lower,upper)

end subroutine

subroutine elsi_set_ntpoly_method_c_wrapper(handle_c,nt_method)&
   bind(C,name="c_elsi_set_ntpoly_method")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: nt_method

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_ntpoly_method(handle_f,nt_method)

end subroutine

subroutine elsi_set_ntpoly_isr_c_wrapper(handle_c,nt_isr)&
   bind(C,name="c_elsi_set_ntpoly_isr")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: nt_isr

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_ntpoly_isr(handle_f,nt_isr)

end subroutine

subroutine elsi_set_ntpoly_tol_c_wrapper(handle_c,nt_tol)&
   bind(C,name="c_elsi_set_ntpoly_tol")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: nt_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_ntpoly_tol(handle_f,nt_tol)

end subroutine

subroutine elsi_set_ntpoly_filter_c_wrapper(handle_c,filter)&
   bind(C,name="c_elsi_set_ntpoly_filter")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: filter

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_ntpoly_filter(handle_f,filter)

end subroutine

subroutine elsi_set_ntpoly_max_iter_c_wrapper(handle_c,max_iter)&
   bind(C,name="c_elsi_set_ntpoly_max_iter")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: max_iter

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_ntpoly_max_iter(handle_f,max_iter)

end subroutine

subroutine elsi_set_mu_broaden_scheme_c_wrapper(handle_c,broaden_scheme)&
   bind(C,name="c_elsi_set_mu_broaden_scheme")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: broaden_scheme

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_broaden_scheme(handle_f,broaden_scheme)

end subroutine

subroutine elsi_set_mu_broaden_width_c_wrapper(handle_c,broaden_width)&
   bind(C,name="c_elsi_set_mu_broaden_width")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: broaden_width

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_broaden_width(handle_f,broaden_width)

end subroutine

subroutine elsi_set_mu_tol_c_wrapper(handle_c,mu_tol)&
   bind(C,name="c_elsi_set_mu_tol")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: mu_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_tol(handle_f,mu_tol)

end subroutine

subroutine elsi_set_mu_spin_degen_c_wrapper(handle_c,spin_degen)&
   bind(C,name="c_elsi_set_mu_spin_degen")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: spin_degen

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_spin_degen(handle_f,spin_degen)

end subroutine

subroutine elsi_set_mu_mp_order_c_wrapper(handle_c,mp_order)&
   bind(C,name="c_elsi_set_mu_mp_order")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: mp_order

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_mp_order(handle_f,mp_order)

end subroutine

subroutine elsi_set_output_log_c_wrapper(handle_c,out_log)&
   bind(C,name="c_elsi_set_output_log")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: out_log

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_output_log(handle_f,out_log)

end subroutine

subroutine elsi_get_pexsi_mu_min_c_wrapper(handle_c,mu_min)&
   bind(C,name="c_elsi_get_pexsi_mu_min")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), intent(out) :: mu_min

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_pexsi_mu_min(handle_f,mu_min)

end subroutine

subroutine elsi_get_pexsi_mu_max_c_wrapper(handle_c,mu_max)&
   bind(C,name="c_elsi_get_pexsi_mu_max")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), intent(out) :: mu_max

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_pexsi_mu_max(handle_f,mu_max)

end subroutine

subroutine elsi_get_initialized_c_wrapper(handle_c,handle_init)&
   bind(C,name="c_elsi_get_initialized")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), intent(out) :: handle_init

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_initialized(handle_f,handle_init)

end subroutine

subroutine elsi_get_version_c_wrapper(major,minor,patch)&
   bind(C,name="c_elsi_get_version")

   implicit none

   integer(kind=c_int), intent(out) :: major
   integer(kind=c_int), intent(out) :: minor
   integer(kind=c_int), intent(out) :: patch

   call elsi_get_version(major,minor,patch)

end subroutine

subroutine elsi_get_datestamp_c_wrapper(datestamp)&
   bind(C,name="c_elsi_get_datestamp")

   implicit none

   integer(kind=c_int), intent(out) :: datestamp

   call elsi_get_datestamp(datestamp)

end subroutine

subroutine elsi_get_n_sing_c_wrapper(handle_c,n_sing)&
   bind(C,name="c_elsi_get_n_sing")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), intent(out) :: n_sing

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_n_sing(handle_f,n_sing)

end subroutine

subroutine elsi_get_mu_c_wrapper(handle_c,mu)&
   bind(C,name="c_elsi_get_mu")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), intent(out) :: mu

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_mu(handle_f,mu)

end subroutine

subroutine elsi_get_edm_real_c_wrapper(handle_c,edm_c)&
   bind(C,name="c_elsi_get_edm_real")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle), pointer :: handle_f
   real(kind=c_double), pointer :: edm_f(:,:)

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(edm_c,edm_f,shape=[lrow,lcol])

   call elsi_get_edm_real(handle_f,edm_f)

end subroutine

subroutine elsi_get_edm_complex_c_wrapper(handle_c,edm_c)&
   bind(C,name="c_elsi_get_edm_complex")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle), pointer :: handle_f
   complex(kind=c_double), pointer :: edm_f(:,:)

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(edm_c,edm_f,shape=[lrow,lcol])

   call elsi_get_edm_complex(handle_f,edm_f)

end subroutine

subroutine elsi_get_edm_real_sparse_c_wrapper(handle_c,edm_c)&
   bind(C,name="c_elsi_get_edm_real_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle), pointer :: handle_f
   real(kind=c_double), pointer :: edm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%bh%nnz_l_sp

   call c_f_pointer(edm_c,edm_f,shape=[nnz_l])

   call elsi_get_edm_real_sparse(handle_f,edm_f)

end subroutine

subroutine elsi_get_edm_complex_sparse_c_wrapper(handle_c,edm_c)&
   bind(C,name="c_elsi_get_edm_complex_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle), pointer :: handle_f
   complex(kind=c_double), pointer :: edm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%bh%nnz_l_sp

   call c_f_pointer(edm_c,edm_f,shape=[nnz_l])

   call elsi_get_edm_complex_sparse(handle_f,edm_f)

end subroutine

subroutine elsi_init_rw_c_wrapper(handle_c,rw_task,parallel_mode,n_basis,&
              n_electron)&
   bind(C,name="c_elsi_init_rw")

   implicit none

   type(c_ptr) :: handle_c
   integer(kind=c_int), value, intent(in) :: rw_task
   integer(kind=c_int), value, intent(in) :: parallel_mode
   integer(kind=c_int), value, intent(in) :: n_basis
   real(kind=c_double), value, intent(in) :: n_electron

   type(elsi_rw_handle), pointer :: handle_f

   allocate(handle_f)

   call elsi_init_rw(handle_f,rw_task,parallel_mode,n_basis,n_electron)

   handle_c = c_loc(handle_f)

end subroutine

subroutine elsi_set_rw_mpi_c_wrapper(handle_c,mpi_comm)&
   bind(C,name="c_elsi_set_rw_mpi")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: mpi_comm

   type(elsi_rw_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_rw_mpi(handle_f,mpi_comm)

end subroutine

subroutine elsi_set_rw_blacs_c_wrapper(handle_c,blacs_ctxt,block_size)&
   bind(C,name="c_elsi_set_rw_blacs")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: blacs_ctxt
   integer(kind=c_int), value, intent(in) :: block_size

   type(elsi_rw_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_rw_blacs(handle_f,blacs_ctxt,block_size)

end subroutine

subroutine elsi_set_rw_csc_c_wrapper(handle_c,nnz_g,nnz_l,n_lcol)&
   bind(C,name="c_elsi_set_rw_csc")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: nnz_g
   integer(kind=c_int), value, intent(in) :: nnz_l
   integer(kind=c_int), value, intent(in) :: n_lcol

   type(elsi_rw_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_rw_csc(handle_f,nnz_g,nnz_l,n_lcol)

end subroutine

subroutine elsi_finalize_rw_c_wrapper(handle_c)&
   bind(C,name="c_elsi_finalize_rw")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c

   type(elsi_rw_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_finalize_rw(handle_f)

end subroutine

subroutine elsi_set_rw_zero_def_c_wrapper(handle_c,zero_def)&
   bind(C,name="c_elsi_set_rw_zero_def")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: zero_def

   type(elsi_rw_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_rw_zero_def(handle_f,zero_def)

end subroutine

subroutine elsi_read_mat_dim_c_wrapper(handle_c,name_c,n_electrons,n_basis,&
              n_lrow,n_lcol)&
   bind(C,name="c_elsi_read_mat_dim")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   character(kind=c_char,len=1) :: name_c(128)
   real(kind=c_double), intent(out) :: n_electrons
   integer(kind=c_int), intent(out) :: n_basis
   integer(kind=c_int), intent(out) :: n_lrow
   integer(kind=c_int), intent(out) :: n_lcol

   type(elsi_rw_handle), pointer :: handle_f

   character(len=:), allocatable :: name_f

   call c_f_pointer(handle_c,handle_f)

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_dim(handle_f,name_f,n_electrons,n_basis,n_lrow,n_lcol)

end subroutine

subroutine elsi_read_mat_dim_sparse_c_wrapper(handle_c,name_c,n_electrons,&
              n_basis,nnz_g,nnz_l,n_lcol)&
   bind(C,name="c_elsi_read_mat_dim_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   character(kind=c_char,len=1) :: name_c(128)
   real(kind=c_double), intent(out) :: n_electrons
   integer(kind=c_int), intent(out) :: n_basis
   integer(kind=c_int), intent(out) :: nnz_g
   integer(kind=c_int), intent(out) :: nnz_l
   integer(kind=c_int), intent(out) :: n_lcol

   type(elsi_rw_handle), pointer :: handle_f

   character(len=:), allocatable :: name_f

   call c_f_pointer(handle_c,handle_f)

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_dim_sparse(handle_f,name_f,n_electrons,n_basis,nnz_g,&
           nnz_l,n_lcol)

end subroutine

subroutine elsi_read_mat_real_c_wrapper(handle_c,name_c,mat_c)&
   bind(C,name="c_elsi_read_mat_real")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   real(kind=c_double), pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_real(handle_f,name_f,mat_f)

end subroutine

subroutine elsi_read_mat_real_sparse_c_wrapper(handle_c,name_c,row_ind_c,&
              col_ptr_c,mat_c)&
   bind(C,name="c_elsi_read_mat_real_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: row_ind_c
   type(c_ptr), value, intent(in) :: col_ptr_c
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   integer(kind=c_int), pointer :: row_ind_f(:)
   integer(kind=c_int), pointer :: col_ptr_f(:)
   real(kind=c_double), pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%bh%nnz_l_sp
   lcol = handle_f%bh%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_real_sparse(handle_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine elsi_write_mat_real_c_wrapper(handle_c,name_c,mat_c)&
   bind(C,name="c_elsi_write_mat_real")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   real(kind=c_double), pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_real(handle_f,name_f,mat_f)

end subroutine

subroutine elsi_write_mat_real_sparse_c_wrapper(handle_c,name_c,row_ind_c,&
              col_ptr_c,mat_c)&
   bind(C,name="c_elsi_write_mat_real_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: row_ind_c
   type(c_ptr), value, intent(in) :: col_ptr_c
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   integer(kind=c_int), pointer :: row_ind_f(:)
   integer(kind=c_int), pointer :: col_ptr_f(:)
   real(kind=c_double), pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%bh%nnz_l_sp
   lcol = handle_f%bh%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_real_sparse(handle_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine elsi_read_mat_complex_c_wrapper(handle_c,name_c,mat_c)&
   bind(C,name="c_elsi_read_mat_complex")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   complex(kind=c_double), pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_complex(handle_f,name_f,mat_f)

end subroutine

subroutine elsi_read_mat_complex_sparse_c_wrapper(handle_c,name_c,row_ind_c,&
              col_ptr_c,mat_c)&
   bind(C,name="c_elsi_read_mat_complex_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: row_ind_c
   type(c_ptr), value, intent(in) :: col_ptr_c
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   integer(kind=c_int), pointer :: row_ind_f(:)
   integer(kind=c_int), pointer :: col_ptr_f(:)
   complex(kind=c_double), pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%bh%nnz_l_sp
   lcol = handle_f%bh%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_complex_sparse(handle_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine elsi_write_mat_complex_c_wrapper(handle_c,name_c,mat_c)&
   bind(C,name="c_elsi_write_mat_complex")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   complex(kind=c_double), pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%bh%n_lrow
   lcol = handle_f%bh%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_complex(handle_f,name_f,mat_f)

end subroutine

subroutine elsi_write_mat_complex_sparse_c_wrapper(handle_c,name_c,row_ind_c,&
              col_ptr_c,mat_c)&
   bind(C,name="c_elsi_write_mat_complex_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: row_ind_c
   type(c_ptr), value, intent(in) :: col_ptr_c
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   integer(kind=c_int), pointer :: row_ind_f(:)
   integer(kind=c_int), pointer :: col_ptr_f(:)
   complex(kind=c_double), pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%bh%nnz_l_sp
   lcol = handle_f%bh%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_complex_sparse(handle_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

end module ELSI_C_INTERFACE
