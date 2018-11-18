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

subroutine elsi_init_c_wrapper(h_c,solver,parallel_mode,matrix_format,n_basis,&
   n_electron,n_state)&
   bind(C,name="c_elsi_init")

   implicit none

   type(c_ptr) :: h_c
   integer(kind=c_int), value, intent(in) :: solver
   integer(kind=c_int), value, intent(in) :: parallel_mode
   integer(kind=c_int), value, intent(in) :: matrix_format
   integer(kind=c_int), value, intent(in) :: n_basis
   real(kind=c_double), value, intent(in) :: n_electron
   integer(kind=c_int), value, intent(in) :: n_state

   type(elsi_handle), pointer :: h_f

   allocate(h_f)

   call elsi_init(h_f,solver,parallel_mode,matrix_format,n_basis,n_electron,&
        n_state)

   h_c = c_loc(h_f)

end subroutine

subroutine elsi_set_mpi_c_wrapper(h_c,mpi_comm)&
   bind(C,name="c_elsi_set_mpi")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: mpi_comm

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mpi(h_f,mpi_comm)

end subroutine

subroutine elsi_set_mpi_global_c_wrapper(h_c,mpi_comm_global)&
   bind(C,name="c_elsi_set_mpi_global")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: mpi_comm_global

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mpi_global(h_f,mpi_comm_global)

end subroutine

subroutine elsi_set_spin_c_wrapper(h_c,n_spin,i_spin)&
   bind(C,name="c_elsi_set_spin")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_spin
   integer(kind=c_int), value, intent(in) :: i_spin

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_spin(h_f,n_spin,i_spin)

end subroutine

subroutine elsi_set_kpoint_c_wrapper(h_c,n_kpt,i_kpt,i_wt)&
   bind(C,name="c_elsi_set_kpoint")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_kpt
   integer(kind=c_int), value, intent(in) :: i_kpt
   real(kind=c_double), value, intent(in) :: i_wt

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_kpoint(h_f,n_kpt,i_kpt,i_wt)

end subroutine

subroutine elsi_set_blacs_c_wrapper(h_c,blacs_ctxt,block_size)&
   bind(C,name="c_elsi_set_blacs")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: blacs_ctxt
   integer(kind=c_int), value, intent(in) :: block_size

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_blacs(h_f,blacs_ctxt,block_size)

end subroutine

subroutine elsi_set_csc_c_wrapper(h_c,nnz_g,nnz_l,n_lcol,row_ind,col_ptr)&
   bind(C,name="c_elsi_set_csc")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: nnz_g
   integer(kind=c_int), value, intent(in) :: nnz_l
   integer(kind=c_int), value, intent(in) :: n_lcol
   integer(kind=c_int), intent(in) :: row_ind(nnz_l)
   integer(kind=c_int), intent(in) :: col_ptr(n_lcol+1)

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_csc(h_f,nnz_g,nnz_l,n_lcol,row_ind,col_ptr)

end subroutine

subroutine elsi_reinit_c_wrapper(h_c)&
   bind(C,name="c_elsi_reinit")

   implicit none

   type(c_ptr), value, intent(in) :: h_c

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_reinit(h_f)

end subroutine

subroutine elsi_finalize_c_wrapper(h_c)&
   bind(C,name="c_elsi_finalize")

   implicit none

   type(c_ptr), value, intent(in) :: h_c

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_finalize(h_f)

end subroutine

subroutine elsi_ev_real_c_wrapper(h_c,ham_c,ovlp_c,eval_c,evec_c)&
   bind(C,name="c_elsi_ev_real")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: ham_f(:,:)
   real(kind=c_double), pointer :: ovlp_f(:,:)
   real(kind=c_double), pointer :: eval_f(:)
   real(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   n_basis = h_f%ph%n_basis
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_real(h_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine elsi_ev_complex_c_wrapper(h_c,ham_c,ovlp_c,eval_c,evec_c)&
   bind(C,name="c_elsi_ev_complex")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: h_f
   complex(kind=c_double_complex), pointer :: ham_f(:,:)
   complex(kind=c_double_complex), pointer :: ovlp_f(:,:)
   real(kind=c_double), pointer :: eval_f(:)
   complex(kind=c_double_complex), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   n_basis = h_f%ph%n_basis
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_complex(h_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine elsi_ev_real_sparse_c_wrapper(h_c,ham_c,ovlp_c,eval_c,evec_c)&
   bind(C,name="c_elsi_ev_real_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: ham_f(:)
   real(kind=c_double), pointer :: ovlp_f(:)
   real(kind=c_double), pointer :: eval_f(:)
   real(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   n_basis = h_f%ph%n_basis
   nnz_l = h_f%bh%nnz_l_sp
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_real_sparse(h_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine elsi_ev_complex_sparse_c_wrapper(h_c,ham_c,ovlp_c,eval_c,evec_c)&
   bind(C,name="c_elsi_ev_complex_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: h_f
   complex(kind=c_double), pointer :: ham_f(:)
   complex(kind=c_double), pointer :: ovlp_f(:)
   real(kind=c_double), pointer :: eval_f(:)
   complex(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   n_basis = h_f%ph%n_basis
   nnz_l = h_f%bh%nnz_l_sp
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_complex_sparse(h_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine elsi_dm_real_c_wrapper(h_c,ham_c,ovlp_c,dm_c,energy)&
   bind(C,name="c_elsi_dm_real")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: ham_f(:,:)
   real(kind=c_double), pointer :: ovlp_f(:,:)
   real(kind=c_double), pointer :: dm_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   n_basis = h_f%ph%n_basis
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_dm_real(h_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine elsi_dm_complex_c_wrapper(h_c,ham_c,ovlp_c,dm_c,energy)&
   bind(C,name="c_elsi_dm_complex")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle), pointer :: h_f
   complex(kind=c_double_complex), pointer :: ham_f(:,:)
   complex(kind=c_double_complex), pointer :: ovlp_f(:,:)
   complex(kind=c_double_complex), pointer :: dm_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   n_basis = h_f%ph%n_basis
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_dm_complex(h_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine elsi_dm_real_sparse_c_wrapper(h_c,ham_c,ovlp_c,dm_c,energy)&
   bind(C,name="c_elsi_dm_real_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: ham_f(:)
   real(kind=c_double), pointer :: ovlp_f(:)
   real(kind=c_double), pointer :: dm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(dm_c,dm_f,shape=[nnz_l])

   call elsi_dm_real_sparse(h_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine elsi_dm_complex_sparse_c_wrapper(h_c,ham_c,ovlp_c,dm_c,energy)&
   bind(C,name="c_elsi_dm_complex_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle), pointer :: h_f
   complex(kind=c_double), pointer :: ham_f(:)
   complex(kind=c_double), pointer :: ovlp_f(:)
   complex(kind=c_double), pointer :: dm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(dm_c,dm_f,shape=[nnz_l])

   call elsi_dm_complex_sparse(h_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine elsi_set_output_c_wrapper(h_c,output)&
   bind(C,name="c_elsi_set_output")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: output

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_output(h_f,output)

end subroutine

subroutine elsi_set_output_log_c_wrapper(h_c,output_log)&
   bind(C,name="c_elsi_set_output_log")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: output_log

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_output_log(h_f,output_log)

end subroutine

subroutine elsi_set_unit_ovlp_c_wrapper(h_c,unit_ovlp)&
   bind(C,name="c_elsi_set_unit_ovlp")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: unit_ovlp

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_unit_ovlp(h_f,unit_ovlp)

end subroutine

subroutine elsi_set_zero_def_c_wrapper(h_c,zero_def)&
   bind(C,name="c_elsi_set_zero_def")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: zero_def

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_zero_def(h_f,zero_def)

end subroutine

subroutine elsi_set_illcond_check_c_wrapper(h_c,illcond_check)&
   bind(C,name="c_elsi_set_illcond_check")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: illcond_check

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_illcond_check(h_f,illcond_check)

end subroutine

subroutine elsi_set_illcond_tol_c_wrapper(h_c,illcond_tol)&
   bind(C,name="c_elsi_set_illcond_tol")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: illcond_tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_illcond_tol(h_f,illcond_tol)

end subroutine

subroutine elsi_set_illcond_abort_c_wrapper(h_c,illcond_abort)&
   bind(C,name="c_elsi_set_illcond_abort")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: illcond_abort

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_illcond_abort(h_f,illcond_abort)

end subroutine

subroutine elsi_set_sing_check_c_wrapper(h_c,sing_check)&
   bind(C,name="c_elsi_set_sing_check")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: sing_check

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_illcond_check(h_f,sing_check)

end subroutine

subroutine elsi_set_sing_tol_c_wrapper(h_c,sing_tol)&
   bind(C,name="c_elsi_set_sing_tol")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: sing_tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_illcond_tol(h_f,sing_tol)

end subroutine

subroutine elsi_set_sing_stop_c_wrapper(h_c,sing_stop)&
   bind(C,name="c_elsi_set_sing_stop")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: sing_stop

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_illcond_abort(h_f,sing_stop)

end subroutine

subroutine elsi_set_elpa_solver_c_wrapper(h_c,solver)&
   bind(C,name="c_elsi_set_elpa_solver")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: solver

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_elpa_solver(h_f,solver)

end subroutine

subroutine elsi_set_elpa_n_single_c_wrapper(h_c,n_single)&
   bind(C,name="c_elsi_set_elpa_n_single")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_single

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_elpa_n_single(h_f,n_single)

end subroutine

subroutine elsi_set_elpa_gpu_c_wrapper(h_c,gpu)&
   bind(C,name="c_elsi_set_elpa_gpu")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: gpu

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_elpa_gpu(h_f,gpu)

end subroutine

subroutine elsi_set_elpa_gpu_kernels_c_wrapper(h_c,gpu_kernels)&
   bind(C,name="c_elsi_set_elpa_gpu_kernels")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: gpu_kernels

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_elpa_gpu_kernels(h_f,gpu_kernels)

end subroutine

subroutine elsi_set_elpa_autotune_c_wrapper(h_c,autotune)&
   bind(C,name="c_elsi_set_elpa_autotune")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: autotune

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_elpa_autotune(h_f,autotune)

end subroutine

subroutine elsi_set_omm_flavor_c_wrapper(h_c,flavor)&
   bind(C,name="c_elsi_set_omm_flavor")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: flavor

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_omm_flavor(h_f,flavor)

end subroutine

subroutine elsi_set_omm_n_elpa_c_wrapper(h_c,n_elpa)&
   bind(C,name="c_elsi_set_omm_n_elpa")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_elpa

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_omm_n_elpa(h_f,n_elpa)

end subroutine

subroutine elsi_set_omm_tol_c_wrapper(h_c,tol)&
   bind(C,name="c_elsi_set_omm_tol")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_omm_tol(h_f,tol)

end subroutine

subroutine elsi_set_pexsi_n_mu_c_wrapper(h_c,n_mu)&
   bind(C,name="c_elsi_set_pexsi_n_mu")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_mu

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_n_mu(h_f,n_mu)

end subroutine

subroutine elsi_set_pexsi_n_pole_c_wrapper(h_c,n_pole)&
   bind(C,name="c_elsi_set_pexsi_n_pole")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_pole

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_n_pole(h_f,n_pole)

end subroutine

subroutine elsi_set_pexsi_np_per_pole_c_wrapper(h_c,np_per_pole)&
   bind(C,name="c_elsi_set_pexsi_np_per_pole")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: np_per_pole

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_np_per_pole(h_f,np_per_pole)

end subroutine

subroutine elsi_set_pexsi_np_symbo_c_wrapper(h_c,np_symbo)&
   bind(C,name="c_elsi_set_pexsi_np_symbo")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: np_symbo

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_np_symbo(h_f,np_symbo)

end subroutine

subroutine elsi_set_pexsi_ordering_c_wrapper(h_c,ordering)&
   bind(C,name="c_elsi_set_pexsi_ordering")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: ordering

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_ordering(h_f,ordering)

end subroutine

subroutine elsi_set_pexsi_temp_c_wrapper(h_c,temp)&
   bind(C,name="c_elsi_set_pexsi_temp")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: temp

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_temp(h_f,temp)

end subroutine

subroutine elsi_set_pexsi_gap_c_wrapper(h_c,gap)&
   bind(C,name="c_elsi_set_pexsi_gap")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: gap

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_gap(h_f,gap)

end subroutine

subroutine elsi_set_pexsi_delta_e_c_wrapper(h_c,delta_e)&
   bind(C,name="c_elsi_set_pexsi_delta_e")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: delta_e

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_delta_e(h_f,delta_e)

end subroutine

subroutine elsi_set_pexsi_mu_min_c_wrapper(h_c,mu_min)&
   bind(C,name="c_elsi_set_pexsi_mu_min")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: mu_min

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_mu_min(h_f,mu_min)

end subroutine

subroutine elsi_set_pexsi_mu_max_c_wrapper(h_c,mu_max)&
   bind(C,name="c_elsi_set_pexsi_mu_max")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: mu_max

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_mu_max(h_f,mu_max)

end subroutine

subroutine elsi_set_pexsi_inertia_tol_c_wrapper(h_c,inertia_tol)&
   bind(C,name="c_elsi_set_pexsi_inertia_tol")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: inertia_tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_inertia_tol(h_f,inertia_tol)

end subroutine

subroutine elsi_set_sips_n_elpa_c_wrapper(h_c,n_elpa)&
   bind(C,name="c_elsi_set_sips_n_elpa")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_elpa

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_n_elpa(h_f,n_elpa)

end subroutine

subroutine elsi_set_sips_n_slice_c_wrapper(h_c,n_slice)&
   bind(C,name="c_elsi_set_sips_n_slice")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_slice

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_n_slice(h_f,n_slice)

end subroutine

subroutine elsi_set_sips_slice_type_c_wrapper(h_c,slice_type)&
   bind(C,name="c_elsi_set_sips_slice_type")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: slice_type

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_slice_type(h_f,slice_type)

end subroutine

subroutine elsi_set_sips_buffer_c_wrapper(h_c,buffer)&
   bind(C,name="c_elsi_set_sips_buffer")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: buffer

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_buffer(h_f,buffer)

end subroutine

subroutine elsi_set_sips_inertia_tol_c_wrapper(h_c,inertia_tol)&
   bind(C,name="c_elsi_set_sips_inertia_tol")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: inertia_tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_inertia_tol(h_f,inertia_tol)

end subroutine

subroutine elsi_set_sips_ev_min_c_wrapper(h_c,ev_min)&
   bind(C,name="c_elsi_set_sips_ev_min")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: ev_min

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_ev_min(h_f,ev_min)

end subroutine

subroutine elsi_set_sips_ev_max_c_wrapper(h_c,ev_max)&
   bind(C,name="c_elsi_set_sips_ev_max")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: ev_max

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_ev_max(h_f,ev_max)

end subroutine

subroutine elsi_set_sips_interval_c_wrapper(h_c,lower,upper)&
   bind(C,name="c_elsi_set_sips_interval")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: lower
   real(kind=c_double), value, intent(in) :: upper

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_ev_min(h_f,lower)
   call elsi_set_sips_ev_max(h_f,upper)

end subroutine

subroutine elsi_set_ntpoly_method_c_wrapper(h_c,method)&
   bind(C,name="c_elsi_set_ntpoly_method")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: method

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_ntpoly_method(h_f,method)

end subroutine

subroutine elsi_set_ntpoly_isr_c_wrapper(h_c,isr)&
   bind(C,name="c_elsi_set_ntpoly_isr")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: isr

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_ntpoly_isr(h_f,isr)

end subroutine

subroutine elsi_set_ntpoly_tol_c_wrapper(h_c,tol)&
   bind(C,name="c_elsi_set_ntpoly_tol")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_ntpoly_tol(h_f,tol)

end subroutine

subroutine elsi_set_ntpoly_filter_c_wrapper(h_c,filter)&
   bind(C,name="c_elsi_set_ntpoly_filter")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: filter

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_ntpoly_filter(h_f,filter)

end subroutine

subroutine elsi_set_ntpoly_max_iter_c_wrapper(h_c,max_iter)&
   bind(C,name="c_elsi_set_ntpoly_max_iter")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: max_iter

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_ntpoly_max_iter(h_f,max_iter)

end subroutine

subroutine elsi_set_mu_broaden_scheme_c_wrapper(h_c,broaden_scheme)&
   bind(C,name="c_elsi_set_mu_broaden_scheme")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: broaden_scheme

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mu_broaden_scheme(h_f,broaden_scheme)

end subroutine

subroutine elsi_set_mu_broaden_width_c_wrapper(h_c,broaden_width)&
   bind(C,name="c_elsi_set_mu_broaden_width")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: broaden_width

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mu_broaden_width(h_f,broaden_width)

end subroutine

subroutine elsi_set_mu_tol_c_wrapper(h_c,tol)&
   bind(C,name="c_elsi_set_mu_tol")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mu_tol(h_f,tol)

end subroutine

subroutine elsi_set_mu_spin_degen_c_wrapper(h_c,spin_degen)&
   bind(C,name="c_elsi_set_mu_spin_degen")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: spin_degen

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mu_spin_degen(h_f,spin_degen)

end subroutine

subroutine elsi_set_mu_mp_order_c_wrapper(h_c,mp_order)&
   bind(C,name="c_elsi_set_mu_mp_order")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: mp_order

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mu_mp_order(h_f,mp_order)

end subroutine

subroutine elsi_get_initialized_c_wrapper(h_c,initialized)&
   bind(C,name="c_elsi_get_initialized")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), intent(out) :: initialized

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_initialized(h_f,initialized)

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

subroutine elsi_get_n_illcond_c_wrapper(h_c,n_illcond)&
   bind(C,name="c_elsi_get_n_illcond")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), intent(out) :: n_illcond

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_n_illcond(h_f,n_illcond)

end subroutine

subroutine elsi_get_n_sing_c_wrapper(h_c,n_sing)&
   bind(C,name="c_elsi_get_n_sing")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), intent(out) :: n_sing

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_n_illcond(h_f,n_sing)

end subroutine

subroutine elsi_get_pexsi_mu_min_c_wrapper(h_c,mu_min)&
   bind(C,name="c_elsi_get_pexsi_mu_min")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), intent(out) :: mu_min

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_pexsi_mu_min(h_f,mu_min)

end subroutine

subroutine elsi_get_pexsi_mu_max_c_wrapper(h_c,mu_max)&
   bind(C,name="c_elsi_get_pexsi_mu_max")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), intent(out) :: mu_max

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_pexsi_mu_max(h_f,mu_max)

end subroutine

subroutine elsi_get_mu_c_wrapper(h_c,mu)&
   bind(C,name="c_elsi_get_mu")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), intent(out) :: mu

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_mu(h_f,mu)

end subroutine

subroutine elsi_get_entropy_c_wrapper(h_c,entropy)&
   bind(C,name="c_elsi_get_entropy")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), intent(out) :: entropy

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_entropy(h_f,entropy)

end subroutine

subroutine elsi_get_edm_real_c_wrapper(h_c,edm_c)&
   bind(C,name="c_elsi_get_edm_real")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: edm_f(:,:)

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(edm_c,edm_f,shape=[lrow,lcol])

   call elsi_get_edm_real(h_f,edm_f)

end subroutine

subroutine elsi_get_edm_complex_c_wrapper(h_c,edm_c)&
   bind(C,name="c_elsi_get_edm_complex")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle), pointer :: h_f
   complex(kind=c_double), pointer :: edm_f(:,:)

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(edm_c,edm_f,shape=[lrow,lcol])

   call elsi_get_edm_complex(h_f,edm_f)

end subroutine

subroutine elsi_get_edm_real_sparse_c_wrapper(h_c,edm_c)&
   bind(C,name="c_elsi_get_edm_real_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: edm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp

   call c_f_pointer(edm_c,edm_f,shape=[nnz_l])

   call elsi_get_edm_real_sparse(h_f,edm_f)

end subroutine

subroutine elsi_get_edm_complex_sparse_c_wrapper(h_c,edm_c)&
   bind(C,name="c_elsi_get_edm_complex_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle), pointer :: h_f
   complex(kind=c_double), pointer :: edm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp

   call c_f_pointer(edm_c,edm_f,shape=[nnz_l])

   call elsi_get_edm_complex_sparse(h_f,edm_f)

end subroutine

subroutine elsi_orthonormalize_ev_real_c_wrapper(h_c,ovlp_c,evec_c)&
   bind(C,name="c_elsi_orthonormalize_ev_real")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: ovlp_f(:,:)
   real(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_orthonormalize_ev_real(h_f,ovlp_f,evec_f)

end subroutine

subroutine elsi_orthonormalize_ev_complex_c_wrapper(h_c,ovlp_c,evec_c)&
   bind(C,name="c_elsi_orthonormalize_ev_complex")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: h_f
   complex(kind=c_double), pointer :: ovlp_f(:,:)
   complex(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_orthonormalize_ev_complex(h_f,ovlp_f,evec_f)

end subroutine

subroutine elsi_extrapolate_dm_real_c_wrapper(h_c,ovlp_c,dm_c)&
   bind(C,name="c_elsi_extrapolate_dm_real")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: ovlp_f(:,:)
   real(kind=c_double), pointer :: dm_f(:,:)

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_extrapolate_dm_real(h_f,ovlp_f,dm_f)

end subroutine

subroutine elsi_extrapolate_dm_complex_c_wrapper(h_c,ovlp_c,dm_c)&
   bind(C,name="c_elsi_extrapolate_dm_complex")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c

   type(elsi_handle), pointer :: h_f
   complex(kind=c_double), pointer :: ovlp_f(:,:)
   complex(kind=c_double), pointer :: dm_f(:,:)

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_extrapolate_dm_complex(h_f,ovlp_f,dm_f)

end subroutine

subroutine elsi_construct_dm_real_c_wrapper(h_c,occ_c,evec_c,dm_c)&
   bind(C,name="c_elsi_construct_dm_real")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: occ_c
   type(c_ptr), value, intent(in) :: evec_c
   type(c_ptr), value, intent(in) :: dm_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: occ_f(:)
   real(kind=c_double), pointer :: evec_f(:,:)
   real(kind=c_double), pointer :: dm_f(:,:)

   integer(kind=c_int) :: n_state
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   n_state = h_f%ph%n_states
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(occ_c,occ_f,shape=[n_state])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_construct_dm_real(h_f,occ_f,evec_f,dm_f)

end subroutine

subroutine elsi_construct_dm_complex_c_wrapper(h_c,occ_c,evec_c,dm_c)&
   bind(C,name="c_elsi_construct_dm_complex")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: occ_c
   type(c_ptr), value, intent(in) :: evec_c
   type(c_ptr), value, intent(in) :: dm_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: occ_f(:)
   complex(kind=c_double), pointer :: evec_f(:,:)
   complex(kind=c_double), pointer :: dm_f(:,:)

   integer(kind=c_int) :: n_state
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   n_state = h_f%ph%n_states
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(occ_c,occ_f,shape=[n_state])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_construct_dm_complex(h_f,occ_f,evec_f,dm_f)

end subroutine

subroutine elsi_construct_edm_real_c_wrapper(h_c,occ_c,eval_c,evec_c,edm_c)&
   bind(C,name="c_elsi_construct_edm_real")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: occ_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: occ_f(:)
   real(kind=c_double), pointer :: eval_f(:)
   real(kind=c_double), pointer :: evec_f(:,:)
   real(kind=c_double), pointer :: edm_f(:,:)

   integer(kind=c_int) :: n_state
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   n_state = h_f%ph%n_states
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(occ_c,occ_f,shape=[n_state])
   call c_f_pointer(eval_c,eval_f,shape=[n_state])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])
   call c_f_pointer(edm_c,edm_f,shape=[lrow,lcol])

   call elsi_construct_edm_real(h_f,occ_f,eval_f,evec_f,edm_f)

end subroutine

subroutine elsi_construct_edm_complex_c_wrapper(h_c,occ_c,eval_c,evec_c,edm_c)&
   bind(C,name="c_elsi_construct_edm_complex")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: occ_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: occ_f(:)
   real(kind=c_double), pointer :: eval_f(:)
   complex(kind=c_double), pointer :: evec_f(:,:)
   complex(kind=c_double), pointer :: edm_f(:,:)

   integer(kind=c_int) :: n_state
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   n_state = h_f%ph%n_states
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(occ_c,occ_f,shape=[n_state])
   call c_f_pointer(eval_c,eval_f,shape=[n_state])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])
   call c_f_pointer(edm_c,edm_f,shape=[lrow,lcol])

   call elsi_construct_edm_complex(h_f,occ_f,eval_f,evec_f,edm_f)

end subroutine

subroutine elsi_init_rw_c_wrapper(h_c,rw_task,parallel_mode,n_basis,n_electron)&
   bind(C,name="c_elsi_init_rw")

   implicit none

   type(c_ptr) :: h_c
   integer(kind=c_int), value, intent(in) :: rw_task
   integer(kind=c_int), value, intent(in) :: parallel_mode
   integer(kind=c_int), value, intent(in) :: n_basis
   real(kind=c_double), value, intent(in) :: n_electron

   type(elsi_rw_handle), pointer :: h_f

   allocate(h_f)

   call elsi_init_rw(h_f,rw_task,parallel_mode,n_basis,n_electron)

   h_c = c_loc(h_f)

end subroutine

subroutine elsi_set_rw_mpi_c_wrapper(h_c,mpi_comm)&
   bind(C,name="c_elsi_set_rw_mpi")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: mpi_comm

   type(elsi_rw_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_rw_mpi(h_f,mpi_comm)

end subroutine

subroutine elsi_set_rw_blacs_c_wrapper(h_c,blacs_ctxt,block_size)&
   bind(C,name="c_elsi_set_rw_blacs")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: blacs_ctxt
   integer(kind=c_int), value, intent(in) :: block_size

   type(elsi_rw_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_rw_blacs(h_f,blacs_ctxt,block_size)

end subroutine

subroutine elsi_set_rw_csc_c_wrapper(h_c,nnz_g,nnz_l,n_lcol)&
   bind(C,name="c_elsi_set_rw_csc")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: nnz_g
   integer(kind=c_int), value, intent(in) :: nnz_l
   integer(kind=c_int), value, intent(in) :: n_lcol

   type(elsi_rw_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_rw_csc(h_f,nnz_g,nnz_l,n_lcol)

end subroutine

subroutine elsi_finalize_rw_c_wrapper(h_c)&
   bind(C,name="c_elsi_finalize_rw")

   implicit none

   type(c_ptr), value, intent(in) :: h_c

   type(elsi_rw_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_finalize_rw(h_f)

end subroutine

subroutine elsi_set_rw_zero_def_c_wrapper(h_c,zero_def)&
   bind(C,name="c_elsi_set_rw_zero_def")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: zero_def

   type(elsi_rw_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_rw_zero_def(h_f,zero_def)

end subroutine

subroutine elsi_read_mat_dim_c_wrapper(h_c,name_c,n_electrons,n_basis,n_lrow,&
   n_lcol)&
   bind(C,name="c_elsi_read_mat_dim")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1) :: name_c(128)
   real(kind=c_double), intent(out) :: n_electrons
   integer(kind=c_int), intent(out) :: n_basis
   integer(kind=c_int), intent(out) :: n_lrow
   integer(kind=c_int), intent(out) :: n_lcol

   type(elsi_rw_handle), pointer :: h_f

   character(len=:), allocatable :: name_f

   call c_f_pointer(h_c,h_f)

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_dim(h_f,name_f,n_electrons,n_basis,n_lrow,n_lcol)

end subroutine

subroutine elsi_read_mat_dim_sparse_c_wrapper(h_c,name_c,n_electrons,n_basis,&
   nnz_g,nnz_l,n_lcol)&
   bind(C,name="c_elsi_read_mat_dim_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1) :: name_c(128)
   real(kind=c_double), intent(out) :: n_electrons
   integer(kind=c_int), intent(out) :: n_basis
   integer(kind=c_int), intent(out) :: nnz_g
   integer(kind=c_int), intent(out) :: nnz_l
   integer(kind=c_int), intent(out) :: n_lcol

   type(elsi_rw_handle), pointer :: h_f

   character(len=:), allocatable :: name_f

   call c_f_pointer(h_c,h_f)

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_dim_sparse(h_f,name_f,n_electrons,n_basis,nnz_g,nnz_l,&
        n_lcol)

end subroutine

subroutine elsi_read_mat_real_c_wrapper(h_c,name_c,mat_c)&
   bind(C,name="c_elsi_read_mat_real")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: h_f
   real(kind=c_double), pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_real(h_f,name_f,mat_f)

end subroutine

subroutine elsi_read_mat_real_sparse_c_wrapper(h_c,name_c,row_ind_c,col_ptr_c,&
   mat_c)&
   bind(C,name="c_elsi_read_mat_real_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: row_ind_c
   type(c_ptr), value, intent(in) :: col_ptr_c
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: h_f
   integer(kind=c_int), pointer :: row_ind_f(:)
   integer(kind=c_int), pointer :: col_ptr_f(:)
   real(kind=c_double), pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp
   lcol = h_f%bh%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_real_sparse(h_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine elsi_write_mat_real_c_wrapper(h_c,name_c,mat_c)&
   bind(C,name="c_elsi_write_mat_real")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: h_f
   real(kind=c_double), pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_real(h_f,name_f,mat_f)

end subroutine

subroutine elsi_write_mat_real_sparse_c_wrapper(h_c,name_c,row_ind_c,col_ptr_c,&
   mat_c)&
   bind(C,name="c_elsi_write_mat_real_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: row_ind_c
   type(c_ptr), value, intent(in) :: col_ptr_c
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: h_f
   integer(kind=c_int), pointer :: row_ind_f(:)
   integer(kind=c_int), pointer :: col_ptr_f(:)
   real(kind=c_double), pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp
   lcol = h_f%bh%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_real_sparse(h_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine elsi_read_mat_complex_c_wrapper(h_c,name_c,mat_c)&
   bind(C,name="c_elsi_read_mat_complex")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: h_f
   complex(kind=c_double), pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_complex(h_f,name_f,mat_f)

end subroutine

subroutine elsi_read_mat_complex_sparse_c_wrapper(h_c,name_c,row_ind_c,&
   col_ptr_c,mat_c)&
   bind(C,name="c_elsi_read_mat_complex_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: row_ind_c
   type(c_ptr), value, intent(in) :: col_ptr_c
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: h_f
   integer(kind=c_int), pointer :: row_ind_f(:)
   integer(kind=c_int), pointer :: col_ptr_f(:)
   complex(kind=c_double), pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp
   lcol = h_f%bh%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_complex_sparse(h_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine elsi_write_mat_complex_c_wrapper(h_c,name_c,mat_c)&
   bind(C,name="c_elsi_write_mat_complex")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: h_f
   complex(kind=c_double), pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_complex(h_f,name_f,mat_f)

end subroutine

subroutine elsi_write_mat_complex_sparse_c_wrapper(h_c,name_c,row_ind_c,&
   col_ptr_c,mat_c)&
   bind(C,name="c_elsi_write_mat_complex_sparse")

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1) :: name_c(128)
   type(c_ptr), value, intent(in) :: row_ind_c
   type(c_ptr), value, intent(in) :: col_ptr_c
   type(c_ptr), value, intent(in) :: mat_c

   type(elsi_rw_handle), pointer :: h_f
   integer(kind=c_int), pointer :: row_ind_f(:)
   integer(kind=c_int), pointer :: col_ptr_f(:)
   complex(kind=c_double), pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp
   lcol = h_f%bh%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_complex_sparse(h_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

end module ELSI_C_INTERFACE
