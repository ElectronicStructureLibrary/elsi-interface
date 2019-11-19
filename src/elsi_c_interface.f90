! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide C interfaces of ELSI.
!!
module ELSI_C_INTERFACE

   use, intrinsic :: ISO_C_BINDING
   use ELSI

   implicit none

   private

contains

!>
!! Convert a C string into a Fortran string. A Fortran string is NOT a character
!! array without a NULL character.
!!
subroutine str_c2f(str_c,str_f)

   implicit none

   character(kind=c_char,len=1), intent(in) :: str_c(*)
   character(len=:), allocatable, intent(out) :: str_f

   integer(kind=c_int) :: str_f_len

   str_f_len = 0

   do
     if(str_c(str_f_len+1) == C_NULL_CHAR) then
        exit
     end if

     str_f_len = str_f_len+1
   end do

   allocate(character(len=str_f_len) :: str_f)

   str_f = transfer(str_c(1:str_f_len),str_f)

end subroutine

subroutine c_elsi_init(h_c,solver,parallel_mode,matrix_format,n_basis,&
   n_electron,n_state) bind(C)

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

subroutine c_elsi_set_mpi(h_c,mpi_comm) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: mpi_comm

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mpi(h_f,mpi_comm)

end subroutine

subroutine c_elsi_set_mpi_global(h_c,mpi_comm_global) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: mpi_comm_global

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mpi_global(h_f,mpi_comm_global)

end subroutine

subroutine c_elsi_set_spin(h_c,n_spin,i_spin) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_spin
   integer(kind=c_int), value, intent(in) :: i_spin

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_spin(h_f,n_spin,i_spin)

end subroutine

subroutine c_elsi_set_kpoint(h_c,n_kpt,i_kpt,i_wt) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_kpt
   integer(kind=c_int), value, intent(in) :: i_kpt
   real(kind=c_double), value, intent(in) :: i_wt

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_kpoint(h_f,n_kpt,i_kpt,i_wt)

end subroutine

subroutine c_elsi_set_blacs(h_c,blacs_ctxt,block_size) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: blacs_ctxt
   integer(kind=c_int), value, intent(in) :: block_size

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_blacs(h_f,blacs_ctxt,block_size)

end subroutine

subroutine c_elsi_set_csc(h_c,nnz_g,nnz_l,n_lcol,row_ind,col_ptr) bind(C)

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

subroutine c_elsi_set_csc_blk(h_c,blk) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: blk

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_csc_blk(h_f,blk)

end subroutine

subroutine c_elsi_set_coo(h_c,nnz_g,nnz_l,row_ind,col_ind) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: nnz_g
   integer(kind=c_int), value, intent(in) :: nnz_l
   integer(kind=c_int), intent(in) :: row_ind(nnz_l)
   integer(kind=c_int), intent(in) :: col_ind(nnz_l)

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_coo(h_f,nnz_g,nnz_l,row_ind,col_ind)

end subroutine

subroutine c_elsi_reinit(h_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_reinit(h_f)

end subroutine

subroutine c_elsi_finalize(h_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_finalize(h_f)

end subroutine

subroutine c_elsi_ev_real(h_c,ham_c,ovlp_c,eval_c,evec_c) bind(C)

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

subroutine c_elsi_ev_complex(h_c,ham_c,ovlp_c,eval_c,evec_c) bind(C)

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

subroutine c_elsi_ev_real_sparse(h_c,ham_c,ovlp_c,eval_c,evec_c) bind(C)

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

subroutine c_elsi_ev_complex_sparse(h_c,ham_c,ovlp_c,eval_c,evec_c) bind(C)

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

subroutine c_elsi_dm_real(h_c,ham_c,ovlp_c,dm_c,energy) bind(C)

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

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_dm_real(h_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine c_elsi_dm_complex(h_c,ham_c,ovlp_c,dm_c,energy) bind(C)

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

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_dm_complex(h_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine c_elsi_dm_real_sparse(h_c,ham_c,ovlp_c,dm_c,energy) bind(C)

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

subroutine c_elsi_dm_complex_sparse(h_c,ham_c,ovlp_c,dm_c,energy) bind(C)

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

subroutine c_elsi_set_input_file(h_c,name_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1), intent(in) :: name_c(*)

   type(elsi_handle), pointer :: h_f

   character(len=:), allocatable :: name_f

   call c_f_pointer(h_c,h_f)

   call str_c2f(name_c,name_f)

   call elsi_set_input_file(h_f,name_f)

end subroutine

subroutine c_elsi_set_output(h_c,output) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: output

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_output(h_f,output)

end subroutine

subroutine c_elsi_set_output_log(h_c,output_log) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: output_log

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_output_log(h_f,output_log)

end subroutine

subroutine c_elsi_set_save_ovlp(h_c,save_ovlp) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: save_ovlp

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_save_ovlp(h_f,save_ovlp)

end subroutine

subroutine c_elsi_set_unit_ovlp(h_c,unit_ovlp) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: unit_ovlp

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_unit_ovlp(h_f,unit_ovlp)

end subroutine

subroutine c_elsi_set_zero_def(h_c,zero_def) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: zero_def

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_zero_def(h_f,zero_def)

end subroutine

subroutine c_elsi_set_illcond_check(h_c,illcond_check) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: illcond_check

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_illcond_check(h_f,illcond_check)

end subroutine

subroutine c_elsi_set_illcond_tol(h_c,illcond_tol) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: illcond_tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_illcond_tol(h_f,illcond_tol)

end subroutine

subroutine c_elsi_set_sing_check(h_c,sing_check) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: sing_check

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_illcond_check(h_f,sing_check)

end subroutine

subroutine c_elsi_set_sing_tol(h_c,sing_tol) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: sing_tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_illcond_tol(h_f,sing_tol)

end subroutine

subroutine c_elsi_set_spin_degeneracy(h_c,spin_degeneracy) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: spin_degeneracy

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_spin_degeneracy(h_f,spin_degeneracy)

end subroutine

subroutine c_elsi_set_energy_gap(h_c,gap) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: gap

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_energy_gap(h_f,gap)

end subroutine

subroutine c_elsi_set_spectrum_width(h_c,width) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: width

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_spectrum_width(h_f,width)

end subroutine

subroutine c_elsi_set_dimensionality(h_c,dimensionality) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: dimensionality

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_dimensionality(h_f,dimensionality)

end subroutine

subroutine c_elsi_set_elpa_solver(h_c,solver) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: solver

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_elpa_solver(h_f,solver)

end subroutine

subroutine c_elsi_set_elpa_n_single(h_c,n_single) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_single

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_elpa_n_single(h_f,n_single)

end subroutine

subroutine c_elsi_set_elpa_gpu(h_c,gpu) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: gpu

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_elpa_gpu(h_f,gpu)

end subroutine

subroutine c_elsi_set_elpa_gpu_kernels(h_c,gpu_kernels) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: gpu_kernels

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_elpa_gpu_kernels(h_f,gpu_kernels)

end subroutine

subroutine c_elsi_set_elpa_autotune(h_c,autotune) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: autotune

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_elpa_autotune(h_f,autotune)

end subroutine

subroutine c_elsi_set_omm_flavor(h_c,flavor) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: flavor

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_omm_flavor(h_f,flavor)

end subroutine

subroutine c_elsi_set_omm_n_elpa(h_c,n_elpa) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_elpa

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_omm_n_elpa(h_f,n_elpa)

end subroutine

subroutine c_elsi_set_omm_tol(h_c,tol) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_omm_tol(h_f,tol)

end subroutine

subroutine c_elsi_set_pexsi_n_mu(h_c,n_mu) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_mu

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_n_mu(h_f,n_mu)

end subroutine

subroutine c_elsi_set_pexsi_n_pole(h_c,n_pole) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_pole

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_n_pole(h_f,n_pole)

end subroutine

subroutine c_elsi_set_pexsi_np_per_pole(h_c,np_per_pole) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: np_per_pole

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_np_per_pole(h_f,np_per_pole)

end subroutine

subroutine c_elsi_set_pexsi_np_symbo(h_c,np_symbo) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: np_symbo

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_np_symbo(h_f,np_symbo)

end subroutine

subroutine c_elsi_set_pexsi_ordering(h_c,ordering) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: ordering

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_ordering(h_f,ordering)

end subroutine

subroutine c_elsi_set_pexsi_temp(h_c,temp) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: temp

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_temp(h_f,temp)

end subroutine

subroutine c_elsi_set_pexsi_gap(h_c,gap) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: gap

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_gap(h_f,gap)

end subroutine

subroutine c_elsi_set_pexsi_delta_e(h_c,delta_e) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: delta_e

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_delta_e(h_f,delta_e)

end subroutine

subroutine c_elsi_set_pexsi_mu_min(h_c,mu_min) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: mu_min

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_mu_min(h_f,mu_min)

end subroutine

subroutine c_elsi_set_pexsi_mu_max(h_c,mu_max) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: mu_max

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_mu_max(h_f,mu_max)

end subroutine

subroutine c_elsi_set_pexsi_inertia_tol(h_c,inertia_tol) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: inertia_tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_pexsi_inertia_tol(h_f,inertia_tol)

end subroutine

subroutine c_elsi_set_eigenexa_method(h_c,method) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: method

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_eigenexa_method(h_f,method)

end subroutine

subroutine c_elsi_set_sips_n_elpa(h_c,n_elpa) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_elpa

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_n_elpa(h_f,n_elpa)

end subroutine

subroutine c_elsi_set_sips_n_slice(h_c,n_slice) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_slice

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_n_slice(h_f,n_slice)

end subroutine

subroutine c_elsi_set_sips_slice_type(h_c,slice_type) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: slice_type

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_slice_type(h_f,slice_type)

end subroutine

subroutine c_elsi_set_sips_buffer(h_c,buffer) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: buffer

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_buffer(h_f,buffer)

end subroutine

subroutine c_elsi_set_sips_inertia_tol(h_c,inertia_tol) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: inertia_tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_inertia_tol(h_f,inertia_tol)

end subroutine

subroutine c_elsi_set_sips_ev_min(h_c,ev_min) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: ev_min

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_ev_min(h_f,ev_min)

end subroutine

subroutine c_elsi_set_sips_ev_max(h_c,ev_max) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: ev_max

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_ev_max(h_f,ev_max)

end subroutine

subroutine c_elsi_set_sips_interval(h_c,lower,upper) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: lower
   real(kind=c_double), value, intent(in) :: upper

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_sips_ev_min(h_f,lower)
   call elsi_set_sips_ev_max(h_f,upper)

end subroutine

subroutine c_elsi_set_ntpoly_method(h_c,method) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: method

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_ntpoly_method(h_f,method)

end subroutine

subroutine c_elsi_set_ntpoly_isr(h_c,isr) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: isr

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_ntpoly_isr(h_f,isr)

end subroutine

subroutine c_elsi_set_ntpoly_tol(h_c,tol) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_ntpoly_tol(h_f,tol)

end subroutine

subroutine c_elsi_set_ntpoly_filter(h_c,filter) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: filter

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_ntpoly_filter(h_f,filter)

end subroutine

subroutine c_elsi_set_ntpoly_max_iter(h_c,max_iter) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: max_iter

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_ntpoly_max_iter(h_f,max_iter)

end subroutine

subroutine c_elsi_set_ntpoly_n_layer(h_c,n_layer) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: n_layer

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_ntpoly_n_layer(h_f,n_layer)

end subroutine

subroutine c_elsi_set_magma_solver(h_c,solver) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: solver

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_magma_solver(h_f,solver)

end subroutine

subroutine c_elsi_set_mu_broaden_scheme(h_c,broaden_scheme) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: broaden_scheme

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mu_broaden_scheme(h_f,broaden_scheme)

end subroutine

subroutine c_elsi_set_mu_broaden_width(h_c,broaden_width) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: broaden_width

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mu_broaden_width(h_f,broaden_width)

end subroutine

subroutine c_elsi_set_mu_tol(h_c,tol) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: tol

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mu_tol(h_f,tol)

end subroutine

subroutine c_elsi_set_mu_mp_order(h_c,mp_order) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: mp_order

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_mu_mp_order(h_f,mp_order)

end subroutine

subroutine c_elsi_get_initialized(h_c,initialized) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), intent(out) :: initialized

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_initialized(h_f,initialized)

end subroutine

subroutine c_elsi_get_version(major,minor,patch) bind(C)

   implicit none

   integer(kind=c_int), intent(out) :: major
   integer(kind=c_int), intent(out) :: minor
   integer(kind=c_int), intent(out) :: patch

   call elsi_get_version(major,minor,patch)

end subroutine

subroutine c_elsi_get_datestamp(datestamp) bind(C)

   implicit none

   integer(kind=c_int), intent(out) :: datestamp

   call elsi_get_datestamp(datestamp)

end subroutine

subroutine c_elsi_get_n_illcond(h_c,n_illcond) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), intent(out) :: n_illcond

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_n_illcond(h_f,n_illcond)

end subroutine

subroutine c_elsi_get_n_sing(h_c,n_sing) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), intent(out) :: n_sing

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_n_illcond(h_f,n_sing)

end subroutine

subroutine c_elsi_get_ovlp_ev_min(h_c,ev_min) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), intent(out) :: ev_min

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_ovlp_ev_min(h_f,ev_min)

end subroutine

subroutine c_elsi_get_ovlp_ev_max(h_c,ev_max) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), intent(out) :: ev_max

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_ovlp_ev_max(h_f,ev_max)

end subroutine

subroutine c_elsi_get_pexsi_mu_min(h_c,mu_min) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), intent(out) :: mu_min

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_pexsi_mu_min(h_f,mu_min)

end subroutine

subroutine c_elsi_get_pexsi_mu_max(h_c,mu_max) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), intent(out) :: mu_max

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_pexsi_mu_max(h_f,mu_max)

end subroutine

subroutine c_elsi_get_mu(h_c,mu) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), intent(out) :: mu

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_mu(h_f,mu)

end subroutine

subroutine c_elsi_get_entropy(h_c,entropy) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), intent(out) :: entropy

   type(elsi_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_get_entropy(h_f,entropy)

end subroutine

subroutine c_elsi_get_edm_real(h_c,edm_c) bind(C)

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

subroutine c_elsi_get_edm_complex(h_c,edm_c) bind(C)

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

subroutine c_elsi_get_edm_real_sparse(h_c,edm_c) bind(C)

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

subroutine c_elsi_get_edm_complex_sparse(h_c,edm_c) bind(C)

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

subroutine c_elsi_orthonormalize_ev_real(h_c,ovlp_c,evec_c) bind(C)

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

subroutine c_elsi_orthonormalize_ev_complex(h_c,ovlp_c,evec_c) bind(C)

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

subroutine c_elsi_orthonormalize_ev_real_sparse(h_c,ovlp_c,evec_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: ovlp_f(:)
   real(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_orthonormalize_ev_real_sparse(h_f,ovlp_f,evec_f)

end subroutine

subroutine c_elsi_orthonormalize_ev_complex_sparse(h_c,ovlp_c,evec_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle), pointer :: h_f
   complex(kind=c_double), pointer :: ovlp_f(:)
   complex(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp
   lrow = h_f%bh%n_lrow
   lcol = h_f%bh%n_lcol

   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_orthonormalize_ev_complex_sparse(h_f,ovlp_f,evec_f)

end subroutine

subroutine c_elsi_extrapolate_dm_real(h_c,ovlp_c,dm_c) bind(C)

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

subroutine c_elsi_extrapolate_dm_complex(h_c,ovlp_c,dm_c) bind(C)

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

subroutine c_elsi_extrapolate_dm_real_sparse(h_c,ovlp_c,dm_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c

   type(elsi_handle), pointer :: h_f
   real(kind=c_double), pointer :: ovlp_f(:)
   real(kind=c_double), pointer :: dm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp

   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(dm_c,dm_f,shape=[nnz_l])

   call elsi_extrapolate_dm_real_sparse(h_f,ovlp_f,dm_f)

end subroutine

subroutine c_elsi_extrapolate_dm_complex_sparse(h_c,ovlp_c,dm_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: dm_c

   type(elsi_handle), pointer :: h_f
   complex(kind=c_double), pointer :: ovlp_f(:)
   complex(kind=c_double), pointer :: dm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(h_c,h_f)

   nnz_l = h_f%bh%nnz_l_sp

   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(dm_c,dm_f,shape=[nnz_l])

   call elsi_extrapolate_dm_complex_sparse(h_f,ovlp_f,dm_f)

end subroutine

subroutine c_elsi_init_rw(h_c,rw_task,parallel_mode,n_basis,n_electron) bind(C)

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

subroutine c_elsi_set_rw_mpi(h_c,mpi_comm) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: mpi_comm

   type(elsi_rw_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_rw_mpi(h_f,mpi_comm)

end subroutine

subroutine c_elsi_set_rw_blacs(h_c,blacs_ctxt,block_size) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: blacs_ctxt
   integer(kind=c_int), value, intent(in) :: block_size

   type(elsi_rw_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_rw_blacs(h_f,blacs_ctxt,block_size)

end subroutine

subroutine c_elsi_set_rw_csc(h_c,nnz_g,nnz_l,n_lcol) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   integer(kind=c_int), value, intent(in) :: nnz_g
   integer(kind=c_int), value, intent(in) :: nnz_l
   integer(kind=c_int), value, intent(in) :: n_lcol

   type(elsi_rw_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_rw_csc(h_f,nnz_g,nnz_l,n_lcol)

end subroutine

subroutine c_elsi_finalize_rw(h_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c

   type(elsi_rw_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_finalize_rw(h_f)

end subroutine

subroutine c_elsi_set_rw_zero_def(h_c,zero_def) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   real(kind=c_double), value, intent(in) :: zero_def

   type(elsi_rw_handle), pointer :: h_f

   call c_f_pointer(h_c,h_f)

   call elsi_set_rw_zero_def(h_f,zero_def)

end subroutine

subroutine c_elsi_read_mat_dim(h_c,name_c,n_electrons,n_basis,n_lrow,n_lcol)&
   bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1), intent(in) :: name_c(*)
   real(kind=c_double), intent(out) :: n_electrons
   integer(kind=c_int), intent(out) :: n_basis
   integer(kind=c_int), intent(out) :: n_lrow
   integer(kind=c_int), intent(out) :: n_lcol

   type(elsi_rw_handle), pointer :: h_f

   character(len=:), allocatable :: name_f

   call c_f_pointer(h_c,h_f)

   call str_c2f(name_c,name_f)

   call elsi_read_mat_dim(h_f,name_f,n_electrons,n_basis,n_lrow,n_lcol)

end subroutine

subroutine c_elsi_read_mat_dim_sparse(h_c,name_c,n_electrons,n_basis,nnz_g,&
   nnz_l,n_lcol) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1), intent(in) :: name_c(*)
   real(kind=c_double), intent(out) :: n_electrons
   integer(kind=c_int), intent(out) :: n_basis
   integer(kind=c_int), intent(out) :: nnz_g
   integer(kind=c_int), intent(out) :: nnz_l
   integer(kind=c_int), intent(out) :: n_lcol

   type(elsi_rw_handle), pointer :: h_f

   character(len=:), allocatable :: name_f

   call c_f_pointer(h_c,h_f)

   call str_c2f(name_c,name_f)

   call elsi_read_mat_dim_sparse(h_f,name_f,n_electrons,n_basis,nnz_g,nnz_l,&
        n_lcol)

end subroutine

subroutine c_elsi_read_mat_real(h_c,name_c,mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1), intent(in) :: name_c(*)
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

   call str_c2f(name_c,name_f)

   call elsi_read_mat_real(h_f,name_f,mat_f)

end subroutine

subroutine c_elsi_read_mat_real_sparse(h_c,name_c,row_ind_c,col_ptr_c,mat_c)&
   bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1), intent(in) :: name_c(*)
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

   call str_c2f(name_c,name_f)

   call elsi_read_mat_real_sparse(h_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine c_elsi_write_mat_real(h_c,name_c,mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1), intent(in) :: name_c(*)
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

   call str_c2f(name_c,name_f)

   call elsi_write_mat_real(h_f,name_f,mat_f)

end subroutine

subroutine c_elsi_write_mat_real_sparse(h_c,name_c,row_ind_c,col_ptr_c,mat_c)&
   bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1), intent(in) :: name_c(*)
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

   call str_c2f(name_c,name_f)

   call elsi_write_mat_real_sparse(h_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine c_elsi_read_mat_complex(h_c,name_c,mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1), intent(in) :: name_c(*)
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

   call str_c2f(name_c,name_f)

   call elsi_read_mat_complex(h_f,name_f,mat_f)

end subroutine

subroutine c_elsi_read_mat_complex_sparse(h_c,name_c,row_ind_c,col_ptr_c,mat_c)&
   bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1), intent(in) :: name_c(*)
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

   call str_c2f(name_c,name_f)

   call elsi_read_mat_complex_sparse(h_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine c_elsi_write_mat_complex(h_c,name_c,mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1), intent(in) :: name_c(*)
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

   call str_c2f(name_c,name_f)

   call elsi_write_mat_complex(h_f,name_f,mat_f)

end subroutine

subroutine c_elsi_write_mat_complex_sparse(h_c,name_c,row_ind_c,col_ptr_c,&
   mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: h_c
   character(kind=c_char,len=1), intent(in) :: name_c(*)
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

   call str_c2f(name_c,name_f)

   call elsi_write_mat_complex_sparse(h_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

end module ELSI_C_INTERFACE
