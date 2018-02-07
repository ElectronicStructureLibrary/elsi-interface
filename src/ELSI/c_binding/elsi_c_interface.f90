! Copyright (c) 2015-2018, the ELSI team. All rights reserved.
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

contains

subroutine c_elsi_init(handle_c,solver,parallel_mode,matrix_format,n_basis,&
              n_electron,n_state) bind(C)

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

   call elsi_init(handle_f,solver,parallel_mode,matrix_format,n_basis,&
           n_electron,n_state)

   handle_c = c_loc(handle_f)

end subroutine

subroutine c_elsi_set_mpi(handle_c,mpi_comm) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: mpi_comm

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mpi(handle_f,mpi_comm)

end subroutine

subroutine c_elsi_set_mpi_global(handle_c,mpi_comm_global) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: mpi_comm_global

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mpi_global(handle_f,mpi_comm_global)

end subroutine

subroutine c_elsi_set_spin(handle_c,n_spin,i_spin) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_spin
   integer(kind=c_int), value, intent(in) :: i_spin

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_spin(handle_f,n_spin,i_spin)

end subroutine

subroutine c_elsi_set_kpoint(handle_c,n_kpt,i_kpt,weight) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_kpt
   integer(kind=c_int), value, intent(in) :: i_kpt
   real(kind=c_double), value, intent(in) :: weight

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_kpoint(handle_f,n_kpt,i_kpt,weight)

end subroutine

subroutine c_elsi_set_blacs(handle_c,blacs_ctxt,block_size) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: blacs_ctxt
   integer(kind=c_int), value, intent(in) :: block_size

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_blacs(handle_f,blacs_ctxt,block_size)

end subroutine

subroutine c_elsi_set_csc(handle_c,nnz_g,nnz_l,n_lcol,row_ind,col_ptr) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: nnz_g
   integer(kind=c_int), value, intent(in) :: nnz_l
   integer(kind=c_int), value, intent(in) :: n_lcol
   integer(kind=c_int),        intent(in) :: row_ind(nnz_l)
   integer(kind=c_int),        intent(in) :: col_ptr(n_lcol+1)

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_csc(handle_f,nnz_g,nnz_l,n_lcol,row_ind,col_ptr)

end subroutine

subroutine c_elsi_finalize(handle_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: handle_c

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_finalize(handle_f)

end subroutine

subroutine c_elsi_ev_real(handle_c,ham_c,ovlp_c,eval_c,evec_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle),   pointer :: handle_f
   real(kind=c_double), pointer :: ham_f(:,:)
   real(kind=c_double), pointer :: ovlp_f(:,:)
   real(kind=c_double), pointer :: eval_f(:)
   real(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%n_basis
   lrow    = handle_f%n_lrow
   lcol    = handle_f%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_real(handle_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine c_elsi_ev_complex(handle_c,ham_c,ovlp_c,eval_c,evec_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle),              pointer :: handle_f
   complex(kind=c_double_complex), pointer :: ham_f(:,:)
   complex(kind=c_double_complex), pointer :: ovlp_f(:,:)
   real(kind=c_double),            pointer :: eval_f(:)
   complex(kind=c_double_complex), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%n_basis
   lrow    = handle_f%n_lrow
   lcol    = handle_f%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_complex(handle_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine c_elsi_ev_real_sparse(handle_c,ham_c,ovlp_c,eval_c,evec_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle),   pointer :: handle_f
   real(kind=c_double), pointer :: ham_f(:)
   real(kind=c_double), pointer :: ovlp_f(:)
   real(kind=c_double), pointer :: eval_f(:)
   real(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%n_basis
   nnz_l   = handle_f%nnz_l_sp
   lrow    = handle_f%n_lrow
   lcol    = handle_f%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_real_sparse(handle_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine c_elsi_ev_complex_sparse(handle_c,ham_c,ovlp_c,eval_c,evec_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: ham_c
   type(c_ptr), value, intent(in) :: ovlp_c
   type(c_ptr), value, intent(in) :: eval_c
   type(c_ptr), value, intent(in) :: evec_c

   type(elsi_handle),      pointer :: handle_f
   complex(kind=c_double), pointer :: ham_f(:)
   complex(kind=c_double), pointer :: ovlp_f(:)
   real(kind=c_double),    pointer :: eval_f(:)
   complex(kind=c_double), pointer :: evec_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%n_basis
   nnz_l   = handle_f%nnz_l_sp
   lrow    = handle_f%n_lrow
   lcol    = handle_f%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(eval_c,eval_f,shape=[n_basis])
   call c_f_pointer(evec_c,evec_f,shape=[lrow,lcol])

   call elsi_ev_complex_sparse(handle_f,ham_f,ovlp_f,eval_f,evec_f)

end subroutine

subroutine c_elsi_dm_real(handle_c,ham_c,ovlp_c,dm_c,energy) bind(C)

   implicit none

   type(c_ptr), value,  intent(in)    :: handle_c
   type(c_ptr), value,  intent(in)    :: ham_c
   type(c_ptr), value,  intent(in)    :: ovlp_c
   type(c_ptr), value,  intent(in)    :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle),   pointer :: handle_f
   real(kind=c_double), pointer :: ham_f(:,:)
   real(kind=c_double), pointer :: ovlp_f(:,:)
   real(kind=c_double), pointer :: dm_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%n_basis
   lrow    = handle_f%n_lrow
   lcol    = handle_f%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_dm_real(handle_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine c_elsi_dm_complex(handle_c,ham_c,ovlp_c,dm_c,energy) bind(C)

   implicit none

   type(c_ptr), value,  intent(in)    :: handle_c
   type(c_ptr), value,  intent(in)    :: ham_c
   type(c_ptr), value,  intent(in)    :: ovlp_c
   type(c_ptr), value,  intent(in)    :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle),              pointer :: handle_f
   complex(kind=c_double_complex), pointer :: ham_f(:,:)
   complex(kind=c_double_complex), pointer :: ovlp_f(:,:)
   complex(kind=c_double_complex), pointer :: dm_f(:,:)

   integer(kind=c_int) :: n_basis
   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   n_basis = handle_f%n_basis
   lrow    = handle_f%n_lrow
   lcol    = handle_f%n_lcol

   call c_f_pointer(ham_c,ham_f,shape=[lrow,lcol])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[lrow,lcol])
   call c_f_pointer(dm_c,dm_f,shape=[lrow,lcol])

   call elsi_dm_complex(handle_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine c_elsi_dm_real_sparse(handle_c,ham_c,ovlp_c,dm_c,energy) bind(C)

   implicit none

   type(c_ptr), value,  intent(in)    :: handle_c
   type(c_ptr), value,  intent(in)    :: ham_c
   type(c_ptr), value,  intent(in)    :: ovlp_c
   type(c_ptr), value,  intent(in)    :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle),   pointer :: handle_f
   real(kind=c_double), pointer :: ham_f(:)
   real(kind=c_double), pointer :: ovlp_f(:)
   real(kind=c_double), pointer :: dm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%nnz_l_sp

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(dm_c,dm_f,shape=[nnz_l])

   call elsi_dm_real_sparse(handle_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine c_elsi_dm_complex_sparse(handle_c,ham_c,ovlp_c,dm_c,energy) bind(C)

   implicit none

   type(c_ptr), value,  intent(in)    :: handle_c
   type(c_ptr), value,  intent(in)    :: ham_c
   type(c_ptr), value,  intent(in)    :: ovlp_c
   type(c_ptr), value,  intent(in)    :: dm_c
   real(kind=c_double), intent(inout) :: energy

   type(elsi_handle),      pointer :: handle_f
   complex(kind=c_double), pointer :: ham_f(:)
   complex(kind=c_double), pointer :: ovlp_f(:)
   complex(kind=c_double), pointer :: dm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%nnz_l_sp

   call c_f_pointer(ham_c,ham_f,shape=[nnz_l])
   call c_f_pointer(ovlp_c,ovlp_f,shape=[nnz_l])
   call c_f_pointer(dm_c,dm_f,shape=[nnz_l])

   call elsi_dm_complex_sparse(handle_f,ham_f,ovlp_f,dm_f,energy)

end subroutine

subroutine c_elsi_set_output(handle_c,out_level) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: out_level

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_output(handle_f,out_level)

end subroutine

subroutine c_elsi_set_unit_ovlp(handle_c,unit_ovlp) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: unit_ovlp

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_unit_ovlp(handle_f,unit_ovlp)

end subroutine

subroutine c_elsi_set_zero_def(handle_c,zero_def) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: zero_def

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_zero_def(handle_f,zero_def)

end subroutine

subroutine c_elsi_set_sing_check(handle_c,sing_check) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: sing_check

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sing_check(handle_f,sing_check)

end subroutine

subroutine c_elsi_set_sing_tol(handle_c,sing_tol) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: sing_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sing_tol(handle_f,sing_tol)

end subroutine

subroutine c_elsi_set_sing_stop(handle_c,sing_stop) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: sing_stop

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sing_stop(handle_f,sing_stop)

end subroutine

subroutine c_elsi_set_uplo(handle_c,uplo) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: uplo

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_uplo(handle_f,uplo)

end subroutine

subroutine c_elsi_set_elpa_solver(handle_c,elpa_solver) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: elpa_solver

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_elpa_solver(handle_f,elpa_solver)

end subroutine

subroutine c_elsi_set_elpa_n_single(handle_c,n_single) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_single

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_elpa_n_single(handle_f,n_single)

end subroutine

subroutine c_elsi_set_omm_flavor(handle_c,omm_flavor) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: omm_flavor

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_omm_flavor(handle_f,omm_flavor)

end subroutine

subroutine c_elsi_set_omm_n_elpa(handle_c,n_elpa) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_elpa

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_omm_n_elpa(handle_f,n_elpa)

end subroutine

subroutine c_elsi_set_omm_tol(handle_c,min_tol) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: min_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_omm_tol(handle_f,min_tol)

end subroutine

subroutine c_elsi_set_omm_ev_shift(handle_c,ev_shift) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: ev_shift

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_omm_ev_shift(handle_f,ev_shift)

end subroutine

subroutine c_elsi_set_pexsi_n_mu(handle_c,n_mu) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_mu

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_n_mu(handle_f,n_mu)

end subroutine

subroutine c_elsi_set_pexsi_n_pole(handle_c,n_pole) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_pole

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_n_pole(handle_f,n_pole)

end subroutine

subroutine c_elsi_set_pexsi_np_per_pole(handle_c,np_per_pole) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: np_per_pole

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_np_per_pole(handle_f,np_per_pole)

end subroutine

subroutine c_elsi_set_pexsi_np_symbo(handle_c,np_symbo) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: np_symbo

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_np_symbo(handle_f,np_symbo)

end subroutine

subroutine c_elsi_set_pexsi_ordering(handle_c,ordering) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: ordering

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_ordering(handle_f,ordering)

end subroutine

subroutine c_elsi_set_pexsi_temp(handle_c,temp) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: temp

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_temp(handle_f,temp)

end subroutine

subroutine c_elsi_set_pexsi_gap(handle_c,gap) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: gap

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_gap(handle_f,gap)

end subroutine

subroutine c_elsi_set_pexsi_delta_e(handle_c,delta_e) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: delta_e

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_delta_e(handle_f,delta_e)

end subroutine

subroutine c_elsi_set_pexsi_mu_min(handle_c,mu_min) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: mu_min

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_mu_min(handle_f,mu_min)

end subroutine

subroutine c_elsi_set_pexsi_mu_max(handle_c,mu_max) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: mu_max

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_mu_max(handle_f,mu_max)

end subroutine

subroutine c_elsi_set_pexsi_inertia_tol(handle_c,inertia_tol) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: inertia_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_pexsi_inertia_tol(handle_f,inertia_tol)

end subroutine

subroutine c_elsi_set_chess_erf_decay(handle_c,decay) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: decay

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_chess_erf_decay(handle_f,decay)

end subroutine

subroutine c_elsi_set_chess_erf_decay_min(handle_c,decay_min) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: decay_min

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_chess_erf_decay_min(handle_f,decay_min)

end subroutine

subroutine c_elsi_set_chess_erf_decay_max(handle_c,decay_max) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: decay_max

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_chess_erf_decay_max(handle_f,decay_max)

end subroutine

subroutine c_elsi_set_chess_ev_ham_min(handle_c,ev_min) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: ev_min

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_chess_ev_ham_min(handle_f,ev_min)

end subroutine

subroutine c_elsi_set_chess_ev_ham_max(handle_c,ev_max) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: ev_max

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_chess_ev_ham_max(handle_f,ev_max)

end subroutine

subroutine c_elsi_set_chess_ev_ovlp_min(handle_c,ev_min) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: ev_min

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_chess_ev_ovlp_min(handle_f,ev_min)

end subroutine

subroutine c_elsi_set_chess_ev_ovlp_max(handle_c,ev_max) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: ev_max

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_chess_ev_ovlp_max(handle_f,ev_max)

end subroutine

subroutine c_elsi_set_sips_n_elpa(handle_c,n_elpa) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_elpa

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_n_elpa(handle_f,n_elpa)

end subroutine

subroutine c_elsi_set_sips_n_slice(handle_c,n_slice) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: n_slice

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_n_slice(handle_f,n_slice)

end subroutine

subroutine c_elsi_set_sips_first_ev(handle_c,first_ev) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: first_ev

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_first_ev(handle_f,first_ev)

end subroutine

subroutine c_elsi_set_sips_slice_type(handle_c,slice_type) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: slice_type

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_slice_type(handle_f,slice_type)

end subroutine

subroutine c_elsi_set_sips_buffer(handle_c,buffer) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: buffer

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_buffer(handle_f,buffer)

end subroutine

subroutine c_elsi_set_sips_inertia_tol(handle_c,inertia_tol) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: inertia_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_inertia_tol(handle_f,inertia_tol)

end subroutine

subroutine c_elsi_set_sips_interval(handle_c,lower,upper) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: lower
   real(kind=c_double), value, intent(in) :: upper

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_sips_interval(handle_f,lower,upper)

end subroutine

subroutine c_elsi_set_dmp_method(handle_c,dmp_method) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: dmp_method

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_dmp_method(handle_f,dmp_method)

end subroutine

subroutine c_elsi_set_dmp_max_step(handle_c,max_step) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: max_step

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_dmp_max_step(handle_f,max_step)

end subroutine

subroutine c_elsi_set_dmp_tol(handle_c,dmp_tol) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: dmp_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_dmp_tol(handle_f,dmp_tol)

end subroutine

subroutine c_elsi_set_mu_broaden_scheme(handle_c,broaden_scheme) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: broaden_scheme

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_broaden_scheme(handle_f,broaden_scheme)

end subroutine

subroutine c_elsi_set_mu_broaden_width(handle_c,broaden_width) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: broaden_width

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_broaden_width(handle_f,broaden_width)

end subroutine

subroutine c_elsi_set_mu_tol(handle_c,mu_tol) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: mu_tol

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_tol(handle_f,mu_tol)

end subroutine

subroutine c_elsi_set_mu_spin_degen(handle_c,spin_degen) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: spin_degen

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_spin_degen(handle_f,spin_degen)

end subroutine

subroutine c_elsi_set_mu_mp_order(handle_c,mp_order) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: mp_order

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_mu_mp_order(handle_f,mp_order)

end subroutine

subroutine c_elsi_get_pexsi_mu_min(handle_c,mu_min) bind(C)

   implicit none

   type(c_ptr), value,  intent(in)  :: handle_c
   real(kind=c_double), intent(out) :: mu_min

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_pexsi_mu_min(handle_f,mu_min)

end subroutine

subroutine c_elsi_get_pexsi_mu_max(handle_c,mu_max) bind(C)

   implicit none

   type(c_ptr), value,  intent(in)  :: handle_c
   real(kind=c_double), intent(out) :: mu_max

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_pexsi_mu_max(handle_f,mu_max)

end subroutine

subroutine c_elsi_get_ovlp_sing(handle_c,ovlp_sing) bind(C)

   implicit none

   type(c_ptr), value,  intent(in)  :: handle_c
   integer(kind=c_int), intent(out) :: ovlp_sing

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_ovlp_sing(handle_f,ovlp_sing)

end subroutine

subroutine c_elsi_get_n_sing(handle_c,n_sing) bind(C)

   implicit none

   type(c_ptr), value,  intent(in)  :: handle_c
   integer(kind=c_int), intent(out) :: n_sing

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_n_sing(handle_f,n_sing)

end subroutine

subroutine c_elsi_get_mu(handle_c,mu) bind(C)

   implicit none

   type(c_ptr), value,  intent(in)  :: handle_c
   real(kind=c_double), intent(out) :: mu

   type(elsi_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_get_mu(handle_f,mu)

end subroutine

subroutine c_elsi_get_edm_real(handle_c,edm_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle),   pointer :: handle_f
   real(kind=c_double), pointer :: edm_f(:,:)

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%n_lrow
   lcol = handle_f%n_lcol

   call c_f_pointer(edm_c,edm_f,shape=[lrow,lcol])

   call elsi_get_edm_real(handle_f,edm_f)

end subroutine

subroutine c_elsi_get_edm_complex(handle_c,edm_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle),      pointer :: handle_f
   complex(kind=c_double), pointer :: edm_f(:,:)

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%n_lrow
   lcol = handle_f%n_lcol

   call c_f_pointer(edm_c,edm_f,shape=[lrow,lcol])

   call elsi_get_edm_complex(handle_f,edm_f)

end subroutine

subroutine c_elsi_get_edm_real_sparse(handle_c,edm_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle),   pointer :: handle_f
   real(kind=c_double), pointer :: edm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%nnz_l_sp

   call c_f_pointer(edm_c,edm_f,shape=[nnz_l])

   call elsi_get_edm_real_sparse(handle_f,edm_f)

end subroutine

subroutine c_elsi_get_edm_complex_sparse(handle_c,edm_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: handle_c
   type(c_ptr), value, intent(in) :: edm_c

   type(elsi_handle),      pointer :: handle_f
   complex(kind=c_double), pointer :: edm_f(:)

   integer(kind=c_int) :: nnz_l

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%nnz_l_sp

   call c_f_pointer(edm_c,edm_f,shape=[nnz_l])

   call elsi_get_edm_complex_sparse(handle_f,edm_f)

end subroutine

subroutine c_elsi_init_rw(handle_c,rw_task,parallel_mode,n_basis,n_electron) &
              bind(C)

   implicit none

   type(c_ptr)                            :: handle_c
   integer(kind=c_int), value, intent(in) :: rw_task
   integer(kind=c_int), value, intent(in) :: parallel_mode
   integer(kind=c_int), value, intent(in) :: n_basis
   real(kind=c_double), value, intent(in) :: n_electron

   type(elsi_rw_handle), pointer :: handle_f

   allocate(handle_f)

   call elsi_init_rw(handle_f,rw_task,parallel_mode,n_basis,n_electron)

   handle_c = c_loc(handle_f)

end subroutine

subroutine c_elsi_set_rw_mpi(handle_c,mpi_comm) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: mpi_comm

   type(elsi_rw_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_rw_mpi(handle_f,mpi_comm)

end subroutine

subroutine c_elsi_set_rw_blacs(handle_c,blacs_ctxt,block_size) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: blacs_ctxt
   integer(kind=c_int), value, intent(in) :: block_size

   type(elsi_rw_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_rw_blacs(handle_f,blacs_ctxt,block_size)

end subroutine

subroutine c_elsi_set_rw_csc(handle_c,nnz_g,nnz_l,n_lcol) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: nnz_g
   integer(kind=c_int), value, intent(in) :: nnz_l
   integer(kind=c_int), value, intent(in) :: n_lcol

   type(elsi_rw_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_rw_csc(handle_f,nnz_g,nnz_l,n_lcol)

end subroutine

subroutine c_elsi_finalize_rw(handle_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in) :: handle_c

   type(elsi_rw_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_finalize_rw(handle_f)

end subroutine

subroutine c_elsi_set_rw_output(handle_c,out_level) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   integer(kind=c_int), value, intent(in) :: out_level

   type(elsi_rw_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_rw_output(handle_f,out_level)

end subroutine

subroutine c_elsi_set_rw_zero_def(handle_c,zero_def) bind(C)

   implicit none

   type(c_ptr),         value, intent(in) :: handle_c
   real(kind=c_double), value, intent(in) :: zero_def

   type(elsi_rw_handle), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call elsi_set_rw_zero_def(handle_f,zero_def)

end subroutine

subroutine c_elsi_read_mat_dim(handle_c,name_c,n_electrons,n_basis,n_lrow,&
              n_lcol) bind(C)

   implicit none

   type(c_ptr), value, intent(in)               :: handle_c
   character(kind=c_char,len=1), dimension(128) :: name_c
   real(c_double),     intent(out)              :: n_electrons
   integer(c_int),     intent(out)              :: n_basis
   integer(c_int),     intent(out)              :: n_lrow
   integer(c_int),     intent(out)              :: n_lcol

   type(elsi_rw_handle), pointer :: handle_f

   character(len=:), allocatable :: name_f

   call c_f_pointer(handle_c,handle_f)

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_dim(handle_f,name_f,n_electrons,n_basis,n_lrow,n_lcol)

end subroutine

subroutine c_elsi_read_mat_dim_sparse(handle_c,name_c,n_electrons,n_basis,&
              nnz_g,nnz_l,n_lcol) bind(C)

   implicit none

   type(c_ptr), value, intent(in)               :: handle_c
   character(kind=c_char,len=1), dimension(128) :: name_c
   real(c_double),     intent(out)              :: n_electrons
   integer(c_int),     intent(out)              :: n_basis
   integer(c_int),     intent(out)              :: nnz_g
   integer(c_int),     intent(out)              :: nnz_l
   integer(c_int),     intent(out)              :: n_lcol

   type(elsi_rw_handle), pointer :: handle_f

   character(len=:), allocatable :: name_f

   call c_f_pointer(handle_c,handle_f)

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_dim_sparse(handle_f,name_f,n_electrons,n_basis,nnz_g,&
           nnz_l,n_lcol)

end subroutine

subroutine c_elsi_read_mat_real(handle_c,name_c,mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in)               :: handle_c
   character(kind=c_char,len=1), dimension(128) :: name_c
   type(c_ptr), value, intent(in)               :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   real(kind=c_double),  pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%n_lrow
   lcol = handle_f%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_real(handle_f,name_f,mat_f)

end subroutine

subroutine c_elsi_read_mat_real_sparse(handle_c,name_c,row_ind_c,col_ptr_c,&
              mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in)               :: handle_c
   character(kind=c_char,len=1), dimension(128) :: name_c
   type(c_ptr), value, intent(in)               :: row_ind_c
   type(c_ptr), value, intent(in)               :: col_ptr_c
   type(c_ptr), value, intent(in)               :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   integer(kind=c_int),  pointer :: row_ind_f(:)
   integer(kind=c_int),  pointer :: col_ptr_f(:)
   real(kind=c_double),  pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%nnz_l_sp
   lcol  = handle_f%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_real_sparse(handle_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine c_elsi_write_mat_real(handle_c,name_c,mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in)               :: handle_c
   character(kind=c_char,len=1), dimension(128) :: name_c
   type(c_ptr), value, intent(in)               :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   real(kind=c_double),  pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%n_lrow
   lcol = handle_f%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_real(handle_f,name_f,mat_f)

end subroutine

subroutine c_elsi_write_mat_real_sparse(handle_c,name_c,row_ind_c,col_ptr_c,&
              mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in)               :: handle_c
   character(kind=c_char,len=1), dimension(128) :: name_c
   type(c_ptr), value, intent(in)               :: row_ind_c
   type(c_ptr), value, intent(in)               :: col_ptr_c
   type(c_ptr), value, intent(in)               :: mat_c

   type(elsi_rw_handle), pointer :: handle_f
   integer(kind=c_int),  pointer :: row_ind_f(:)
   integer(kind=c_int),  pointer :: col_ptr_f(:)
   real(kind=c_double),  pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%nnz_l_sp
   lcol  = handle_f%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_real_sparse(handle_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine c_elsi_read_mat_complex(handle_c,name_c,mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in)               :: handle_c
   character(kind=c_char,len=1), dimension(128) :: name_c
   type(c_ptr), value, intent(in)               :: mat_c

   type(elsi_rw_handle),   pointer :: handle_f
   complex(kind=c_double), pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%n_lrow
   lcol = handle_f%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_complex(handle_f,name_f,mat_f)

end subroutine

subroutine c_elsi_read_mat_complex_sparse(handle_c,name_c,row_ind_c,col_ptr_c,&
              mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in)               :: handle_c
   character(kind=c_char,len=1), dimension(128) :: name_c
   type(c_ptr), value, intent(in)               :: row_ind_c
   type(c_ptr), value, intent(in)               :: col_ptr_c
   type(c_ptr), value, intent(in)               :: mat_c

   type(elsi_rw_handle),   pointer :: handle_f
   integer(kind=c_int),    pointer :: row_ind_f(:)
   integer(kind=c_int),    pointer :: col_ptr_f(:)
   complex(kind=c_double), pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%nnz_l_sp
   lcol  = handle_f%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_read_mat_complex_sparse(handle_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

subroutine c_elsi_write_mat_complex(handle_c,name_c,mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in)               :: handle_c
   character(kind=c_char,len=1), dimension(128) :: name_c
   type(c_ptr), value, intent(in)               :: mat_c

   type(elsi_rw_handle),   pointer :: handle_f
   complex(kind=c_double), pointer :: mat_f(:,:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: lrow
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   lrow = handle_f%n_lrow
   lcol = handle_f%n_lcol

   call c_f_pointer(mat_c,mat_f,shape=[lrow,lcol])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_complex(handle_f,name_f,mat_f)

end subroutine

subroutine c_elsi_write_mat_complex_sparse(handle_c,name_c,row_ind_c,col_ptr_c,&
              mat_c) bind(C)

   implicit none

   type(c_ptr), value, intent(in)               :: handle_c
   character(kind=c_char,len=1), dimension(128) :: name_c
   type(c_ptr), value, intent(in)               :: row_ind_c
   type(c_ptr), value, intent(in)               :: col_ptr_c
   type(c_ptr), value, intent(in)               :: mat_c

   type(elsi_rw_handle),   pointer :: handle_f
   integer(kind=c_int),    pointer :: row_ind_f(:)
   integer(kind=c_int),    pointer :: col_ptr_f(:)
   complex(kind=c_double), pointer :: mat_f(:)

   character(len=:), allocatable :: name_f

   integer(kind=c_int) :: nnz_l
   integer(kind=c_int) :: lcol

   call c_f_pointer(handle_c,handle_f)

   nnz_l = handle_f%nnz_l_sp
   lcol  = handle_f%n_lcol_sp

   call c_f_pointer(row_ind_c,row_ind_f,shape=[nnz_l])
   call c_f_pointer(col_ptr_c,col_ptr_f,shape=[lcol+1])
   call c_f_pointer(mat_c,mat_f,shape=[nnz_l])

   name_f = c_string_to_f_string(name_c)

   call elsi_write_mat_complex_sparse(handle_f,name_f,row_ind_f,col_ptr_f,mat_f)

end subroutine

end module ELSI_C_INTERFACE
