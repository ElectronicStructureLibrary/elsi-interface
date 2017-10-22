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
!! This module sets up work arrays in ELSI.
!!
module ELSI_MATRICES

   use ELSI_CONSTANTS, only: BLACS_DENSE,OMM_SOLVER,MULTI_PROC,FULL_MAT,UT_MAT,&
                             LT_MAT
   use ELSI_DATATYPE
   use ELSI_MALLOC
   use ELSI_PRECISION, only: i4,r8
   use MATRIXSWITCH,   only: m_register_pdbc

   implicit none

   private

   public :: elsi_set_ham
   public :: elsi_set_ovlp
   public :: elsi_set_eval
   public :: elsi_set_evec
   public :: elsi_set_dm
   public :: elsi_set_row_ind
   public :: elsi_set_col_ptr
   public :: elsi_set_sparse_ham
   public :: elsi_set_sparse_ovlp
   public :: elsi_set_sparse_dm
   public :: elsi_set_full_mat

   interface elsi_set_ham
      module procedure elsi_set_real_ham,&
                       elsi_set_complex_ham
   end interface

   interface elsi_set_sparse_ham
      module procedure elsi_set_sparse_real_ham,&
                       elsi_set_sparse_complex_ham
   end interface

   interface elsi_set_ovlp
      module procedure elsi_set_real_ovlp,&
                       elsi_set_complex_ovlp
   end interface

   interface elsi_set_sparse_ovlp
      module procedure elsi_set_sparse_real_ovlp,&
                       elsi_set_sparse_complex_ovlp
   end interface

   interface elsi_set_evec
      module procedure elsi_set_real_evec,&
                       elsi_set_complex_evec
   end interface

   interface elsi_set_dm
      module procedure elsi_set_real_dm,&
                       elsi_set_complex_dm
   end interface

   interface elsi_set_sparse_dm
      module procedure elsi_set_sparse_real_dm,&
                       elsi_set_sparse_complex_dm
   end interface

   interface elsi_set_full_mat
      module procedure elsi_set_full_mat_real,&
                       elsi_set_full_mat_complex
   end interface

contains

!>
!! This routine sets the real hamiltonian matrix.
!!
subroutine elsi_set_real_ham(e_h,h_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                         !< Handle
   real(kind=r8),     intent(inout), target :: h_in(e_h%n_lrow,e_h%n_lcol) !< Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_real_ham"

   if(e_h%matrix_format == BLACS_DENSE) then
      call elsi_set_full_mat(e_h,h_in)
   endif

   if(e_h%solver == OMM_SOLVER) then
      call m_register_pdbc(e_h%ham_omm,h_in,e_h%sc_desc)
   else
      e_h%ham_real => h_in
   endif

end subroutine

!>
!! This routine sets the complex hamiltonian matrix.
!!
subroutine elsi_set_complex_ham(e_h,h_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                         !< Handle
   complex(kind=r8),  intent(inout), target :: h_in(e_h%n_lrow,e_h%n_lcol) !< Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_complex_ham"

   if(e_h%matrix_format == BLACS_DENSE) then
      call elsi_set_full_mat(e_h,h_in)
   endif

   if(e_h%solver == OMM_SOLVER) then
      call m_register_pdbc(e_h%ham_omm,h_in,e_h%sc_desc)
   else
      e_h%ham_cmplx => h_in
   endif

end subroutine

!>
!! This routine sets the sparse real Hamiltonian matrix.
!!
subroutine elsi_set_sparse_real_ham(e_h,h_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                !< Handle
   real(kind=r8),     intent(inout), target :: h_in(e_h%nnz_l_sp) !< Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_sparse_real_ham"

   e_h%ham_real_ccs => h_in

end subroutine

!>
!! This routine sets the sparse complex Hamiltonian matrix.
!!
subroutine elsi_set_sparse_complex_ham(e_h,h_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                !< Handle
   complex(kind=r8),  intent(inout), target :: h_in(e_h%nnz_l_sp) !< Hamiltonian matirx

   character*40, parameter :: caller = "elsi_set_sparse_complex_ham"

   e_h%ham_cmplx_ccs => h_in

end subroutine

!>
!! This routine sets the real overlap matrix.
!!
subroutine elsi_set_real_ovlp(e_h,s_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                         !< Handle
   real(kind=r8),     intent(inout), target :: s_in(e_h%n_lrow,e_h%n_lcol) !< Overlap matrix

   character*40, parameter :: caller = "elsi_set_real_ovlp"

   if(.not. e_h%ovlp_is_unit) then
      if(e_h%matrix_format == BLACS_DENSE .and. e_h%n_elsi_calls == 1) then
         call elsi_set_full_mat(e_h,s_in)
      endif

      if(e_h%solver == OMM_SOLVER) then
         call m_register_pdbc(e_h%ovlp_omm,s_in,e_h%sc_desc)
      else
         e_h%ovlp_real => s_in
      endif
   endif

end subroutine

!>
!! This routine sets the complex ovlp matrix.
!!
subroutine elsi_set_complex_ovlp(e_h,s_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                         !< Handle
   complex(kind=r8),  intent(inout), target :: s_in(e_h%n_lrow,e_h%n_lcol) !< Overlap matrix

   character*40, parameter :: caller = "elsi_set_complex_ovlp"

   if(.not. e_h%ovlp_is_unit) then
      if(e_h%matrix_format == BLACS_DENSE .and. e_h%n_elsi_calls == 1) then
         call elsi_set_full_mat(e_h,s_in)
      endif

      if(e_h%solver == OMM_SOLVER) then
         call m_register_pdbc(e_h%ovlp_omm,s_in,e_h%sc_desc)
      else
         e_h%ovlp_cmplx => s_in
      endif
   endif

end subroutine

!>
!! This routine sets the sparse real overlap matrix.
!!
subroutine elsi_set_sparse_real_ovlp(e_h,s_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                !< Handle
   real(kind=r8),     intent(inout), target :: s_in(e_h%nnz_l_sp) !< Overlap matrix

   character*40, parameter :: caller = "elsi_set_sparse_real_ovlp"

   if(.not. e_h%ovlp_is_unit) then
      e_h%ovlp_real_ccs => s_in
   endif

end subroutine

!>
!! This routine sets the sparse complex overlap matrix.
!!
subroutine elsi_set_sparse_complex_ovlp(e_h,s_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                !< Handle
   complex(kind=r8),  intent(inout), target :: s_in(e_h%nnz_l_sp) !< Overlap matrix

   character*40, parameter :: caller = "elsi_set_sparse_complex_ovlp"

   if(.not. e_h%ovlp_is_unit) then
      e_h%ovlp_cmplx_ccs => s_in
   endif

end subroutine

!>
!! This routine sets the eigenvalues.
!!
subroutine elsi_set_eval(e_h,eval_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                  !< Handle
   real(kind=r8),     intent(inout), target :: eval_in(e_h%n_basis) !< Eigenvalues

   character*40, parameter :: caller = "elsi_set_eval"

   e_h%eval => eval_in

end subroutine

!>
!! This routine sets the real eigenvectors.
!!
subroutine elsi_set_real_evec(e_h,evec_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                            !< Handle
   real(kind=r8),     intent(inout), target :: evec_in(e_h%n_lrow,e_h%n_lcol) !< Eigenvectors

   character*40, parameter :: caller = "elsi_set_real_evec"

   e_h%evec_real => evec_in

end subroutine

!>
!! This routine sets the complex eigenvectors.
!!
subroutine elsi_set_complex_evec(e_h,evec_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                            !< Handle
   complex(kind=r8),  intent(inout), target :: evec_in(e_h%n_lrow,e_h%n_lcol) !< Eigenvectors

   character*40, parameter :: caller = "elsi_set_complex_evec"

   e_h%evec_cmplx => evec_in

end subroutine

!>
!! This routine sets the real density matrix.
!!
subroutine elsi_set_real_dm(e_h,d_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                         !< Handle
   real(kind=r8),     intent(inout), target :: d_in(e_h%n_lrow,e_h%n_lcol) !< Density matrix

   character*40, parameter :: caller = "elsi_set_real_dm"

   if(e_h%solver == OMM_SOLVER) then
      call m_register_pdbc(e_h%dm_omm,d_in,e_h%sc_desc)
   else
      e_h%dm_real => d_in
   endif

end subroutine

!>
!! This routine sets the complex density matrix.
!!
subroutine elsi_set_complex_dm(e_h,d_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                         !< Handle
   complex(kind=r8),  intent(inout), target :: d_in(e_h%n_lrow,e_h%n_lcol) !< Density matrix

   character*40, parameter :: caller = "elsi_set_complex_dm"

   if(e_h%solver == OMM_SOLVER) then
      call m_register_pdbc(e_h%dm_omm,d_in,e_h%sc_desc)
   else
      e_h%dm_cmplx => d_in
   endif

end subroutine

!>
!! This routine sets the sparse real density matrix.
!!
subroutine elsi_set_sparse_real_dm(e_h,d_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                !< Handle
   real(kind=r8),     intent(inout), target :: d_in(e_h%nnz_l_sp) !< Density matrix

   character*40, parameter :: caller = "elsi_set_sparse_real_dm"

   e_h%dm_real_ccs => d_in

end subroutine

!>
!! This routine sets the sparse complex density matrix.
!!
subroutine elsi_set_sparse_complex_dm(e_h,d_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                !< Handle
   complex(kind=r8),  intent(inout), target :: d_in(e_h%nnz_l_sp) !< Density matrix

   character*40, parameter :: caller = "elsi_set_sparse_complex_dm"

   e_h%dm_cmplx_ccs => d_in

end subroutine

!>
!! This routine sets the row index.
!!
subroutine elsi_set_row_ind(e_h,row_ind_in)

   implicit none

   type(elsi_handle), intent(inout)        :: e_h                      !< Handle
   integer(kind=i4),  intent(in),   target :: row_ind_in(e_h%nnz_l_sp) !< Row index

   character*40, parameter :: caller = "elsi_set_row_ind"

   e_h%row_ind_ccs => row_ind_in

end subroutine

!>
!! This routine sets the column pointer.
!!
subroutine elsi_set_col_ptr(e_h,col_ptr_in)

   implicit none

   type(elsi_handle), intent(inout)        :: e_h                         !< Handle
   integer(kind=i4),  intent(in),   target :: col_ptr_in(e_h%n_lcol_sp+1) !< Column pointer

   character*40, parameter :: caller = "elsi_set_col_ptr"

   e_h%col_ptr_ccs => col_ptr_in

end subroutine

!>
!! This routine sets a full matrix from a (upper or lower) triangular matrix.
!! The size of matrix should be the same as the Hamiltonian matrix.
!!
subroutine elsi_set_full_mat_real(e_h,mat)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                        !< Handle
   real(kind=r8),     intent(inout) :: mat(e_h%n_lrow,e_h%n_lcol) !< Matrix

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   real(kind=r8), allocatable :: tmp_real(:,:)

   character*40, parameter :: caller = "elsi_set_full_mat_real"

   if(e_h%uplo /= FULL_MAT .and. e_h%parallel_mode == MULTI_PROC) then
      call elsi_allocate(e_h,tmp_real,e_h%n_lrow,e_h%n_lcol,"tmp_real",caller)

      call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,mat,1,1,e_h%sc_desc,0.0_r8,&
              tmp_real,1,1,e_h%sc_desc)

      if(e_h%uplo == UT_MAT) then ! Upper triangular
         do i_col = 1,e_h%n_basis-1
            if(e_h%loc_col(i_col) == 0) cycle

            do i_row = i_col+1,e_h%n_basis
               if(e_h%loc_row(i_row) > 0) then
                  mat(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                     tmp_real(e_h%loc_row(i_row),e_h%loc_col(i_col))
               endif
            enddo
         enddo
      elseif(e_h%uplo == LT_MAT) then ! Lower triangular
         do i_col = 2,e_h%n_basis
            if(e_h%loc_col(i_col) == 0) cycle

            do i_row = 1,i_col-1
               if(e_h%loc_row(i_row) > 0) then
                  mat(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                     tmp_real(e_h%loc_row(i_row),e_h%loc_col(i_col))
               endif
            enddo
         enddo
      endif

      call elsi_deallocate(e_h,tmp_real,"tmp_real")
   endif

end subroutine

!>
!! This routine sets a full matrix from a (upper or lower) triangular matrix.
!! The size of matrix should be the same as the Hamiltonian matrix.
!!
subroutine elsi_set_full_mat_complex(e_h,mat)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                        !< Handle
   complex(kind=r8),  intent(inout) :: mat(e_h%n_lrow,e_h%n_lcol) !< Matrix

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character*40, parameter :: caller = "elsi_set_full_mat_complex"

   if(e_h%uplo /= FULL_MAT .and. e_h%parallel_mode == MULTI_PROC) then
      call elsi_allocate(e_h,tmp_cmplx,e_h%n_lrow,e_h%n_lcol,"tmp_cmplx",caller)

      call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),mat,1,1,e_h%sc_desc,&
              (0.0_r8,0.0_r8),tmp_cmplx,1,1,e_h%sc_desc)

      if(e_h%uplo == UT_MAT) then ! Upper triangular
         do i_col = 1,e_h%n_basis-1
            if(e_h%loc_col(i_col) == 0) cycle

            do i_row = i_col+1,e_h%n_basis
               if(e_h%loc_row(i_row) > 0) then
                  mat(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                     tmp_cmplx(e_h%loc_row(i_row),e_h%loc_col(i_col))
               endif
            enddo
         enddo
      elseif(e_h%uplo == LT_MAT) then ! Lower triangular
         do i_col = 2,e_h%n_basis
            if(e_h%loc_col(i_col) == 0) cycle

            do i_row = 1,i_col-1
               if(e_h%loc_row(i_row) > 0) then
                  mat(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                     tmp_cmplx(e_h%loc_row(i_row),e_h%loc_col(i_col))
               endif
            enddo
         enddo
      endif

      call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")

      ! Make diagonal real
      do i_col = 1,e_h%n_basis
         if(e_h%loc_col(i_col) == 0 .or. e_h%loc_row(i_col) == 0) cycle

         mat(e_h%loc_row(i_col),e_h%loc_col(i_col)) = &
            real(mat(e_h%loc_row(i_col),e_h%loc_col(i_col)),kind=r8)
      enddo
   endif

end subroutine

end module ELSI_MATRICES
