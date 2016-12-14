!Copyright (c) 2016, ELSI consortium
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in the
!   documentation and/or other materials provided with the distribution.
! * Neither the name of the ELSI project nor the
!   names of its contributors may be used to endorse or promote products
!   derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! This module is only compiled into the code if the actual libOMM is not
! availabel. The purpose is to make the (in that case unused) libOMM parts
! of the code compile.

module MatrixSwitch

   implicit none

   public :: m_register_pdbc
   public :: m_add
   public :: m_allocate
   public :: m_deallocate
   public :: ms_scalapack_setup
   public :: omm

   integer, parameter :: dp=selected_real_kind(15,300)

   type psp_matrix_spm
      character(3) :: str_type
      logical :: is_initialized=.false.
      logical :: is_serial
      logical :: is_real
      logical :: is_square
      integer :: nnz
      integer :: glb_dim1
      integer :: glb_dim2
      integer :: loc_dim1
      integer :: loc_dim2
      integer, allocatable :: row_ind(:)
      integer, allocatable :: col_ind(:)
      integer, allocatable :: col_ptr(:)
      integer :: desc(9)
      real(dp), allocatable :: dval(:)
      complex(dp), allocatable :: zval(:)
   end type psp_matrix_spm

   type matrix
      character(3) :: str_type ! label identifying the storage format
      logical :: is_initialized=.false.
      logical :: is_serial
      logical :: is_real
      logical :: is_square
      logical :: is_sparse
      integer :: dim1
      integer :: dim2
      integer, pointer :: iaux1(:) => null()
      integer, pointer :: iaux2(:) => null()
      integer, pointer :: iaux3(:) => null()
      integer, pointer :: iaux4(:) => null()
      real(dp), pointer :: dval(:,:) => null()
      complex(dp), pointer :: zval(:,:) => null()
      type(psp_matrix_spm) :: spm
   end type matrix

interface m_register_pdbc
   module procedure m_register_pddbc,m_register_pzdbc
end interface

interface m_add
   module procedure m_dadd,m_zadd
end interface

contains

subroutine m_allocate(m_name,i,j,label)

   implicit none

   type(matrix) :: m_name
   integer :: i
   integer :: j
   character(5), optional :: label

   write(*,"(A)") " A libOMM stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine m_register_pddbc(m_name,A,desc)

   implicit none

   type(matrix) :: m_name
   real(dp), target :: A(:,:)
   integer, target :: desc(9)

   write(*,"(A)") " A libOMM stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine m_register_pzdbc(m_name,A,desc)

   implicit none

   type(matrix) :: m_name 
   complex(dp), target :: A(:,:)
   integer, target :: desc(9) 

   write(*,"(A)") " A libOMM stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop
 
end subroutine

subroutine m_deallocate(m_name)

   implicit none

   type(matrix) :: m_name

   write(*,"(A)") " A libOMM stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine m_dadd(A,opA,C,alpha,beta,label)

   implicit none

   type(matrix) :: A
   character(1) :: opA
   type(matrix) :: C
   real(dp) :: alpha
   real(dp) :: beta
   character(3), optional :: label

   write(*,"(A)") " A libOMM stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine m_zadd(A,opA,C,alpha,beta,label)

   implicit none

   type(matrix) :: A
   character(1) :: opA
   type(matrix) :: C
   complex(dp) :: alpha
   complex(dp) :: beta
   character(3), optional :: label

   write(*,"(A)") " A libOMM stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine omm(m,n,H,S,new_S,e_min,D_min,calc_ED,eta,C_min,init_C,&
               T,scale_T,flavour,np,ip,cg_tol,long_out,dealloc,&
               m_storage,m_operation,mpi_rank)

   implicit none

   integer :: m
   integer :: n
   type(matrix) :: H
   type(matrix) :: S
   logical :: new_S
   real(dp) :: e_min
   type(matrix) :: D_min
   logical :: calc_ED
   real(dp) :: eta
   type(matrix) :: C_min
   logical :: init_C
   type(matrix) :: T
   real(dp) :: scale_T
   integer :: flavour
   integer :: np
   integer :: ip
   real(dp) :: cg_tol
   logical :: long_out
   logical :: dealloc
   character(5) :: m_storage
   character(3) :: m_operation
   integer :: mpi_rank

   write(*,"(A)") " A libOMM stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine ms_scalapack_setup(mpi_comm,nprow,order,bs_def,bs_list,icontxt)

   implicit none

   integer :: mpi_comm
   integer :: nprow
   character(1) :: order
   integer :: bs_def
   integer, optional :: bs_list(:)
   integer, optional :: icontxt

   write(*,"(A)") " A libOMM stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

end module
