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
!! This module contains utilities for ELSI, including array allocations and
!! printings for debugging.
!!
module ELSI_UTILS

   use ELSI_CONSTANTS
   use ELSI_DATATYPE
   use ELSI_PRECISION
   use FOE_BASE, only: foe_data_deallocate
   use F_PPEXSI_INTERFACE
   use MATRIXSWITCH, only: m_register_pdbc,m_deallocate
   use M_QETSC
   use SPARSEMATRIX_BASE, only: deallocate_sparse_matrix,deallocate_matrices

   implicit none

   include "mpif.h"

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
   public :: elsi_statement_print
   public :: elsi_allocate
   public :: elsi_deallocate
   public :: elsi_cleanup
   public :: elsi_stop
   public :: elsi_check
   public :: elsi_check_handle
   public :: elsi_get_global_row
   public :: elsi_get_global_col
   public :: elsi_get_local_nnz
   public :: elsi_set_full_mat
   public :: elsi_init_timers
   public :: elsi_get_time
   public :: elsi_final_print

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

   interface elsi_allocate
      module procedure elsi_allocate_integer4_1d,&
                       elsi_allocate_integer4_2d,&
                       elsi_allocate_integer4_3d,&
                       elsi_allocate_integer8_1d,&
                       elsi_allocate_real8_1d,&
                       elsi_allocate_real8_2d,&
                       elsi_allocate_real8_3d,&
                       elsi_allocate_complex16_1d,&
                       elsi_allocate_complex16_2d,&
                       elsi_allocate_complex16_3d
   end interface

   interface elsi_deallocate
      module procedure elsi_deallocate_integer4_1d,&
                       elsi_deallocate_integer4_2d,&
                       elsi_deallocate_integer4_3d,&
                       elsi_deallocate_integer8_1d,&
                       elsi_deallocate_real8_1d,&
                       elsi_deallocate_real8_2d,&
                       elsi_deallocate_real8_3d,&
                       elsi_deallocate_complex16_1d,&
                       elsi_deallocate_complex16_2d,&
                       elsi_deallocate_complex16_3d
   end interface

   interface elsi_get_local_nnz
      module procedure elsi_get_local_nnz_real,&
                       elsi_get_local_nnz_complex
   end interface

   interface elsi_set_full_mat
      module procedure elsi_set_full_mat_real,&
                       elsi_set_full_mat_complex
   end interface

contains

!>
!! This routine prints a statement.
!!
subroutine elsi_statement_print(info_str,e_h)

   implicit none

   character(len=*),  intent(in) :: info_str !< Message to print
   type(elsi_handle), intent(in) :: e_h      !< Handle

   if(print_info) then
      if(e_h%myid_all == 0) then
         write(*,"(A)") trim(info_str)
      endif
   endif

end subroutine

!>
!! This routine allocates a 1D array with real(kind=r8).
!!
subroutine elsi_allocate_real8_1d(e_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h       !< Handle
   real(kind=r8),     intent(inout), allocatable :: array(:)  !< Data
   integer(kind=i4),  intent(in)                 :: dim1      !< Size
   character(len=*),  intent(in)                 :: arrayname !< Name
   character(len=*),  intent(in)                 :: caller    !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info_str

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*8

      write(info_str,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info_str,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 1D array with integer(kind=i4).
!!
subroutine elsi_allocate_integer4_1d(e_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h       !< Handle
   integer(kind=i4),  intent(inout), allocatable :: array(:)  !< Data
   integer(kind=i4),  intent(in)                 :: dim1      !< Size
   character(len=*),  intent(in)                 :: arrayname !< Name
   character(len=*),  intent(in)                 :: caller    !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info_str

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*4

      write(info_str,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info_str,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 1D array with integer(kind=i8).
!!
subroutine elsi_allocate_integer8_1d(e_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h       !< Handle
   integer(kind=i8),  intent(inout), allocatable :: array(:)  !< Data
   integer(kind=i4),  intent(in)                 :: dim1      !< Size
   character(len=*),  intent(in)                 :: arrayname !< Name
   character(len=*),  intent(in)                 :: caller    !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info_str

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*8

      write(info_str,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info_str,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 1D array with complex(kind=r8).
!!
subroutine elsi_allocate_complex16_1d(e_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h       !< Handle
   complex(kind=r8),  intent(inout), allocatable :: array(:)  !< Data
   integer(kind=i4),  intent(in)                 :: dim1      !< Size
   character(len=*),  intent(in)                 :: arrayname !< Name
   character(len=*),  intent(in)                 :: caller    !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info_str

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*16

      write(info_str,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info_str,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = (0.0_r8,0.0_r8)

end subroutine

!>
!! This routine allocates a 2D array of real(kind=r8).
!!
subroutine elsi_allocate_real8_2d(e_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h        !< Handle
   real(kind=r8),     intent(inout), allocatable :: array(:,:) !< Data
   integer(kind=i4),  intent(in)                 :: dim1       !< Size
   integer(kind=i4),  intent(in)                 :: dim2       !< Size
   character(len=*),  intent(in)                 :: arrayname  !< Name
   character(len=*),  intent(in)                 :: caller     !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info_str

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*8

      write(info_str,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(info_str,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 2D array of integer(kind=i4).
!!
subroutine elsi_allocate_integer4_2d(e_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h        !< Handle
   integer(kind=i4),  intent(inout), allocatable :: array(:,:) !< Data
   integer(kind=i4),  intent(in)                 :: dim1       !< Size
   integer(kind=i4),  intent(in)                 :: dim2       !< Size
   character(len=*),  intent(in)                 :: arrayname  !< Name
   character(len=*),  intent(in)                 :: caller     !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info_str

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*4

      write(info_str,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(info_str,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 2D array of complex(kind=r8).
!!
subroutine elsi_allocate_complex16_2d(e_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h        !< Handle
   complex(kind=r8),  intent(inout), allocatable :: array(:,:) !< Data
   integer(kind=i4),  intent(in)                 :: dim1       !< Size
   integer(kind=i4),  intent(in)                 :: dim2       !< Size
   character(len=*),  intent(in)                 :: arrayname  !< Name
   character(len=*),  intent(in)                 :: caller     !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info_str

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*16

      write(info_str,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(info_str,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = (0.0_r8,0.0_r8)

end subroutine

!>
!! This routine allocates a 3D array of real(kind=r8).
!!
subroutine elsi_allocate_real8_3d(e_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h          !< Handle
   real(kind=r8),     intent(inout), allocatable :: array(:,:,:) !< Data
   integer(kind=i4),  intent(in)                 :: dim1         !< Size
   integer(kind=i4),  intent(in)                 :: dim2         !< Size
   integer(kind=i4),  intent(in)                 :: dim3         !< Size
   character(len=*),  intent(in)                 :: arrayname    !< Name
   character(len=*),  intent(in)                 :: caller       !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info_str

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*dim3*8

      write(info_str,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   allocate(array(dim1,dim2,dim3),stat=error)

   if(error > 0) then
      write(info_str,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 3D array of integer(kind=i4).
!!
subroutine elsi_allocate_integer4_3d(e_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h          !< Handle
   integer(kind=i4),  intent(inout), allocatable :: array(:,:,:) !< Data
   integer(kind=i4),  intent(in)                 :: dim1         !< Size
   integer(kind=i4),  intent(in)                 :: dim2         !< Size
   integer(kind=i4),  intent(in)                 :: dim3         !< Size
   character(len=*),  intent(in)                 :: arrayname    !< Name
   character(len=*),  intent(in)                 :: caller       !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info_str

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*dim3*4

      write(info_str,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   allocate(array(dim1,dim2,dim3),stat=error)

   if(error > 0) then
      write(info_str,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 3D array of complex(kind=r8).
!!
subroutine elsi_allocate_complex16_3d(e_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h          !< Handle
   complex(kind=r8),  intent(inout), allocatable :: array(:,:,:) !< Data
   integer(kind=i4),  intent(in)                 :: dim1         !< Size
   integer(kind=i4),  intent(in)                 :: dim2         !< Size
   integer(kind=i4),  intent(in)                 :: dim3         !< Size
   character(len=*),  intent(in)                 :: arrayname    !< Name
   character(len=*),  intent(in)                 :: caller       !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info_str

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*dim3*16

      write(info_str,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   allocate(array(dim1,dim2,dim3),stat=error)

   if(error > 0) then
      write(info_str,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = (0.0_r8,0.0_r8)

end subroutine

!>
!! This routine deallocates a 1D array with real(kind=r8).
!!
subroutine elsi_deallocate_real8_1d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h       !< Handle
   real(kind=r8),     intent(inout), allocatable :: array(:)  !< Data
   character(len=*),  intent(in)                 :: arrayname !< Name

   character*200 :: info_str

   if(print_mem) then
      write(info_str,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 1D array with integer(kind=i4).
!!
subroutine elsi_deallocate_integer4_1d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h       !< Handle
   integer(kind=i4),  intent(inout), allocatable :: array(:)  !< Data
   character(len=*),  intent(in)                 :: arrayname !< Name

   character*200 :: info_str

   if(print_mem) then
      write(info_str,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 1D array with integer(kind=i8).
!!
subroutine elsi_deallocate_integer8_1d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h       !< Handle
   integer(kind=i8),  intent(inout), allocatable :: array(:)  !< Data
   character(len=*),  intent(in)                 :: arrayname !< Name

   character*200 :: info_str

   if(print_mem) then
      write(info_str,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 1D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex16_1d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h       !< Handle
   complex(kind=r8),  intent(inout), allocatable :: array(:)  !< Data
   character(len=*),  intent(in)                 :: arrayname !< Name

   character*200 :: info_str

   if(print_mem) then
      write(info_str,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 2D array with real(kind=r8).
!!
subroutine elsi_deallocate_real8_2d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h        !< Handle
   real(kind=r8),     intent(inout), allocatable :: array(:,:) !< Data
   character(len=*),  intent(in)                 :: arrayname  !< Name

   character*200 :: info_str

   if(print_mem) then
      write(info_str,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 2D array with integer(kind=i4).
!!
subroutine elsi_deallocate_integer4_2d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h        !< Handle
   integer(kind=i4),  intent(inout), allocatable :: array(:,:) !< Data
   character(len=*),  intent(in)                 :: arrayname  !< Name

   character*200 :: info_str

   if(print_mem) then
      write(info_str,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 2D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex16_2d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h        !< Handle
   complex(kind=r8),  intent(inout), allocatable :: array(:,:) !< Data
   character(len=*),  intent(in)                 :: arrayname  !< Name

   character*200 :: info_str

   if(print_mem) then
      write(info_str,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 3D array with real(kind=r8).
!!
subroutine elsi_deallocate_real8_3d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h          !< Handle
   real(kind=r8),     intent(inout), allocatable :: array(:,:,:) !< Data
   character(len=*),  intent(in)                 :: arrayname    !< Name

   character*200 :: info_str

   if(print_mem) then
      write(info_str,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 3D array with integer(kind=i4).
!!
subroutine elsi_deallocate_integer4_3d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h          !< Handle
   integer(kind=i4),  intent(inout), allocatable :: array(:,:,:) !< Data
   character(len=*),  intent(in)                 :: arrayname    !< Name

   character*200 :: info_str

   if(print_mem) then
      write(info_str,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 3D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex16_3d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h          !< Handle
   complex(kind=r8),  intent(inout), allocatable :: array(:,:,:) !< Data
   character(len=*),  intent(in)                 :: arrayname    !< Name

   character*200 :: info_str

   if(print_mem) then
      write(info_str,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info_str,e_h)
   endif

   deallocate(array)

end subroutine

!>
!! Clean shutdown in case of errors.
!!
subroutine elsi_stop(info,e_h,caller)

   implicit none

   character(len=*),  intent(in) :: info   !< Error message to print
   type(elsi_handle), intent(in) :: e_h    !< Handle
   character(len=*),  intent(in) :: caller !< The subroutine in trouble

   character(len=4096) :: info_str
   integer :: i_task
   integer :: mpierr

   if(e_h%global_mpi_ready) then
      do i_task = 0,e_h%n_procs_all-1
         if(e_h%myid_all == i_task) then
            write(info_str,"(' Error! MPI task',I5,' in ',A,': ',A)")&
               e_h%myid_all,trim(caller),trim(info)
            write(*,"(A)") trim(info_str)
            exit
         endif
         call MPI_Barrier(e_h%mpi_comm_all,mpierr)
      enddo
   elseif(e_h%mpi_ready) then
      do i_task = 0,e_h%n_procs-1
         if(e_h%myid == i_task) then
            write(info_str,"(' Error! MPI task',I5,' in ',A,': ',A)")&
               e_h%myid,trim(caller),trim(info)
            write(*,"(A)") trim(info_str)
            exit
         endif
         call MPI_Barrier(e_h%mpi_comm,mpierr)
      enddo
   else
      write(info_str,"(' Error!',A,': ',A)") trim(caller),trim(info)
      write(*,"(A)") trim(info_str)
   endif

   stop

end subroutine

!>
!! This routine sets the real hamiltonian matrix.
!!
subroutine elsi_set_real_ham(e_h,h_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                             !< Handle
   real(kind=r8),     intent(inout), target :: h_in(e_h%n_l_rows,e_h%n_l_cols) !< Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_real_ham"

   if(e_h%matrix_format == BLACS_DENSE) then
      call elsi_set_full_mat(e_h,h_in)
   endif

   if(e_h%solver == LIBOMM) then
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

   type(elsi_handle), intent(inout)         :: e_h                             !< Handle
   complex(kind=r8),  intent(inout), target :: h_in(e_h%n_l_rows,e_h%n_l_cols) !< Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_complex_ham"

   if(e_h%matrix_format == BLACS_DENSE) then
      call elsi_set_full_mat(e_h,h_in)
   endif

   if(e_h%solver == LIBOMM) then
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

   type(elsi_handle), intent(inout)         :: e_h                             !< Handle
   real(kind=r8),     intent(inout), target :: s_in(e_h%n_l_rows,e_h%n_l_cols) !< Overlap matrix

   character*40, parameter :: caller = "elsi_set_real_ovlp"

   if(.not. e_h%ovlp_is_unit) then
      if(e_h%matrix_format == BLACS_DENSE .and. e_h%n_elsi_calls == 1) then
         call elsi_set_full_mat(e_h,s_in)
      endif

      if(e_h%solver == LIBOMM) then
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

   type(elsi_handle), intent(inout)         :: e_h                             !< Handle
   complex(kind=r8),  intent(inout), target :: s_in(e_h%n_l_rows,e_h%n_l_cols) !< Overlap matrix

   character*40, parameter :: caller = "elsi_set_complex_ovlp"

   if(.not. e_h%ovlp_is_unit) then
      if(e_h%matrix_format == BLACS_DENSE .and. e_h%n_elsi_calls == 1) then
         call elsi_set_full_mat(e_h,s_in)
      endif

      if(e_h%solver == LIBOMM) then
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

   type(elsi_handle), intent(inout)         :: e_h                                !< Handle
   real(kind=r8),     intent(inout), target :: evec_in(e_h%n_l_rows,e_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_set_real_evec"

   e_h%evec_real => evec_in

end subroutine

!>
!! This routine sets the complex eigenvectors.
!!
subroutine elsi_set_complex_evec(e_h,evec_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                                !< Handle
   complex(kind=r8),  intent(inout), target :: evec_in(e_h%n_l_rows,e_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_set_complex_evec"

   e_h%evec_cmplx => evec_in

end subroutine

!>
!! This routine sets the real density matrix.
!!
subroutine elsi_set_real_dm(e_h,d_in)

   implicit none

   type(elsi_handle), intent(inout)         :: e_h                             !< Handle
   real(kind=r8),     intent(inout), target :: d_in(e_h%n_l_rows,e_h%n_l_cols) !< Density matrix

   character*40, parameter :: caller = "elsi_set_real_dm"

   if(e_h%solver == LIBOMM) then
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

   type(elsi_handle), intent(inout)         :: e_h                             !< Handle
   complex(kind=r8),  intent(inout), target :: d_in(e_h%n_l_rows,e_h%n_l_cols) !< Density matrix

   character*40, parameter :: caller = "elsi_set_complex_dm"

   if(e_h%solver == LIBOMM) then
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

   type(elsi_handle), intent(inout)        :: e_h                           !< Handle
   integer(kind=i4),  intent(in),   target :: col_ptr_in(e_h%n_l_cols_sp+1) !< Column pointer

   character*40, parameter :: caller = "elsi_set_col_ptr"

   e_h%col_ptr_ccs => col_ptr_in

end subroutine

!>
!! This routine frees memory.
!!
subroutine elsi_cleanup(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: ierr

   character*40, parameter :: caller = "elsi_cleanup"

   ! Nullify pointers
   if(associated(e_h%ham_real))        nullify(e_h%ham_real)
   if(associated(e_h%ham_cmplx))       nullify(e_h%ham_cmplx)
   if(associated(e_h%ovlp_real))       nullify(e_h%ovlp_real)
   if(associated(e_h%ovlp_cmplx))      nullify(e_h%ovlp_cmplx)
   if(associated(e_h%eval))            nullify(e_h%eval)
   if(associated(e_h%evec_real))       nullify(e_h%evec_real)
   if(associated(e_h%evec_cmplx))      nullify(e_h%evec_cmplx)
   if(associated(e_h%dm_real))         nullify(e_h%dm_real)
   if(associated(e_h%dm_cmplx))        nullify(e_h%dm_cmplx)
   if(associated(e_h%ham_real_ccs))    nullify(e_h%ham_real_ccs)
   if(associated(e_h%ham_cmplx_ccs))   nullify(e_h%ham_cmplx_ccs)
   if(associated(e_h%ovlp_real_ccs))   nullify(e_h%ovlp_real_ccs)
   if(associated(e_h%ovlp_cmplx_ccs))  nullify(e_h%ovlp_cmplx_ccs)
   if(associated(e_h%dm_real_ccs))     nullify(e_h%dm_real_ccs)
   if(associated(e_h%dm_cmplx_ccs))    nullify(e_h%dm_cmplx_ccs)
   if(associated(e_h%row_ind_ccs))     nullify(e_h%row_ind_ccs)
   if(associated(e_h%col_ptr_ccs))     nullify(e_h%col_ptr_ccs)

   ! ELPA
   if(allocated(e_h%ham_real_elpa))    call elsi_deallocate(e_h,e_h%ham_real_elpa,"ham_real_elpa")
   if(allocated(e_h%ham_cmplx_elpa))   call elsi_deallocate(e_h,e_h%ham_cmplx_elpa,"ham_cmplx_elpa")
   if(allocated(e_h%ovlp_real_elpa))   call elsi_deallocate(e_h,e_h%ovlp_real_elpa,"ovlp_real_elpa")
   if(allocated(e_h%ovlp_cmplx_elpa))  call elsi_deallocate(e_h,e_h%ovlp_cmplx_elpa,"ovlp_cmplx_elpa")
   if(allocated(e_h%eval_elpa))        call elsi_deallocate(e_h,e_h%eval_elpa,"eval_elpa")
   if(allocated(e_h%evec_real_elpa))   call elsi_deallocate(e_h,e_h%evec_real_elpa,"evec_real_elpa")
   if(allocated(e_h%evec_cmplx_elpa))  call elsi_deallocate(e_h,e_h%evec_cmplx_elpa,"evec_cmplx_elpa")
   if(allocated(e_h%dm_real_elpa))     call elsi_deallocate(e_h,e_h%dm_real_elpa,"dm_real_elpa")
   if(allocated(e_h%dm_cmplx_elpa))    call elsi_deallocate(e_h,e_h%dm_cmplx_elpa,"dm_cmplx_elpa")
   if(allocated(e_h%occ_num))          call elsi_deallocate(e_h,e_h%occ_num,"occ_num")
   if(allocated(e_h%k_weight))         call elsi_deallocate(e_h,e_h%k_weight,"k_weight")
   if(allocated(e_h%eval_all))         call elsi_deallocate(e_h,e_h%eval_all,"eval_all")

   ! libOMM
   if(e_h%ham_omm%is_initialized)      call m_deallocate(e_h%ham_omm)
   if(e_h%ovlp_omm%is_initialized)     call m_deallocate(e_h%ovlp_omm)
   if(e_h%dm_omm%is_initialized)       call m_deallocate(e_h%dm_omm)
   if(e_h%coeff%is_initialized)        call m_deallocate(e_h%coeff)
   if(e_h%tdm_omm%is_initialized)      call m_deallocate(e_h%tdm_omm)
   if(allocated(e_h%ovlp_real_omm))    call elsi_deallocate(e_h,e_h%ovlp_real_omm,"ovlp_real_omm")
   if(allocated(e_h%ovlp_cmplx_omm))   call elsi_deallocate(e_h,e_h%ovlp_cmplx_omm,"ovlp_cmplx_omm")

   ! PEXSI
   if(allocated(e_h%ham_real_pexsi))   call elsi_deallocate(e_h,e_h%ham_real_pexsi,"ham_real_pexsi")
   if(allocated(e_h%ham_cmplx_pexsi))  call elsi_deallocate(e_h,e_h%ham_cmplx_pexsi,"ham_cmplx_pexsi")
   if(allocated(e_h%ovlp_real_pexsi))  call elsi_deallocate(e_h,e_h%ovlp_real_pexsi,"ovlp_real_pexsi")
   if(allocated(e_h%ovlp_cmplx_pexsi)) call elsi_deallocate(e_h,e_h%ovlp_cmplx_pexsi,"ovlp_cmplx_pexsi")
   if(allocated(e_h%dm_real_pexsi))    call elsi_deallocate(e_h,e_h%dm_real_pexsi,"dm_real_pexsi")
   if(allocated(e_h%dm_cmplx_pexsi))   call elsi_deallocate(e_h,e_h%dm_cmplx_pexsi,"dm_cmplx_pexsi")
   if(allocated(e_h%row_ind_pexsi))    call elsi_deallocate(e_h,e_h%row_ind_pexsi,"row_ind_pexsi")
   if(allocated(e_h%col_ptr_pexsi))    call elsi_deallocate(e_h,e_h%col_ptr_pexsi,"col_ptr_pexsi")
   if(allocated(e_h%ne_vec))           call elsi_deallocate(e_h,e_h%ne_vec,"ne_vec")

   ! CheSS
   if(allocated(e_h%ham_real_chess))   call elsi_deallocate(e_h,e_h%ham_real_chess,"ham_real_chess")
   if(allocated(e_h%ovlp_real_chess))  call elsi_deallocate(e_h,e_h%ovlp_real_chess,"ovlp_real_chess")
   if(allocated(e_h%row_ind_chess))    call elsi_deallocate(e_h,e_h%row_ind_chess,"row_ind_chess")
   if(allocated(e_h%col_ptr_chess))    call elsi_deallocate(e_h,e_h%col_ptr_chess,"col_ptr_chess")
   if(allocated(e_h%row_ind_buf))      call elsi_deallocate(e_h,e_h%row_ind_buf,"row_ind_buf")
   if(allocated(e_h%col_ptr_buf))      call elsi_deallocate(e_h,e_h%col_ptr_buf,"col_ptr_buf")

   ! SIPs
   if(allocated(e_h%ham_real_sips))    call elsi_deallocate(e_h,e_h%ham_real_sips,"ham_real_sips")
   if(allocated(e_h%ham_cmplx_sips))   call elsi_deallocate(e_h,e_h%ham_cmplx_sips,"ham_cmplx_sips")
   if(allocated(e_h%ovlp_real_sips))   call elsi_deallocate(e_h,e_h%ovlp_real_sips,"ovlp_real_sips")
   if(allocated(e_h%ovlp_cmplx_sips))  call elsi_deallocate(e_h,e_h%ovlp_cmplx_sips,"ovlp_cmplx_sips")
   if(allocated(e_h%dm_real_sips))     call elsi_deallocate(e_h,e_h%dm_real_sips,"dm_real_sips")
   if(allocated(e_h%dm_cmplx_sips))    call elsi_deallocate(e_h,e_h%dm_cmplx_sips,"dm_cmplx_sips")
   if(allocated(e_h%row_ind_sips))     call elsi_deallocate(e_h,e_h%row_ind_sips,"row_ind_sips")
   if(allocated(e_h%col_ptr_sips))     call elsi_deallocate(e_h,e_h%col_ptr_sips,"col_ptr_sips")
   if(allocated(e_h%slices))           call elsi_deallocate(e_h,e_h%slices,"slices")

   if(allocated(e_h%loc_row))          call elsi_deallocate(e_h,e_h%loc_row,"loc_row")
   if(allocated(e_h%loc_col))          call elsi_deallocate(e_h,e_h%loc_col,"loc_col")

   ! Finalize PEXSI
   if(e_h%pexsi_started) then
      call f_ppexsi_plan_finalize(e_h%pexsi_plan,ierr)
      call MPI_Comm_free(e_h%comm_among_pole,ierr)
      call MPI_Comm_free(e_h%comm_in_pole,ierr)
      call MPI_Comm_free(e_h%comm_among_point,ierr)
      call MPI_Comm_free(e_h%comm_in_point,ierr)
   endif

   ! Finalize CheSS
   if(e_h%chess_started) then
      call deallocate_sparse_matrix(e_h%sparse_mat(1))
      call deallocate_sparse_matrix(e_h%sparse_mat(2))
      call foe_data_deallocate(e_h%ice_obj)
      call foe_data_deallocate(e_h%foe_obj)
      call deallocate_matrices(e_h%ham_chess)
      call deallocate_matrices(e_h%ovlp_chess)
      call deallocate_matrices(e_h%dm_chess)
      call deallocate_matrices(e_h%edm_chess)
      call deallocate_matrices(e_h%ovlp_inv_sqrt(1))
      call f_lib_finalize()
   endif

   ! Finalize SIPs
   if(e_h%sips_started) then
      call clean_qetsc()
   endif

   ! Reset e_h
   call elsi_reset_handle(e_h)

end subroutine

!>
!! This routine resets an ELSI handle.
!!
subroutine elsi_reset_handle(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*40, parameter :: caller = "elsi_reset_handle"

   e_h%handle_ready     = .false.
   e_h%solver           = UNSET
   e_h%matrix_data_type = UNSET
   e_h%matrix_format    = UNSET
   e_h%uplo             = FULL_MAT
   e_h%parallel_mode    = UNSET
   e_h%n_elsi_calls     = 0
   e_h%n_b_rows         = UNSET
   e_h%n_b_cols         = UNSET
   e_h%n_p_rows         = UNSET
   e_h%n_p_cols         = UNSET
   e_h%n_l_rows         = UNSET
   e_h%n_l_cols         = UNSET
   e_h%myid             = UNSET
   e_h%myid_all         = UNSET
   e_h%n_procs          = UNSET
   e_h%n_procs_all      = UNSET
   e_h%mpi_comm         = UNSET
   e_h%mpi_comm_all     = UNSET
   e_h%mpi_comm_row     = UNSET
   e_h%mpi_comm_col     = UNSET
   e_h%my_p_row         = UNSET
   e_h%my_p_col         = UNSET
   e_h%mpi_ready        = .false.
   e_h%global_mpi_ready = .false.
   e_h%blacs_ctxt       = UNSET
   e_h%sc_desc          = UNSET
   e_h%blacs_ready      = .false.
   e_h%nnz_g            = UNSET
   e_h%nnz_l            = UNSET
   e_h%nnz_l_sp         = UNSET
   e_h%n_l_cols_sp      = UNSET
   e_h%zero_threshold   = 1.0e-13_r8
   e_h%sparsity_ready   = .false.
   e_h%ovlp_is_unit     = .false.
   e_h%ovlp_is_sing     = .false.
   e_h%no_sing_check    = .false.
   e_h%sing_tol         = 1.0e-5_r8
   e_h%stop_sing        = .false.
   e_h%n_nonsing        = UNSET
   e_h%n_electrons      = 0.0_r8
   e_h%mu               = 0.0_r8
   e_h%n_basis          = UNSET
   e_h%n_spins          = 1
   e_h%n_kpts           = 1
   e_h%n_states         = UNSET
   e_h%n_states_solve   = UNSET
   e_h%i_spin           = 1
   e_h%i_kpt            = 1
   e_h%i_weight         = 1.0_r8
   e_h%energy_hdm       = 0.0_r8
   e_h%energy_sedm      = 0.0_r8
   e_h%free_energy      = 0.0_r8
   e_h%broaden_scheme   = 0
   e_h%broaden_width    = 1.0e-2_r8
   e_h%occ_tolerance    = 1.0e-13_r8
   e_h%max_mu_steps     = 100
   e_h%spin_degen       = 0.0_r8
   e_h%spin_is_set      = .false.
   e_h%mu_ready         = .false.
   e_h%edm_ready_real   = .false.
   e_h%edm_ready_cmplx  = .false.
   e_h%elpa_solver      = UNSET
   e_h%elpa_output      = .false.
   e_h%n_states_omm     = UNSET
   e_h%n_elpa_steps     = UNSET
   e_h%new_overlap      = .true.
   e_h%coeff_ready      = .false.
   e_h%omm_flavor       = UNSET
   e_h%scale_kinetic    = 0.0_r8
   e_h%calc_ed          = .false.
   e_h%eta              = 0.0_r8
   e_h%min_tol          = 1.0e-12_r8
   e_h%omm_output       = .false.
   e_h%do_dealloc       = .false.
   e_h%use_psp          = .false.
   e_h%n_p_per_pole     = UNSET
   e_h%n_p_per_point    = UNSET
   e_h%my_p_row_pexsi   = UNSET
   e_h%my_p_col_pexsi   = UNSET
   e_h%n_p_rows_pexsi   = UNSET
   e_h%n_p_cols_pexsi   = UNSET
   e_h%my_point         = UNSET
   e_h%myid_point       = UNSET
   e_h%comm_among_pole  = UNSET
   e_h%comm_in_pole     = UNSET
   e_h%comm_among_point = UNSET
   e_h%comm_in_point    = UNSET
   e_h%ne_pexsi         = 0.0_r8
   e_h%pexsi_started    = .false.
   e_h%erf_decay        = 0.0_r8
   e_h%erf_decay_min    = 0.0_r8
   e_h%erf_decay_max    = 0.0_r8
   e_h%ev_ham_min       = 0.0_r8
   e_h%ev_ham_max       = 0.0_r8
   e_h%ev_ovlp_min      = 0.0_r8
   e_h%ev_ovlp_max      = 0.0_r8
   e_h%beta             = 0.0_r8
   e_h%chess_started    = .false.
   e_h%n_p_per_slice    = UNSET
   e_h%n_inertia_steps  = UNSET
   e_h%slicing_method   = UNSET
   e_h%inertia_option   = UNSET
   e_h%unbound          = UNSET
   e_h%n_slices         = UNSET
   e_h%interval         = 0.0_r8
   e_h%slice_buffer     = 0.0_r8
   e_h%sips_started     = .false.
   e_h%clock_rate       = UNSET

end subroutine

!>
!! This routine guarantees that there are no mutually conflicting parameters.
!!
subroutine elsi_check(e_h,caller)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   character(len=*),  intent(in)    :: caller !< Caller

   ! General check of solver, parallel mode, data type, matrix format
   if(e_h%solver < 0 .or. e_h%solver >= N_SOLVERS) then
      call elsi_stop(" An unsupported solver has been chosen. Please choose"//&
              " ELPA, LIBOMM, or PEXSI solver. Exiting...",e_h,caller)
   endif

   if(e_h%parallel_mode < 0 .or. e_h%parallel_mode >= N_PARALLEL_MODES) then
      call elsi_stop(" An unsupported parallel mode has been chosen. Please"//&
              " choose either SINGLE_PROC or MULTI_PROC. Exiting...",e_h,caller)
   endif

   if(e_h%matrix_data_type < 0 .or. &
      e_h%matrix_data_type >= N_MATRIX_DATA_TYPES) then
      call elsi_stop(" An unsupported matirx data type has been chosen."//&
              " Please choose either REAL_VALUES or COMPLEX_VALUES."//&
              " Exiting...",e_h,caller)
   endif

   if(e_h%matrix_format < 0 .or. &
      e_h%matrix_format >= N_MATRIX_STORAGE_FORMATS) then
      call elsi_stop(" An unsupported matirx storage format has been set."//&
              " Please choose either BLACS_DENSE or PEXSI_CSC. Exiting...",e_h,&
              caller)
   endif

   if(e_h%uplo /= FULL_MAT) then
      if(e_h%matrix_format /= BLACS_DENSE .or. &
         e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(" Upper/lower triangular input matrix only supported"//&
                 " with BLACS_DENSE matrix storage format and MULTI_PROC"//&
                 " parallel mode. Exiting...",e_h,caller)
      endif

      if(e_h%uplo /= UT_MAT .and. e_h%uplo /= LT_MAT) then
         call elsi_stop(" Invalid choice of uplo. Exiting...",e_h,caller)
      endif
   endif

   ! Spin and k-point
   if(e_h%n_spins*e_h%n_kpts > 1) then
      if(.not. e_h%global_mpi_ready) then
         call elsi_stop(" Spin/k-point calculations require a global MPI"//&
                 " communicator. Exiting...",e_h,caller)
      endif
   endif

   if(.not. e_h%global_mpi_ready) then
      e_h%mpi_comm_all = e_h%mpi_comm
      e_h%n_procs_all  = e_h%n_procs
      e_h%myid_all     = e_h%myid
   endif

   ! Specific check for each solver
   if(e_h%solver == AUTO) then
      call elsi_stop(" AUTO not yet available. Please choose ELPA, LIBOMM,"//&
               " or PEXSI solver. Exiting...",e_h,caller)
   elseif(e_h%solver == ELPA) then
      if(e_h%parallel_mode == MULTI_PROC) then
         if(.not. e_h%mpi_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI being set"//&
                    " up before calling the solver. Exiting...",e_h,caller)
         endif

         if(.not. e_h%blacs_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires BLACS being"//&
                    " set up before calling the solver. Exiting...",e_h,caller)
         endif
      endif

      if(e_h%matrix_format == PEXSI_CSC) then
         if(.not. e_h%sparsity_ready) then
            call elsi_stop(" The PEXSI_CSC format has been chosen. Please"//&
                    " set the sparsity pattern before calling the solver."//&
                    " Exiting...",e_h,caller)
         endif
      endif
   elseif(e_h%solver == LIBOMM) then
      if(e_h%parallel_mode == MULTI_PROC) then
         if(.not. e_h%mpi_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI being set"//&
                    " up before calling the solver. Exiting...",e_h,caller)
         endif

         if(.not. e_h%blacs_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires BLACS being"//&
                    " set up before calling the solver. Exiting...",e_h,caller)
         endif
      else
         call elsi_stop(" libOMM has been chosen as the solver. Please"//&
                 " choose MULTI_PROC parallel mode. Exiting...",e_h,caller)
      endif

      if(e_h%matrix_format == PEXSI_CSC) then
         if(.not. e_h%sparsity_ready) then
            call elsi_stop(" The PEXSI_CSC format has been chosen. Please"//&
                    " set the sparsity pattern before calling the solver."//&
                    " Exiting...",e_h,caller)
         endif
      endif
   elseif(e_h%solver == PEXSI) then
      if(e_h%parallel_mode == MULTI_PROC) then
         if(.not. e_h%mpi_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI being set"//&
                    " up before calling the solver. Exiting...",e_h,caller)
         endif
      else
         call elsi_stop(" PEXSI has been chosen as the solver. Please choose"//&
                 " MULTI_PROC parallel mode. Exiting...",e_h,caller)
      endif

      if(e_h%matrix_format == BLACS_DENSE) then
         if(.not. e_h%blacs_ready) then
            call elsi_stop(" The BLACS_DENSE format has been chosen. Please"//&
                    " set up BLACS before calling the solver. Exiting...",e_h,&
                    caller)
         endif
      else
         if(.not. e_h%sparsity_ready) then
            call elsi_stop(" The PEXSI_CSC format has been chosen. Please"//&
                    " set the sparsity pattern before calling the solver."//&
                    " Exiting...",e_h,caller)
         endif
      endif

      if(e_h%n_basis < e_h%n_p_per_pole) then
         call elsi_stop(" PEXSI has been chosen as the solver. The matrix"//&
                 " size is too small for this number of MPI tasks. Exiting...",&
                 e_h,caller)
      endif

      if(e_h%n_p_per_pole == UNSET) then
         if(mod(e_h%n_procs,e_h%pexsi_options%numPole*&
            e_h%pexsi_options%nPoints) /= 0) then
            call elsi_stop("  The total number of MPI tasks must be a"//&
                    " multiple of the number of MPI tasks per pole times the"//&
                    " number of mu points. Exiting...",e_h,caller)
         endif
      else
         if(mod(e_h%n_procs,e_h%n_p_per_pole*&
            e_h%pexsi_options%nPoints) /= 0) then
            call elsi_stop("  The total number of MPI tasks must be a"//&
                    " multiple of the number of MPI tasks per pole times the"//&
                    " number of mu points. Exiting...",e_h,caller)
         endif

         if(e_h%n_p_per_pole*e_h%pexsi_options%numPole*&
            e_h%pexsi_options%nPoints < e_h%n_procs) then
            call elsi_stop("  Specified number of MPI tasks per pole is too"//&
                    " small for the total number of MPI tasks. Exiting...",e_h,&
                    caller)
         endif
      endif
   elseif(e_h%solver == CHESS) then
      call elsi_statement_print("  ATTENTION! CheSS is EXPERIMENTAL.",e_h)

      if(e_h%n_basis < e_h%n_procs) then
         call elsi_stop(" CheSS has been chosen as the solver. The matrix"//&
                 " size is too small for this number of MPI tasks. Exiting...",&
                 e_h,caller)
      endif

      if(e_h%parallel_mode == MULTI_PROC) then
         if(.not. e_h%mpi_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI being set"//&
                    " up before calling the solver. Exiting...",e_h,caller)
         endif
      else
         call elsi_stop(" CheSS has been chosen as the solver. Please choose"//&
                 " MULTI_PROC parallel mode. Exiting...",e_h,caller)
      endif

      if(e_h%matrix_format == BLACS_DENSE) then
         if(.not. e_h%blacs_ready) then
            call elsi_stop(" The BLACS_DENSE format has been chosen. Please"//&
                    " set up BLACS before calling the solver. Exiting...",e_h,&
                    caller)
         endif

         if(e_h%ovlp_is_unit) then
            call elsi_stop(" CheSS with an identity overlap matrix not yet"//&
                    " available. Exiting...",e_h,caller)
         endif
      else
         call elsi_stop(" CheSS has been chosen as the solver. Please choose"//&
                 " BLACS_DENSE matrix format. Exiting...",e_h,caller)
      endif
   elseif(e_h%solver == SIPS) then
      call elsi_statement_print("  ATTENTION! SIPs is EXPERIMENTAL.",e_h)

      if(e_h%n_basis < e_h%n_procs) then
         call elsi_stop(" SIPs has been chosen as the solver. The matrix"//&
                 " size is too small for this number of MPI tasks. Exiting...",&
                 e_h,caller)
      endif

      if(e_h%parallel_mode == MULTI_PROC) then
         if(.not. e_h%mpi_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI being set"//&
                    " up before calling the solver. Exiting...",e_h,caller)
         endif
      else
         call elsi_stop(" SIPs has been chosen as the solver. Please choose"//&
                 " MULTI_PROC parallel mode. Exiting...",e_h,caller)
      endif

      if(e_h%matrix_format == BLACS_DENSE) then
         if(.not. e_h%blacs_ready) then
            call elsi_stop(" The BLACS_DENSE format has been chosen. Please"//&
                    " set up BLACS before calling the solver. Exiting...",e_h,&
                    caller)
         endif
      else
         if(.not. e_h%sparsity_ready) then
            call elsi_stop(" The PEXSI_CSC format has been chosen. Please"//&
                    " set the sparsity pattern before calling the solver."//&
                    " Exiting...",e_h,caller)
         endif
      endif
   else
      call elsi_stop(" No supported solver has been chosen. Please choose"//&
              " ELPA, LIBOMM, or PEXSI solver. Exiting...",e_h,caller)
   endif

end subroutine

!>
!! This routine checks whether the input handle has been properly
!! initialized by the user or not.
!!
!! Do NOT call this subroutine outside the (main) ELSI module.
!!
subroutine elsi_check_handle(e_h,caller)

   implicit none

   type(elsi_handle), intent(in) :: e_h    !< Handle
   character(len=*),  intent(in) :: caller !< Caller

   if(.not. e_h%handle_ready) then
      call elsi_stop(" Invalid handle! An ELSI handle must be properly"//&
              " initialized with 'elsi_init' before being used. Exiting...",&
              e_h,caller)
   endif

end subroutine

!>
!! This routine computes the global row index based on the local row index.
!!
subroutine elsi_get_global_row(e_h,global_idx,local_idx)

   implicit none

   type(elsi_handle), intent(in)  :: e_h        !< Handle
   integer(kind=i4),  intent(in)  :: local_idx  !< Local index
   integer(kind=i4),  intent(out) :: global_idx !< Global index

   integer(kind=i4) :: block
   integer(kind=i4) :: idx

   character*40, parameter :: caller = "elsi_get_global_row"

   block = (local_idx-1)/e_h%n_b_rows
   idx = local_idx-block*e_h%n_b_rows

   global_idx = e_h%my_p_row*e_h%n_b_rows+block*e_h%n_b_rows*e_h%n_p_rows+idx

end subroutine

!>
!! This routine computes the global column index based on the local column index.
!!
subroutine elsi_get_global_col(e_h,global_idx,local_idx)

   implicit none

   type(elsi_handle), intent(in)  :: e_h        !< Handle
   integer(kind=i4),  intent(in)  :: local_idx  !< Local index
   integer(kind=i4),  intent(out) :: global_idx !< Global index

   integer(kind=i4) :: block
   integer(kind=i4) :: idx

   character*40, parameter :: caller = "elsi_get_global_col"

   block = (local_idx-1)/e_h%n_b_cols
   idx = local_idx-block*e_h%n_b_cols

   global_idx = e_h%my_p_col*e_h%n_b_cols+block*e_h%n_b_cols*e_h%n_p_cols+idx

end subroutine

!>
!! This routine counts the local number of non_zero elements.
!!
subroutine elsi_get_local_nnz_real(e_h,matrix,n_rows,n_cols,nnz)

   implicit none

   type(elsi_handle), intent(in)  :: e_h                   !< Handle
   real(kind=r8),     intent(in)  :: matrix(n_rows,n_cols) !< Local matrix
   integer(kind=i4),  intent(in)  :: n_rows                !< Local rows
   integer(kind=i4),  intent(in)  :: n_cols                !< Local cols
   integer(kind=i4),  intent(out) :: nnz                   !< Number of non-zero

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col

   character*40, parameter :: caller = "elsi_get_local_nnz_real"

   nnz = 0

   do i_col = 1,n_cols
      do i_row = 1,n_rows
         if(abs(matrix(i_row,i_col)) > e_h%zero_threshold) then
            nnz = nnz+1
         endif
      enddo
   enddo

end subroutine

!>
!! This routine counts the local number of non_zero elements.
!!
subroutine elsi_get_local_nnz_complex(e_h,matrix,n_rows,n_cols,nnz)

   implicit none

   type(elsi_handle), intent(in)  :: e_h                   !< Handle
   complex(kind=r8),  intent(in)  :: matrix(n_rows,n_cols) !< Local matrix
   integer(kind=i4),  intent(in)  :: n_rows                !< Local rows
   integer(kind=i4),  intent(in)  :: n_cols                !< Local cols
   integer(kind=i4),  intent(out) :: nnz                   !< Number of non-zero

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col

   character*40, parameter :: caller = "elsi_get_local_nnz_complex"

   nnz = 0

   do i_col = 1,n_cols
      do i_row = 1,n_rows
         if(abs(matrix(i_row,i_col)) > e_h%zero_threshold) then
            nnz = nnz+1
         endif
      enddo
   enddo

end subroutine

!>
!! This routine sets a full matrix from a (upper or lower) triangular matrix.
!! The size of matrix should be the same as the Hamiltonian matrix.
!!
subroutine elsi_set_full_mat_real(e_h,mat)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                            !< Handle
   real(kind=r8),     intent(inout) :: mat(e_h%n_l_rows,e_h%n_l_cols) !< Matrix

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   real(kind=r8), allocatable :: tmp_real(:,:)

   character*40, parameter :: caller = "elsi_set_full_mat_real"

   if(e_h%uplo /= FULL_MAT .and. e_h%parallel_mode == MULTI_PROC) then
      call elsi_allocate(e_h,tmp_real,e_h%n_l_rows,e_h%n_l_cols,"tmp_real",&
              caller)

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

   type(elsi_handle), intent(inout) :: e_h                            !< Handle
   complex(kind=r8),  intent(inout) :: mat(e_h%n_l_rows,e_h%n_l_cols) !< Matrix

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character*40, parameter :: caller = "elsi_set_full_mat_complex"

   if(e_h%uplo /= FULL_MAT .and. e_h%parallel_mode == MULTI_PROC) then
      call elsi_allocate(e_h,tmp_cmplx,e_h%n_l_rows,e_h%n_l_cols,"tmp_cmplx",&
              caller)

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
            dble(mat(e_h%loc_row(i_col),e_h%loc_col(i_col)))
      enddo
   endif

end subroutine

!>
!! This routine initializes the timer.
!!
subroutine elsi_init_timers(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: initial_time
   integer(kind=i4) :: clock_max

   character*40, parameter :: caller = "elsi_init_timers"

   call system_clock(initial_time,e_h%clock_rate,clock_max)

end subroutine

!>
!! This routine gets the current wallclock time.
!!
subroutine elsi_get_time(e_h,wtime)

   implicit none

   type(elsi_handle), intent(in)  :: e_h   !< Handle
   real(kind=r8),     intent(out) :: wtime !< Time

   integer(kind=i4) :: tics

   character*40, parameter :: caller = "elsi_get_time"

   call system_clock(tics)

   wtime = 1.0_r8*tics/e_h%clock_rate

end subroutine

!>
!! This routine prints a final output.
!!
subroutine elsi_final_print(e_h)

   implicit none

   type(elsi_handle), intent(in) :: e_h !< Handle

   real(kind=r8) :: sparsity

   character*40, parameter :: caller = "elsi_final_print"

   if(print_info) then
      if(e_h%myid_all == 0) then
         write(*,"('  |------------------------------------------')")
         write(*,"('  | Final ELSI Output')")
         write(*,"('  |------------------------------------------')")

         write(*,"('  | Number of basis functions :',I13)") e_h%n_basis
         if(e_h%solver == PEXSI .or. e_h%solver == SIPS) then
            write(*,"('  | Number of nonzeros        :',I13)") e_h%nnz_g

            sparsity = 1.0_r8-(1.0_r8*e_h%nnz_g/e_h%n_basis/e_h%n_basis)
            write(*,"('  | Sparsity                  :',F13.3)") sparsity
         endif
         write(*,"('  | Number of electrons       :',F13.1)") e_h%n_electrons
         write(*,"('  | Number of states          :',I13)") e_h%n_states
         write(*,"('  | Number of spins           :',I13)") e_h%n_spins
         write(*,"('  | Number of k-points        :',I13)") e_h%n_kpts

         if(e_h%solver == ELPA) then
            write(*,"('  | Solver                    :',A13)") "ELPA"
         elseif(e_h%solver == LIBOMM) then
            write(*,"('  | Solver                    :',A13)") "libOMM"
         elseif(e_h%solver == PEXSI) then
            write(*,"('  | Solver                    :',A13)") "PEXSI"
         elseif(e_h%solver == CHESS) then
            write(*,"('  | Solver                    :',A13)") "CheSS"
         elseif(e_h%solver == SIPS) then
            write(*,"('  | Solver                    :',A13)") "SIPs"
         endif

         if(e_h%parallel_mode == MULTI_PROC) then
            write(*,"('  | Parallel mode             :',A13)") "MULTI_PROC"
         elseif(e_h%parallel_mode == SINGLE_PROC) then
            write(*,"('  | Parallel mode             :',A13)") "SINGLE_PROC"
         endif

         if(e_h%matrix_format == BLACS_DENSE) then
            write(*,"('  | Matrix format             :',A13)") "BLACS_DENSE"
         elseif(e_h%matrix_format == PEXSI_CSC) then
            write(*,"('  | Matrix format             :',A13)") "PEXSI_CSC"
         endif

         write(*,"('  |------------------------------------------')")
         write(*,"('  | ELSI Project (c)  elsi-interchange.org')")
         write(*,"('  |------------------------------------------')")
      endif
   endif

end subroutine

end module ELSI_UTILS
