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

   use ELSI_CONSTANTS, only: AUTO,ELPA,LIBOMM,PEXSI,CHESS,SIPS,BLACS_DENSE,&
                             PEXSI_CSC,MULTI_PROC,SINGLE_PROC,FULL_MAT,UT_MAT,&
                             LT_MAT,N_SOLVERS,N_MATRIX_DATA_TYPES,&
                             N_MATRIX_STORAGE_FORMATS,N_PARALLEL_MODES,UNSET
   use ELSI_DATATYPE, only: elsi_handle,print_info,print_mem
   use ELSI_PRECISION, only: r8,i4
   use FOE_BASE, only: foe_data_deallocate
   use F_PPEXSI_INTERFACE
   use MATRIXSWITCH, only: matrix,m_register_pdbc,m_deallocate
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
      module procedure elsi_allocate_integer_1d,&
                       elsi_allocate_integer_2d,&
                       elsi_allocate_integer_3d,&
                       elsi_allocate_real_1d,&
                       elsi_allocate_real_2d,&
                       elsi_allocate_real_3d,&
                       elsi_allocate_complex_1d,&
                       elsi_allocate_complex_2d,&
                       elsi_allocate_complex_3d
   end interface

   interface elsi_deallocate
      module procedure elsi_deallocate_integer_1d,&
                       elsi_deallocate_integer_2d,&
                       elsi_deallocate_integer_3d,&
                       elsi_deallocate_real_1d,&
                       elsi_deallocate_real_2d,&
                       elsi_deallocate_real_3d,&
                       elsi_deallocate_complex_1d,&
                       elsi_deallocate_complex_2d,&
                       elsi_deallocate_complex_3d
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
subroutine elsi_statement_print(info,elsi_h)

   implicit none

   character(len=*),  intent(in) :: info   !< Message to print
   type(elsi_handle), intent(in) :: elsi_h !< Handle

   if(print_info) then
      if(elsi_h%myid_all == 0) then
         write(*,'(A)') trim(info)
      endif
   endif

end subroutine

!>
!! This routine allocates a 1D array with real(kind=r8).
!!
subroutine elsi_allocate_real_1d(elsi_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle),          intent(in)    :: elsi_h    !< Handle
   real(kind=r8), allocatable, intent(inout) :: array(:)  !< Data
   integer(kind=i4),           intent(in)    :: dim1      !< Size
   character(len=*),           intent(in)    :: arrayname !< Name
   character(len=*),           intent(in)    :: caller    !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*8

      write(info,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info,elsi_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 1D array with integer.
!!
subroutine elsi_allocate_integer_1d(elsi_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h    !< Handle
   integer(kind=i4), allocatable, intent(inout) :: array(:)  !< Data
   integer(kind=i4),              intent(in)    :: dim1      !< Size
   character(len=*),              intent(in)    :: arrayname !< Name
   character(len=*),              intent(in)    :: caller    !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*4

      write(info,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info,elsi_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 1D array with complex(kind=r8).
!!
subroutine elsi_allocate_complex_1d(elsi_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h    !< Handle
   complex(kind=r8), allocatable, intent(inout) :: array(:)  !< Data
   integer(kind=i4),              intent(in)    :: dim1      !< Size
   character(len=*),              intent(in)    :: arrayname !< Name
   character(len=*),              intent(in)    :: caller    !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*16

      write(info,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info,elsi_h,caller)
   endif

   array = cmplx(0.0_r8,0.0_r8)

end subroutine

!>
!! This routine allocates a 2D array of real(kind=r8).
!!
subroutine elsi_allocate_real_2d(elsi_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle),          intent(in)    :: elsi_h     !< Handle
   real(kind=r8), allocatable, intent(inout) :: array(:,:) !< Data
   integer(kind=i4),           intent(in)    :: dim1       !< Size
   integer(kind=i4),           intent(in)    :: dim2       !< Size
   character(len=*),           intent(in)    :: arrayname  !< Name
   character(len=*),           intent(in)    :: caller     !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*8

      write(info,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(info,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info,elsi_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 2D array of integer.
!!
subroutine elsi_allocate_integer_2d(elsi_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h     !< Handle
   integer(kind=i4), allocatable, intent(inout) :: array(:,:) !< Data
   integer(kind=i4),              intent(in)    :: dim1       !< Size
   integer(kind=i4),              intent(in)    :: dim2       !< Size
   character(len=*),              intent(in)    :: arrayname  !< Name
   character(len=*),              intent(in)    :: caller     !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*4

      write(info,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(info,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info,elsi_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 2D array of complex(kind=r8).
!!
subroutine elsi_allocate_complex_2d(elsi_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h     !< Handle
   complex(kind=r8), allocatable, intent(inout) :: array(:,:) !< Data
   integer(kind=i4),              intent(in)    :: dim1       !< Size
   integer(kind=i4),              intent(in)    :: dim2       !< Size
   character(len=*),              intent(in)    :: arrayname  !< Name
   character(len=*),              intent(in)    :: caller     !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*16

      write(info,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(info,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info,elsi_h,caller)
   endif

   array = cmplx(0.0_r8,0.0_r8)

end subroutine

!>
!! This routine allocates a 3D array of real(kind=r8).
!!
subroutine elsi_allocate_real_3d(elsi_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   type(elsi_handle),          intent(in)    :: elsi_h       !< Handle
   real(kind=r8), allocatable, intent(inout) :: array(:,:,:) !< Data
   integer(kind=i4),           intent(in)    :: dim1         !< Size
   integer(kind=i4),           intent(in)    :: dim2         !< Size
   integer(kind=i4),           intent(in)    :: dim3         !< Size
   character(len=*),           intent(in)    :: arrayname    !< Name
   character(len=*),           intent(in)    :: caller       !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*dim3*8

      write(info,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   allocate(array(dim1,dim2,dim3),stat=error)

   if(error > 0) then
      write(info,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info,elsi_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 3D array of integer.
!!
subroutine elsi_allocate_integer_3d(elsi_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h       !< Handle
   integer(kind=i4), allocatable, intent(inout) :: array(:,:,:) !< Data
   integer(kind=i4),              intent(in)    :: dim1         !< Size
   integer(kind=i4),              intent(in)    :: dim2         !< Size
   integer(kind=i4),              intent(in)    :: dim3         !< Size
   character(len=*),              intent(in)    :: arrayname    !< Name
   character(len=*),              intent(in)    :: caller       !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*dim3*4

      write(info,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   allocate(array(dim1,dim2,dim3),stat=error)

   if(error > 0) then
      write(info,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info,elsi_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 3D array of complex(kind=r8).
!!
subroutine elsi_allocate_complex_3d(elsi_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h       !< Handle
   complex(kind=r8), allocatable, intent(inout) :: array(:,:,:) !< Data
   integer(kind=i4),              intent(in)    :: dim1         !< Size
   integer(kind=i4),              intent(in)    :: dim2         !< Size
   integer(kind=i4),              intent(in)    :: dim3         !< Size
   character(len=*),              intent(in)    :: arrayname    !< Name
   character(len=*),              intent(in)    :: caller       !< Caller

   real(kind=r8) :: arraysize
   integer       :: error
   character*200 :: info

   if(print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*dim3*16

      write(info,"(A,F12.3,A,A)") "    Allocating ",arraysize," MB for ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   allocate(array(dim1,dim2,dim3),stat=error)

   if(error > 0) then
      write(info,"(A,A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info,elsi_h,caller)
   endif

   array = cmplx(0.0_r8,0.0_r8)

end subroutine

!>
!! This routine deallocates a 1D array with real(kind=r8).
!!
subroutine elsi_deallocate_real_1d(elsi_h,array,arrayname)

   implicit none

   type(elsi_handle),          intent(in)    :: elsi_h    !< Handle
   real(kind=r8), allocatable, intent(inout) :: array(:)  !< Data
   character(len=*),           intent(in)    :: arrayname !< Name

   character*200 :: info

   if(print_mem) then
      write(info,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 1D array with integer.
!!
subroutine elsi_deallocate_integer_1d(elsi_h,array,arrayname)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h    !< Handle
   integer(kind=i4), allocatable, intent(inout) :: array(:)  !< Data
   character(len=*),              intent(in)    :: arrayname !< Name

   character*200 :: info

   if(print_mem) then
      write(info,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 1D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex_1d(elsi_h,array,arrayname)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h    !< Handle
   complex(kind=r8), allocatable, intent(inout) :: array(:)  !< Data
   character(len=*),              intent(in)    :: arrayname !< Name

   character*200 :: info

   if(print_mem) then
      write(info,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 2D array with real(kind=r8).
!!
subroutine elsi_deallocate_real_2d(elsi_h,array,arrayname)

   implicit none

   type(elsi_handle),          intent(in)    :: elsi_h     !< Handle
   real(kind=r8), allocatable, intent(inout) :: array(:,:) !< Data
   character(len=*),           intent(in)    :: arrayname  !< Name

   character*200 :: info

   if(print_mem) then
      write(info,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 2D array with integer.
!!
subroutine elsi_deallocate_integer_2d(elsi_h,array,arrayname)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h     !< Handle
   integer(kind=i4), allocatable, intent(inout) :: array(:,:) !< Data
   character(len=*),              intent(in)    :: arrayname  !< Name

   character*200 :: info

   if(print_mem) then
      write(info,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 2D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex_2d(elsi_h,array,arrayname)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h     !< Handle
   complex(kind=r8), allocatable, intent(inout) :: array(:,:) !< Data
   character(len=*),              intent(in)    :: arrayname  !< Name

   character*200 :: info

   if(print_mem) then
      write(info,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 3D array with real(kind=r8).
!!
subroutine elsi_deallocate_real_3d(elsi_h,array,arrayname)

   implicit none

   type(elsi_handle),          intent(in)    :: elsi_h       !< Handle
   real(kind=r8), allocatable, intent(inout) :: array(:,:,:) !< Data
   character(len=*),           intent(in)    :: arrayname    !< Name

   character*200 :: info

   if(print_mem) then
      write(info,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 3D array with integer.
!!
subroutine elsi_deallocate_integer_3d(elsi_h,array,arrayname)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h       !< Handle
   integer(kind=i4), allocatable, intent(inout) :: array(:,:,:) !< Data
   character(len=*),              intent(in)    :: arrayname    !< Name

   character*200 :: info

   if(print_mem) then
      write(info,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 3D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex_3d(elsi_h,array,arrayname)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h       !< Handle
   complex(kind=r8), allocatable, intent(inout) :: array(:,:,:) !< Data
   character(len=*),              intent(in)    :: arrayname    !< Name

   character*200 :: info

   if(print_mem) then
      write(info,"(A,A)") "    Deallocating ",trim(arrayname)
      call elsi_statement_print(info,elsi_h)
   endif

   deallocate(array)

end subroutine

!>
!! Clean shutdown in case of errors.
!!
subroutine elsi_stop(info,elsi_h,caller)

   implicit none

   character(len=*),  intent(in) :: info   !< Error message to print
   type(elsi_handle), intent(in) :: elsi_h !< Handle
   character(len=*),  intent(in) :: caller !< The subroutine in trouble

   character(len=4096) :: string_info
   integer :: i_task
   integer :: mpierr

   if(elsi_h%global_mpi_ready) then
      do i_task = 0,elsi_h%n_procs_all-1
         if(elsi_h%myid_all == i_task) then
            write(string_info,"(1X,' Error! MPI task',I5,' in ',A,': ',A)")&
               elsi_h%myid_all,trim(caller),trim(info)
            write(*,'(A)') trim(string_info)
            exit
         endif
         call MPI_Barrier(elsi_h%mpi_comm_all,mpierr)
      enddo
   elseif(elsi_h%mpi_ready) then
      do i_task = 0,elsi_h%n_procs-1
         if(elsi_h%myid == i_task) then
            write(string_info,"(1X,' Error! MPI task',I5,' in ',A,': ',A)")&
               elsi_h%myid,trim(caller),trim(info)
            write(*,'(A)') trim(string_info)
            exit
         endif
         call MPI_Barrier(elsi_h%mpi_comm,mpierr)
      enddo
   else
      write(string_info,"(1X,' Error!',A,': ',A)") trim(caller),trim(info)
      write(*,'(A)') trim(string_info)
   endif

   stop

end subroutine

!>
!! This routine sets the real hamiltonian matrix.
!!
subroutine elsi_set_real_ham(elsi_h,h_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     target        :: h_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_real_ham"

   if(elsi_h%solver == ELPA) then
      if(elsi_h%matrix_format == BLACS_DENSE) then
         call elsi_set_full_mat(elsi_h,h_in)
      endif

      elsi_h%ham_real => h_in
   elseif(elsi_h%solver == LIBOMM) then
      if(elsi_h%matrix_format == BLACS_DENSE) then
         call elsi_set_full_mat(elsi_h,h_in)
      endif

      call m_register_pdbc(elsi_h%ham_omm,h_in,elsi_h%sc_desc)
   endif

end subroutine

!>
!! This routine sets the complex hamiltonian matrix.
!!
subroutine elsi_set_complex_ham(elsi_h,h_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   complex(kind=r8),  target        :: h_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_complex_ham"

   if(elsi_h%solver == ELPA) then
      if(elsi_h%matrix_format == BLACS_DENSE) then
         call elsi_set_full_mat(elsi_h,h_in)
      endif

      elsi_h%ham_complex => h_in
   elseif(elsi_h%solver == LIBOMM) then
      if(elsi_h%matrix_format == BLACS_DENSE) then
         call elsi_set_full_mat(elsi_h,h_in)
      endif

      call m_register_pdbc(elsi_h%ham_omm,h_in,elsi_h%sc_desc)
   endif

end subroutine

!>
!! This routine sets the sparse real Hamiltonian matrix.
!!
subroutine elsi_set_sparse_real_ham(elsi_h,h_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                !< Handle
   real(kind=r8),     target        :: h_in(elsi_h%nnz_l_sp) !< Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_sparse_real_ham"

   if((elsi_h%solver == PEXSI) .or. (elsi_h%solver == SIPS)) then
      elsi_h%ham_real_ccs => h_in
   endif

end subroutine

!>
!! This routine sets the sparse complex Hamiltonian matrix.
!!
subroutine elsi_set_sparse_complex_ham(elsi_h,h_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                !< Handle
   complex(kind=r8),  target        :: h_in(elsi_h%nnz_l_sp) !< Hamiltonian matirx

   character*40, parameter :: caller = "elsi_set_sparse_complex_ham"

   if((elsi_h%solver == PEXSI) .or. (elsi_h%solver == SIPS)) then
      elsi_h%ham_complex_ccs => h_in
   endif

end subroutine

!>
!! This routine sets the real overlap matrix.
!!
subroutine elsi_set_real_ovlp(elsi_h,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     target        :: s_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap matrix

   character*40, parameter :: caller = "elsi_set_real_ovlp"

   if(.not. elsi_h%ovlp_is_unit) then
      if(elsi_h%solver == ELPA) then
         if((elsi_h%matrix_format == BLACS_DENSE) .and. &
            (elsi_h%n_elsi_calls == 1)) then
            call elsi_set_full_mat(elsi_h,s_in)
         endif

         elsi_h%ovlp_real => s_in
      elseif(elsi_h%solver == LIBOMM) then
         if((elsi_h%matrix_format == BLACS_DENSE) .and. &
            (elsi_h%n_elsi_calls == 1)) then
            call elsi_set_full_mat(elsi_h,s_in)
         endif

         call m_register_pdbc(elsi_h%ovlp_omm,s_in,elsi_h%sc_desc)
      endif
   endif

end subroutine

!>
!! This routine sets the complex ovlp matrix.
!!
subroutine elsi_set_complex_ovlp(elsi_h,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   complex(kind=r8),  target        :: s_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap matrix

   character*40, parameter :: caller = "elsi_set_complex_ovlp"

   if(.not. elsi_h%ovlp_is_unit) then
      if(elsi_h%solver == ELPA) then
         if((elsi_h%matrix_format == BLACS_DENSE) .and. &
            (elsi_h%n_elsi_calls == 1)) then
            call elsi_set_full_mat(elsi_h,s_in)
         endif

         elsi_h%ovlp_complex => s_in
      elseif(elsi_h%solver == LIBOMM) then
         if((elsi_h%matrix_format == BLACS_DENSE) .and. &
            (elsi_h%n_elsi_calls == 1)) then
            call elsi_set_full_mat(elsi_h,s_in)
         endif

         call m_register_pdbc(elsi_h%ovlp_omm,s_in,elsi_h%sc_desc)
      endif
   endif
end subroutine

!>
!! This routine sets the sparse real overlap matrix.
!!
subroutine elsi_set_sparse_real_ovlp(elsi_h,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                !< Handle
   real(kind=r8),     target        :: s_in(elsi_h%nnz_l_sp) !< Overlap matrix

   character*40, parameter :: caller = "elsi_set_sparse_real_ovlp"

   if(.not. elsi_h%ovlp_is_unit) then
      if((elsi_h%solver == PEXSI) .or. (elsi_h%solver == SIPS)) then
         elsi_h%ovlp_real_ccs => s_in
      endif
   endif

end subroutine

!>
!! This routine sets the sparse complex overlap matrix.
!!
subroutine elsi_set_sparse_complex_ovlp(elsi_h,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                !< Handle
   complex(kind=r8),  target        :: s_in(elsi_h%nnz_l_sp) !< Overlap matrix

   character*40, parameter :: caller = "elsi_set_sparse_complex_ovlp"

   if(.not. elsi_h%ovlp_is_unit) then
      if((elsi_h%solver == PEXSI) .or. (elsi_h%solver == SIPS)) then
         elsi_h%ovlp_complex_ccs => s_in
      endif
   endif

end subroutine

!>
!! This routine sets the eigenvalues.
!!
subroutine elsi_set_eval(elsi_h,e_val_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                   !< Handle
   real(kind=r8),     target        :: e_val_in(elsi_h%n_basis) !< Eigenvalues

   character*40, parameter :: caller = "elsi_set_eval"

   if((elsi_h%solver == ELPA) .or. (elsi_h%solver == SIPS)) then
      elsi_h%eval => e_val_in
   endif

end subroutine

!>
!! This routine sets the real eigenvectors.
!!
subroutine elsi_set_real_evec(elsi_h,e_vec_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                    !< Handle
   real(kind=r8),     target        :: e_vec_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_set_real_evec"

   if((elsi_h%solver == ELPA) .or. (elsi_h%solver == SIPS)) then
      elsi_h%evec_real => e_vec_in
   endif

end subroutine

!>
!! This routine sets the complex eigenvectors.
!!
subroutine elsi_set_complex_evec(elsi_h,e_vec_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                    !< Handle
   complex(kind=r8),  target        :: e_vec_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_set_complex_evec"

   if((elsi_h%solver == ELPA) .or. (elsi_h%solver == SIPS)) then
      elsi_h%evec_complex => e_vec_in
   endif

end subroutine

!>
!! This routine sets the real density matrix.
!!
subroutine elsi_set_real_dm(elsi_h,d_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     target        :: d_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix

   character*40, parameter :: caller = "elsi_set_real_dm"

   if(elsi_h%solver == ELPA) then
      elsi_h%dm_real => d_in
   elseif(elsi_h%solver == LIBOMM) then
      call m_register_pdbc(elsi_h%dm_omm,d_in,elsi_h%sc_desc)
   endif

end subroutine

!>
!! This routine sets the complex density matrix.
!!
subroutine elsi_set_complex_dm(elsi_h,d_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   complex(kind=r8),  target        :: d_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix

   character*40, parameter :: caller = "elsi_set_complex_dm"

   if(elsi_h%solver == ELPA) then
      elsi_h%dm_complex => d_in
   elseif(elsi_h%solver == LIBOMM) then
      call m_register_pdbc(elsi_h%dm_omm,d_in,elsi_h%sc_desc)
   endif

end subroutine

!>
!! This routine sets the sparse real density matrix.
!!
subroutine elsi_set_sparse_real_dm(elsi_h,d_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                !< Handle
   real(kind=r8),     target        :: d_in(elsi_h%nnz_l_sp) !< Density matrix

   character*40, parameter :: caller = "elsi_set_sparse_real_dm"

   if(elsi_h%solver == PEXSI) then
      elsi_h%dm_real_ccs => d_in
   endif

end subroutine

!>
!! This routine sets the sparse complex density matrix.
!!
subroutine elsi_set_sparse_complex_dm(elsi_h,d_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                !< Handle
   complex(kind=r8),  target        :: d_in(elsi_h%nnz_l_sp) !< Density matrix

   character*40, parameter :: caller = "elsi_set_sparse_complex_dm"

   if(elsi_h%solver == PEXSI) then
      elsi_h%dm_complex_ccs => d_in
   endif

end subroutine

!>
!! This routine sets the row index.
!!
subroutine elsi_set_row_ind(elsi_h,row_ind_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                      !< Handle
   integer(kind=i4),  target        :: row_ind_in(elsi_h%nnz_l_sp) !< Row index

   character*40, parameter :: caller = "elsi_set_row_ind"

   elsi_h%row_ind_ccs => row_ind_in

end subroutine

!>
!! This routine sets the column pointer.
!!
subroutine elsi_set_col_ptr(elsi_h,col_ptr_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                           !< Handle
   integer(kind=i4),  target        :: col_ptr_in(elsi_h%n_l_cols_sp+1) !< Column pointer

   character*40, parameter :: caller = "elsi_set_col_ptr"

   elsi_h%col_ptr_ccs => col_ptr_in

end subroutine

!>
!! This routine frees memory.
!!
subroutine elsi_cleanup(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: ierr

   character*40, parameter :: caller = "elsi_cleanup"

   ! Nullify pointers
   if(associated(elsi_h%ham_real))          nullify(elsi_h%ham_real)
   if(associated(elsi_h%ham_complex))       nullify(elsi_h%ham_complex)
   if(associated(elsi_h%ovlp_real))         nullify(elsi_h%ovlp_real)
   if(associated(elsi_h%ovlp_complex))      nullify(elsi_h%ovlp_complex)
   if(associated(elsi_h%eval))              nullify(elsi_h%eval)
   if(associated(elsi_h%evec_real))         nullify(elsi_h%evec_real)
   if(associated(elsi_h%evec_complex))      nullify(elsi_h%evec_complex)
   if(associated(elsi_h%dm_real))           nullify(elsi_h%dm_real)
   if(associated(elsi_h%dm_complex))        nullify(elsi_h%dm_complex)
   if(associated(elsi_h%ham_real_ccs))      nullify(elsi_h%ham_real_ccs)
   if(associated(elsi_h%ham_complex_ccs))   nullify(elsi_h%ham_complex_ccs)
   if(associated(elsi_h%ovlp_real_ccs))     nullify(elsi_h%ovlp_real_ccs)
   if(associated(elsi_h%ovlp_complex_ccs))  nullify(elsi_h%ovlp_complex_ccs)
   if(associated(elsi_h%dm_real_ccs))       nullify(elsi_h%dm_real_ccs)
   if(associated(elsi_h%dm_complex_ccs))    nullify(elsi_h%dm_complex_ccs)
   if(associated(elsi_h%row_ind_ccs))       nullify(elsi_h%row_ind_ccs)
   if(associated(elsi_h%col_ptr_ccs))       nullify(elsi_h%col_ptr_ccs)

   ! ELPA
   if(allocated(elsi_h%ham_real_elpa))      call elsi_deallocate(elsi_h,elsi_h%ham_real_elpa,"ham_real_elpa")
   if(allocated(elsi_h%ham_complex_elpa))   call elsi_deallocate(elsi_h,elsi_h%ham_complex_elpa,"ham_complex_elpa")
   if(allocated(elsi_h%ovlp_real_elpa))     call elsi_deallocate(elsi_h,elsi_h%ovlp_real_elpa,"ovlp_real_elpa")
   if(allocated(elsi_h%ovlp_complex_elpa))  call elsi_deallocate(elsi_h,elsi_h%ovlp_complex_elpa,"ovlp_complex_elpa")
   if(allocated(elsi_h%eval_elpa))          call elsi_deallocate(elsi_h,elsi_h%eval_elpa,"eval_elpa")
   if(allocated(elsi_h%evec_real_elpa))     call elsi_deallocate(elsi_h,elsi_h%evec_real_elpa,"evec_real_elpa")
   if(allocated(elsi_h%evec_complex_elpa))  call elsi_deallocate(elsi_h,elsi_h%evec_complex_elpa,"evec_complex_elpa")
   if(allocated(elsi_h%dm_real_elpa))       call elsi_deallocate(elsi_h,elsi_h%dm_real_elpa,"dm_real_elpa")
   if(allocated(elsi_h%dm_complex_elpa))    call elsi_deallocate(elsi_h,elsi_h%dm_complex_elpa,"dm_complex_elpa")
   if(allocated(elsi_h%occ_num))            call elsi_deallocate(elsi_h,elsi_h%occ_num,"occ_num")
   if(allocated(elsi_h%k_weight))           call elsi_deallocate(elsi_h,elsi_h%k_weight,"k_weight")
   if(allocated(elsi_h%eval_all))           call elsi_deallocate(elsi_h,elsi_h%eval_all,"eval_all")

   ! libOMM
   if(elsi_h%ham_omm%is_initialized)        call m_deallocate(elsi_h%ham_omm)
   if(elsi_h%ovlp_omm%is_initialized)       call m_deallocate(elsi_h%ovlp_omm)
   if(elsi_h%dm_omm%is_initialized)         call m_deallocate(elsi_h%dm_omm)
   if(elsi_h%coeff_omm%is_initialized)      call m_deallocate(elsi_h%coeff_omm)
   if(elsi_h%tdm_omm%is_initialized)        call m_deallocate(elsi_h%tdm_omm)
   if(allocated(elsi_h%ovlp_real_omm))      call elsi_deallocate(elsi_h,elsi_h%ovlp_real_omm,"ovlp_real_omm")
   if(allocated(elsi_h%ovlp_complex_omm))   call elsi_deallocate(elsi_h,elsi_h%ovlp_complex_omm,"ovlp_complex_omm")

   ! PEXSI
   if(allocated(elsi_h%ham_real_pexsi))     call elsi_deallocate(elsi_h,elsi_h%ham_real_pexsi,"ham_real_pexsi")
   if(allocated(elsi_h%ham_complex_pexsi))  call elsi_deallocate(elsi_h,elsi_h%ham_complex_pexsi,"ham_complex_pexsi")
   if(allocated(elsi_h%ovlp_real_pexsi))    call elsi_deallocate(elsi_h,elsi_h%ovlp_real_pexsi,"ovlp_real_pexsi")
   if(allocated(elsi_h%ovlp_complex_pexsi)) call elsi_deallocate(elsi_h,elsi_h%ovlp_complex_pexsi,"ovlp_complex_pexsi")
   if(allocated(elsi_h%dm_real_pexsi))      call elsi_deallocate(elsi_h,elsi_h%dm_real_pexsi,"dm_real_pexsi")
   if(allocated(elsi_h%dm_complex_pexsi))   call elsi_deallocate(elsi_h,elsi_h%dm_complex_pexsi,"dm_complex_pexsi")
   if(allocated(elsi_h%row_ind_pexsi))      call elsi_deallocate(elsi_h,elsi_h%row_ind_pexsi,"row_ind_pexsi")
   if(allocated(elsi_h%col_ptr_pexsi))      call elsi_deallocate(elsi_h,elsi_h%col_ptr_pexsi,"col_ptr_pexsi")
   if(allocated(elsi_h%ne_vec))             call elsi_deallocate(elsi_h,elsi_h%ne_vec,"ne_vec")

   ! CheSS
   if(allocated(elsi_h%ham_real_chess))     call elsi_deallocate(elsi_h,elsi_h%ham_real_chess,"ham_real_chess")
   if(allocated(elsi_h%ovlp_real_chess))    call elsi_deallocate(elsi_h,elsi_h%ovlp_real_chess,"ovlp_real_chess")
   if(allocated(elsi_h%row_ind_chess))      call elsi_deallocate(elsi_h,elsi_h%row_ind_chess,"row_ind_chess")
   if(allocated(elsi_h%col_ptr_chess))      call elsi_deallocate(elsi_h,elsi_h%col_ptr_chess,"col_ptr_chess")
   if(allocated(elsi_h%row_ind_buffer))     call elsi_deallocate(elsi_h,elsi_h%row_ind_buffer,"row_ind_buffer")
   if(allocated(elsi_h%col_ptr_buffer))     call elsi_deallocate(elsi_h,elsi_h%col_ptr_buffer,"col_ptr_buffer")

   ! SIPs
   if(allocated(elsi_h%ham_real_sips))      call elsi_deallocate(elsi_h,elsi_h%ham_real_sips,"ham_real_sips")
   if(allocated(elsi_h%ovlp_real_sips))     call elsi_deallocate(elsi_h,elsi_h%ovlp_real_sips,"ovlp_real_sips")
   if(allocated(elsi_h%dm_real_sips))       call elsi_deallocate(elsi_h,elsi_h%dm_real_sips,"dm_real_sips")
   if(allocated(elsi_h%row_ind_sips))       call elsi_deallocate(elsi_h,elsi_h%row_ind_sips,"row_ind_sips")
   if(allocated(elsi_h%col_ptr_sips))       call elsi_deallocate(elsi_h,elsi_h%col_ptr_sips,"col_ptr_sips")
   if(allocated(elsi_h%slices))             call elsi_deallocate(elsi_h,elsi_h%slices,"slices")

   if(allocated(elsi_h%local_row))          call elsi_deallocate(elsi_h,elsi_h%local_row,"local_row")
   if(allocated(elsi_h%local_col))          call elsi_deallocate(elsi_h,elsi_h%local_col,"local_col")

   ! Finalize PEXSI
   if(elsi_h%pexsi_started) then
      call f_ppexsi_plan_finalize(elsi_h%pexsi_plan,ierr)
      call MPI_Comm_free(elsi_h%comm_among_pole,ierr)
      call MPI_Comm_free(elsi_h%comm_in_pole,ierr)
      call MPI_Comm_free(elsi_h%comm_among_point,ierr)
      call MPI_Comm_free(elsi_h%comm_in_point,ierr)
   endif

   ! Finalize CheSS
   if(elsi_h%chess_started) then
      call deallocate_sparse_matrix(elsi_h%sparse_mat(1))
      call deallocate_sparse_matrix(elsi_h%sparse_mat(2))
      call foe_data_deallocate(elsi_h%ice_obj)
      call foe_data_deallocate(elsi_h%foe_obj)
      call deallocate_matrices(elsi_h%ham_chess)
      call deallocate_matrices(elsi_h%ovlp_chess)
      call deallocate_matrices(elsi_h%dm_chess)
      call deallocate_matrices(elsi_h%edm_chess)
      call deallocate_matrices(elsi_h%ovlp_inv_sqrt(1))
      call f_lib_finalize()
   endif

   ! Finalize SIPs
   if(elsi_h%sips_started) then
      call clean_qetsc()
   endif

   ! Reset elsi_h
   call elsi_reset_handle(elsi_h)

end subroutine

!>
!! This routine resets an ELSI handle.
!!
subroutine elsi_reset_handle(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_reset_handle"

   elsi_h%handle_ready      = .false.
   elsi_h%solver            = UNSET
   elsi_h%matrix_data_type  = UNSET
   elsi_h%matrix_format     = UNSET
   elsi_h%uplo              = FULL_MAT
   elsi_h%parallel_mode     = UNSET
   elsi_h%n_elsi_calls      = 0
   elsi_h%n_b_rows          = UNSET
   elsi_h%n_b_cols          = UNSET
   elsi_h%n_p_rows          = UNSET
   elsi_h%n_p_cols          = UNSET
   elsi_h%n_l_rows          = UNSET
   elsi_h%n_l_cols          = UNSET
   elsi_h%myid              = UNSET
   elsi_h%myid_all          = UNSET
   elsi_h%n_procs           = UNSET
   elsi_h%n_procs_all       = UNSET
   elsi_h%mpi_comm          = UNSET
   elsi_h%mpi_comm_all      = UNSET
   elsi_h%mpi_comm_row      = UNSET
   elsi_h%mpi_comm_col      = UNSET
   elsi_h%my_p_row          = UNSET
   elsi_h%my_p_col          = UNSET
   elsi_h%mpi_ready         = .false.
   elsi_h%global_mpi_ready  = .false.
   elsi_h%blacs_ctxt        = UNSET
   elsi_h%sc_desc           = UNSET
   elsi_h%blacs_ready       = .false.
   elsi_h%nnz_g             = UNSET
   elsi_h%nnz_l             = UNSET
   elsi_h%nnz_l_sp          = UNSET
   elsi_h%n_l_cols_sp       = UNSET
   elsi_h%zero_threshold    = 1.0e-15_r8
   elsi_h%sparsity_ready    = .false.
   elsi_h%ovlp_is_unit      = .false.
   elsi_h%ovlp_is_sing      = .false.
   elsi_h%no_sing_check     = .false.
   elsi_h%sing_tol          = 1.0e-5_r8
   elsi_h%stop_sing         = .false.
   elsi_h%n_nonsing         = UNSET
   elsi_h%n_electrons       = 0.0_r8
   elsi_h%mu                = 0.0_r8
   elsi_h%n_basis           = UNSET
   elsi_h%n_spins           = 1
   elsi_h%n_kpts            = 1
   elsi_h%n_states          = UNSET
   elsi_h%n_states_solve    = UNSET
   elsi_h%i_spin            = 1
   elsi_h%i_kpt             = 1
   elsi_h%i_weight          = 1.0_r8
   elsi_h%energy_hdm        = 0.0_r8
   elsi_h%energy_sedm       = 0.0_r8
   elsi_h%free_energy       = 0.0_r8
   elsi_h%broaden_scheme    = 0
   elsi_h%broaden_width     = 1.0e-2_r8
   elsi_h%occ_tolerance     = 1.0e-13_r8
   elsi_h%max_mu_steps      = 100
   elsi_h%spin_degen        = 0.0_r8
   elsi_h%spin_is_set       = .false.
   elsi_h%mu_ready          = .false.
   elsi_h%edm_ready_real    = .false.
   elsi_h%edm_ready_complex = .false.
   elsi_h%elpa_solver       = UNSET
   elsi_h%elpa_output       = .false.
   elsi_h%n_states_omm      = UNSET
   elsi_h%n_elpa_steps      = UNSET
   elsi_h%new_overlap       = .true.
   elsi_h%coeff_ready       = .false.
   elsi_h%omm_flavor        = UNSET
   elsi_h%scale_kinetic     = 0.0_r8
   elsi_h%calc_ed           = .false.
   elsi_h%eta               = 0.0_r8
   elsi_h%min_tol           = 1.0e-12_r8
   elsi_h%omm_output        = .false.
   elsi_h%do_dealloc        = .false.
   elsi_h%use_psp           = .false.
   elsi_h%n_p_per_pole      = UNSET
   elsi_h%n_p_per_point     = UNSET
   elsi_h%my_p_row_pexsi    = UNSET
   elsi_h%my_p_col_pexsi    = UNSET
   elsi_h%n_p_rows_pexsi    = UNSET
   elsi_h%n_p_cols_pexsi    = UNSET
   elsi_h%my_point          = UNSET
   elsi_h%myid_point        = UNSET
   elsi_h%comm_among_pole   = UNSET
   elsi_h%comm_in_pole      = UNSET
   elsi_h%comm_among_point  = UNSET
   elsi_h%comm_in_point     = UNSET
   elsi_h%ne_pexsi          = 0.0_r8
   elsi_h%pexsi_started     = .false.
   elsi_h%erf_decay         = 0.0_r8
   elsi_h%erf_decay_min     = 0.0_r8
   elsi_h%erf_decay_max     = 0.0_r8
   elsi_h%ev_ham_min        = 0.0_r8
   elsi_h%ev_ham_max        = 0.0_r8
   elsi_h%ev_ovlp_min       = 0.0_r8
   elsi_h%ev_ovlp_max       = 0.0_r8
   elsi_h%beta              = 0.0_r8
   elsi_h%chess_started     = .false.
   elsi_h%n_p_per_slice     = UNSET
   elsi_h%n_inertia_steps   = UNSET
   elsi_h%slicing_method    = UNSET
   elsi_h%inertia_option    = UNSET
   elsi_h%unbound           = UNSET
   elsi_h%n_slices          = UNSET
   elsi_h%interval          = 0.0_r8
   elsi_h%slice_buffer      = 0.0_r8
   elsi_h%sips_started      = .false.
   elsi_h%clock_rate        = UNSET

end subroutine

!>
!! This routine guarantees that there are no mutually conflicting parameters.
!!
subroutine elsi_check(elsi_h,caller)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   character(len=*),  intent(in)    :: caller !< Caller

   ! General check of solver, parallel mode, matrix data type, matrix storage format
   if((elsi_h%solver < 0) .or. (elsi_h%solver >= N_SOLVERS)) then
      call elsi_stop(" An unsupported solver has been chosen."//&
              " Please choose ELPA, LIBOMM, or PEXSI solver."//&
              " Exiting...",elsi_h,caller)
   endif

   if((elsi_h%parallel_mode < 0) .or. &
      (elsi_h%parallel_mode >= N_PARALLEL_MODES)) then
      call elsi_stop(" An unsupported parallel mode has been chosen."//&
              " Please choose either SINGLE_PROC or MULTI_PROC."//&
              " Exiting...",elsi_h,caller)
   endif

   if((elsi_h%matrix_data_type < 0) .or. &
      (elsi_h%matrix_data_type >= N_MATRIX_DATA_TYPES)) then
      call elsi_stop(" An unsupported matirx data type has been chosen."//&
              " Please choose either REAL_VALUES or COMPLEX_VALUES."//&
              " Exiting...",elsi_h,caller)
   endif

   if((elsi_h%matrix_format < 0) .or. &
      (elsi_h%matrix_format >= N_MATRIX_STORAGE_FORMATS)) then
      call elsi_stop(" An unsupported matirx storage format has been"//&
              " set. Please choose either BLACS_DENSE or PEXSI_CSC."//&
              " Exiting...",elsi_h,caller)
   endif

   if(elsi_h%uplo /= FULL_MAT) then
      if((elsi_h%matrix_format /= BLACS_DENSE) .or. &
         (elsi_h%parallel_mode /= MULTI_PROC)) then
         call elsi_stop(" Upper/lower triangular input matrix only"//&
                 " supported with BLACS_DENSE matrix storage format and"//&
                 " MULTI_PROC parallel mode. Exiting...",elsi_h,caller)
      endif

      if((elsi_h%uplo /= UT_MAT) .and. (elsi_h%uplo /= LT_MAT)) then
         call elsi_stop(" Invalid choice of uplo. Exiting...",elsi_h,caller)
      endif
   endif

   ! Spin and k-point
   if(elsi_h%n_spins*elsi_h%n_kpts > 1) then
      if(.not. elsi_h%global_mpi_ready) then
         call elsi_stop(" Spin/k-point calculations require a global"//&
                 " MPI communicator. Exiting...",elsi_h,caller)
      endif
   endif

   if(.not. elsi_h%global_mpi_ready) then
      elsi_h%mpi_comm_all = elsi_h%mpi_comm
      elsi_h%n_procs_all  = elsi_h%n_procs
      elsi_h%myid_all     = elsi_h%myid
   endif

   ! Specific check for each solver
   if(elsi_h%solver == AUTO) then
      call elsi_stop(" AUTO not yet available."//&
               " Please choose ELPA, LIBOMM, or PEXSI solver."//&
               " Exiting...",elsi_h,caller)
   elseif(elsi_h%solver == ELPA) then
      if(elsi_h%parallel_mode == MULTI_PROC) then
         if(.not. elsi_h%mpi_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI"//&
                    " being set up before calling the solver."//&
                    " Exiting...",elsi_h,caller)
         endif

         if(.not. elsi_h%blacs_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires BLACS"//&
                    " being set up before calling the solver."//&
                    " Exiting...",elsi_h,caller)
         endif
      endif

      if(elsi_h%matrix_format == PEXSI_CSC) then
         if(.not. elsi_h%sparsity_ready) then
            call elsi_stop(" The PEXSI_CSC format has been chosen."//&
                    " Please set the sparsity pattern before"//&
                    " calling the solver. Exiting...",elsi_h,caller)
         endif
      endif
   elseif(elsi_h%solver == LIBOMM) then
      if(elsi_h%parallel_mode == MULTI_PROC) then
         if(.not. elsi_h%mpi_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI"//&
                    " being set up before calling the solver."//&
                    " Exiting...",elsi_h,caller)
         endif

         if(.not. elsi_h%blacs_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires BLACS"//&
                    " being set up before calling the solver."//&
                    " Exiting...",elsi_h,caller)
         endif
      else
         call elsi_stop(" libOMM has been chosen as the solver."//&
                 " Please choose MULTI_PROC parallel mode."//&
                 " Exiting...",elsi_h,caller)
      endif

      if(elsi_h%matrix_format == PEXSI_CSC) then
         if(.not. elsi_h%sparsity_ready) then
            call elsi_stop(" The PEXSI_CSC format has been chosen."//&
                    " Please set the sparsity pattern before"//&
                    " calling the solver. Exiting...",elsi_h,caller)
         endif
      endif
   elseif(elsi_h%solver == PEXSI) then
      if(elsi_h%parallel_mode == MULTI_PROC) then
         if(.not. elsi_h%mpi_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI"//&
                    " being set up before calling the solver."//&
                    " Exiting...",elsi_h,caller)
         endif
      else
         call elsi_stop(" PEXSI has been chosen as the solver."//&
                 " Please choose MULTI_PROC parallel mode."//&
                 " Exiting...",elsi_h,caller)
      endif

      if(elsi_h%matrix_format == BLACS_DENSE) then
         if(.not. elsi_h%blacs_ready) then
            call elsi_stop(" The BLACS_DENSE format has been chosen."//&
                    " Please set up BLACS before calling the"//&
                    " solver. Exiting...",elsi_h,caller)
         endif

         if(elsi_h%ovlp_is_unit) then
            call elsi_stop(" PEXSI with BLACS_DENSE format and an"//&
                    " identity overlap matrix not yet available."//&
                    " Exiting...",elsi_h,caller)
         endif
      else
         if(.not. elsi_h%sparsity_ready) then
            call elsi_stop(" The PEXSI_CSC format has been chosen."//&
                    " Please set the sparsity pattern before"//&
                    " calling the solver. Exiting...",elsi_h,caller)
         endif
      endif

      if(elsi_h%n_basis < elsi_h%n_p_per_pole) then
         call elsi_stop(" PEXSI has been chosen as the solver."//&
                 " The size of matrix is too small for this"//&
                 " number of MPI tasks. Exiting...",elsi_h,caller)
      endif

      if(elsi_h%n_p_per_pole == UNSET) then
         if(mod(elsi_h%n_procs,elsi_h%pexsi_options%numPole*&
            elsi_h%pexsi_options%nPoints) /= 0) then
            call elsi_stop("  The total number of MPI tasks must be"//&
                    " a multiple of the number of MPI tasks per pole"//&
                    " times the number of mu points. Exiting...",&
                    elsi_h,caller)
         endif
      else
         if(mod(elsi_h%n_procs,elsi_h%n_p_per_pole*&
            elsi_h%pexsi_options%nPoints) /= 0) then
            call elsi_stop("  The total number of MPI tasks must be"//&
                    " a multiple of the number of MPI tasks per pole"//&
                    " times the number of mu points. Exiting...",&
                    elsi_h,caller)
         endif

         if(elsi_h%n_p_per_pole*elsi_h%pexsi_options%numPole*&
            elsi_h%pexsi_options%nPoints < elsi_h%n_procs) then
            call elsi_stop("  Specified number of MPI tasks per pole"//&
                    " is too small for the total number of MPI tasks."//&
                    " Exiting...",elsi_h,caller)
         endif
      endif
   elseif(elsi_h%solver == CHESS) then
      call elsi_statement_print("  ATTENTION! CheSS is EXPERIMENTAL.",elsi_h)

      if(elsi_h%n_basis < elsi_h%n_procs) then
         call elsi_stop(" CheSS has been chosen as the solver."//&
                 " The size of matrix is too small for this"//&
                 " number of MPI tasks. Exiting...",elsi_h,caller)
      endif

      if(elsi_h%parallel_mode == MULTI_PROC) then
         if(.not. elsi_h%mpi_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI"//&
                    " being set up before calling the solver."//&
                    " Exiting...",elsi_h,caller)
         endif
      else
         call elsi_stop(" CheSS has been chosen as the solver."//&
                 " Please choose MULTI_PROC parallel mode."//&
                 " Exiting...",elsi_h,caller)
      endif

      if(elsi_h%matrix_format == BLACS_DENSE) then
         if(.not. elsi_h%blacs_ready) then
            call elsi_stop(" The BLACS_DENSE format has been chosen."//&
                    " Please set up BLACS before calling the solver."//&
                    " Exiting...",elsi_h,caller)
         endif

         if(elsi_h%ovlp_is_unit) then
            call elsi_stop(" CheSS with BLACS_DENSE format and an"//&
                    " identity overlap matrix not yet available."//&
                    " Exiting...",elsi_h,caller)
         endif
      else
         call elsi_stop(" CheSS has been chosen as the solver."//&
                 " Please choose BLACS_DENSE matrix format."//&
                 " Exiting...",elsi_h,caller)
      endif
   elseif(elsi_h%solver == SIPS) then
      call elsi_statement_print("  ATTENTION! SIPs is EXPERIMENTAL.",elsi_h)

      if(elsi_h%n_basis < elsi_h%n_procs) then
         call elsi_stop(" SIPs has been chosen as the solver."//&
                 " The size of matrix is too small for this"//&
                 " number of MPI tasks. Exiting...",elsi_h,caller)
      endif

      if(elsi_h%parallel_mode == MULTI_PROC) then
         if(.not. elsi_h%mpi_ready) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI"//&
                    " being set up before calling the solver."//&
                    " Exiting...",elsi_h,caller)
         endif
      else
         call elsi_stop(" SIPs has been chosen as the solver."//&
                 " Please choose MULTI_PROC parallel mode."//&
                 " Exiting...",elsi_h,caller)
      endif

      if(elsi_h%matrix_format == BLACS_DENSE) then
         if(.not. elsi_h%blacs_ready) then
            call elsi_stop(" The BLACS_DENSE format has been chosen."//&
                    " Please set up BLACS before calling the solver."//&
                    " Exiting...",elsi_h,caller)
         endif

         if(elsi_h%ovlp_is_unit) then
            call elsi_stop(" SIPs with BLACS_DENSE format and an"//&
                    " identity overlap matrix not yet available."//&
                    " Exiting...",elsi_h,caller)
         endif
      else
         if(.not. elsi_h%sparsity_ready) then
            call elsi_stop(" The PEXSI_CSC format has been chosen."//&
                    " Please set the sparsity pattern before"//&
                    " calling the solver. Exiting...",elsi_h,caller)
         endif
      endif
   else
      call elsi_stop(" No supported solver has been chosen."//&
              " Please choose ELPA, LIBOMM, or PEXSI solver."//&
              " Exiting...",elsi_h,caller)
   endif

end subroutine

!>
!! This routine checks whether the input elsi_h has been properly
!! initialized by the user or not.
!!
!! Do NOT call this subroutine outside the (main) ELSI module.
!!
subroutine elsi_check_handle(elsi_h,caller)

   implicit none

   type(elsi_handle), intent(in) :: elsi_h !< Handle
   character(len=*),  intent(in) :: caller !< Caller

   if(.not. elsi_h%handle_ready) then
      call elsi_stop(" Invalid handle! An ELSI handle must be properly"//&
              " initialized with 'elsi_init' before being used."//&
              " Exiting...",elsi_h,caller)
   endif

end subroutine

!>
!! This routine computes the global row index based on the local row index.
!!
subroutine elsi_get_global_row(elsi_h,global_idx,local_idx)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h     !< Handle
   integer(kind=i4),  intent(in)  :: local_idx  !< Local index
   integer(kind=i4),  intent(out) :: global_idx !< Global index

   integer(kind=i4) :: block
   integer(kind=i4) :: idx

   character*40, parameter :: caller = "elsi_get_global_row"

   block = (local_idx-1)/elsi_h%n_b_rows
   idx = local_idx-block*elsi_h%n_b_rows

   global_idx = elsi_h%my_p_row*elsi_h%n_b_rows+&
                   block*elsi_h%n_b_rows*elsi_h%n_p_rows+idx

end subroutine

!>
!! This routine computes the global column index based on the local column index.
!!
subroutine elsi_get_global_col(elsi_h,global_idx,local_idx)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h     !< Handle
   integer(kind=i4),  intent(in)  :: local_idx  !< Local index
   integer(kind=i4),  intent(out) :: global_idx !< Global index

   integer(kind=i4) :: block
   integer(kind=i4) :: idx

   character*40, parameter :: caller = "elsi_get_global_col"

   block = (local_idx-1)/elsi_h%n_b_cols
   idx = local_idx-block*elsi_h%n_b_cols

   global_idx = elsi_h%my_p_col*elsi_h%n_b_cols+&
                   block*elsi_h%n_b_cols*elsi_h%n_p_cols+idx

end subroutine

!>
!! This routine counts the local number of non_zero elements.
!!
subroutine elsi_get_local_nnz_real(elsi_h,matrix,n_rows,n_cols,nnz)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h                !< Handle
   real(kind=r8),     intent(in)  :: matrix(n_rows,n_cols) !< Local matrix
   integer(kind=i4),  intent(in)  :: n_rows                !< Local rows
   integer(kind=i4),  intent(in)  :: n_cols                !< Local cols
   integer(kind=i4),  intent(out) :: nnz                   !< Number of non-zero

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_nz

   character*40, parameter :: caller = "elsi_get_local_nnz_real"

   nnz = 0

   do i_col = 1,n_cols
      do i_row = 1,n_rows
         if(abs(matrix(i_row,i_col)) > elsi_h%zero_threshold) then
            nnz = nnz+1
         endif
      enddo
   enddo

end subroutine

!>
!! This routine counts the local number of non_zero elements.
!!
subroutine elsi_get_local_nnz_complex(elsi_h,matrix,n_rows,n_cols,nnz)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h                !< Handle
   complex(kind=r8),  intent(in)  :: matrix(n_rows,n_cols) !< Local matrix
   integer(kind=i4),  intent(in)  :: n_rows                !< Local rows
   integer(kind=i4),  intent(in)  :: n_cols                !< Local cols
   integer(kind=i4),  intent(out) :: nnz                   !< Number of non-zero

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_nz

   character*40, parameter :: caller = "elsi_get_local_nnz_complex"

   nnz = 0

   do i_col = 1,n_cols
      do i_row = 1,n_rows
         if(abs(matrix(i_row,i_col)) > elsi_h%zero_threshold) then
            nnz = nnz+1
         endif
      enddo
   enddo

end subroutine

!>
!! This routine sets a full matrix from a (upper or lower) triangular matrix.
!! The size of matrix should be the same as the Hamiltonian matrix.
!!
subroutine elsi_set_full_mat_real(elsi_h,mat)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                               !< Handle
   real(kind=r8),     intent(inout) :: mat(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Matrix

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   real(kind=r8), allocatable :: tmp_real(:,:)

   character*40, parameter :: caller = "elsi_set_full_mat_real"

   if((elsi_h%uplo /= FULL_MAT) .and. (elsi_h%parallel_mode == MULTI_PROC)) then
      call elsi_allocate(elsi_h,tmp_real,elsi_h%n_l_rows,elsi_h%n_l_cols,"tmp_real",caller)

      call pdtran(elsi_h%n_basis,elsi_h%n_basis,1.0_r8,mat,1,1,elsi_h%sc_desc,&
              0.0_r8,tmp_real,1,1,elsi_h%sc_desc)

      if(elsi_h%uplo == UT_MAT) then ! Upper triangular
         do i_col = 1,elsi_h%n_basis-1
            if(elsi_h%local_col(i_col) == 0) cycle

            do i_row = i_col+1,elsi_h%n_basis
               if(elsi_h%local_row(i_row) > 0) then
                  mat(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                     tmp_real(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
               endif
            enddo
         enddo
      elseif(elsi_h%uplo == LT_MAT) then ! Lower triangular
         do i_col = 2,elsi_h%n_basis
            if(elsi_h%local_col(i_col) == 0) cycle

            do i_row = 1,i_col-1
               if(elsi_h%local_row(i_row) > 0) then
                  mat(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                     tmp_real(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
               endif
            enddo
         enddo
      endif

      call elsi_deallocate(elsi_h,tmp_real,"tmp_real")
   endif

end subroutine

!>
!! This routine sets a full matrix from a (upper or lower) triangular matrix.
!! The size of matrix should be the same as the Hamiltonian matrix.
!!
subroutine elsi_set_full_mat_complex(elsi_h,mat)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                               !< Handle
   complex(kind=r8),  intent(inout) :: mat(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Matrix

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   complex(kind=r8), allocatable :: tmp_complex(:,:)

   character*40, parameter :: caller = "elsi_set_full_mat_complex"

   if((elsi_h%uplo /= FULL_MAT) .and. (elsi_h%parallel_mode == MULTI_PROC)) then
      call elsi_allocate(elsi_h,tmp_complex,elsi_h%n_l_rows,elsi_h%n_l_cols,&
              "tmp_complex",caller)

      call pztranc(elsi_h%n_basis,elsi_h%n_basis,(1.0_r8,0.0_r8),mat,1,1,&
              elsi_h%sc_desc,(0.0_r8,0.0_r8),tmp_complex,1,1,elsi_h%sc_desc)

      if(elsi_h%uplo == UT_MAT) then ! Upper triangular
         do i_col = 1,elsi_h%n_basis-1
            if(elsi_h%local_col(i_col) == 0) cycle

            do i_row = i_col+1,elsi_h%n_basis
               if(elsi_h%local_row(i_row) > 0) then
                  mat(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                     tmp_complex(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
               endif
            enddo
         enddo
      elseif(elsi_h%uplo == LT_MAT) then ! Lower triangular
         do i_col = 2,elsi_h%n_basis
            if(elsi_h%local_col(i_col) == 0) cycle

            do i_row = 1,i_col-1
               if(elsi_h%local_row(i_row) > 0) then
                  mat(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                     tmp_complex(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
               endif
            enddo
         enddo
      endif

      call elsi_deallocate(elsi_h,tmp_complex,"tmp_complex")

      ! Make diagonal real
      do i_col = 1,elsi_h%n_basis
         if((elsi_h%local_col(i_col) == 0) .or. (elsi_h%local_row(i_col) == 0)) cycle

         mat(elsi_h%local_row(i_col),elsi_h%local_col(i_col)) = &
            dble(mat(elsi_h%local_row(i_col),elsi_h%local_col(i_col)))
      enddo
   endif

end subroutine

!>
!! This routine initializes the timer.
!!
subroutine elsi_init_timers(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: initial_time
   integer(kind=i4) :: clock_max

   character*40, parameter :: caller = "elsi_init_timers"

   call system_clock(initial_time,elsi_h%clock_rate,clock_max)

end subroutine

!>
!! This routine gets the current wallclock time.
!!
subroutine elsi_get_time(elsi_h,wtime)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h !< Handle
   real(kind=r8),     intent(out) :: wtime !< Time

   integer(kind=i4) :: tics

   character*40, parameter :: caller = "elsi_get_time"

   call system_clock(tics)

   wtime = 1.0_r8*tics/elsi_h%clock_rate

end subroutine

!>
!! This routine prints a final output.
!!
subroutine elsi_final_print(elsi_h)

   implicit none

   type(elsi_handle), intent(in) :: elsi_h !< Handle

   real(kind=r8)    :: sparsity
   integer(kind=i4) :: i_proc

   character*40, parameter :: caller = "elsi_final_print"

   if(print_info) then
      if(elsi_h%myid_all == 0) then
         write(*,"('  |------------------------------------------')")
         write(*,"('  | Final ELSI Output')")
         write(*,"('  |------------------------------------------')")

         write(*,"('  | Number of basis functions :',I13)") elsi_h%n_basis
         if((elsi_h%solver == PEXSI) .or. (elsi_h%solver == SIPS)) then
            write(*,"('  | Number of nonzeros        :',I13)") elsi_h%nnz_g
            sparsity = 1.0_r8-(1.0_r8*elsi_h%nnz_g/elsi_h%n_basis)/elsi_h%n_basis
            write(*,"('  | Sparsity                  :',F13.3)") sparsity
         endif
         write(*,"('  | Number of electrons       :',F13.1)") elsi_h%n_electrons
         write(*,"('  | Number of states          :',I13)") elsi_h%n_states
         write(*,"('  | Number of spins           :',I13)") elsi_h%n_spins
         write(*,"('  | Number of k-points        :',I13)") elsi_h%n_kpts

         if(elsi_h%solver == ELPA) then
            write(*,"('  | Solver                    :',A13)") "ELPA"
         elseif(elsi_h%solver == LIBOMM) then
            write(*,"('  | Solver                    :',A13)") "libOMM"
         elseif(elsi_h%solver == PEXSI) then
            write(*,"('  | Solver                    :',A13)") "PEXSI"
         elseif(elsi_h%solver == CHESS) then
            write(*,"('  | Solver                    :',A13)") "CheSS"
         elseif(elsi_h%solver == SIPS) then
            write(*,"('  | Solver                    :',A13)") "SIPs"
         endif

         if(elsi_h%parallel_mode == MULTI_PROC) then
            write(*,"('  | Parallel mode             :',A13)") "MULTI_PROC"
         elseif(elsi_h%parallel_mode == SINGLE_PROC) then
            write(*,"('  | Parallel mode             :',A13)") "SINGLE_PROC"
         endif

         if(elsi_h%matrix_format == BLACS_DENSE) then
            write(*,"('  | Matrix format             :',A13)") "BLACS_DENSE"
         elseif(elsi_h%matrix_format == PEXSI_CSC) then
            write(*,"('  | Matrix format             :',A13)") "PEXSI_CSC"
         endif

         write(*,"('  |------------------------------------------')")
         write(*,"('  | ELSI Project (c)  elsi-interchange.org')")
         write(*,"('  |------------------------------------------')")
      endif
   endif

end subroutine

end module ELSI_UTILS
