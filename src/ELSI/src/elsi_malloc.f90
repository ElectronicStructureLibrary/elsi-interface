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
!! This module contains customized allocate and deallocate.
!!
module ELSI_MALLOC

   use ELSI_DATATYPE,  only: elsi_handle
   use ELSI_IO,        only: elsi_say
   use ELSI_MPI,       only: elsi_stop
   use ELSI_PRECISION, only: i4,i8,r4,r8

   implicit none

   private

   public :: elsi_allocate
   public :: elsi_deallocate

   interface elsi_allocate
      module procedure elsi_allocate_integer4_1d,&
                       elsi_allocate_integer4_2d,&
                       elsi_allocate_integer4_3d,&
                       elsi_allocate_integer8_1d,&
                       elsi_allocate_real8_1d,&
                       elsi_allocate_real8_2d,&
                       elsi_allocate_real8_3d,&
                       elsi_allocate_real4_1d,&
                       elsi_allocate_real4_2d,&
                       elsi_allocate_complex16_1d,&
                       elsi_allocate_complex16_2d,&
                       elsi_allocate_complex16_3d,&
                       elsi_allocate_complex8_2d
   end interface

   interface elsi_deallocate
      module procedure elsi_deallocate_integer4_1d,&
                       elsi_deallocate_integer4_2d,&
                       elsi_deallocate_integer4_3d,&
                       elsi_deallocate_integer8_1d,&
                       elsi_deallocate_real8_1d,&
                       elsi_deallocate_real8_2d,&
                       elsi_deallocate_real8_3d,&
                       elsi_deallocate_real4_1d,&
                       elsi_deallocate_real4_2d,&
                       elsi_deallocate_complex16_1d,&
                       elsi_deallocate_complex16_2d,&
                       elsi_deallocate_complex16_3d,&
                       elsi_deallocate_complex8_2d
   end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ELSI_ALLOCATE MODULE PROCEDURES !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! This routine allocates a 1D array with real(kind=r4).
!!
subroutine elsi_allocate_real4_1d(e_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   real(kind=r4),     intent(inout), allocatable :: array(:)
   integer(kind=i4),  intent(in)                 :: dim1
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*4

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0.0_r4

end subroutine

!>
!! This routine allocates a 1D array with real(kind=r8).
!!
subroutine elsi_allocate_real8_1d(e_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   real(kind=r8),     intent(inout), allocatable :: array(:)
   integer(kind=i4),  intent(in)                 :: dim1
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*8

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 1D array with integer(kind=i4).
!!
subroutine elsi_allocate_integer4_1d(e_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   integer(kind=i4),  intent(inout), allocatable :: array(:)
   integer(kind=i4),  intent(in)                 :: dim1
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*4

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 1D array with integer(kind=i8).
!!
subroutine elsi_allocate_integer8_1d(e_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   integer(kind=i8),  intent(inout), allocatable :: array(:)
   integer(kind=i4),  intent(in)                 :: dim1
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*8

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 1D array with complex(kind=r8).
!!
subroutine elsi_allocate_complex16_1d(e_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   complex(kind=r8),  intent(inout), allocatable :: array(:)
   integer(kind=i4),  intent(in)                 :: dim1
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*16

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = (0.0_r8,0.0_r8)

end subroutine

!>
!! This routine allocates a 2D array of real(kind=r4).
!!
subroutine elsi_allocate_real4_2d(e_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   real(kind=r4),     intent(inout), allocatable :: array(:,:)
   integer(kind=i4),  intent(in)                 :: dim1
   integer(kind=i4),  intent(in)                 :: dim2
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*4

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0.0_r4

end subroutine

!>
!! This routine allocates a 2D array of real(kind=r8).
!!
subroutine elsi_allocate_real8_2d(e_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   real(kind=r8),     intent(inout), allocatable :: array(:,:)
   integer(kind=i4),  intent(in)                 :: dim1
   integer(kind=i4),  intent(in)                 :: dim2
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*8

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 2D array of integer(kind=i4).
!!
subroutine elsi_allocate_integer4_2d(e_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   integer(kind=i4),  intent(inout), allocatable :: array(:,:)
   integer(kind=i4),  intent(in)                 :: dim1
   integer(kind=i4),  intent(in)                 :: dim2
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*4

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 2D array of complex(kind=r4).
!!
subroutine elsi_allocate_complex8_2d(e_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   complex(kind=r4),  intent(inout), allocatable :: array(:,:)
   integer(kind=i4),  intent(in)                 :: dim1
   integer(kind=i4),  intent(in)                 :: dim2
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*8

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = (0.0_r4,0.0_r4)

end subroutine

!>
!! This routine allocates a 2D array of complex(kind=r8).
!!
subroutine elsi_allocate_complex16_2d(e_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   complex(kind=r8),  intent(inout), allocatable :: array(:,:)
   integer(kind=i4),  intent(in)                 :: dim1
   integer(kind=i4),  intent(in)                 :: dim2
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*16

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = (0.0_r8,0.0_r8)

end subroutine

!>
!! This routine allocates a 3D array of real(kind=r8).
!!
subroutine elsi_allocate_real8_3d(e_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   real(kind=r8),     intent(inout), allocatable :: array(:,:,:)
   integer(kind=i4),  intent(in)                 :: dim1
   integer(kind=i4),  intent(in)                 :: dim2
   integer(kind=i4),  intent(in)                 :: dim3
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*dim3*8

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1,dim2,dim3),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 3D array of integer(kind=i4).
!!
subroutine elsi_allocate_integer4_3d(e_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   integer(kind=i4),  intent(inout), allocatable :: array(:,:,:)
   integer(kind=i4),  intent(in)                 :: dim1
   integer(kind=i4),  intent(in)                 :: dim2
   integer(kind=i4),  intent(in)                 :: dim3
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*dim3*4

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1,dim2,dim3),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 3D array of complex(kind=r8).
!!
subroutine elsi_allocate_complex16_3d(e_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   complex(kind=r8),  intent(inout), allocatable :: array(:,:,:)
   integer(kind=i4),  intent(in)                 :: dim1
   integer(kind=i4),  intent(in)                 :: dim2
   integer(kind=i4),  intent(in)                 :: dim3
   character(len=*),  intent(in)                 :: arrayname
   character(len=*),  intent(in)                 :: caller

   real(kind=r8)    :: arraysize
   integer(kind=i4) :: error
   character*200    :: info_str

   if(e_h%print_mem) then
      arraysize = 1.0e-6_r8*dim1*dim2*dim3*16

      write(info_str,"(A,F10.3,2A)") "    Allocating ",arraysize," MB for ",&
         trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   allocate(array(dim1,dim2,dim3),stat=error)

   if(error > 0) then
      write(info_str,"(2A)") " Error in allocating ",trim(arrayname)
      call elsi_stop(info_str,e_h,caller)
   endif

   array = (0.0_r8,0.0_r8)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ELSI_DEALLOCATE MODULE PROCEDURES !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! This routine deallocates a 1D array with real(kind=r4).
!!
subroutine elsi_deallocate_real4_1d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   real(kind=r4),     intent(inout), allocatable :: array(:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 1D array with real(kind=r8).
!!
subroutine elsi_deallocate_real8_1d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   real(kind=r8),     intent(inout), allocatable :: array(:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 1D array with integer(kind=i4).
!!
subroutine elsi_deallocate_integer4_1d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   integer(kind=i4),  intent(inout), allocatable :: array(:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 1D array with integer(kind=i8).
!!
subroutine elsi_deallocate_integer8_1d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   integer(kind=i8),  intent(inout), allocatable :: array(:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 1D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex16_1d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   complex(kind=r8),  intent(inout), allocatable :: array(:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 2D array with real(kind=r4).
!!
subroutine elsi_deallocate_real4_2d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   real(kind=r4),     intent(inout), allocatable :: array(:,:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 2D array with real(kind=r8).
!!
subroutine elsi_deallocate_real8_2d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   real(kind=r8),     intent(inout), allocatable :: array(:,:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 2D array with integer(kind=i4).
!!
subroutine elsi_deallocate_integer4_2d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   integer(kind=i4),  intent(inout), allocatable :: array(:,:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 2D array with complex(kind=r4).
!!
subroutine elsi_deallocate_complex8_2d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   complex(kind=r4),  intent(inout), allocatable :: array(:,:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!>
!! This routine deallocates a 2D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex16_2d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   complex(kind=r8),  intent(inout), allocatable :: array(:,:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 3D array with real(kind=r8).
!!
subroutine elsi_deallocate_real8_3d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   real(kind=r8),     intent(inout), allocatable :: array(:,:,:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 3D array with integer(kind=i4).
!!
subroutine elsi_deallocate_integer4_3d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   integer(kind=i4),  intent(inout), allocatable :: array(:,:,:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

!>
!! This routine deallocates a 3D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex16_3d(e_h,array,arrayname)

   implicit none

   type(elsi_handle), intent(in)                 :: e_h
   complex(kind=r8),  intent(inout), allocatable :: array(:,:,:)
   character(len=*),  intent(in)                 :: arrayname

   character*200 :: info_str

   if(e_h%print_mem) then
      write(info_str,"(2A)") "    Deallocating ",trim(arrayname)
      call elsi_say(e_h,info_str)
   endif

   deallocate(array)

end subroutine

end module ELSI_MALLOC
