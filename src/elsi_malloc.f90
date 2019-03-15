! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide wrapped allocate and deallocate routines with basic error handling
!! and memory usage output.
!!
module ELSI_MALLOC

   use ELSI_DATATYPE, only: elsi_basic_t
   use ELSI_MPI, only: elsi_stop
   use ELSI_OUTPUT, only: elsi_say
   use ELSI_PRECISION, only: i4,i8,r4,r8

   implicit none

   private

   public :: elsi_allocate
   public :: elsi_deallocate

   interface elsi_allocate
      module procedure elsi_allocate_integer4_1d
      module procedure elsi_allocate_integer4_2d
      module procedure elsi_allocate_integer4_3d
      module procedure elsi_allocate_integer8_1d
      module procedure elsi_allocate_real8_1d
      module procedure elsi_allocate_real8_2d
      module procedure elsi_allocate_real8_3d
      module procedure elsi_allocate_real4_1d
      module procedure elsi_allocate_real4_2d
      module procedure elsi_allocate_complex16_1d
      module procedure elsi_allocate_complex16_2d
      module procedure elsi_allocate_complex16_3d
      module procedure elsi_allocate_complex8_2d
   end interface

   interface elsi_deallocate
      module procedure elsi_deallocate_integer4_1d
      module procedure elsi_deallocate_integer4_2d
      module procedure elsi_deallocate_integer4_3d
      module procedure elsi_deallocate_integer8_1d
      module procedure elsi_deallocate_real8_1d
      module procedure elsi_deallocate_real8_2d
      module procedure elsi_deallocate_real8_3d
      module procedure elsi_deallocate_real4_1d
      module procedure elsi_deallocate_real4_2d
      module procedure elsi_deallocate_complex16_1d
      module procedure elsi_deallocate_complex16_2d
      module procedure elsi_deallocate_complex16_3d
      module procedure elsi_deallocate_complex8_2d
   end interface

contains

!>
!! Allocate a 1D array with real(kind=r4).
!!
subroutine elsi_allocate_real4_1d(bh,array,dim_1,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   real(kind=r4), intent(inout), allocatable :: array(:)
   integer(kind=i4), intent(in) :: dim_1
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*4

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = 0.0_r4

end subroutine

!>
!! Allocate a 1D array with real(kind=r8).
!!
subroutine elsi_allocate_real8_1d(bh,array,dim_1,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout), allocatable :: array(:)
   integer(kind=i4), intent(in) :: dim_1
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*8

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = 0.0_r8

end subroutine

!>
!! Allocate a 1D array with integer(kind=i4).
!!
subroutine elsi_allocate_integer4_1d(bh,array,dim_1,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(inout), allocatable :: array(:)
   integer(kind=i4), intent(in) :: dim_1
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*4

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = 0

end subroutine

!>
!! Allocate a 1D array with integer(kind=i8).
!!
subroutine elsi_allocate_integer8_1d(bh,array,dim_1,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i8), intent(inout), allocatable :: array(:)
   integer(kind=i4), intent(in) :: dim_1
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*8

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = 0

end subroutine

!>
!! Allocate a 1D array with complex(kind=r8).
!!
subroutine elsi_allocate_complex16_1d(bh,array,dim_1,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout), allocatable :: array(:)
   integer(kind=i4), intent(in) :: dim_1
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*16

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = (0.0_r8,0.0_r8)

end subroutine

!>
!! Allocate a 2D array of real(kind=r4).
!!
subroutine elsi_allocate_real4_2d(bh,array,dim_1,dim_2,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   real(kind=r4), intent(inout), allocatable :: array(:,:)
   integer(kind=i4), intent(in) :: dim_1
   integer(kind=i4), intent(in) :: dim_2
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*dim_2*4

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1,dim_2),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = 0.0_r4

end subroutine

!>
!! Allocate a 2D array of real(kind=r8).
!!
subroutine elsi_allocate_real8_2d(bh,array,dim_1,dim_2,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout), allocatable :: array(:,:)
   integer(kind=i4), intent(in) :: dim_1
   integer(kind=i4), intent(in) :: dim_2
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*dim_2*8

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1,dim_2),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = 0.0_r8

end subroutine

!>
!! Allocate a 2D array of integer(kind=i4).
!!
subroutine elsi_allocate_integer4_2d(bh,array,dim_1,dim_2,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(inout), allocatable :: array(:,:)
   integer(kind=i4), intent(in) :: dim_1
   integer(kind=i4), intent(in) :: dim_2
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*dim_2*4

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1,dim_2),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = 0

end subroutine

!>
!! Allocate a 2D array of complex(kind=r4).
!!
subroutine elsi_allocate_complex8_2d(bh,array,dim_1,dim_2,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r4), intent(inout), allocatable :: array(:,:)
   integer(kind=i4), intent(in) :: dim_1
   integer(kind=i4), intent(in) :: dim_2
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*dim_2*8

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1,dim_2),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = (0.0_r4,0.0_r4)

end subroutine

!>
!! Allocate a 2D array of complex(kind=r8).
!!
subroutine elsi_allocate_complex16_2d(bh,array,dim_1,dim_2,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout), allocatable :: array(:,:)
   integer(kind=i4), intent(in) :: dim_1
   integer(kind=i4), intent(in) :: dim_2
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*dim_2*16

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1,dim_2),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = (0.0_r8,0.0_r8)

end subroutine

!>
!! Allocate a 3D array of real(kind=r8).
!!
subroutine elsi_allocate_real8_3d(bh,array,dim_1,dim_2,dim_3,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout), allocatable :: array(:,:,:)
   integer(kind=i4), intent(in) :: dim_1
   integer(kind=i4), intent(in) :: dim_2
   integer(kind=i4), intent(in) :: dim_3
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*dim_2*dim_3*8

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1,dim_2,dim_3),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = 0.0_r8

end subroutine

!>
!! Allocate a 3D array of integer(kind=i4).
!!
subroutine elsi_allocate_integer4_3d(bh,array,dim_1,dim_2,dim_3,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(inout), allocatable :: array(:,:,:)
   integer(kind=i4), intent(in) :: dim_1
   integer(kind=i4), intent(in) :: dim_2
   integer(kind=i4), intent(in) :: dim_3
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*dim_2*dim_3*4

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1,dim_2,dim_3),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = 0

end subroutine

!>
!! Allocate a 3D array of complex(kind=r8).
!!
subroutine elsi_allocate_complex16_3d(bh,array,dim_1,dim_2,dim_3,label,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout), allocatable :: array(:,:,:)
   integer(kind=i4), intent(in) :: dim_1
   integer(kind=i4), intent(in) :: dim_2
   integer(kind=i4), intent(in) :: dim_3
   character(len=*), intent(in) :: label
   character(len=*), intent(in) :: caller

   real(kind=r8) :: mem
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   if(bh%print_info > 2) then
      mem = 1.0e-6_r8*dim_1*dim_2*dim_3*16

      write(msg,"(4X,A,F10.3,2A)") "Allocating ",mem," MB for ",trim(label)
      call elsi_say(bh,msg)
   end if

   allocate(array(dim_1,dim_2,dim_3),stat=ierr)

   if(ierr > 0) then
      write(msg,"(2A)") "Error in allocating ",trim(label)
      call elsi_stop(bh,msg,caller)
   end if

   array = (0.0_r8,0.0_r8)

end subroutine

!>
!! Deallocate a 1D array with real(kind=r4).
!!
subroutine elsi_deallocate_real4_1d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   real(kind=r4), intent(inout), allocatable :: array(:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!! Deallocate a 1D array with real(kind=r8).
!!
subroutine elsi_deallocate_real8_1d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout), allocatable :: array(:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!! Deallocate a 1D array with integer(kind=i4).
!!
subroutine elsi_deallocate_integer4_1d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(inout), allocatable :: array(:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!! Deallocate a 1D array with integer(kind=i8).
!!
subroutine elsi_deallocate_integer8_1d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i8), intent(inout), allocatable :: array(:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!! Deallocate a 1D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex16_1d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout), allocatable :: array(:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!! Deallocate a 2D array with real(kind=r4).
!!
subroutine elsi_deallocate_real4_2d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   real(kind=r4), intent(inout), allocatable :: array(:,:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!! Deallocate a 2D array with real(kind=r8).
!!
subroutine elsi_deallocate_real8_2d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout), allocatable :: array(:,:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!! Deallocate a 2D array with integer(kind=i4).
!!
subroutine elsi_deallocate_integer4_2d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(inout), allocatable :: array(:,:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!! Deallocate a 2D array with complex(kind=r4).
!!
subroutine elsi_deallocate_complex8_2d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r4), intent(inout), allocatable :: array(:,:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!>
!! Deallocate a 2D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex16_2d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout), allocatable :: array(:,:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!! Deallocate a 3D array with real(kind=r8).
!!
subroutine elsi_deallocate_real8_3d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout), allocatable :: array(:,:,:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!! Deallocate a 3D array with integer(kind=i4).
!!
subroutine elsi_deallocate_integer4_3d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(inout), allocatable :: array(:,:,:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

!>
!! Deallocate a 3D array with complex(kind=r8).
!!
subroutine elsi_deallocate_complex16_3d(bh,array,label)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout), allocatable :: array(:,:,:)
   character(len=*), intent(in) :: label

   character(len=200) :: msg

   if(bh%print_info > 2) then
      write(msg,"(4X,2A)") "Deallocating ",trim(label)
      call elsi_say(bh,msg)
   end if

   deallocate(array)

end subroutine

end module ELSI_MALLOC
