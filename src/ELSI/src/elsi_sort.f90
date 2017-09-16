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
!! This module provides implementations of the heapsort algorithm.
!!
module ELSI_SORT

   use ELSI_PRECISION, only: r8,i4,i8

   implicit none

   private

   public :: elsi_heapsort

   interface swap
      module procedure swap_i8,&
                       swap_i4,&
                       swap_r8,&
                       swap_c16
   end interface

   interface downheap
      module procedure downheap_real_v1,&
                       downheap_real_v2,&
                       downheap_real_v3,&
                       downheap_real_v4,&
                       downheap_complex_v1,&
                       downheap_complex_v2,&
                       downheap_complex_v3,&
                       downheap_complex_v4
   end interface

   interface elsi_heapsort
      module procedure heapsort_real_v1,&
                       heapsort_real_v2,&
                       heapsort_real_v3,&
                       heapsort_real_v4,&
                       heapsort_complex_v1,&
                       heapsort_complex_v2,&
                       heapsort_complex_v3,&
                       heapsort_complex_v4
   end interface

contains

!>
!! This routine sorts an integer(kind=i8) array by heapsort, moves two
!! real(kind=r8) arrays and two integer(kind=i4) arrays accordingly.
!!
subroutine heapsort_real_v1(length,a_i8,b_r8,c_r8,d_i4,e_i4)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length) !< i8 array to be sorted
   real(kind=r8)   , intent(inout) :: b_r8(length) !< r8 array to be moved
   real(kind=r8)   , intent(inout) :: c_r8(length) !< r8 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length) !< i4 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i8,b_r8,c_r8,d_i4,e_i4,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i8,1,i)
      call swap(length,b_r8,1,i)
      call swap(length,c_r8,1,i)
      call swap(length,d_i4,1,i)
      call swap(length,e_i4,1,i)

      i = i-1

      call downheap(length,a_i8,b_r8,c_r8,d_i4,e_i4,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_real_v1(length,a_i8,b_r8,c_r8,d_i4,e_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length) !< i8 array to be sorted
   real(kind=r8)   , intent(inout) :: b_r8(length) !< r8 array to be moved
   real(kind=r8)   , intent(inout) :: c_r8(length) !< r8 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(in)    :: top          !< Top of heap
   integer(kind=i4), intent(in)    :: bottom       !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i8(w+1) > a_i8(w)) then
            w = w+1
         endif
      endif

      if(a_i8(v) >= a_i8(w)) then
         return
      else
         call swap(length,a_i8,v,w)
         call swap(length,b_r8,v,w)
         call swap(length,c_r8,v,w)
         call swap(length,d_i4,v,w)
         call swap(length,e_i4,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine sorts an integer(kind=i8) array by heapsort, moves one
!! real(kind=r8) array and two integer(kind=i4) arrays accordingly.
!!
subroutine heapsort_real_v2(length,a_i8,b_r8,c_i4,d_i4)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length) !< i8 array to be sorted
   real(kind=r8)   , intent(inout) :: b_r8(length) !< r8 array to be moved
   integer(kind=i4), intent(inout) :: c_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length) !< i4 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i8,b_r8,c_i4,d_i4,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i8,1,i)
      call swap(length,b_r8,1,i)
      call swap(length,c_i4,1,i)
      call swap(length,d_i4,1,i)

      i = i-1

      call downheap(length,a_i8,b_r8,c_i4,d_i4,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_real_v2(length,a_i8,b_r8,c_i4,d_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length) !< i8 array to be sorted
   real(kind=r8)   , intent(inout) :: b_r8(length) !< r8 array to be moved
   integer(kind=i4), intent(inout) :: c_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(in)    :: top          !< Top of heap
   integer(kind=i4), intent(in)    :: bottom       !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i8(w+1) > a_i8(w)) then
            w = w+1
         endif
      endif

      if(a_i8(v) >= a_i8(w)) then
         return
      else
         call swap(length,a_i8,v,w)
         call swap(length,b_r8,v,w)
         call swap(length,c_i4,v,w)
         call swap(length,d_i4,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine sorts an integer(kind=i4) array by heapsort, moves two
!! real(kind=r8) arrays and two integer(kind=i4) arrays accordingly.
!!
subroutine heapsort_real_v3(length,a_i4,b_r8,c_r8,d_i4,e_i4)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i4), intent(inout) :: a_i4(length) !< i4 array to be sorted
   real(kind=r8)   , intent(inout) :: b_r8(length) !< r8 array to be moved
   real(kind=r8)   , intent(inout) :: c_r8(length) !< r8 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length) !< i4 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i4,b_r8,c_r8,d_i4,e_i4,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i4,1,i)
      call swap(length,b_r8,1,i)
      call swap(length,c_r8,1,i)
      call swap(length,d_i4,1,i)
      call swap(length,e_i4,1,i)

      i = i-1

      call downheap(length,a_i4,b_r8,c_r8,d_i4,e_i4,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_real_v3(length,a_i4,b_r8,c_r8,d_i4,e_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i4), intent(inout) :: a_i4(length) !< i4 array to be sorted
   real(kind=r8)   , intent(inout) :: b_r8(length) !< r8 array to be moved
   real(kind=r8)   , intent(inout) :: c_r8(length) !< r8 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(in)    :: top          !< Top of heap
   integer(kind=i4), intent(in)    :: bottom       !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i4(w+1) > a_i4(w)) then
            w = w+1
         endif
      endif

      if(a_i4(v) >= a_i4(w)) then
         return
      else
         call swap(length,a_i4,v,w)
         call swap(length,b_r8,v,w)
         call swap(length,c_r8,v,w)
         call swap(length,d_i4,v,w)
         call swap(length,e_i4,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine sorts an integer(kind=i4) array by heapsort, moves one
!! real(kind=r8) arrays and two integer(kind=i4) arrays accordingly.
!!
subroutine heapsort_real_v4(length,a_i4,b_r8,c_i4,d_i4)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i4), intent(inout) :: a_i4(length) !< i4 array to be sorted
   real(kind=r8)   , intent(inout) :: b_r8(length) !< r8 array to be moved
   integer(kind=i4), intent(inout) :: c_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length) !< i4 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i4,b_r8,c_i4,d_i4,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i4,1,i)
      call swap(length,b_r8,1,i)
      call swap(length,c_i4,1,i)
      call swap(length,d_i4,1,i)

      i = i-1

      call downheap(length,a_i4,b_r8,c_i4,d_i4,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_real_v4(length,a_i4,b_r8,c_i4,d_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i4), intent(inout) :: a_i4(length) !< i4 array to be sorted
   real(kind=r8)   , intent(inout) :: b_r8(length) !< r8 array to be moved
   integer(kind=i4), intent(inout) :: c_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(in)    :: top          !< Top of heap
   integer(kind=i4), intent(in)    :: bottom       !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i4(w+1) > a_i4(w)) then
            w = w+1
         endif
      endif

      if(a_i4(v) >= a_i4(w)) then
         return
      else
         call swap(length,a_i4,v,w)
         call swap(length,b_r8,v,w)
         call swap(length,c_i4,v,w)
         call swap(length,d_i4,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine sorts an integer(kind=i8) array by heapsort, moves two
!! complex(kind=r8) arrays and two integer(kind=i4) arrays accordingly.
!!
subroutine heapsort_complex_v1(length,a_i8,b_c16,c_c16,d_i4,e_i4)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length)  !< i8 array to be sorted
   complex(kind=r8), intent(inout) :: b_c16(length) !< c16 array to be moved
   complex(kind=r8), intent(inout) :: c_c16(length) !< c16 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length)  !< i4 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i8,b_c16,c_c16,d_i4,e_i4,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i8,1,i)
      call swap(length,b_c16,1,i)
      call swap(length,c_c16,1,i)
      call swap(length,d_i4,1,i)
      call swap(length,e_i4,1,i)

      i = i-1

      call downheap(length,a_i8,b_c16,c_c16,d_i4,e_i4,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_complex_v1(length,a_i8,b_c16,c_c16,d_i4,e_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length)  !< i8 array to be sorted
   complex(kind=r8), intent(inout) :: b_c16(length) !< c16 array to be moved
   complex(kind=r8), intent(inout) :: c_c16(length) !< c16 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(in)    :: top           !< Top of heap
   integer(kind=i4), intent(in)    :: bottom        !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i8(w+1) > a_i8(w)) then
            w = w+1
         endif
      endif

      if(a_i8(v) >= a_i8(w)) then
         return
      else
         call swap(length,a_i8,v,w)
         call swap(length,b_c16,v,w)
         call swap(length,c_c16,v,w)
         call swap(length,d_i4,v,w)
         call swap(length,e_i4,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine sorts an integer(kind=i8) array by heapsort, moves one
!! complex(kind=r8) array and two integer(kind=i4) arrays accordingly.
!!
subroutine heapsort_complex_v2(length,a_i8,b_c16,c_i4,d_i4)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length)  !< i8 array to be sorted
   complex(kind=r8), intent(inout) :: b_c16(length) !< c16 array to be moved
   integer(kind=i4), intent(inout) :: c_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length)  !< i4 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i8,b_c16,c_i4,d_i4,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i8,1,i)
      call swap(length,b_c16,1,i)
      call swap(length,c_i4,1,i)
      call swap(length,d_i4,1,i)

      i = i-1

      call downheap(length,a_i8,b_c16,c_i4,d_i4,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_complex_v2(length,a_i8,b_c16,c_i4,d_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length)  !< i8 array to be sorted
   complex(kind=r8), intent(inout) :: b_c16(length) !< c16 array to be moved
   integer(kind=i4), intent(inout) :: c_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(in)    :: top           !< Top of heap
   integer(kind=i4), intent(in)    :: bottom        !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i8(w+1) > a_i8(w)) then
            w = w+1
         endif
      endif

      if(a_i8(v) >= a_i8(w)) then
         return
      else
         call swap(length,a_i8,v,w)
         call swap(length,b_c16,v,w)
         call swap(length,c_i4,v,w)
         call swap(length,d_i4,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine sorts an integer(kind=i4) array by heapsort, moves two
!! complex(kind=r8) arrays and two integer(kind=i4) arrays accordingly.
!!
subroutine heapsort_complex_v3(length,a_i4,b_c16,c_c16,d_i4,e_i4)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i4), intent(inout) :: a_i4(length)  !< i4 array to be sorted
   complex(kind=r8), intent(inout) :: b_c16(length) !< c16 array to be moved
   complex(kind=r8), intent(inout) :: c_c16(length) !< c16 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length)  !< i4 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i4,b_c16,c_c16,d_i4,e_i4,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i4,1,i)
      call swap(length,b_c16,1,i)
      call swap(length,c_c16,1,i)
      call swap(length,d_i4,1,i)
      call swap(length,e_i4,1,i)

      i = i-1

      call downheap(length,a_i4,b_c16,c_c16,d_i4,e_i4,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_complex_v3(length,a_i4,b_c16,c_c16,d_i4,e_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i4), intent(inout) :: a_i4(length)  !< i4 array to be sorted
   complex(kind=r8), intent(inout) :: b_c16(length) !< c16 array to be moved
   complex(kind=r8), intent(inout) :: c_c16(length) !< c16 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(in)    :: top           !< Top of heap
   integer(kind=i4), intent(in)    :: bottom        !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i4(w+1) > a_i4(w)) then
            w = w+1
         endif
      endif

      if(a_i4(v) >= a_i4(w)) then
         return
      else
         call swap(length,a_i4,v,w)
         call swap(length,b_c16,v,w)
         call swap(length,c_c16,v,w)
         call swap(length,d_i4,v,w)
         call swap(length,e_i4,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine sorts an integer(kind=i4) array by heapsort, moves one
!! complex(kind=r8) array and two integer(kind=i4) arrays accordingly.
!!
subroutine heapsort_complex_v4(length,a_i4,b_c16,c_i4,d_i4)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i4), intent(inout) :: a_i4(length)  !< i4 array to be sorted
   complex(kind=r8), intent(inout) :: b_c16(length) !< c16 array to be moved
   integer(kind=i4), intent(inout) :: c_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length)  !< i4 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i4,b_c16,c_i4,d_i4,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i4,1,i)
      call swap(length,b_c16,1,i)
      call swap(length,c_i4,1,i)
      call swap(length,d_i4,1,i)

      i = i-1

      call downheap(length,a_i4,b_c16,c_i4,d_i4,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_complex_v4(length,a_i4,b_c16,c_i4,d_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i4), intent(inout) :: a_i4(length)  !< i4 array to be sorted
   complex(kind=r8), intent(inout) :: b_c16(length) !< c16 array to be moved
   integer(kind=i4), intent(inout) :: c_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(in)    :: top           !< Top of heap
   integer(kind=i4), intent(in)    :: bottom        !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i4(w+1) > a_i4(w)) then
            w = w+1
         endif
      endif

      if(a_i4(v) >= a_i4(w)) then
         return
      else
         call swap(length,a_i4,v,w)
         call swap(length,b_c16,v,w)
         call swap(length,c_i4,v,w)
         call swap(length,d_i4,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine swaps two numbers in an integer(kind=i8) array.
!!
subroutine swap_i8(length,array,i,j)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i8), intent(inout) :: array(length) !< i8 array
   integer(kind=i4), intent(in)    :: i             !< Index to be swapped
   integer(kind=i4), intent(in)    :: j             !< Index to be swapped

   integer(kind=i8) :: tmp

   tmp      = array(i)
   array(i) = array(j)
   array(j) = tmp

end subroutine

!>
!! This routine swaps two numbers in an integer(kind=i4) array.
!!
subroutine swap_i4(length,array,i,j)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i4), intent(inout) :: array(length) !< i4 array
   integer(kind=i4), intent(in)    :: i             !< Index to be swapped
   integer(kind=i4), intent(in)    :: j             !< Index to be swapped

   integer(kind=i4) :: tmp

   tmp      = array(i)
   array(i) = array(j)
   array(j) = tmp

end subroutine

!>
!! This routine swaps two numbers in a real(kind=r8) array.
!!
subroutine swap_r8(length,array,i,j)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   real(kind=r8),    intent(inout) :: array(length) !< r8 array
   integer(kind=i4), intent(in)    :: i             !< Index to be swapped
   integer(kind=i4), intent(in)    :: j             !< Index to be swapped

   real(kind=r8) :: tmp

   tmp      = array(i)
   array(i) = array(j)
   array(j) = tmp

end subroutine

!>
!! This routine swaps two numbers in a complex(kind=r8) array.
!!
subroutine swap_c16(length,array,i,j)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   complex(kind=r8), intent(inout) :: array(length) !< c16 array
   integer(kind=i4), intent(in)    :: i             !< Index to be swapped
   integer(kind=i4), intent(in)    :: j             !< Index to be swapped

   complex(kind=r8) :: tmp

   tmp      = array(i)
   array(i) = array(j)
   array(j) = tmp

end subroutine

end module ELSI_SORT
