! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide implementations of the heapsort algorithm.
!!
module ELSI_SORT

   use ELSI_PRECISION, only: r8,i4,i8

   implicit none

   private

   public :: elsi_heapsort
   public :: elsi_permute
   public :: elsi_unpermute

   interface swap
      module procedure swap_i8
      module procedure swap_i4
      module procedure swap_r8
   end interface

   interface downheap
      module procedure downheap_i4
      module procedure downheap_i8
      module procedure downheap_r8
   end interface

   interface elsi_heapsort
      module procedure heapsort_i4
      module procedure heapsort_i8
      module procedure heapsort_r8
   end interface

   interface elsi_permute
      module procedure permute_i4
      module procedure permute_r8
      module procedure permute_c16
   end interface

   interface elsi_unpermute
      module procedure unpermute_r8
      module procedure unpermute_c16
   end interface

contains

!>
!! Sort a real(kind=r8) array by heapsort, returns an array of permutation.
!!
subroutine heapsort_r8(length,array,perm)

   implicit none

   integer(kind=i4), intent(in) :: length
   real(kind=r8), intent(inout) :: array(length)
   integer(kind=i4), intent(out) :: perm(length)

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   do i = 1,length
      perm(i) = i
   end do

   top = length/2

   do i = top,1,-1
      call downheap(length,array,perm,i,length)
   end do

   i = length

   do while(i > 1)
      call swap(length,array,1,i)
      call swap(length,perm,1,i)

      i = i-1

      call downheap(length,array,perm,1,i)
   end do

end subroutine

!>
!! Restore the max-heap structure.
!!
subroutine downheap_r8(length,a_r8,b_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in) :: length
   real(kind=r8), intent(inout) :: a_r8(length)
   integer(kind=i4), intent(inout) :: b_i4(length)
   integer(kind=i4), intent(in) :: top
   integer(kind=i4), intent(in) :: bottom

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while(w <= bottom)
      if(w+1 <= bottom) then
         if(a_r8(w+1) > a_r8(w)) then
            w = w+1
         end if
      end if

      if(a_r8(v) >= a_r8(w)) then
         return
      else
         call swap(length,a_r8,v,w)
         call swap(length,b_i4,v,w)

         v = w
         w = 2*v
      end if
   end do

end subroutine

!>
!! Sort an integer(kind=i8) array by heapsort, returns an array of permutation.
!!
subroutine heapsort_i8(length,array,perm)

   implicit none

   integer(kind=i4), intent(in) :: length
   integer(kind=i8), intent(inout) :: array(length)
   integer(kind=i4), intent(out) :: perm(length)

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   do i = 1,length
      perm(i) = i
   end do

   top = length/2

   do i = top,1,-1
      call downheap(length,array,perm,i,length)
   end do

   i = length

   do while(i > 1)
      call swap(length,array,1,i)
      call swap(length,perm,1,i)

      i = i-1

      call downheap(length,array,perm,1,i)
   end do

end subroutine

!>
!! Restore the max-heap structure.
!!
subroutine downheap_i8(length,a_i8,b_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in) :: length
   integer(kind=i8), intent(inout) :: a_i8(length)
   integer(kind=i4), intent(inout) :: b_i4(length)
   integer(kind=i4), intent(in) :: top
   integer(kind=i4), intent(in) :: bottom

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while(w <= bottom)
      if(w+1 <= bottom) then
         if(a_i8(w+1) > a_i8(w)) then
            w = w+1
         end if
      end if

      if(a_i8(v) >= a_i8(w)) then
         return
      else
         call swap(length,a_i8,v,w)
         call swap(length,b_i4,v,w)

         v = w
         w = 2*v
      end if
   end do

end subroutine

!>
!! Sort an integer(kind=i4) array by heapsort, returns an array of permutation.
!!
subroutine heapsort_i4(length,array,perm)

   implicit none

   integer(kind=i4), intent(in) :: length
   integer(kind=i4), intent(inout) :: array(length)
   integer(kind=i4), intent(out) :: perm(length)

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   do i = 1,length
      perm(i) = i
   end do

   top = length/2

   do i = top,1,-1
      call downheap(length,array,perm,i,length)
   end do

   i = length

   do while(i > 1)
      call swap(length,array,1,i)
      call swap(length,perm,1,i)

      i = i-1

      call downheap(length,array,perm,1,i)
   end do

end subroutine

!>
!! Restore the max-heap structure.
!!
subroutine downheap_i4(length,a_i4,b_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in) :: length
   integer(kind=i4), intent(inout) :: a_i4(length)
   integer(kind=i4), intent(inout) :: b_i4(length)
   integer(kind=i4), intent(in) :: top
   integer(kind=i4), intent(in) :: bottom

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while(w <= bottom)
      if(w+1 <= bottom) then
         if(a_i4(w+1) > a_i4(w)) then
            w = w+1
         end if
      end if

      if(a_i4(v) >= a_i4(w)) then
         return
      else
         call swap(length,a_i4,v,w)
         call swap(length,b_i4,v,w)

         v = w
         w = 2*v
      end if
   end do

end subroutine

!>
!! Permute an integer(kind=i4) array.
!!
subroutine permute_i4(length,perm,array)

   implicit none

   integer(kind=i4), intent(in) :: length
   integer(kind=i4), intent(in) :: perm(length)
   integer(kind=i4), intent(inout) :: array(length)

   integer(kind=i4) :: i
   integer(kind=i4) :: work(length)

   do i = 1,length
      work(i) = array(perm(i))
   end do

   array = work

end subroutine

!>
!! Permute a real(kind=r8) array.
!!
subroutine permute_r8(length,perm,array)

   implicit none

   integer(kind=i4), intent(in) :: length
   integer(kind=i4), intent(in) :: perm(length)
   real(kind=r8), intent(inout) :: array(length)

   real(kind=r8) :: work(length)
   integer(kind=i4) :: i

   do i = 1,length
      work(i) = array(perm(i))
   end do

   array = work

end subroutine

!>
!! Permute a complex(kind=r8) array.
!!
subroutine permute_c16(length,perm,array)

   implicit none

   integer(kind=i4), intent(in) :: length
   integer(kind=i4), intent(in) :: perm(length)
   complex(kind=r8), intent(inout) :: array(length)

   complex(kind=r8) :: work(length)
   integer(kind=i4) :: i

   do i = 1,length
      work(i) = array(perm(i))
   end do

   array = work

end subroutine

!>
!! Unpermute a real(kind=r8) array.
!!
subroutine unpermute_r8(length,perm,array)

   implicit none

   integer(kind=i4), intent(in) :: length
   integer(kind=i4), intent(in) :: perm(length)
   real(kind=r8), intent(inout) :: array(length)

   real(kind=r8) :: work(length)
   integer(kind=i4) :: i

   do i = 1,length
      work(perm(i)) = array(i)
   end do

   array = work

end subroutine

!>
!! Unpermute a complex(kind=r8) array.
!!
subroutine unpermute_c16(length,perm,array)

   implicit none

   integer(kind=i4), intent(in) :: length
   integer(kind=i4), intent(in) :: perm(length)
   complex(kind=r8), intent(inout) :: array(length)

   complex(kind=r8) :: work(length)
   integer(kind=i4) :: i

   do i = 1,length
      work(perm(i)) = array(i)
   end do

   array = work

end subroutine

!>
!! Swap two numbers in an integer(kind=i8) array.
!!
subroutine swap_i8(length,array,i,j)

   implicit none

   integer(kind=i4), intent(in) :: length
   integer(kind=i8), intent(inout) :: array(length)
   integer(kind=i4), intent(in) :: i
   integer(kind=i4), intent(in) :: j

   integer(kind=i8) :: tmp

   tmp = array(i)
   array(i) = array(j)
   array(j) = tmp

end subroutine

!>
!! Swap two numbers in an integer(kind=i4) array.
!!
subroutine swap_i4(length,array,i,j)

   implicit none

   integer(kind=i4), intent(in) :: length
   integer(kind=i4), intent(inout) :: array(length)
   integer(kind=i4), intent(in) :: i
   integer(kind=i4), intent(in) :: j

   integer(kind=i4) :: tmp

   tmp = array(i)
   array(i) = array(j)
   array(j) = tmp

end subroutine

!>
!! Swap two numbers in a real(kind=r8) array.
!!
subroutine swap_r8(length,array,i,j)

   implicit none

   integer(kind=i4), intent(in) :: length
   real(kind=r8), intent(inout) :: array(length)
   integer(kind=i4), intent(in) :: i
   integer(kind=i4), intent(in) :: j

   real(kind=r8) :: tmp

   tmp = array(i)
   array(i) = array(j)
   array(j) = tmp

end subroutine

end module ELSI_SORT
