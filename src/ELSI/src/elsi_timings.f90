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
!! This module contains a collection of utilities related to timings in ELSI.
!!
module ELSI_TIMINGS

   use ELSI_CONSTANTS, only: TIMING_STRING_LEN,UNSET,UNSET_STRING
   use ELSI_DATATYPE,  only: elsi_handle,elsi_timings_handle
   use ELSI_IO,        only: elsi_say
   use ELSI_PRECISION, only: i4,r8

   implicit none

   private

   integer(kind=i4), parameter :: STARTING_SIZE_TIMINGS = 1

   public :: elsi_init_timer
   public :: elsi_get_time
   public :: elsi_init_timings
   public :: elsi_add_timing
   public :: elsi_print_timings
   public :: elsi_finalize_timings

contains

!>
!! This routine initializes the timer.
!!
subroutine elsi_init_timer(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   integer(kind=i4) :: initial_time
   integer(kind=i4) :: clock_max

   character(len=40), parameter :: caller = "elsi_init_timer"

   if(e_h%handle_ready) then
      e_h%handle_changed = .true.
   endif

   call system_clock(initial_time,e_h%clock_rate,clock_max)

end subroutine

!>
!! This routine gets the current wallclock time.
!!
subroutine elsi_get_time(e_h,wtime)

   implicit none

   type(elsi_handle), intent(in)  :: e_h
   real(kind=r8),     intent(out) :: wtime

   integer(kind=i4) :: tics

   character(len=40), parameter :: caller = "elsi_get_time"

   call system_clock(tics)

   wtime = 1.0_r8*tics/e_h%clock_rate

end subroutine

!>
!! This routine initializes a timing handle.
!!
subroutine elsi_init_timings(t_h,set_label)

   implicit none

   type(elsi_timings_handle),  intent(inout)        :: t_h
   character(len=*),           intent(in), optional :: set_label

   character(len=40), parameter :: caller = "elsi_init_timings"

   t_h%size_timings = STARTING_SIZE_TIMINGS
   t_h%n_timings    = 0
   t_h%user_tag     = UNSET_STRING

   if(present(set_label)) then
      t_h%set_label = set_label
   else
      t_h%set_label = UNSET_STRING
   endif

   allocate(t_h%times(t_h%size_timings))
   allocate(t_h%elsi_tags(t_h%size_timings))
   allocate(t_h%user_tags(t_h%size_timings))

end subroutine

!>
!! This routine resizes arrays in the timing handle by a factor of two.
!!
subroutine elsi_resize_timing_arrays(t_h)

   implicit none

   type(elsi_timings_handle), intent(inout) :: t_h

   integer :: i_timing
   real(kind=r8),                    allocatable :: tmp_times(:)
   character(len=TIMING_STRING_LEN), allocatable :: tmp_elsi_tags(:)
   character(len=TIMING_STRING_LEN), allocatable :: tmp_user_tags(:)

   character(len=40), parameter :: caller = "elsi_resize_timing_arrays"

   allocate(tmp_times(t_h%size_timings))
   allocate(tmp_elsi_tags(t_h%size_timings))
   allocate(tmp_user_tags(t_h%size_timings))

   ! Save old timings
   do i_timing = 1,t_h%size_timings
      tmp_times(i_timing)     = t_h%times(i_timing)
      tmp_elsi_tags(i_timing) = t_h%elsi_tags(i_timing)
      tmp_user_tags(i_timing) = t_h%user_tags(i_timing)
   enddo

   ! Resize arrays by factor of two
   deallocate(t_h%times)
   deallocate(t_h%elsi_tags)
   deallocate(t_h%user_tags)

   t_h%size_timings = 2*t_h%size_timings

   allocate(t_h%times(t_h%size_timings))
   allocate(t_h%elsi_tags(t_h%size_timings))
   allocate(t_h%user_tags(t_h%size_timings))

   ! Fill in arrays with old timings
   do i_timing = 1,t_h%size_timings/2
      t_h%times(i_timing)     = tmp_times(i_timing)
      t_h%elsi_tags(i_timing) = tmp_elsi_tags(i_timing)
      t_h%user_tags(i_timing) = tmp_user_tags(i_timing)
   enddo

   deallocate(tmp_times)
   deallocate(tmp_elsi_tags)
   deallocate(tmp_user_tags)

end subroutine

!>
!! This routine adds provided timing to current handle.
!!
subroutine elsi_add_timing(t_h,time,elsi_tag,user_tag_in,iter_in)

   implicit none

   type(elsi_timings_handle), intent(inout)        :: t_h
   real(kind=r8),             intent(in)           :: time
   character(len=*),          intent(in)           :: elsi_tag
   character(len=*),          intent(in), optional :: user_tag_in
   integer(kind=i4),          intent(in), optional :: iter_in

   character(len=TIMING_STRING_LEN) :: user_tag
   integer(kind=i4)                 :: iter

   character(len=40), parameter :: caller = "elsi_add_timings"

   if(present(iter_in)) then
      iter = iter_in
   else
      t_h%n_timings = t_h%n_timings+1
      iter          = t_h%n_timings
   endif

   ! If more timings than maximum array size, resize arrays
   do while (iter >= t_h%size_timings)
     call elsi_resize_timing_arrays(t_h)
   enddo

   if(present(user_tag_in)) then
      user_tag = user_tag_in
   else
      user_tag = t_h%user_tag
   endif

   t_h%times(iter)     = time
   t_h%elsi_tags(iter) = elsi_tag
   t_h%user_tags(iter) = user_tag

end subroutine

!>
!! This routine prints out timings collected so far for current handle.
!!
subroutine elsi_print_timings(e_h,t_h)

   implicit none

   type(elsi_handle),         intent(in) :: e_h
   type(elsi_timings_handle), intent(in) :: t_h

   character(len=200) :: info_str
   integer(kind=i4)   :: iter

   character(len=40), parameter :: caller = "elsi_print_timings"

   write(info_str,"(A,A)") "Timing Set:        ",t_h%set_label
   call elsi_say(e_h,info_str)
   write(info_str,"(A,I4)") "Number of timings: ",t_h%n_timings
   call elsi_say(e_h,info_str)
   call elsi_say(e_h,"   #  system_clock [s]  elsi_tag             user_tag            ")

   do iter = 1,t_h%n_timings
      write(info_str,"(I4,1X,F17.3,2X,A,1X,A)") iter,t_h%times(iter),&
         t_h%elsi_tags(iter),t_h%user_tags(iter)
      call elsi_say(e_h,info_str)
   enddo

end subroutine

!>
!! This routine deallocates the timings handle.
!!
subroutine elsi_finalize_timings(t_h)

   implicit none

   type(elsi_timings_handle), intent(inout) :: t_h

   character(len=40), parameter :: caller = "elsi_finalize_timings"

   if(allocated(t_h%times)) then
      deallocate(t_h%times)
   endif
   if(allocated(t_h%elsi_tags)) then
      deallocate(t_h%elsi_tags)
   endif
   if(allocated(t_h%user_tags)) then
      deallocate(t_h%user_tags)
   endif

   t_h%n_timings    = 0
   t_h%size_timings = UNSET
   t_h%user_tag     = UNSET_STRING
   t_h%set_label    = UNSET_STRING

end subroutine

end module ELSI_TIMINGS
