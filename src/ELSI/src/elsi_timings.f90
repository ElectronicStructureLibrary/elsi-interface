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

! WPH: Since this is in progress, I'm putting my todo list directly in the code:
! - Create subroutines to resize arrays on-the-fly, should STARTING_SIZE_TIMINGS be
!   too small.
! - Add error checking to elsi_add_timing to make sure input strings aren't too long.
! - Do we time the solver itself, or the ELSI wrapper around the subroutine?
!   The latter makes sense for the decision layer, and the former for solver
!   development.  Both are useful.  Should do both.
! - Add ELSI interface subroutines for solver timings, i.e. elsi_print_solver_timings
!   as a wrapper around elsi_print_timings(e_h,e_h%solver_timings).  The calling code
!   should not have access to any ELSI timing handle, it should always go through an
!   interface subroutine defined on the ELSI handle.

! NOTES:
! - I've made the choice to not use elsi_malloc, because it considerably increases the
!   code density for little payoff; the allocated arrays here will be on the order of
!   kB, if that

!>
!! This module contains a collection of utilities related to timings in ELSI.
!!
module ELSI_TIMINGS

   use ELSI_CONSTANTS, only: TIMING_STRING_LEN,UNSET,UNSET_STRING
   use ELSI_DATATYPE,  only: elsi_handle,elsi_timings_handle
   use ELSI_PRECISION, only: i4,r8
   use ELSI_UTILS,     only: elsi_say

   implicit none

   private

   integer(kind=i4), parameter :: STARTING_SIZE_TIMINGS = 1024

   ! Global timing subroutines (i.e. those not attached to any timing handle)
   public  :: elsi_init_timer
   public  :: elsi_get_time
   ! ELSI Interface subroutines (for modifying a timing handle via the ELSI
   ! interface, but without exposing the timing API to the calling code)
   public  :: elsi_set_solver_timing_tag
   ! Timing handle subroutines
   public  :: elsi_init_timings
   public  :: elsi_add_timing
   public  :: elsi_print_timings
   public  :: elsi_finalize_timings
!   ! Error checking subroutines (private to this module)
!   private :: check_string_length

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GLOBAL TIMING SUBROUTINES !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! This routine initializes the timer.
!!
subroutine elsi_init_timer(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: initial_time
   integer(kind=i4) :: clock_max

   character*40, parameter :: caller = "elsi_init_timer"

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ELSI INTERFACE SUBROUTINES !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! This routine sets the next user_tag for the solver timings
!!
subroutine elsi_set_solver_timing_tag(e_h,user_tag)

   implicit none

   type(elsi_handle), intent(inout) :: e_h   !< Handle
   character(len=*),  intent(in)    :: user_tag

   character*40, parameter :: caller = "elsi_set_solver_timing_tag"

   e_h%solver_timings%next_user_tag = user_tag

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TIMING HANDLE SUBROUTINES !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! This routine initializes a timing handle.
!!
subroutine elsi_init_timings(t_h,set_label)

   implicit none

   type(elsi_timings_handle),  intent(inout) :: t_h       !< Handle
   character(len=*), optional, intent(in)    :: set_label 

   character*40, parameter :: caller = "elsi_init_timings"

   t_h%size_timings  = STARTING_SIZE_TIMINGS
   t_h%n_timings     = 0
   t_h%next_user_tag = UNSET_STRING
   if(present(set_label)) then
      t_h%set_label = set_label
   else
      t_h%set_label = UNSET_STRING
   endif

   allocate( t_h%times(t_h%size_timings) )
   allocate( t_h%elsi_tags(t_h%size_timings) )
   allocate( t_h%user_tags(t_h%size_timings) )

end subroutine

!>
!! This routine adds provided timing to current handle.
!!
subroutine elsi_add_timing(t_h,time,elsi_tag,user_tag_in,iter_in)

   implicit none

   type(elsi_timings_handle),           intent(inout) :: t_h   !< Handle
   real(kind=r8),                       intent(in)    :: time
   character(len=*),                    intent(in)    :: elsi_tag
   character(len=*),          optional, intent(in)    :: user_tag_in
   integer(kind=i4),          optional, intent(in)    :: iter_in

   character(len=TIMING_STRING_LEN) :: user_tag
   integer(kind=i4)                 :: iter

   character*40, parameter :: caller = "elsi_add_timings"
   
   ! By default, increment the number of iterations and update
   ! the next iteration
   ! If the iteration is manually specified, the number of 
   ! iterations will NOT be updated.
   if(present(iter_in)) then
      iter = iter_in
   else
      t_h%n_timings = t_h%n_timings + 1
      iter = t_h%n_timings 
   endif

   ! Temporary check to make sure we don't have more timings than
   ! the array size.  Arrays should be resized on-the-fly instead. 
   if(iter.le.t_h%size_timings) then
      if(present(user_tag_in)) then
         user_tag = user_tag_in
      else
         user_tag           = t_h%next_user_tag
         t_h%next_user_tag  = UNSET_STRING
      endif
  
      t_h%times(iter)     = time
      t_h%elsi_tags(iter) = elsi_tag
      t_h%user_tags(iter) = user_tag
   else
      t_h%n_timings = t_h%n_timings - 1
   !   call elsi_resize_timing_arrays(t_h)
   endif

end subroutine

!>
!! This routine prints out timings collected so far for current handle.
!!
subroutine elsi_print_timings(e_h,t_h)

   implicit none

   type(elsi_handle),          intent(in) :: e_h       !< Handle
   type(elsi_timings_handle),  intent(in) :: t_h !< Handle

   character*200    :: info_str
   integer          :: iter

   character*40, parameter :: caller = "elsi_print_timings"

   write(info_str,"(2X,A,A)")  "|   Timing Set:        ", t_h%set_label
   call elsi_say(e_h,info_str)
   write(info_str,"(2X,A,I4)") "|   Number of timings: ", t_h%n_timings
   call elsi_say(e_h,info_str)

   call elsi_say(e_h,"  |      #  system_clock [s]  elsi_tag             user_tag            ")
   do iter = 1, t_h%n_timings
      write(info_str,"(2X,A,I4,1X,F17.3,2X,A,1X,A)") &
           "|   ", &
           iter, &
           t_h%times(iter), &
           t_h%elsi_tags(iter), &
           t_h%user_tags(iter)
      call elsi_say(e_h,info_str)
   end do

end subroutine


!>
!! This routine deallocates the timings handle.
!!
subroutine elsi_finalize_timings(t_h)

   implicit none

   type(elsi_timings_handle), intent(inout)  :: t_h   !< Handle

   character*40, parameter :: caller = "elsi_finalize_timings"

   if(allocated(t_h%times)) then
      deallocate(t_h%times)
   endif
   if(allocated(t_h%elsi_tags)) then
      deallocate(t_h%elsi_tags)
   endif
   if(allocated(t_h%user_tags)) then
      deallocate(t_h%user_tags)
   endif
   t_h%n_timings     = 0
   t_h%size_timings  = UNSET
   t_h%next_user_tag = UNSET_STRING
   t_h%set_label     = UNSET_STRING

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ERROR CHECKING SUBROUTINES !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! WPH:  I don't like the following because it would need a dependence on e_h
!       due to MPI_Abort in subroutines which really shouldn't depend on e_h
!subroutine check_string_length(e_h,string_in,my_caller)
!
!   implicit none
!
!   type(elsi_handle), intent(in) :: e_h       !< Handle
!   character(len=*),  intent(in) :: string_in
!   character(len=*),  intent(in) :: my_caller
!
!   character*200 :: info_str
!
!   character*40, parameter :: caller = "check_string_length"
!
!   if (len(string_in).gt.TIMING_STRING_LEN) then
!      write(info_str,"(A,A,A,I4)") "String ", string_in, " exceeds the &
!           &maximum string length of ", TIMING_STRING_LEN
!      call elsi_stop(info_str,e_h,my_caller)
!   endif
!
!end subroutine

end module ELSI_TIMINGS
