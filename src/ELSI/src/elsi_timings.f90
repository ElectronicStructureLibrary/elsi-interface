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
! - Unify naming scheme.  These are all solver timings.  We may want other timings in
!   the future, so should be specific.
!   - Abstract out timings into their own datatype?
! - Create subroutines to resize arrays on-the-fly, should STARTING_MAX_N_TIMINGS be
!   too small.
! - Extend the elsi_{,de}allocate subroutines to arrays used here and replace bare
!   {,de}allocate statements.
! - Do we time the solver itself, or the ELSI wrapper around the subroutine?
!   The latter makes sense for the decision layer, and the former for solver
!   development.  Both are useful.  Should do both.
! - Possibly distinguish between DM and EV as its own tag?
! - Eliminate hard-coded format string in elsi_print_timings
! - Move pre-existing timing subroutines in elsi_utils here.  Only once this module is
!   stable-sh, though.

!>
!! This module contains a collection of utilities related to tracking timings in ELSI.
!!
module ELSI_TIMINGS

   use ELSI_CONSTANTS, only: TIMER_STRING_LEN,UNSET
   use ELSI_DATATYPE,  only: elsi_handle
   use ELSI_PRECISION, only: i4, r8
   use ELSI_UTILS,     only: elsi_say

   implicit none

   private

   integer(kind=i4), parameter :: STARTING_MAX_N_TIMINGS = 1024

   public :: elsi_init_timings
   public :: elsi_add_timing
   public :: elsi_print_timings
   public :: elsi_finalize_timings

contains

!>
!! This routine initializes timings.
!!
subroutine elsi_init_timings(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*40, parameter :: caller = "elsi_init_timings"

   e_h%max_n_timings = STARTING_MAX_N_TIMINGS
   e_h%n_timings     = 0
   allocate( e_h%solver_timings_times(e_h%max_n_timings) )
   allocate( e_h%solver_timings_solvers(e_h%max_n_timings) )
   allocate( e_h%solver_timings_tags(e_h%max_n_timings) )

end subroutine

!>
!! This routine adds timings to the list of timings
!!
subroutine elsi_add_timing(e_h,time,tag_in,solver_in,iter_in)

   implicit none

   type(elsi_handle),                      intent(inout) :: e_h   !< Handle
   real(kind=r8),                             intent(in) :: time
   character(len=TIMER_STRING_LEN), optional, intent(in) :: tag_in
   integer(kind=i4),                optional, intent(in) :: solver_in
   integer(kind=i4),                optional, intent(in) :: iter_in

   integer(kind=i4)                :: iter
   integer(kind=i4)                :: solver
   character(len=TIMER_STRING_LEN) :: tag

   character*40, parameter :: caller = "elsi_add_timings"

   ! By default, use UNSET for the tag
   if(present(tag_in)) then
      tag = tag_in
   else
      tag = "UNSET"
   endif

   ! By default, use the handle's solver for the solver
   if(present(solver_in)) then
      solver = solver_in
   else
      solver = e_h%solver
   endif

   ! By default, increment the number of iterations and update
   ! the next iteration
   ! If the iteration is manually specified, the number of 
   ! iterations will NOT be updated.
   if(present(iter_in)) then
      iter = iter_in
   else
      e_h%n_timings = e_h%n_timings + 1
      iter = e_h%n_timings 
   endif

   if(iter.le.e_h%max_n_timings) then
      e_h%solver_timings_times(iter)   = time
      e_h%solver_timings_solvers(iter) = solver
      e_h%solver_timings_tags(iter)    = tag
  !else
  !   call elsi_resize_timing_arrays(e_h)
   endif

end subroutine

!>
!! This routine prints out timings collected so far
!!
subroutine elsi_print_timings(e_h)

   implicit none

   type(elsi_handle), intent(in) :: e_h !< Handle

   character*200 :: info_str
   integer       :: iter

   character*40, parameter :: caller = "elsi_print_timings"

   call elsi_say(e_h,"  |------------------------------------------")
   call elsi_say(e_h,"  | ELSI Timings                             ")
   call elsi_say(e_h,"  |------------------------------------------")

   do iter = 1, e_h%n_timings
      write(info_str,"(I7,1X,F13.3,1X,I2,1X,A100)") iter, e_h%solver_timings_times(iter), &
           e_h%solver_timings_solvers(iter), e_h%solver_timings_tags(iter)
      call elsi_say(e_h,info_str)
   end do

end subroutine


!>
!! This routine deallocates timings.
!!
subroutine elsi_finalize_timings(e_h)

   implicit none

   type(elsi_handle), intent(inout)  :: e_h   !< Handle

   character*40, parameter :: caller = "elsi_finalize_timings"

   if(allocated(e_h%solver_timings_times)) then
      deallocate(e_h%solver_timings_times)
   endif
   if(allocated(e_h%solver_timings_solvers)) then
      deallocate(e_h%solver_timings_solvers)
   endif
   if(allocated(e_h%solver_timings_tags)) then
      deallocate(e_h%solver_timings_tags)
   endif
   e_h%n_timings     = 0
   e_h%max_n_timings = UNSET

end subroutine

end module ELSI_TIMINGS
