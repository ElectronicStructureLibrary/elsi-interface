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

!> ELSI Interchange
!! This module provides timers for ELSI.
!!

module ELSI_TIMERS

  use iso_c_binding
  use ELSI_DIMENSIONS

  implicit none
  private

  !> Timers for the total time
  logical:: total_time = .false.
  real*8 :: walltime_total
  real*8 :: walltime_total_start

  !> Timers for solving the eigenvalue problem
  logical:: solve_time = .false.
  real*8 :: walltime_solve_evp
  real*8 :: walltime_solve_evp_start

  !> system clock specifics
  integer :: clock_rate
  integer :: clock_max

  !> public subroutines
  public :: elsi_init_timers
  public :: elsi_print_setup
  public :: elsi_print_timers
  public :: elsi_start_total_time
  public :: elsi_stop_total_time
  public :: elsi_start_solve_evp_time
  public :: elsi_stop_solve_evp_time

contains

!>
!!  This routine sets all timers to zero.
!!
subroutine elsi_init_timers()

   implicit none

   integer :: initial_time

   walltime_total           = 0d0
   walltime_total_start     = 0d0
   walltime_solve_evp       = 0d0
   walltime_solve_evp_start = 0d0

   call system_clock (initial_time, clock_rate, clock_max)

end subroutine

!>
!!  This routine prints the setup.
!!
subroutine elsi_print_setup()

   implicit none

   if(myid == 0) then
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| Initial ELSI Output:                                  ')")
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| ELSI worked on a eigenvalue problem with dimensions   ')")
      write(*,"('| Rank                           : ',I13)") n_g_rank
      if(method == ELPA) then
         write(*,"('| Method:                        : ',A13)") "ELPA"
      endif
      if(method == LIBOMM) then
         write(*,"('| Method:                        : ',A13)") "libOMM"
      endif
      if(method == PEXSI) then
         write(*,"('| Method:                        : ',A13)") "PEXSI"
      endif
      if(method == CHESS) then
         write(*,"('| Method:                        : ',A13)") "CheSS"
      endif
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| Parallel Distribution:                                ')")
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| Process grid                   : ',I5,' x ',I5)")&
            n_p_rows, n_p_cols
      write(*,"('| Local matrix size process 0    : ',I5,' x ',I5)")&
            n_l_rows, n_l_cols
      write(*,"('| Local matrix blocking process 0: ',I5,' x ',I5)")&
            n_b_rows, n_b_cols
      write(*,"('|-------------------------------------------------------')")
   endif

end subroutine

!>
!!  This routine gets the current wallclock time.
!!
subroutine elsi_get_time(wtime)

   implicit none
   !> Wallclock time in seconds
   real*8, intent(out) :: wtime

   integer :: tics

   call system_clock(tics)

   wtime = 1d0 * tics / clock_rate

end subroutine

!>
!!  This routine starts the total timer.
!!
subroutine elsi_start_total_time()

   implicit none

   walltime_total           = 0d0
   walltime_total_start     = 0d0

   total_time = .true.
   call elsi_get_time(walltime_total_start)

end subroutine

!>
!!  This routine stops the total timer.
!!
subroutine elsi_stop_total_time()

   implicit none

   real*8 :: stop_time
   
   call elsi_get_time(stop_time)
   walltime_total = walltime_total &
     + stop_time - walltime_total_start

   if(myid == 0) then
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| Final ELSI Output:                                    ')")
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| ELSI worked on a eigenvalue problem with dimensions   ')")
      write(*,"('| Rank                           : ',I13)") n_g_rank
      write(*,"('| Number of States               : ',I13)") n_states
      if (method == ELPA) then
         write(*,"('| Method:                        : ',A13)") "ELPA"
      end if
      if(method == LIBOMM) then
         write(*,"('| Method:                        : ',A13)") "libOMM"
      endif
      if(method == PEXSI) then
         write(*,"('| Method:                        : ',A13)") "PEXSI"
      endif
      if(method == CHESS) then
         write(*,"('| Method:                        : ',A13)") "CheSS"
      endif

      if(total_time) then
         write(*,"('|-------------------------------------------------------')")
         write(*,"('| Total Time                     :',F13.3,' s')")&
               walltime_total
      endif
      write(*,"('|-------------------------------------------------------')")
   endif

end subroutine

!>
!!  This routine starts the solve evp timer.
!!
subroutine elsi_start_solve_evp_time()

   implicit none

   walltime_solve_evp       = 0d0
   walltime_solve_evp_start = 0d0

   solve_time = .true.
   call elsi_get_time(walltime_solve_evp_start)

end subroutine

!>
!!  This routine stops the solve evp timer.
!!
subroutine elsi_stop_solve_evp_time()

   implicit none

   real*8 :: stop_time

   call elsi_get_time(stop_time)
   walltime_solve_evp = walltime_solve_evp &
     + stop_time - walltime_solve_evp_start

   if(myid == 0) then
      if(solve_time) then
         write(*,"('|-------------------------------------------------------')")
         write(*,"('| Solving the Eigenvalue problem :',F13.3,' s')")&
               walltime_solve_evp
         write(*,"('|-------------------------------------------------------')")
      endif
   endif

end subroutine

end module ELSI_TIMERS
