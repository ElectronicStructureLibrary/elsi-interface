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

!>
!! This module provides timers for ELSI.
!!

module ELSI_TIMERS

   use iso_c_binding
   use ELSI_DIMENSIONS

   implicit none
   private

   real*8  :: walltime_solve_evp
   real*8  :: walltime_solve_evp_start

   real*8  :: walltime_2dbc_to_1dccs
   real*8  :: walltime_2dbc_to_1dccs_start
  
   real*8  :: walltime_1dccs_to_2dbc
   real*8  :: walltime_1dccs_to_2dbc_start

   integer :: clock_rate
   integer :: clock_max

   public :: elsi_init_timers
   public :: elsi_print_timers
   public :: elsi_start_solve_evp_time
   public :: elsi_stop_solve_evp_time
   public :: elsi_start_2dbc_to_1dccs_time
   public :: elsi_stop_2dbc_to_1dccs_time
   public :: elsi_start_1dccs_to_2dbc_time
   public :: elsi_stop_1dccs_to_2dbc_time

contains

!>
!! This routine sets all timers to zero.
!!
subroutine elsi_init_timers()

   implicit none

   integer :: initial_time

   walltime_solve_evp       = 0d0
   walltime_solve_evp_start = 0d0

   walltime_2dbc_to_1dccs       = 0d0
   walltime_2dbc_to_1dccs_start = 0d0

   walltime_1dccs_to_2dbc       = 0d0
   walltime_1dccs_to_2dbc_start = 0d0

   call system_clock(initial_time, clock_rate, clock_max)

end subroutine

!>
!! This routine prints the timing results.
!!
subroutine elsi_print_timers()

   implicit none

   integer :: i_proc

   if(myid == 0) then
      write(*,"('  |-------------------------------------------------------')")
      write(*,"('  | Final ELSI Output:                                    ')")
      write(*,"('  |-------------------------------------------------------')")
      write(*,"('  | Eigenvalue problem size             : ',I13)") n_g_size
      if(method == PEXSI) then
         write(*,"('  | Non zero elements                   : ',I13)") &
               nnz_g
         write(*,"('  | Sparsity                            : ',F13.3)") &
               1d0-(1d0*nnz_g/n_g_size)/n_g_size
      endif
      write(*,"('  | Number of electrons                 : ',F13.1)") n_electrons
      write(*,"('  | Number of states                    : ',I13)") n_states
      if(method == ELPA) then
         write(*,"('  | Method                              : ',A13)") "ELPA"
      endif
      if(method == LIBOMM) then
         write(*,"('  | Method                              : ',A13)") "libOMM"
      endif
      if(method == PEXSI) then
         write(*,"('  | Method                              : ',A13)") "PEXSI"
      endif
      write(*,"('  |-------------------------------------------------------')")
      write(*,"('  | ELSI Project (c)  elsi-interchange.org                ')")
      write(*,"('  |-------------------------------------------------------')")
   endif

end subroutine

!>
!! This routine gets the current wallclock time.
!!
subroutine elsi_get_time(wtime)

   implicit none
   real*8, intent(out) :: wtime
 
   integer :: tics

   call system_clock(tics)

   wtime = 1d0 * tics / clock_rate

end subroutine

!>
!! This routine starts solve_evp timer.
!!
subroutine elsi_start_solve_evp_time()

   implicit none
   
   call elsi_get_time(walltime_solve_evp_start)

end subroutine

!>
!! This routine ends solve_evp timer.
!!
subroutine elsi_stop_solve_evp_time()

   implicit none

   real*8 :: stop_time
   
   call elsi_get_time(stop_time)
   walltime_solve_evp = stop_time - walltime_solve_evp_start

   if(myid == 0) write(*,"('  | ELSI solver finished in ',F13.3,' s')") &
                       walltime_solve_evp

end subroutine

!>
!! This routine starts 2dbc_to_1dccs timer.
!!
subroutine elsi_start_2dbc_to_1dccs_time()

   implicit none
   
   call elsi_get_time(walltime_2dbc_to_1dccs_start)

end subroutine

!>
!! This routine ends 2dbc_to_1dccs timer.
!!
subroutine elsi_stop_2dbc_to_1dccs_time()

   implicit none

   real*8 :: stop_time
   
   call elsi_get_time(stop_time)
   walltime_2dbc_to_1dccs = stop_time - walltime_2dbc_to_1dccs_start

   if(myid == 0) write(*,"('  | ELSI matrix conversion done in ',F13.3,' s')") &
                       walltime_2dbc_to_1dccs

end subroutine

!>
!! This routine starts 1dccs_to_2dbc timer.
!!
subroutine elsi_start_1dccs_to_2dbc_time()

   implicit none

   call elsi_get_time(walltime_1dccs_to_2dbc_start)

end subroutine

!>
!! This routine ends 1dccs_to_2dbc timer.
!!
subroutine elsi_stop_1dccs_to_2dbc_time()

   implicit none

   real*8 :: stop_time

   call elsi_get_time(stop_time)
   walltime_1dccs_to_2dbc = stop_time - walltime_1dccs_to_2dbc_start

   if(myid == 0) write(*,"('  | ELSI matrix conversion done in ',F13.3,' s')") &
                       walltime_1dccs_to_2dbc

end subroutine

end module ELSI_TIMERS
