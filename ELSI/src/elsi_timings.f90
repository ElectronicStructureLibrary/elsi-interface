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
!! This module is the actual ELSI Interface, providing functions for setting up
!! and solving or circumventing an eigenvalue problem using ELPA, OMM, or PEXSI
!!

module ELSI_TIMERS

  use iso_c_binding
  use ELSI_DIMENSIONS

  implicit none
  private

  !> Timers for the file reading of the eigenvalue problem
  logical:: total_time = .false.
  real*8 :: walltime_total  
  real*8 :: walltime_total_start  

  !> Timers for the file reading of the eigenvalue problem
  logical:: read_time = .false.
  real*8 :: walltime_read_evp  
  real*8 :: walltime_read_evp_start  
  
  !> Timers for the file writing of the eigenvalue problem
  logical:: write_time = .false.
  real*8 :: walltime_write_evp  
  real*8 :: walltime_write_evp_start  

  !> Timers for solving the eigenvalue problem
  logical:: solve_time = .false.
  real*8 :: walltime_solve_evp  
  real*8 :: walltime_solve_evp_start

  !> Timers for solving the eigenvalue problem
  logical:: dist_pexsi_time = .false.
  real*8 :: walltime_dist_pexsi  
  real*8 :: walltime_dist_pexsi_start   
  
  !> system clock specifics
  integer :: clock_rate
  integer :: clock_max

  !> public subroutines
  public :: elsi_initialize_timers 
  public :: elsi_print_setup 
  public :: elsi_print_timers 
  public :: elsi_start_total_time 
  public :: elsi_stop_total_time
  public :: elsi_start_read_evp_time
  public :: elsi_stop_read_evp_time
  public :: elsi_start_write_evp_time
  public :: elsi_stop_write_evp_time
  public :: elsi_start_solve_evp_time
  public :: elsi_stop_solve_evp_time
  public :: elsi_start_dist_pexsi_time
  public :: elsi_stop_dist_pexsi_time



contains

!>
!!  This routine sets all timers to zero
!!
subroutine elsi_initialize_timers()

   implicit none

   integer :: initial_time

   !> Timers for the file reading of the eigenvalue problem
   walltime_total        = 0d0
   walltime_total_start  = 0d0

   !> Timers for the file reading of the eigenvalue problem
   walltime_read_evp        = 0d0 
   walltime_read_evp_start  = 0d0
   
   !> Timers for the file writing of the eigenvalue problem
   walltime_write_evp        = 0d0 
   walltime_write_evp_start  = 0d0

   !> Timers for solving the eigenvalue problem
   walltime_solve_evp        = 0d0
   walltime_solve_evp_start  = 0d0

    !> Timers for communicating the eigenvalue problem to PEXSI style
   walltime_dist_pexsi        = 0d0
   walltime_dist_pexsi_start  = 0d0

   call system_clock (initial_time, clock_rate, clock_max)

end subroutine

!>
!!  This routine prints the timing results
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
      if(method == OMM_DENSE) then 
         write(*,"('| Method:                        : ',A13)") "OMM_DENSE"
      endif
      if(method == PEXSI) then 
         write(*,"('| Method:                        : ',A13)") "PEXSI"
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
!!  This routine prints the timing results
!!
subroutine elsi_print_timers()

   implicit none

   integer :: i_proc

   if(myid == 0) then
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| Final ELSI Output:                                    ')")
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| ELSI worked on a eigenvalue problem with dimensions   ')")
      write(*,"('| Rank                           : ',I13)") n_g_rank
      write(*,"('| Non zero elements              : ',I13)") n_g_nonzero
      write(*,"('| Sparcity                       : ',F13.3)") &
            (1d0 * n_g_nonzero / n_g_rank) / n_g_rank
      write(*,"('| Number of Electrons            : ',F13.3)") &
            n_electrons
      write(*,"('| Number of States               : ',I13)") n_eigenvectors   
      if (method == ELPA) then 
      write(*,"('| Method:                        : ',A13)") "ELPA"
      end if
      if(method == OMM_DENSE) then 
         write(*,"('| Method:                        : ',A13)") "OMM_DENSE"
      endif
      if(method == PEXSI) then 
         write(*,"('| Method:                        : ',A13)") "PEXSI"
      endif
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| Parallel Distribution:                                ')")
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| Process grid                   : ',I5,' x ',I5)") &
            n_p_rows, n_p_cols
    endif

    do i_proc = 0, n_procs - 1
       if(i_proc == myid) then
          write(*,"('| Local matrix size process     ',I5,' : ',I5,' x ',I5)")&
                myid, n_l_rows, n_l_cols
          write(*,"('| Local matrix blocking process ',I5,' : ',I5,' x ',I5)")&
                myid, n_b_rows, n_b_cols
       endif
       call MPI_Barrier(mpi_comm_global,mpierr)
    enddo

    if(myid == 0) then
       write(*,"('|-------------------------------------------------------')")
       write(*,"('| Timings:                                              ')")
       write(*,"('|-------------------------------------------------------')")
       if(read_time) then
          write(*,"('| Reading the Eigenvalue problem :',F13.3,' s')")&
                walltime_read_evp  
       endif
       if(write_time) then
          write(*,"('| Writing the Eigenvalue problem :',F13.3,' s')")&
                walltime_write_evp
       endif
       if(dist_pexsi_time) then
          write(*,"('| Distributing to PEXSI          :',F13.3,' s')")&
                walltime_dist_pexsi
       endif
       if(solve_time) then
          write(*,"('| Solving the Eigenvalue problem :',F13.3,' s')")&
                walltime_solve_evp
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
!!  This routine gets the current wallclock time
!!
subroutine elsi_get_time(wtime)

   implicit none
   !> The Wallclocktime in seconds
   real*8, intent(out) :: wtime
   
   integer :: tics

   call system_clock(tics)

   wtime = 1d0 * tics / clock_rate

end subroutine

!>
!!  This routine starts the total timer
!!
subroutine elsi_start_total_time()

   implicit none
  
   total_time = .true. 
   call elsi_get_time(walltime_total_start)

end subroutine

!>
!!  This routine stops the total timer
!!
subroutine elsi_stop_total_time()

   implicit none

   real*8 :: stop_time
   
   call elsi_get_time(stop_time)
   walltime_total = walltime_total &
     + stop_time - walltime_total_start

end subroutine

!>
!!  This routine starts the read evp timer
!!
subroutine elsi_start_read_evp_time()

   implicit none
   
   read_time = .true. 
   call elsi_get_time(walltime_read_evp_start)

end subroutine

!>
!!  This routine stops the read evp timer
!!
subroutine elsi_stop_read_evp_time()

   implicit none

   real*8 :: stop_time
   
   call elsi_get_time(stop_time)
   walltime_read_evp = walltime_read_evp &
     + stop_time - walltime_read_evp_start

end subroutine

!>
!!  This routine starts the read evp timer
!!
subroutine elsi_start_write_evp_time()

   implicit none
   
   write_time = .true. 
   call elsi_get_time(walltime_write_evp_start)

end subroutine

!>
!!  This routine stops the read evp timer
!!
subroutine elsi_stop_write_evp_time()

   implicit none

   real*8 :: stop_time
   
   call elsi_get_time(stop_time)
   walltime_write_evp = walltime_write_evp &
     + stop_time - walltime_write_evp_start

end subroutine

!>
!!  This routine starts the read evp timer
!!
subroutine elsi_start_solve_evp_time()

   implicit none
   
   solve_time = .true. 
   call elsi_get_time(walltime_solve_evp_start)

end subroutine

!>
!!  This routine stops the read evp timer
!!
subroutine elsi_stop_solve_evp_time()

   implicit none

   real*8 :: stop_time
   
   call elsi_get_time(stop_time)
   walltime_solve_evp = walltime_solve_evp &
     + stop_time - walltime_solve_evp_start

end subroutine

!>
!!  This routine starts the dist pexsi timer
!!
subroutine elsi_start_dist_pexsi_time()

   implicit none
   
   dist_pexsi_time = .true. 
   call elsi_get_time(walltime_dist_pexsi_start)

end subroutine

!>
!!  This routine stops the dist pexsi timer
!!
subroutine elsi_stop_dist_pexsi_time()

   implicit none

   real*8 :: stop_time
   
   call elsi_get_time(stop_time)
   walltime_dist_pexsi = walltime_dist_pexsi &
     + stop_time - walltime_dist_pexsi_start

end subroutine

end module ELSI_TIMERS
