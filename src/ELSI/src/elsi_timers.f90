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
!! This module provides timers for ELSI.
!!
module ELSI_TIMERS

   use iso_c_binding
   use ELSI_CONSTANTS, only: ELPA,LIBOMM,PEXSI,CHESS,SIPS
   use ELSI_DIMENSIONS, only: elsi_handle,print_info
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS

   implicit none
   private

   integer(kind=i4) :: clock_rate
   integer(kind=i4) :: clock_max

   public :: elsi_init_timers
   public :: elsi_final_print
   public :: elsi_start_generalized_evp_time
   public :: elsi_stop_generalized_evp_time
   public :: elsi_start_redistribution_time
   public :: elsi_stop_redistribution_time
   public :: elsi_start_transform_evp_time
   public :: elsi_stop_transform_evp_time
   public :: elsi_start_back_transform_ev_time
   public :: elsi_stop_back_transform_ev_time
   public :: elsi_start_singularity_check_time
   public :: elsi_stop_singularity_check_time
   public :: elsi_start_standard_evp_time
   public :: elsi_stop_standard_evp_time
   public :: elsi_start_density_matrix_time
   public :: elsi_stop_density_matrix_time
   public :: elsi_start_cholesky_time
   public :: elsi_stop_cholesky_time
   public :: elsi_start_inertia_time
   public :: elsi_stop_inertia_time

contains

!>
!! This routine sets all timers to zero.
!!
subroutine elsi_init_timers(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   integer(kind=i4)        :: initial_time
   character*40, parameter :: caller = "elsi_init_timers"

   elsi_h%t_generalized_evp         = 0.0_r8
   elsi_h%t_generalized_evp_start   = 0.0_r8
   elsi_h%t_redistribution          = 0.0_r8
   elsi_h%t_redistribution_start    = 0.0_r8
   elsi_h%t_transform_evp           = 0.0_r8
   elsi_h%t_transform_evp_start     = 0.0_r8
   elsi_h%t_back_transform_ev       = 0.0_r8
   elsi_h%t_back_transform_ev_start = 0.0_r8
   elsi_h%t_singularity_check       = 0.0_r8
   elsi_h%t_singularity_check_start = 0.0_r8
   elsi_h%t_standard_evp            = 0.0_r8
   elsi_h%t_standard_evp_start      = 0.0_r8
   elsi_h%t_density_matrix          = 0.0_r8
   elsi_h%t_density_matrix_start    = 0.0_r8
   elsi_h%t_cholesky                = 0.0_r8
   elsi_h%t_cholesky_start          = 0.0_r8
   elsi_h%t_inertia                 = 0.0_r8
   elsi_h%t_inertia_start           = 0.0_r8

   call system_clock(initial_time,clock_rate,clock_max)

end subroutine

!>
!! This routine prints a final output.
!!
subroutine elsi_final_print(elsi_h)

   implicit none

   type(elsi_handle), intent(in) :: elsi_h

   real(kind=r8)           :: sparsity
   integer(kind=i4)        :: i_proc
   character*40, parameter :: caller = "elsi_final_print"

   if(print_info) then
      if(elsi_h%myid == 0) then
         write(*,"('  |-----------------------------------------')")
         write(*,"('  | Final ELSI Output:')")
         write(*,"('  |-----------------------------------------')")
         write(*,"('  | Eigenvalue problem size : ',I13)") elsi_h%n_g_size
         if((elsi_h%solver == PEXSI) .or. (elsi_h%solver == SIPS)) then
            write(*,"('  | Non zero elements       : ',I13)") elsi_h%nnz_g
            sparsity = 1.0_r8-(1.0_r8*elsi_h%nnz_g/elsi_h%n_g_size)/elsi_h%n_g_size
            write(*,"('  | Sparsity                : ',F13.3)") sparsity
         endif
         write(*,"('  | Number of electrons     : ',F13.1)") elsi_h%n_electrons
         write(*,"('  | Number of states        : ',I13)") elsi_h%n_states
         if(elsi_h%solver == ELPA) then
            write(*,"('  | Method                  : ',A13)") "ELPA"
         elseif(elsi_h%solver == LIBOMM) then
            write(*,"('  | Method                  : ',A13)") "libOMM"
         elseif(elsi_h%solver == PEXSI) then
            write(*,"('  | Method                  : ',A13)") "PEXSI"
         elseif(elsi_h%solver == CHESS) then
            write(*,"('  | Method                  : ',A13)") "CheSS"
         elseif(elsi_h%solver == SIPS) then
            write(*,"('  | Method                  : ',A13)") "SIPs"
         endif
         write(*,"('  |-----------------------------------------')")
         write(*,"('  | ELSI Project (c)  elsi-interchange.org')")
         write(*,"('  |-----------------------------------------')")
      endif
   endif

end subroutine

!>
!! This routine gets the current wallclock time.
!!
subroutine elsi_get_time(wtime)

   implicit none

   real(kind=r8), intent(out) :: wtime

   integer(kind=i4)        :: tics
   character*40, parameter :: caller = "elsi_get_time"

   call system_clock(tics)

   wtime = 1.0_r8*tics/clock_rate

end subroutine

!>
!! This routine starts the generalized_evp timer.
!!
subroutine elsi_start_generalized_evp_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   character*40, parameter :: caller = "elsi_start_generalized_evp_time"

   call elsi_get_time(elsi_h%t_generalized_evp_start)

end subroutine

!>
!! This routine ends the generalized_evp timer.
!!
subroutine elsi_stop_generalized_evp_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   real(kind=r8)           :: stop_time
   character*200           :: info_str
   character*40, parameter :: caller = "elsi_stop_generalized_evp_time"

   call elsi_get_time(stop_time)
   elsi_h%t_generalized_evp = stop_time-elsi_h%t_generalized_evp_start

   write(info_str,"('  Finished solving generalized eigenproblem')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") elsi_h%t_generalized_evp
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine starts the density_matrix timer.
!!
subroutine elsi_start_density_matrix_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   character*40, parameter :: caller = "elsi_start_density_matrix_time"

   call elsi_get_time(elsi_h%t_density_matrix_start)

end subroutine

!>
!! This routine ends the density_matrix timer.
!!
subroutine elsi_stop_density_matrix_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   real(kind=r8)           :: stop_time
   character*200           :: info_str
   character*40, parameter :: caller = "elsi_stop_density_matrix_time"

   call elsi_get_time(stop_time)
   elsi_h%t_density_matrix = stop_time-elsi_h%t_density_matrix_start

   write(info_str,"('  Finished density matrix calculation')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") elsi_h%t_density_matrix
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine starts the redistribution timer.
!!
subroutine elsi_start_redistribution_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   character*40, parameter :: caller = "elsi_start_redistribution_time"
   
   call elsi_get_time(elsi_h%t_redistribution_start)

end subroutine

!>
!! This routine ends the redistribution timer.
!!
subroutine elsi_stop_redistribution_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   real(kind=r8)           :: stop_time
   character*200           :: info_str
   character*40, parameter :: caller = "elsi_stop_redistribution_time"

   call elsi_get_time(stop_time)
   elsi_h%t_redistribution = stop_time-elsi_h%t_redistribution_start

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") elsi_h%t_redistribution
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine starts the transform_evp timer.
!!
subroutine elsi_start_transform_evp_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   character*40, parameter :: caller = "elsi_start_transform_evp_time"

   call elsi_get_time(elsi_h%t_transform_evp_start)

end subroutine

!>
!! This routine ends the transform_evp timer.
!!
subroutine elsi_stop_transform_evp_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   real(kind=r8)           :: stop_time
   character*200           :: info_str
   character*40, parameter :: caller = "elsi_stop_transform_evp_time"

   call elsi_get_time(stop_time)
   elsi_h%t_transform_evp = stop_time-elsi_h%t_transform_evp_start

   write(info_str,"('  Finished transformation to standard eigenproblem')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") elsi_h%t_transform_evp
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine starts the back_transform_ev timer.
!!
subroutine elsi_start_back_transform_ev_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   character*40, parameter :: caller = "elsi_start_back_transform_ev_time"

   call elsi_get_time(elsi_h%t_back_transform_ev_start)

end subroutine

!>
!! This routine ends the back_transform_ev timer.
!!
subroutine elsi_stop_back_transform_ev_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   real(kind=r8)           :: stop_time
   character*200           :: info_str
   character*40, parameter :: caller = "elsi_stop_back_transform_ev_time"

   call elsi_get_time(stop_time)
   elsi_h%t_back_transform_ev = stop_time-elsi_h%t_back_transform_ev_start

   write(info_str,"('  Finished back-transformation of eigenvectors')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") elsi_h%t_back_transform_ev
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine starts the singularity_check timer.
!!
subroutine elsi_start_singularity_check_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   character*40, parameter :: caller = "elsi_start_singularity_check_time"

   call elsi_get_time(elsi_h%t_singularity_check_start)

end subroutine

!>
!! This routine ends the singularity_check timer.
!!
subroutine elsi_stop_singularity_check_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   real(kind=r8)           :: stop_time
   character*200           :: info_str
   character*40, parameter :: caller = "elsi_stop_singularity_check_time"

   call elsi_get_time(stop_time)
   elsi_h%t_singularity_check = stop_time-elsi_h%t_singularity_check_start

   write(info_str,"('  Finished singularity check of overlap matrix')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") elsi_h%t_singularity_check
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine starts the standard_evp timer.
!!
subroutine elsi_start_standard_evp_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   character*40, parameter :: caller = "elsi_start_standard_evp_time"

   call elsi_get_time(elsi_h%t_standard_evp_start)

end subroutine

!>
!! This routine ends the standard_evp timer.
!!
subroutine elsi_stop_standard_evp_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   real(kind=r8)           :: stop_time
   character*200           :: info_str
   character*40, parameter :: caller = "elsi_stop_standard_evp_time"

   call elsi_get_time(stop_time)
   elsi_h%t_standard_evp = stop_time-elsi_h%t_standard_evp_start

   write(info_str,"('  Finished solving standard eigenproblem')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") elsi_h%t_standard_evp
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine starts the cholesky timer.
!!
subroutine elsi_start_cholesky_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   character*40, parameter :: caller = "elsi_start_cholesky_time"

   call elsi_get_time(elsi_h%t_cholesky_start)

end subroutine

!>
!! This routine ends the cholesky timer.
!!
subroutine elsi_stop_cholesky_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   real(kind=r8)           :: stop_time
   character*200           :: info_str
   character*40, parameter :: caller = "elsi_stop_cholesky_time"

   call elsi_get_time(stop_time)
   elsi_h%t_cholesky = stop_time-elsi_h%t_cholesky_start

   write(info_str,"('  Finished Cholesky decomposition')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") elsi_h%t_cholesky
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine starts the inertia timer.
!!
subroutine elsi_start_inertia_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   character*40, parameter :: caller = "elsi_start_inertia_time"

   call elsi_get_time(elsi_h%t_inertia_start)

end subroutine

!>
!! This routine ends the inertia timer.
!!
subroutine elsi_stop_inertia_time(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   real(kind=r8)           :: stop_time
   character*200           :: info_str
   character*40, parameter :: caller = "elsi_stop_inertia_time"

   call elsi_get_time(stop_time)
   elsi_h%t_inertia = stop_time-elsi_h%t_inertia_start

   write(info_str,"('  Finished inertia counting')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") elsi_h%t_inertia
   call elsi_statement_print(info_str,elsi_h)

end subroutine

end module ELSI_TIMERS
