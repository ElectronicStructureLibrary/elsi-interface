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

   use elsi_precision, only : dp
   use iso_c_binding
   use ELSI_DIMENSIONS
   use ELSI_UTILS

   implicit none
   private

   real(kind=dp) :: t_generalized_evp
   real(kind=dp) :: t_generalized_evp_start
   real(kind=dp) :: t_redistribution
   real(kind=dp) :: t_redistribution_start
   real(kind=dp) :: t_transform_evp
   real(kind=dp) :: t_transform_evp_start
   real(kind=dp) :: t_back_transform_ev
   real(kind=dp) :: t_back_transform_ev_start
   real(kind=dp) :: t_singularity_check
   real(kind=dp) :: t_singularity_check_start
   real(kind=dp) :: t_standard_evp
   real(kind=dp) :: t_standard_evp_start
   real(kind=dp) :: t_density_matrix
   real(kind=dp) :: t_density_matrix_start
   real(kind=dp) :: t_cholesky
   real(kind=dp) :: t_cholesky_start

   integer :: clock_rate
   integer :: clock_max

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

contains

!>
!! This routine sets all timers to zero.
!!
subroutine elsi_init_timers()

   implicit none

   integer :: initial_time

   t_generalized_evp         = 0.0_dp
   t_generalized_evp_start   = 0.0_dp
   t_redistribution          = 0.0_dp
   t_redistribution_start    = 0.0_dp
   t_transform_evp           = 0.0_dp
   t_transform_evp_start     = 0.0_dp
   t_back_transform_ev       = 0.0_dp
   t_back_transform_ev_start = 0.0_dp
   t_singularity_check       = 0.0_dp
   t_singularity_check_start = 0.0_dp
   t_standard_evp            = 0.0_dp
   t_standard_evp_start      = 0.0_dp
   t_density_matrix          = 0.0_dp
   t_density_matrix_start    = 0.0_dp
   t_cholesky                = 0.0_dp
   t_cholesky_start          = 0.0_dp

   call system_clock(initial_time,clock_rate,clock_max)

end subroutine

!>
!! This routine prints a final output.
!!
subroutine elsi_final_print()

   implicit none

   real(kind=dp)  :: sparsity
   integer :: i_proc

   if(print_info) then
      if(myid == 0) then
         write(*,"('  |-----------------------------------------')")
         write(*,"('  | Final ELSI Output:')")
         write(*,"('  |-----------------------------------------')")
         write(*,"('  | Eigenvalue problem size : ',I13)") n_g_size
         if(method == PEXSI) then
            write(*,"('  | Non zero elements       : ',I13)") nnz_g
            sparsity = 1.0_dp-(1.0_dp*nnz_g/n_g_size)/n_g_size
            write(*,"('  | Sparsity                : ',F13.3)") sparsity
         endif
         write(*,"('  | Number of electrons     : ',F13.1)") n_electrons
         write(*,"('  | Number of states        : ',I13)") n_states
         if(method == ELPA) then
            write(*,"('  | Method                  : ',A13)") "ELPA"
         elseif(method == LIBOMM) then
            write(*,"('  | Method                  : ',A13)") "libOMM"
         elseif(method == PEXSI) then
            write(*,"('  | Method                  : ',A13)") "PEXSI"
         elseif(method == CHESS) then
            write(*,"('  | Method                  : ',A13)") "CheSS"
         elseif(method == SIPS) then
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
   real(kind=dp), intent(out) :: wtime
 
   integer :: tics

   call system_clock(tics)

   wtime = 1.0_dp*tics/clock_rate

end subroutine

!>
!! This routine starts generalized_evp timer.
!!
subroutine elsi_start_generalized_evp_time()

   implicit none
   
   call elsi_get_time(t_generalized_evp_start)

end subroutine

!>
!! This routine ends generalized_evp timer.
!!
subroutine elsi_stop_generalized_evp_time()

   implicit none

   real(kind=dp) :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_generalized_evp = stop_time-t_generalized_evp_start

   if(method == SIPS) then
      write(info_str,"('  Finished solving generalized eigenproblem')")
      call elsi_statement_print(info_str)
      write(info_str,"('  | Time :',F10.3,' s')") t_generalized_evp
      call elsi_statement_print(info_str)
   endif

end subroutine

!>
!! This routine starts density_matrix timer.
!!
subroutine elsi_start_density_matrix_time()

   implicit none

   call elsi_get_time(t_density_matrix_start)

end subroutine

!>
!! This routine ends density_matrix timer.
!!
subroutine elsi_stop_density_matrix_time()

   implicit none

   real(kind=dp) :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_density_matrix = stop_time-t_density_matrix_start

   write(info_str,"('  Finished density matrix calculation')")
   call elsi_statement_print(info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t_density_matrix
   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine starts redistribution timer.
!!
subroutine elsi_start_redistribution_time()

   implicit none
   
   call elsi_get_time(t_redistribution_start)

end subroutine

!>
!! This routine ends redistribution timer.
!!
subroutine elsi_stop_redistribution_time()

   implicit none

   real(kind=dp) :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_redistribution = stop_time-t_redistribution_start

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_statement_print(info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t_redistribution
   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine starts transform_evp timer.
!!
subroutine elsi_start_transform_evp_time()

   implicit none

   call elsi_get_time(t_transform_evp_start)

end subroutine

!>
!! This routine ends transform_evp timer.
!!
subroutine elsi_stop_transform_evp_time()

   implicit none

   real(kind=dp) :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_transform_evp = stop_time-t_transform_evp_start

   write(info_str,"('  Finished transformation to standard eigenproblem')")
   call elsi_statement_print(info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t_transform_evp
   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine starts back_transform_ev timer.
!!
subroutine elsi_start_back_transform_ev_time()

   implicit none

   call elsi_get_time(t_back_transform_ev_start)

end subroutine

!>
!! This routine ends back_transform_ev timer.
!!
subroutine elsi_stop_back_transform_ev_time()

   implicit none

   real(kind=dp) :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_back_transform_ev = stop_time-t_back_transform_ev_start

   write(info_str,"('  Finished back-transformation of eigenvectors')")
   call elsi_statement_print(info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t_back_transform_ev
   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine starts singularity_check timer.
!!
subroutine elsi_start_singularity_check_time()

   implicit none

   call elsi_get_time(t_singularity_check_start)

end subroutine

!>
!! This routine ends singularity_check timer.
!!
subroutine elsi_stop_singularity_check_time()

   implicit none

   real(kind=dp) :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_singularity_check = stop_time-t_singularity_check_start

   write(info_str,"('  Finished singularity check of overlap matrix')")
   call elsi_statement_print(info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t_singularity_check
   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine starts standard_evp timer.
!!
subroutine elsi_start_standard_evp_time()

   implicit none

   call elsi_get_time(t_standard_evp_start)

end subroutine

!>
!! This routine ends standard_evp timer.
!!
subroutine elsi_stop_standard_evp_time()

   implicit none

   real(kind=dp) :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_standard_evp = stop_time-t_standard_evp_start

   write(info_str,"('  Finished solving standard eigenproblem')")
   call elsi_statement_print(info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t_standard_evp
   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine starts cholesky timer.
!!
subroutine elsi_start_cholesky_time()

   implicit none

   call elsi_get_time(t_cholesky_start)

end subroutine

!>
!! This routine ends cholesky timer.
!!
subroutine elsi_stop_cholesky_time()

   implicit none

   real(kind=dp) :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_cholesky = stop_time-t_cholesky_start

   write(info_str,"('  Finished Cholesky decomposition')")
   call elsi_statement_print(info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t_cholesky
   call elsi_statement_print(info_str)

end subroutine

end module ELSI_TIMERS
