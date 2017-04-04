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
   use ELSI_DIMENSIONS
   use ELSI_UTILS

   implicit none
   private

   real*8 :: t_generalized_evp
   real*8 :: t_generalized_evp_start
   real*8 :: t_blacs_to_pexsi
   real*8 :: t_blacs_to_pexsi_start
   real*8 :: t_pexsi_to_blacs
   real*8 :: t_pexsi_to_blacs_start
   real*8 :: t_blacs_to_sips
   real*8 :: t_blacs_to_sips_start
   real*8 :: t_sips_to_blacs
   real*8 :: t_sips_to_blacs_start
   real*8 :: t_transform_evp
   real*8 :: t_transform_evp_start
   real*8 :: t_back_transform_ev
   real*8 :: t_back_transform_ev_start
   real*8 :: t_singularity_check
   real*8 :: t_singularity_check_start
   real*8 :: t_standard_evp
   real*8 :: t_standard_evp_start
   real*8 :: t_density_matrix
   real*8 :: t_density_matrix_start
   real*8 :: t_cholesky
   real*8 :: t_cholesky_start

   integer :: clock_rate
   integer :: clock_max

   public :: elsi_init_timers
   public :: elsi_final_print
   public :: elsi_start_generalized_evp_time
   public :: elsi_stop_generalized_evp_time
   public :: elsi_start_blacs_to_pexsi_time
   public :: elsi_stop_blacs_to_pexsi_time
   public :: elsi_start_pexsi_to_blacs_time
   public :: elsi_stop_pexsi_to_blacs_time
   public :: elsi_start_blacs_to_sips_time
   public :: elsi_stop_blacs_to_sips_time
   public :: elsi_start_sips_to_blacs_time
   public :: elsi_stop_sips_to_blacs_time
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

   t_generalized_evp         = 0.0d0
   t_generalized_evp_start   = 0.0d0
   t_blacs_to_pexsi          = 0.0d0
   t_blacs_to_pexsi_start    = 0.0d0
   t_pexsi_to_blacs          = 0.0d0
   t_pexsi_to_blacs_start    = 0.0d0
   t_blacs_to_sips           = 0.0d0
   t_blacs_to_sips_start     = 0.0d0
   t_sips_to_blacs           = 0.0d0
   t_sips_to_blacs_start     = 0.0d0
   t_transform_evp           = 0.0d0
   t_transform_evp_start     = 0.0d0
   t_back_transform_ev       = 0.0d0
   t_back_transform_ev_start = 0.0d0
   t_singularity_check       = 0.0d0
   t_singularity_check_start = 0.0d0
   t_standard_evp            = 0.0d0
   t_standard_evp_start      = 0.0d0
   t_density_matrix          = 0.0d0
   t_density_matrix_start    = 0.0d0
   t_cholesky                = 0.0d0
   t_cholesky_start          = 0.0d0

   call system_clock(initial_time,clock_rate,clock_max)

end subroutine

!>
!! This routine prints a final output.
!!
subroutine elsi_final_print()

   implicit none

   integer :: i_proc

   if(print_info) then
      if(myid == 0) then
         write(*,"('  |-----------------------------------------')")
         write(*,"('  | Final ELSI Output:')")
         write(*,"('  |-----------------------------------------')")
         write(*,"('  | Eigenvalue problem size : ',I13)") n_g_size
         if(method == PEXSI) then
            write(*,"('  | Non zero elements       : ',I13)") &
                  nnz_g
            write(*,"('  | Sparsity                : ',F13.3)") &
                  1.0d0-(1.0d0*nnz_g/n_g_size)/n_g_size
         endif
         write(*,"('  | Number of electrons     : ',F13.1)") n_electrons
         write(*,"('  | Number of states        : ',I13)") n_states
         if(method == ELPA) then
            write(*,"('  | Method                  : ',A13)") "ELPA"
         endif
         if(method == LIBOMM) then
            write(*,"('  | Method                  : ',A13)") "libOMM"
         endif
         if(method == PEXSI) then
            write(*,"('  | Method                  : ',A13)") "PEXSI"
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
   real*8, intent(out) :: wtime
 
   integer :: tics

   call system_clock(tics)

   wtime = 1.0d0*tics/clock_rate

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

   real*8 :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_generalized_evp = stop_time-t_generalized_evp_start

   if(method == SIPS) then
      write(info_str,"('  | Solution to generalized eigenproblem (SIPs) :',F10.3,' s')")&
         t_generalized_evp
   endif

   call elsi_statement_print(info_str)

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

   real*8 :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_density_matrix = stop_time-t_density_matrix_start

   if(method == ELPA) then
      write(info_str,"('  | Density matrix calculation                  :',F10.3,' s')")&
         t_density_matrix
   elseif(method == LIBOMM) then
      write(info_str,"('  | Density matrix calculation (libOMM)         :',F10.3,' s')")&
         t_density_matrix
   elseif(method == PEXSI) then
      write(info_str,"('  | Density matrix calculation (PEXSI)          :',F10.3,' s')")&
         t_density_matrix
   elseif(method == CHESS) then
      write(info_str,"('  | Density matrix calculation (CheSS)          :',F10.3,' s')")&
         t_density_matrix
   endif

   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine starts blacs_to_pexsi timer.
!!
subroutine elsi_start_blacs_to_pexsi_time()

   implicit none
   
   call elsi_get_time(t_blacs_to_pexsi_start)

end subroutine

!>
!! This routine ends blacs_to_pexsi timer.
!!
subroutine elsi_stop_blacs_to_pexsi_time()

   implicit none

   real*8 :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_blacs_to_pexsi = stop_time-t_blacs_to_pexsi_start

   write(info_str,"('  | Matrix redistribution                       :',F10.3,' s')")&
      t_blacs_to_pexsi
   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine starts pexsi_to_blacs timer.
!!
subroutine elsi_start_pexsi_to_blacs_time()

   implicit none

   call elsi_get_time(t_pexsi_to_blacs_start)

end subroutine

!>
!! This routine ends pexsi_to_blacs timer.
!!
subroutine elsi_stop_pexsi_to_blacs_time()

   implicit none

   real*8 :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_pexsi_to_blacs = stop_time-t_pexsi_to_blacs_start

   write(info_str,"('  | Matrix redistribution                       :',F10.3,' s')")&
      t_pexsi_to_blacs
   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine starts blacs_to_sips timer.
!!
subroutine elsi_start_blacs_to_sips_time()

   implicit none

   call elsi_get_time(t_blacs_to_sips_start)

end subroutine

!>
!! This routine ends blacs_to_sips timer.
!!
subroutine elsi_stop_blacs_to_sips_time()

   implicit none

   real*8 :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_blacs_to_sips = stop_time-t_blacs_to_sips_start

   write(info_str,"('  | Matrix redistribution                       :',F10.3,' s')")&
      t_blacs_to_sips
   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine starts sips_to_blacs timer.
!!
subroutine elsi_start_sips_to_blacs_time()

   implicit none

   call elsi_get_time(t_sips_to_blacs_start)

end subroutine

!>
!! This routine ends sips_to_blacs timer.
!!
subroutine elsi_stop_sips_to_blacs_time()

   implicit none

   real*8 :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_sips_to_blacs = stop_time-t_sips_to_blacs_start

   write(info_str,"('  | Matrix redistribution                       :',F10.3,' s')")&
      t_sips_to_blacs
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

   real*8 :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_transform_evp = stop_time-t_transform_evp_start

   write(info_str,"('  | Transformation to standard eigenproblem     :',F10.3,' s')")&
      t_transform_evp
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

   real*8 :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_back_transform_ev = stop_time-t_back_transform_ev_start

   write(info_str,"('  | Transformation to original eigenvectors     :',F10.3,' s')")&
      t_back_transform_ev
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

   real*8 :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_singularity_check = stop_time-t_singularity_check_start

   write(info_str,"('  | Singularity check of overlap matrix         :',F10.3,' s')")&
      t_singularity_check
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

   real*8 :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_standard_evp = stop_time-t_standard_evp_start

   write(info_str,"('  | Solution to standard eigenproblem (ELPA)    :',F10.3,' s')")&
      t_standard_evp
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

   real*8 :: stop_time
   character*200 :: info_str

   call elsi_get_time(stop_time)
   t_cholesky = stop_time-t_cholesky_start

   write(info_str,"('  | Cholesky decomposition                      :',F10.3,' s')")&
      t_cholesky
   call elsi_statement_print(info_str)

end subroutine

end module ELSI_TIMERS
