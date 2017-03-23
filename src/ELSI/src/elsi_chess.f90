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
!! This module provides interfaces to CheSS.
!!

module ELSI_CHESS

   use ELSI_DIMENSIONS
   use ELSI_TIMERS
   use ELSI_UTILS

   implicit none
   private

   public :: elsi_solve_evp_chess
   public :: elsi_set_chess_default_options
   public :: elsi_print_chess_options

contains

!=========================
! ELSI routines for CheSS
!=========================

!>
!! This routine interfaces to CheSS.
!!
subroutine elsi_solve_evp_chess()

   implicit none
   include "mpif.h"

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_solve_evp_chess"

   call elsi_start_solve_evp_time()

   ! Solve the eigenvalue problem
   call elsi_statement_print("  Starting CheSS density matrix solver")

   call MPI_Barrier(mpi_comm_global,mpierr)
   call elsi_stop_solve_evp_time()

end subroutine

!>
!! Set CheSS variables to ELSI default.
!!
subroutine elsi_set_chess_default_options()

   implicit none

end subroutine

!>
!! Print CheSS settings.
!!
subroutine elsi_print_chess_options()

   implicit none

   character*200 :: info_str

end subroutine

end module ELSI_CHESS
