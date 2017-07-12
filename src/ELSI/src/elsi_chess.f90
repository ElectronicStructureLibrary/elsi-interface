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

   use ISO_C_BINDING
   use ELSI_CONSTANTS, only: REAL_VALUES
   use ELSI_DIMENSIONS, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4
   use ELSI_TIMERS
   use ELSI_UTILS
   use SPARSEMATRIX_BASE
   use SPARSEMATRIX_HIGHLEVEL, only: matrix_fermi_operator_expansion
   use FOE_BASE, only: foe_data


   implicit none
   private

   public :: elsi_solve_evp_chess
   public :: elsi_compute_edm_chess
   public :: elsi_set_chess_default
   public :: elsi_print_chess_options

contains

!>
!! This routine interfaces to CheSS.
!!
subroutine elsi_solve_evp_chess(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_solve_evp_chess"

   call MPI_Barrier(elsi_h%mpi_comm,mpierr)
   call elsi_stop_density_matrix_time(elsi_h)

end subroutine

!> 
!! This routine computes the energy-weighted density matrix.
!! 
subroutine elsi_compute_edm_chess(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_compute_edm_chess"

end subroutine

!> 
!! This routine sets default CheSS parameters.
!! 
subroutine elsi_set_chess_default(elsi_h)
   
   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_set_chess_default"

end subroutine

!>
!! This routine prints CheSS settings.
!!          
subroutine elsi_print_chess_options(elsi_h)

   implicit none

   type(elsi_handle), intent(in) :: elsi_h !< Handle

   character*200 :: info_str

   character*40, parameter :: caller = "elsi_print_chess_options"

   write(info_str,"(A)") "  CheSS settings:"
   call elsi_statement_print(info_str,elsi_h)

end subroutine

end module ELSI_CHESS
