! Copyright (c) 2015-2018, the ELSI team. All rights reserved.
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
!! This module contains constants used in ELSI.
!!
module ELSI_CONSTANTS

   use ELSI_PRECISION, only: r8,i4

   implicit none

   character(len=8), parameter :: RELEASE_DATE = "20170527"

   real(kind=r8), parameter :: SQRT_PI        = 1.7724538509055160273_r8
   real(kind=r8), parameter :: INVERT_SQRT_PI = 0.5641895835477563275_r8

   integer(kind=i4), parameter :: UNSET            = -910910
   integer(kind=i4), parameter :: N_SOLVERS        = 7
   integer(kind=i4), parameter :: N_MATRIX_FORMATS = 2
   integer(kind=i4), parameter :: N_PARALLEL_MODES = 2

   ! Method names
   integer(kind=i4), parameter :: AUTO         = 0
   integer(kind=i4), parameter :: ELPA_SOLVER  = 1
   integer(kind=i4), parameter :: OMM_SOLVER   = 2
   integer(kind=i4), parameter :: PEXSI_SOLVER = 3
   integer(kind=i4), parameter :: CHESS_SOLVER = 4
   integer(kind=i4), parameter :: SIPS_SOLVER  = 5
   integer(kind=i4), parameter :: DMP_SOLVER   = 6

   ! Real or complex data
   integer(kind=i4), parameter :: REAL_VALUES    = 0
   integer(kind=i4), parameter :: COMPLEX_VALUES = 1

   ! Output for solver subroutines (density matrix or eigenvectors)
   integer(kind=i4), parameter :: OUTPUT_EV = 0
   integer(kind=i4), parameter :: OUTPUT_DM = 1

   ! Matrix formats
   integer(kind=i4), parameter :: BLACS_DENSE = 0
   integer(kind=i4), parameter :: PEXSI_CSC   = 1

   ! Triangular matrix
   integer(kind=i4), parameter :: FULL_MAT = 0
   integer(kind=i4), parameter :: UT_MAT   = 1
   integer(kind=i4), parameter :: LT_MAT   = 2

   ! Parallel modes
   integer(kind=i4), parameter :: SINGLE_PROC = 0
   integer(kind=i4), parameter :: MULTI_PROC  = 1

   ! Broadening schemes
   integer(kind=i4), parameter :: GAUSSIAN            = 0
   integer(kind=i4), parameter :: FERMI               = 1
   integer(kind=i4), parameter :: METHFESSEL_PAXTON_0 = 2
   integer(kind=i4), parameter :: METHFESSEL_PAXTON_1 = 3
   integer(kind=i4), parameter :: CUBIC               = 4

   ! Density matrix purification methods
   integer(kind=i4), parameter :: TRACE_CORRECTING = 0
   integer(kind=i4), parameter :: CANONICAL        = 1

   ! Matrix reading and writing
   integer(kind=i4), parameter :: HEADER_SIZE  = 16
   integer(kind=i4), parameter :: FILE_VERSION = 170915

   ! Reading and writing tasks
   integer(kind=i4), parameter :: READ_FILE  = 0
   integer(kind=i4), parameter :: WRITE_FILE = 1

   ! Constants for ELSI file IO
   character(len=*), parameter :: UNSET_STRING  = "N/A"
   integer(kind=i4), parameter :: FILE_NAME_LEN = 80
   integer(kind=i4), parameter :: HUMAN_READ    = 0
   integer(kind=i4), parameter :: JSON          = 1

   ! Control placement of commas in JSON records
   integer(kind=i4), parameter :: COMMA_AFTER  = 0
   integer(kind=i4), parameter :: COMMA_BEFORE = 1
   integer(kind=i4), parameter :: NO_COMMA     = 2

   ! Constants for ELSI timings
   integer(kind=i4),             parameter :: TIMING_STRING_LEN           = 20
   integer(kind=i4),             parameter :: SOLVER_TIMINGS_UNIT_DEFAULT = 66
   character(len=FILE_NAME_LEN), parameter :: SOLVER_TIMINGS_FILE_DEFAULT = "elsi_solver_timings.json"
   integer(kind=i4),             parameter :: DATETIME_LEN                = 29

end module ELSI_CONSTANTS
