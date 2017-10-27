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
!! This module contains constants used in ELSI.
!!
module ELSI_CONSTANTS

   use ELSI_PRECISION, only: r8,i4

   implicit none

   real(kind=r8), parameter :: INVERT_SQRT_PI = 0.5641895835477563275_r8 ! 1/sqrt(pi)

   integer(kind=i4), parameter :: UNSET                = -910910
   integer(kind=i4), parameter :: N_SOLVERS            = 6
   integer(kind=i4), parameter :: N_MATRIX_FORMATS     = 2
   integer(kind=i4), parameter :: N_PARALLEL_MODES     = 2
   integer(kind=i4), parameter :: N_BROADENING_SCHEMES = 4

   ! Method names
   integer(kind=i4), parameter :: AUTO   = 0
   integer(kind=i4), parameter :: ELPA_SOLVER  = 1
   integer(kind=i4), parameter :: OMM_SOLVER   = 2
   integer(kind=i4), parameter :: PEXSI_SOLVER = 3
   integer(kind=i4), parameter :: CHESS_SOLVER = 4
   integer(kind=i4), parameter :: SIPS_SOLVER  = 5

   ! Real or complex data
   integer(kind=i4), parameter :: REAL_VALUES    = 0
   integer(kind=i4), parameter :: COMPLEX_VALUES = 1

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

   ! Matrix reading and writing
   integer(kind=i4), parameter :: HEADER_SIZE  = 16
   integer(kind=i4), parameter :: FILE_VERSION = 170915

   ! Reading and writing
   integer(kind=i4), parameter :: READ_FILE  = 0
   integer(kind=i4), parameter :: WRITE_FILE = 1

   ! File formats
   integer(kind=i4), parameter :: DENSE_FILE = 0
   integer(kind=i4), parameter :: CSC_FILE   = 1

end module ELSI_CONSTANTS
