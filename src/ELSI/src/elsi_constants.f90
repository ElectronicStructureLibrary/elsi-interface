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

   use, intrinsic :: ISO_C_BINDING
   use ELSI_PRECISION, only: r8,i4

   implicit none

   real(kind=r8), parameter :: INVERT_SQRT_PI = 0.5641895835477563275_r8 !< Constant: 1/sqrt(pi)

   integer(kind=i4), parameter :: N_SOLVERS                = 6
   integer(kind=i4), parameter :: N_MATRIX_DATA_TYPES      = 2
   integer(kind=i4), parameter :: N_MATRIX_STORAGE_FORMATS = 2
   integer(kind=i4), parameter :: N_PARALLEL_MODES         = 2
   integer(kind=i4), parameter :: N_BROADENING_SCHEMES     = 4
   integer(kind=i4), parameter :: UNSET                    = -910910
   integer(kind=i4), parameter :: HEADER_SIZE              = 12

   ! Method names
   enum, bind(C)
      enumerator :: AUTO, ELPA, LIBOMM, PEXSI, CHESS, SIPS
   end enum

   ! Real or complex data
   enum, bind(C)
      enumerator :: REAL_VALUES, COMPLEX_VALUES
   end enum

   ! Storage formats
   enum, bind(C)
      enumerator :: BLACS_DENSE, PEXSI_CSC
   end enum

   ! Triangular matrix
   enum, bind(C)
      enumerator :: FULL_MAT, UT_MAT, LT_MAT
   end enum

   ! Parallel modes
   enum, bind(C)
      enumerator :: SINGLE_PROC, MULTI_PROC
   end enum

   ! Broadening schemes
   enum, bind(C)
      enumerator :: GAUSSIAN, FERMI, METHFESSEL_PAXTON_0, METHFESSEL_PAXTON_1
   end enum

end module ELSI_CONSTANTS
