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
!! This module contains constant values and choices for runtime parameters to be reused throughout ELSI.
!!

module ELSI_CONSTANTS

   use ELSI_PRECISION, only : dp

   real(kind=dp), parameter :: INVERT_SQRT_PI = 0.564189583547756_dp !< Constant: 1/sqrt(pi)

! ========= ALIAS =========

   integer, parameter :: UNSET = -1

   integer, parameter :: N_SOLVERS = 6
   !> Solvers (AUTO=0, ELPA=1, LIBOMM=2, PEXSI=3, CHESS=4, SIPS=5)
   enum, bind( C )
      enumerator :: AUTO, ELPA, LIBOMM, PEXSI, CHESS, SIPS
   end enum

   integer, parameter :: N_MATRIX_DATA_TYPES = 2
   !> Matrix data types, i.e. real or complex (REAL_VALUES=0, COMPLEX_VALUES=1)
   enum, bind( C )
      enumerator :: REAL_VALUES, COMPLEX_VALUES
   end enum

   integer, parameter :: N_MATRIX_STORAGE_FORMATS = 2
   !> Matrix storage format (BLACS_DENSE=0, PEXSI_CSC=1)
   enum, bind( C )
      enumerator :: BLACS_DENSE, PEXSI_CSC
   end enum

   integer, parameter :: N_PARALLEL_MODES = 2
   !> Parallel modes (SINGLE_PROC=0, MULTI_PROC=1)
   enum, bind( C )
      enumerator :: SINGLE_PROC, MULTI_PROC
   end enum

   integer, parameter :: N_BROADENING_SCHEMES = 4
   !> Broadening schemes (used if ELPA is chosen to compute density matrix)
   enum, bind( C )
      enumerator :: GAUSSIAN, FERMI, METHFESSEL_PAXTON_0, METHFESSEL_PAXTON_1
   end enum

end module ELSI_CONSTANTS
