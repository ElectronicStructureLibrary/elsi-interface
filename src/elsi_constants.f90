! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains constants used in ELSI.
!!
module ELSI_CONSTANTS

   use ELSI_PRECISION, only: r8,i4

   implicit none

   real(kind=r8), parameter :: SQRT_PI        = 1.7724538509055160273_r8
   real(kind=r8), parameter :: INVERT_SQRT_PI = 0.5641895835477563275_r8

   integer(kind=i4), parameter :: UNSET            = -910910
   integer(kind=i4), parameter :: N_SOLVERS        = 7
   integer(kind=i4), parameter :: N_MATRIX_FORMATS = 3
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
   integer(kind=i4), parameter :: REAL_DATA  = 0
   integer(kind=i4), parameter :: CMPLX_DATA = 1

   ! Matrix formats
   integer(kind=i4), parameter :: BLACS_DENSE = 0
   integer(kind=i4), parameter :: PEXSI_CSC   = 1
   integer(kind=i4), parameter :: SIESTA_CSC  = 2

   ! Triangular matrix
   integer(kind=i4), parameter :: FULL_MAT = 0
   integer(kind=i4), parameter :: UT_MAT   = 1
   integer(kind=i4), parameter :: LT_MAT   = 2

   ! Parallel modes
   integer(kind=i4), parameter :: SINGLE_PROC = 0
   integer(kind=i4), parameter :: MULTI_PROC  = 1

   ! Broadening schemes
   integer(kind=i4), parameter :: GAUSSIAN          = 0
   integer(kind=i4), parameter :: FERMI             = 1
   integer(kind=i4), parameter :: METHFESSEL_PAXTON = 2
   integer(kind=i4), parameter :: CUBIC             = 3
   integer(kind=i4), parameter :: COLD              = 4

   ! Density matrix purification methods
   integer(kind=i4), parameter :: TRACE_CORRECTING = 0
   integer(kind=i4), parameter :: CANONICAL        = 1

   ! Matrix reading and writing
   integer(kind=i4), parameter :: HEADER_SIZE  = 16
   integer(kind=i4), parameter :: FILE_VERSION = 170915
   integer(kind=i4), parameter :: READ_FILE  = 0
   integer(kind=i4), parameter :: WRITE_FILE = 1

   ! IO
   integer(kind=i4), parameter :: STR_LEN      = 20
   integer(kind=i4), parameter :: DATETIME_LEN = 29
   integer(kind=i4), parameter :: UUID_LEN     = 36 ! RFC 4122 format

end module ELSI_CONSTANTS
