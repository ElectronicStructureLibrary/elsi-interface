! Copyright (c) 2015-2020, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Contain constants used in ELSI.
!!
module ELSI_CONSTANT

   use ELSI_PRECISION, only: r8,i4

   implicit none

   real(kind=r8), parameter :: SQRT_PI = 1.7724538509055160273_r8
   real(kind=r8), parameter :: INVERT_SQRT_PI = 0.5641895835477563275_r8

   integer(kind=i4), parameter :: UNSET = -910910
   integer(kind=i4), parameter :: N_SOLVERS = 9
   integer(kind=i4), parameter :: N_MATRIX_FORMATS = 4
   integer(kind=i4), parameter :: N_PARALLEL_MODES = 2

   ! Solver
   integer(kind=i4), parameter :: AUTO_SOLVER = 0
   integer(kind=i4), parameter :: ELPA_SOLVER = 1
   integer(kind=i4), parameter :: OMM_SOLVER = 2
   integer(kind=i4), parameter :: PEXSI_SOLVER = 3
   integer(kind=i4), parameter :: EIGENEXA_SOLVER = 4
   integer(kind=i4), parameter :: SIPS_SOLVER = 5
   integer(kind=i4), parameter :: NTPOLY_SOLVER = 6
   integer(kind=i4), parameter :: MAGMA_SOLVER = 7
   integer(kind=i4), parameter :: BSEPACK_SOLVER = 8

   ! Real or complex data
   integer(kind=i4), parameter :: REAL_DATA = 0
   integer(kind=i4), parameter :: CMPLX_DATA = 1

   ! Matrix format
   integer(kind=i4), parameter :: BLACS_DENSE = 0
   integer(kind=i4), parameter :: PEXSI_CSC = 1
   integer(kind=i4), parameter :: SIESTA_CSC = 2
   integer(kind=i4), parameter :: GENERIC_COO = 3

   ! Matrix conversion mask
   integer(kind=i4), parameter :: MASK_HS = 0
   integer(kind=i4), parameter :: MASK_H = 1
   integer(kind=i4), parameter :: MASK_S = 2

   ! Triangular matrix
   integer(kind=i4), parameter :: FULL_MAT = 0
   integer(kind=i4), parameter :: UT_MAT = 1
   integer(kind=i4), parameter :: LT_MAT = 2

   ! Parallelization mode
   integer(kind=i4), parameter :: SINGLE_PROC = 0
   integer(kind=i4), parameter :: MULTI_PROC = 1

   ! Broadening
   integer(kind=i4), parameter :: GAUSSIAN = 0
   integer(kind=i4), parameter :: FERMI = 1
   integer(kind=i4), parameter :: METHFESSEL_PAXTON = 2
   integer(kind=i4), parameter :: CUBIC = 3
   integer(kind=i4), parameter :: COLD = 4

   ! Density matrix
   integer(kind=i4), parameter :: GET_DM = 0
   integer(kind=i4), parameter :: GET_EDM = 1
   integer(kind=i4), parameter :: GET_FDM = 2

   ! Density matrix purification
   integer(kind=i4), parameter :: NTPOLY_PM = 0
   integer(kind=i4), parameter :: NTPOLY_TRS2 = 1
   integer(kind=i4), parameter :: NTPOLY_TRS4 = 2
   integer(kind=i4), parameter :: NTPOLY_HPCP = 3

   ! Density matrix extrapolation
   integer(kind=i4), parameter :: EXTRA_FACTOR = 0
   integer(kind=i4), parameter :: EXTRA_TRS2 = 1

   ! Matrix reading and writing
   integer(kind=i4), parameter :: HEADER_SIZE = 16
   integer(kind=i4), parameter :: FILE_VERSION = 170915
   integer(kind=i4), parameter :: READ_FILE = 0
   integer(kind=i4), parameter :: WRITE_FILE = 1

end module ELSI_CONSTANT
