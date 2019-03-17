! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Define precision constants used across ELSI.
!!
module ELSI_PRECISION

   use, intrinsic :: ISO_C_BINDING, only: c_float,c_double,c_int32_t,c_int64_t

   implicit none

   integer, parameter :: r4 = c_float
   integer, parameter :: r8 = c_double
   integer, parameter :: i4 = c_int32_t
   integer, parameter :: i8 = c_int64_t

end module ELSI_PRECISION
