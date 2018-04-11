! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the root directory.

!> Definitions for precisions for primitive values (ISO_C_BINDING used for
!! portability)
module FJSON_PRECISION

   use, intrinsic :: ISO_C_BINDING, only: c_float, c_double, c_int32_t, &
                                          c_int64_t

   implicit none

   !> Precision for single-precision floating point numbers
   integer, parameter :: r4 = c_float

   !> Precision for double-precision floating point numbers
   integer, parameter :: r8 = c_double

   !> Precision for integer numbers
   integer, parameter :: i4 = c_int32_t

   !> Precision for long integer numbers
   integer, parameter :: i8 = c_int64_t

end module FJSON_PRECISION
