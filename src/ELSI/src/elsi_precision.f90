module ELSI_PRECISION

  use iso_c_binding, only: c_float, c_double, c_int32_t, c_int64_t

  implicit none

  integer, parameter :: sp = c_float
  integer, parameter :: dp = c_double
  integer, parameter :: i4 = c_int32_t
  integer, parameter :: i8 = c_int64_t

end module ELSI_PRECISION
