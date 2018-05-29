! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the root directory.

!> Various helper subroutines that are not specific to FortJSON and could be
!! reused for other purposes.
module FJSON_UTILITIES

   use FJSON_PRECISION, only: i4, i8, r4, r8

   implicit none

   private

   public :: fjson_resize_array
   public :: fjson_convert_var_to_value
   public :: fjson_append_string
   public :: fjson_truncate_string
   public :: fjson_get_datetime_rfc3339

   interface fjson_resize_array
      module procedure fjson_resize_array_1D_i4
   end interface

   interface fjson_convert_var_to_value
      module procedure fjson_convert_i4_to_value
      module procedure fjson_convert_i8_to_value
      module procedure fjson_convert_r4_to_value
      module procedure fjson_convert_r8_to_value
      module procedure fjson_convert_str_to_value
      module procedure fjson_convert_bool_to_value
      module procedure fjson_convert_null_to_value
   end interface

contains

   !> This routine resizes an 1D array of 4-byte integers.  If the new size is
   !! zero or negative, the array will be deallocated.  If the new size is
   !! smaller than the older size, the array will be truncated.  If the new size
   !! is larger than the older size, the additional elements will be initialized
   !! to zero.
   subroutine fjson_resize_array_1D_i4(array, new_size)

      implicit none

      !> Array to be resized
      integer(kind=i4), allocatable, intent(inout) :: array(:)

      !> New size for array
      integer(kind=i4), intent(in) :: new_size

      integer(kind=i4), allocatable :: old_array(:)
      integer(kind=i4) :: old_size

      character(len=*), parameter :: caller = "fjson_resize_array_1D_i4"

      if (.not.allocated(array)) then
         allocate(array(new_size))
         array = 0_i4
         return
      end if

      if (new_size <= 0) then
         deallocate(array)
         return
      end if

      old_size = SIZE(array)
      if (new_size == old_size) return

      allocate(old_array(old_size))
      old_array = array
      deallocate(array)
      allocate(array(new_size))

      if (old_size < new_size) then
         array(1:old_size) = old_array(1:old_size)
         array(old_size+1:new_size) = 0_i4
      else ! new_size < old_size
         array(1:new_size) = old_array(1:new_size)
      end if

      deallocate(old_array)

   end subroutine

   !> Converts an integer(kind=i4) to a properly-formatted JSON value string
   function fjson_convert_i4_to_value(var) result(val_string)

      implicit none

      !> Variable to be converted
      integer(kind=i4), intent(in) :: var

      !> JSON string
      character(len=:), allocatable :: val_string

      character(len=40) :: int_string ! Number of digits for unsigned 8-byte int

      character(len=*), parameter :: caller = "fjson_convert_i4_to_value"

      write(int_string, "(I20)") var
      val_string = trim(adjustl(int_string))

   end function

   !> Converts an integer(kind=i8) to a properly-formatted JSON value string
   function fjson_convert_i8_to_value(var) result(val_string)

      implicit none

      !> Variable to be converted
      integer(kind=i8), intent(in) :: var

      !> JSON string
      character(len=:), allocatable :: val_string

      character(len=40) :: int_string ! Number of digits for unsigned 8-byte int

      character(len=*), parameter :: caller = "fjson_convert_i8_to_value"

      write(int_string, "(I30)") var
      val_string = trim(adjustl(int_string))

   end function

   !> Converts a real(kind=r4) to a properly-formatted JSON value string
   function fjson_convert_r4_to_value(var) result(val_string)

      implicit none

      !> Variable to be converted
      real(kind=r4), intent(in) :: var

      !> JSON string
      character(len=:), allocatable :: val_string

      character(len=40) :: real_string ! Num. characters in format specifier

      character(len=*), parameter :: caller = "fjson_convert_r4_to_value"

      write(real_string, "(E20.8)") var
      val_string = trim(adjustl(real_string))

   end function

   !> Converts a real(kind=r8) to a properly-formatted JSON value string
   function fjson_convert_r8_to_value(var) result(val_string)

      implicit none

      !> Variable to be converted
      real(kind=r8), intent(in) :: var

      !> JSON string
      character(len=:), allocatable :: val_string

      character(len=40) :: real_string ! Num. characters in format specifier

      character(len=*), parameter :: caller = "fjson_convert_r8_to_value"

      ! WPH (23 April 2018): Exponent field width of 3 needed to circumvent
      !     a bug in ifort 14, where a negative exponent of -100 or lower will
      !     "clip" into the E specifier, causing a malformed string.
      write(real_string, "(E30.15E3)") var
      val_string = trim(adjustl(real_string))

   end function

   !> Converts a string to a properly-formatted JSON value string
   function fjson_convert_str_to_value(var) result(val_string)

      implicit none

      !> Variable to be converted
      character(len=*), intent(in) :: var

      !> JSON string
      character(len=:), allocatable :: val_string

      character(len=:), allocatable :: str_string

      character(len=*), parameter :: caller = "fjson_convert_str_to_value"

      str_string = '"' // var // '"'
      val_string = trim(adjustl(str_string))

   end function

   !> Converts a boolean to a properly-formatted JSON value string
   !> The literal name tokens are generated (true, false), not the strings
   !! ("true", "false")
   function fjson_convert_bool_to_value(var) result(val_string)

      implicit none

      !> Variable to be converted
      logical, intent(in) :: var

      !> JSON string
      character(len=:), allocatable :: val_string

      character(len=:), allocatable :: bool_string

      character(len=*), parameter :: caller = "fjson_convert_bool_to_value"

      if (var) then
         bool_string = "true"
      else
         bool_string = "false"
      end if
      val_string = trim(adjustl(bool_string))

   end function

   !> Return a null literal name token as a properly-formatted JSON value string
   function fjson_convert_null_to_value() result(val_string)

      implicit none

      !> JSON string
      character(len=:), allocatable :: val_string

      character(len=:), allocatable :: null_string

      character(len=*), parameter :: caller = "fjson_convert_null_to_value"

      null_string = "null"
      val_string = trim(adjustl(null_string))

   end function

   !> Appends a string to the end of a target string.  If the target string is
   !! unallocated, it will be allocated.
   subroutine fjson_append_string(target_str, append_str)

      implicit none

      !> Target string
      character(len=:), allocatable, intent(inout) :: target_str

      !> String to append to end of target string
      character(len=*), intent(in) :: append_str

      character(len=:), allocatable :: tmp_string

      character(len=*), parameter :: caller = "fjson_append_string"

      if (allocated(target_str)) then
         tmp_string = target_str // append_str
      else
         tmp_string = append_str
      end if

      if (allocated(target_str)) then
         deallocate(target_str)
      end if

      target_str = tmp_string

      deallocate(tmp_string)

   end subroutine

   !> Removes a set number of characters from the end of a string.  If the
   !! number of characters equals or exceeds the length of the string, the
   !! string will be deallocated.
   subroutine fjson_truncate_string(target_str, n_chars_to_remove)

      implicit none

      !> String to be truncated
      character(len=:), allocatable, intent(inout) :: target_str

      !> Number of characters to remove from end of string
      integer(kind=i4), intent(in) :: n_chars_to_remove

      integer(kind=i4) :: size_new_str
      character(len=:), allocatable :: tmp_str

      character(len=*), parameter :: caller = "fjson_truncate_string"

      if (allocated(target_str)) then
         size_new_str = len(target_str) - n_chars_to_remove
      else
         return
      end if

      if (size_new_str < 1) then
         deallocate(target_str)
         return
      end if

      tmp_str = target_str(1:size_new_str)

      deallocate(target_str)

      target_str = tmp_str

      deallocate(tmp_str)

   end subroutine

   !> Returns the current datetime, formatted as an RFC3339 string with time
   !! zone offset.
   subroutine fjson_get_datetime_rfc3339(dt)

      use FJSON_CONSTANTS, only: DATETIME_LEN

      implicit none

      !> Current datetime, formatted as an RFD3339 string with time zone offset.
      character(len=DATETIME_LEN), intent(out) :: dt

      integer(kind=i4) :: datetime(8)
      integer(kind=i4) :: tmp_int
      character(len=4) :: year
      character(len=2) :: month
      character(len=2) :: day
      character(len=2) :: hour
      character(len=2) :: minute
      character(len=2) :: second
      character(len=3) :: millisecond
      character(len=1) :: timezone_sign
      character(len=2) :: timezone_hour
      character(len=2) :: timezone_min

      character(len=*), parameter :: caller = "fjson_get_datetime_rfc3339"

      call date_and_time(values=datetime)

      ! Get year
      if (datetime(1) < 10) then
         write(year, "(A3,I1)") "000", datetime(1)
      else if (datetime(1) < 100) then
         write(year, "(A2,I2)") "00", datetime(1)
      else if (datetime(1) < 1000) then
         write(year, "(A1,I3)") "0", datetime(1)
      else
         write(year,"(I4)") datetime(1)
      end if

      ! Get month
      if (datetime(2) < 10) then
         write(month, "(A1,I1)") "0", datetime(2)
      else
         write(month, "(I2)") datetime(2)
      end if

      ! Get day
      if (datetime(3) < 10) then
         write(day, "(A1,I1)") "0", datetime(3)
      else
         write(day, "(I2)") datetime(3)
      end if

      ! Get hour
      if (datetime(5) < 10) then
         write(hour, "(A1,I1)") "0", datetime(5)
      else
         write(hour, "(I2)") datetime(5)
      end if

      ! Get minute
      if (datetime(6) < 10) then
         write(minute, "(A1,I1)") "0", datetime(6)
      else
         write(minute, "(I2)") datetime(6)
      end if

      ! Get second
      if (datetime(7) < 10) then
         write(second, "(A1,I1)") "0", datetime(7)
      else
         write(second, "(I2)") datetime(7)
      end if

      ! Get millisecond
      if (datetime(8) < 10) then
         write(millisecond, "(A2,I1)") "00", datetime(8)
      else if (datetime(8) < 100) then
         write(millisecond, "(A1,I2)") "0", datetime(8)
      else
         write(millisecond, "(I3)") datetime(8)
      end if

      ! Get time zone sign (ahead or behind UTC)
      if (datetime(4) < 0) then
         timezone_sign = "-"
         datetime(4)   = -1*datetime(4)
      else
         timezone_sign = "+"
      end if

      ! Get timezone minutes
      tmp_int = mod(datetime(4), 60)
      if (tmp_int < 10) then
         write(timezone_min, "(A1,I1)") "0", tmp_int
      else
         write(timezone_min, "(I2)") tmp_int
      end if

      ! Get timezone hours
      tmp_int = datetime(4)/60
      if (tmp_int < 10) then
         write(timezone_hour, "(A1,I1)") "0", tmp_int
      else
         write(timezone_hour, "(I2)") tmp_int
      end if

      write(dt, "(A4,A1,A2,A1,A2,A1,A2,A1,A2,A1,A2,A1,A3,A1,A2,A1,A2)") &
         year, "-", month, "-", day, "T", hour, ":", minute, ":", second, ".", &
         millisecond, timezone_sign, timezone_hour, ":", timezone_min

   end subroutine

end module FJSON_UTILITIES
