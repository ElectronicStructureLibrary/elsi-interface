! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the root directory.

!> The meat-and-potatoes module of FortJSON; provides subroutines to write to
!! files in a standard JSON format (ECMA-404)
module FJSON_RW

   use FJSON_DATATYPE,  only: fjson_handle
   use FJSON_PRECISION, only: r4, r8, i4, i8

   implicit none

   private

   private :: fjson_write
   ! Subroutines involving strings, numbers, true/false, and null
   public  :: fjson_write_value
   private :: fjson_write_value_generic
   public  :: fjson_write_name_value
   private :: fjson_write_name_value_generic
   ! Subroutines involving objects
   public :: fjson_start_object
   public :: fjson_start_name_object
   public :: fjson_finish_object
   ! Subroutines involving arrays
   public :: fjson_start_array
   public :: fjson_start_name_array
   public :: fjson_finish_array

   ! TODO: Move these into an FJSON_ERROR (or similar) module
   public :: fjson_raise_error
   public :: fjson_is_handle_init
   public :: fjson_error_message
   public :: fjson_get_last_error

   interface fjson_write_value
      module procedure fjson_write_value_i4, &
                       fjson_write_value_i8, &
                       fjson_write_value_r4, &
                       fjson_write_value_r8, &
                       fjson_write_value_str, &
                       fjson_write_value_bool, &
                       fjson_write_value_null, &
                       fjson_write_value_array_i4, &
                       fjson_write_value_array_i8, &
                       fjson_write_value_array_r4, &
                       fjson_write_value_array_r8, &
                       fjson_write_value_array_bool
   end interface

   interface fjson_write_name_value
      module procedure fjson_write_name_value_i4, &
                       fjson_write_name_value_i8, &
                       fjson_write_name_value_r4, &
                       fjson_write_name_value_r8, &
                       fjson_write_name_value_str, &
                       fjson_write_name_value_bool, &
                       fjson_write_name_value_null, &
                       fjson_write_name_value_array_i4, &
                       fjson_write_name_value_array_i8, &
                       fjson_write_name_value_array_r4, &
                       fjson_write_name_value_array_r8, &
                       fjson_write_name_value_array_bool
   end interface

contains

   !> Writes a message to the indicated unit.  All writing to file should be
   !! done through this subroutine, to minimize code rewrite for non-standard
   !! IO infrastructures
   subroutine fjson_write(fj_h, message, should_advance)

      use FJSON_CONSTANTS, only: FILE_WRITE_ERROR

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> Message to be written
      character(len=*), intent(in) :: message

      !> Should there be a line break afterwards?  (default:. true..)
      logical, optional, intent(in) :: should_advance

      character(len=*), parameter :: caller = "fjson_write"

      if (.not.fjson_is_handle_init(fj_h)) return

      ! Make sure both the file and unit are opened
      if (.not.allocated(fj_h%file_name)) then
         call fjson_raise_error(fj_h, FILE_WRITE_ERROR)
         return
      end if

      ! WPH, 10 May 2018: I've uncovered a strange error in Intel 18.0.2
      ! on my local development cluster in which the INQUIRE built-in will
      ! claim that a file or unit has not been opened unless it has been
      ! previously written to within the present scope.  This occurs even if
      ! the file has been written to previously during program execution.
      ! This bug does not occur with Intel 14.0.1, GNU 4.9.0, or PGI 17.10.
      ! Accordingly, I've left the check in here but commented it out.
      !call fjson_check_file_io(fj_h)
      !if (fjson_get_last_error(fj_h) == FILE_WRITE_ERROR) return

      if (.not.present(should_advance)) then
         write(fj_h%print_unit,"(A)") message
      else if (should_advance) then
         write(fj_h%print_unit,"(A)") message
      else
         write(fj_h%print_unit,"(A)", advance="no") message
      end if

   end subroutine

   !===========================================================================!
   !=                             NAME/VALUE PAIRS                            =!
   !===========================================================================!

   !> Generic handler for name/value pairs
   subroutine fjson_write_name_value_generic(fj_h, label, setting)

      use FJSON_CONSTANTS, only: JSON_OBJECT, NAME_VALUE_ERROR

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value in name/value pair
      character(len=*), intent(in) :: setting

      character(len=:), allocatable :: message

      character(len=*), parameter :: caller = "fjson_write_name_value_str"

      if (.not.fjson_is_handle_init(fj_h)) return

      if (fj_h%scope(fj_h%n_entries_scope) /= JSON_OBJECT) then
         call fjson_raise_error(fj_h, NAME_VALUE_ERROR)
      end if

      if (.not.fj_h%first_entry) call fjson_write(fj_h, ",")

      if (allocated(fj_h%prefix)) then
         message = fj_h%prefix // '"' // trim(adjustl(label)) // '": ' // &
            trim(adjustl(setting))
      else
         message = '"' // trim(adjustl(label)) // '": ' // &
            trim(adjustl(setting))
      end if
      call fjson_write(fj_h, message, should_advance=.false.)

      fj_h%first_entry = .false.

   end subroutine

   !> Writes out a name/value pair, where the value is an integer(kind=i4)
   subroutine fjson_write_name_value_i4(fj_h, label, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value in name/value pair
      integer(kind=i4), intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_name_value_i4"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_name_value_generic(fj_h, label, val_string)

   end subroutine

   !> Writes out a name/value pair, where the value is an integer(kind=i8)
   subroutine fjson_write_name_value_i8(fj_h, label, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value in name/value pair
      integer(kind=i8), intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_name_value_i8"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_name_value_generic(fj_h, label, val_string)

   end subroutine

   !> Writes out a name/value pair, where the value is a single-precision
   !! floating-point number
   subroutine fjson_write_name_value_r4(fj_h, label, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value in name/value pair
      real(kind=r4), intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_name_value_r4"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_name_value_generic(fj_h, label, val_string)

   end subroutine

   !> Writes out a name/value pair, where the value is a double-precision
   !! floating-point number
   subroutine fjson_write_name_value_r8(fj_h, label, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value in name/value pair
      real(kind=r8), intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_name_value_r8"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_name_value_generic(fj_h, label, val_string)

   end subroutine

   !> Writes out a name/value pair, where the value is a string
   subroutine fjson_write_name_value_str(fj_h, label, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value in name/value pair
      character(len=*), intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_name_value_str"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_name_value_generic(fj_h, label, val_string)

   end subroutine

   !> Writes out a name/value pair, where the value is boolean.  The resulting
   !! JSON output will be a literal name token (true, false), not a string
   !! ("true", "false")
   subroutine fjson_write_name_value_bool(fj_h, label, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value in name/value pair
      logical, intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_name_value_bool"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_name_value_generic(fj_h, label, val_string)

   end subroutine

   !> Writes out a name/value pair, where the value is literal name token
   !! null (not the string "null")
   subroutine fjson_write_name_value_null(fj_h, label)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_name_value_null"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value()
      call fjson_write_name_value_generic(fj_h, label, val_string)

   end subroutine

   !> Writes out a name/value pair, where the value is an array of
   !! integer(kind=i4)
   subroutine fjson_write_name_value_array_i4(fj_h, label, setting, len_arr)

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value
      integer(kind=i4), intent(in) :: setting(*)

      !> number of elements to print
      integer(kind=i4), intent(in) :: len_arr

      integer(kind=i4) :: num

      character(len=*), parameter :: caller = "fjson_write_value_array_i4"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_start_name_array(fj_h, label)
      do num = 1, len_arr
         call fjson_write_value(fj_h, setting(num))
      end do
      call fjson_finish_array(fj_h)

   end subroutine

   !> Writes out a name/value pair, where the value is an array of
   !! integer(kind=i8)
   subroutine fjson_write_name_value_array_i8(fj_h, label, setting, len_arr)

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value in name/value pair
      integer(kind=i8), intent(in) :: setting(*)

      !> number of elements to print
      integer(kind=i4), intent(in) :: len_arr

      integer(kind=i4) :: num

      character(len=*), parameter :: caller = "fjson_write_value_array_i8"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_start_name_array(fj_h, label)
      do num = 1, len_arr
         call fjson_write_value(fj_h, setting(num))
      end do
      call fjson_finish_array(fj_h)

   end subroutine

   !> Writes out a name/value pair , where the value is an array of
   !! single-precision floating-point number
   subroutine fjson_write_name_value_array_r4(fj_h, label, setting, len_arr)

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value in name/value pair
      real(kind=r4), intent(in) :: setting(*)

      !> number of elements to print
      integer(kind=i4), intent(in) :: len_arr

      integer(kind=i4) :: num

      character(len=*), parameter :: caller = "fjson_write_value_array_r4"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_start_name_array(fj_h, label)
      do num = 1, len_arr
         call fjson_write_value(fj_h, setting(num))
      end do
      call fjson_finish_array(fj_h)

   end subroutine

   !> Writes out a value, where the value is an array of double-precision
   !! floating-point number
   subroutine fjson_write_name_value_array_r8(fj_h, label, setting, len_arr)

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value in name/value pair
      real(kind=r8), intent(in) :: setting(*)

      !> number of elements to print
      integer(kind=i4), intent(in) :: len_arr

      integer(kind=i4) :: num

      character(len=*), parameter :: caller = "fjson_write_value_array_r8"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_start_name_array(fj_h, label)
      do num = 1, len_arr
         call fjson_write_value(fj_h, setting(num))
      end do
      call fjson_finish_array(fj_h)

   end subroutine

   !> Writes out a value, where the value is an array of logical variables
   subroutine fjson_write_name_value_array_bool(fj_h, label, setting, len_arr)

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      !> value in name/value pair
      logical, intent(in) :: setting(*)

      !> number of elements to print
      integer(kind=i4), intent(in) :: len_arr

      integer(kind=i4) :: num

      character(len=*), parameter :: caller = "fjson_write_value_array_bool"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_start_name_array(fj_h, label)
      do num = 1, len_arr
         call fjson_write_value(fj_h, setting(num))
      end do
      call fjson_finish_array(fj_h)

   end subroutine

   !===========================================================================!
   !=                                  VALUES                                 =!
   !===========================================================================!

   !> Generic handler for values
   subroutine fjson_write_value_generic(fj_h, setting)

      use FJSON_CONSTANTS, only: JSON_ROOT, JSON_OBJECT, FINISHED_ERROR, &
                                 VALUE_ERROR

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      character(len=*), intent(in) :: setting

      character(len=:), allocatable :: message

      character(len=*), parameter :: caller = "fjson_write_value_generic"

      if (.not.fjson_is_handle_init(fj_h)) return

      if (fj_h%scope(fj_h%n_entries_scope) == JSON_ROOT .and. &
          fj_h%finished) then
         call fjson_raise_error(fj_h, FINISHED_ERROR)
      end if

      if (fj_h%scope(fj_h%n_entries_scope) == JSON_OBJECT) then
         call fjson_raise_error(fj_h, VALUE_ERROR)
      end if

      if (.not.fj_h%first_entry) call fjson_write(fj_h, ",")

      if (allocated(fj_h%prefix)) then
         message = fj_h%prefix // trim(adjustl(setting))
      else
         message =  trim(adjustl(setting))
      end if
      call fjson_write(fj_h, message, should_advance=.false.)

      fj_h%first_entry = .false.

      ! Only one JSON value is allowed as the root entry of the text
      if (fj_h%n_entries_scope == 1) fj_h%finished = .true.

   end subroutine

   !> Writes out a value, where the value is an integer(kind=i4)
   subroutine fjson_write_value_i4(fj_h, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      integer(kind=i4), intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_value_i4"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_value_generic(fj_h, val_string)

   end subroutine

   !> Writes out a value, where the value is an integer(kind=i8)
   subroutine fjson_write_value_i8(fj_h, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      integer(kind=i8), intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_value_i8"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_value_generic(fj_h, val_string)

   end subroutine

   !> Writes out a value, where the value is a single-precision
   !! floating-point number
   subroutine fjson_write_value_r4(fj_h, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      real(kind=r4), intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_value_r4"

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_value_generic(fj_h, val_string)

   end subroutine

   !> Writes out a value, where the value is a double-precision
   !! floating-point number
   subroutine fjson_write_value_r8(fj_h, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      real(kind=r8), intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_value_r8"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_value_generic(fj_h, val_string)

   end subroutine

   !> Writes out a value, where the value is a string
   subroutine fjson_write_value_str(fj_h, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      character(len=*), intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_value_str"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_value_generic(fj_h, val_string)

   end subroutine

   !> Writes out a value, where the value is boolean.  The resulting
   !! JSON output will be a literal name token (true, false), not a string
   !! ("true", "false")
   subroutine fjson_write_value_bool(fj_h, setting)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      logical, intent(in) :: setting

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_value_bool"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value(setting)
      call fjson_write_value_generic(fj_h, val_string)

   end subroutine

   !> Writes out a value, where the value is the literal name token null (not
   !! the string "null")
   subroutine fjson_write_value_null(fj_h)

      use FJSON_UTILITIES, only : fjson_convert_var_to_value

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      character(len=:), allocatable :: val_string

      character(len=*), parameter :: caller = "fjson_write_value_null"

      if (.not.fjson_is_handle_init(fj_h)) return

      val_string = fjson_convert_var_to_value()
      call fjson_write_value_generic(fj_h, val_string)

   end subroutine

   !> Writes out a value, where the value is an array of integer(kind=i4)
   subroutine fjson_write_value_array_i4(fj_h, setting, len_arr)

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      integer(kind=i4), intent(in) :: setting(*)

      !> number of elements to print
      integer(kind=i4), intent(in) :: len_arr

      integer(kind=i4) :: num

      character(len=*), parameter :: caller = "fjson_write_value_array_i4"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_start_array(fj_h)
      do num = 1, len_arr
         call fjson_write_value(fj_h, setting(num))
      end do
      call fjson_finish_array(fj_h)

   end subroutine

   !> Writes out a value, where the value is an array of integer(kind=i8)
   subroutine fjson_write_value_array_i8(fj_h, setting, len_arr)

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      integer(kind=i8), intent(in) :: setting(*)

      !> number of elements to print
      integer(kind=i4), intent(in) :: len_arr

      integer(kind=i4) :: num

      character(len=*), parameter :: caller = "fjson_write_value_array_i8"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_start_array(fj_h)
      do num = 1, len_arr
         call fjson_write_value(fj_h, setting(num))
      end do
      call fjson_finish_array(fj_h)

   end subroutine

   !> Writes out a value, where the value is an array of single-precision
   !! floating-point number
   subroutine fjson_write_value_array_r4(fj_h, setting, len_arr)

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      real(kind=r4), intent(in) :: setting(*)

      !> number of elements to print
      integer(kind=i4), intent(in) :: len_arr

      integer(kind=i4) :: num

      character(len=*), parameter :: caller = "fjson_write_value_array_r4"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_start_array(fj_h)
      do num = 1, len_arr
         call fjson_write_value(fj_h, setting(num))
      end do
      call fjson_finish_array(fj_h)

   end subroutine

   !> Writes out a value, where the value is an array of double-precision
   !! floating-point number
   subroutine fjson_write_value_array_r8(fj_h, setting, len_arr)

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      real(kind=r8), intent(in) :: setting(*)

      !> number of elements to print
      integer(kind=i4), intent(in) :: len_arr

      integer(kind=i4) :: num

      character(len=*), parameter :: caller = "fjson_write_value_array_r8"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_start_array(fj_h)
      do num = 1, len_arr
         call fjson_write_value(fj_h, setting(num))
      end do
      call fjson_finish_array(fj_h)

   end subroutine

   !> Writes out a value, where the value is an array of logical variables
   subroutine fjson_write_value_array_bool(fj_h, setting, len_arr)

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> value
      logical, intent(in) :: setting(*)

      !> number of elements to print
      integer(kind=i4), intent(in) :: len_arr

      integer(kind=i4) :: num

      character(len=*), parameter :: caller = "fjson_write_value_array_bool"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_start_array(fj_h)
      do num = 1, len_arr
         call fjson_write_value(fj_h, setting(num))
      end do
      call fjson_finish_array(fj_h)

   end subroutine

   !===========================================================================!
   !=                                 OBJECTS                                 =!
   !===========================================================================!

   !> Starts a new object as a value
   subroutine fjson_start_object(fj_h)

      use FJSON_CONSTANTS, only: JSON_OBJECT
      use FJSON_UTILITIES, only: fjson_append_string, fjson_resize_array

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      character(len=*), parameter :: caller = "fjson_start_object"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_write_value_generic(fj_h, "{")
      ! Override default behavior and add a line break
      fj_h%first_entry = .true.
      call fjson_write(fj_h, "")

      call fjson_append_string(fj_h%prefix, "  ")

      if (fj_h%n_entries_scope == size(fj_h%scope)) then
         call fjson_resize_array(fj_h%scope, 2*size(fj_h%scope))
      end if
      fj_h%n_entries_scope = fj_h%n_entries_scope + 1
      fj_h%scope(fj_h%n_entries_scope) = JSON_OBJECT

   end subroutine

   !> Starts a new object as a name/value pair
   subroutine fjson_start_name_object(fj_h, label)

      use FJSON_CONSTANTS, only: JSON_OBJECT
      use FJSON_UTILITIES, only: fjson_append_string, fjson_resize_array

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      character(len=*), parameter :: caller = "fjson_start_name_object"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_write_name_value_generic(fj_h, label, "{")
      ! Override default behavior and add a line break
      fj_h%first_entry = .true.
      call fjson_write(fj_h, "")

      call fjson_append_string(fj_h%prefix, "  ")

      if (fj_h%n_entries_scope == size(fj_h%scope)) then
         call fjson_resize_array(fj_h%scope, 2*size(fj_h%scope))
      end if
      fj_h%n_entries_scope = fj_h%n_entries_scope + 1
      fj_h%scope(fj_h%n_entries_scope) = JSON_OBJECT

   end subroutine

   !> Finishes the current object, whether it was created as a value or a
   !! name/value pair.
   subroutine fjson_finish_object(fj_h)

      use FJSON_CONSTANTS, only: JSON_NULL, JSON_OBJECT, OBJECT_ERROR
      use FJSON_UTILITIES, only: fjson_truncate_string

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      character(len=:), allocatable :: message

      character(len=*), parameter :: caller = "fjson_finish_object"

      if (.not.fjson_is_handle_init(fj_h)) return

      if (fj_h%scope(fj_h%n_entries_scope) /= JSON_OBJECT) then
         call fjson_raise_error(fj_h, OBJECT_ERROR)
      end if

      if (fj_h%first_entry) then
         ! The object is empty, which is allowed by the standard
         fj_h%first_entry = .false.
      else
         call fjson_write(fj_h, "")
      end if

      call fjson_truncate_string(fj_h%prefix, 2)

      if (allocated(fj_h%prefix)) then
         message = fj_h%prefix // '}'
      else
         message = '}'
      end if
      call fjson_write(fj_h, message, should_advance=.false.)

      fj_h%scope(fj_h%n_entries_scope) = JSON_NULL
      fj_h%n_entries_scope = fj_h%n_entries_scope - 1
      fj_h%first_entry = .false.

      ! Only one JSON value is allowed as the root entry of the text
      if (fj_h%n_entries_scope == 1) fj_h%finished = .true.

   end subroutine

   !===========================================================================!
   !=                                 ARRAYS                                  =!
   !===========================================================================!

   !> Starts a new array as a value
   subroutine fjson_start_array(fj_h)

      use FJSON_CONSTANTS, only: JSON_ARRAY
      use FJSON_UTILITIES, only: fjson_append_string, fjson_resize_array

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      character(len=*), parameter :: caller = "fjson_start_array"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_write_value_generic(fj_h, "[")
      ! Override default behavior and add a line break
      fj_h%first_entry = .true.
      call fjson_write(fj_h, "")

      call fjson_append_string(fj_h%prefix, "  ")

      if (fj_h%n_entries_scope == size(fj_h%scope)) then
         call fjson_resize_array(fj_h%scope, 2*size(fj_h%scope))
      end if
      fj_h%n_entries_scope = fj_h%n_entries_scope + 1
      fj_h%scope(fj_h%n_entries_scope) = JSON_ARRAY

   end subroutine

   !> Starts a new array as a name/value pair
   subroutine fjson_start_name_array(fj_h, label)

      use FJSON_CONSTANTS, only: JSON_ARRAY
      use FJSON_UTILITIES, only: fjson_append_string, fjson_resize_array

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> name in name/value pair
      character(len=*), intent(in) :: label

      character(len=*), parameter :: caller = "fjson_start_name_array"

      if (.not.fjson_is_handle_init(fj_h)) return

      call fjson_write_name_value_generic(fj_h, label, "[")
      ! Override default behavior and add a line break
      fj_h%first_entry = .true.
      call fjson_write(fj_h, "")

      call fjson_append_string(fj_h%prefix, "  ")

      if (fj_h%n_entries_scope == size(fj_h%scope)) then
         call fjson_resize_array(fj_h%scope, 2*size(fj_h%scope))
      end if
      fj_h%n_entries_scope = fj_h%n_entries_scope + 1
      fj_h%scope(fj_h%n_entries_scope) = JSON_ARRAY

   end subroutine

   !> Finishes the current array, whether it was created as a value or a
   !! name/value pair.
   subroutine fjson_finish_array(fj_h)

      use FJSON_CONSTANTS, only: JSON_NULL, JSON_ARRAY, ARRAY_ERROR
      use FJSON_UTILITIES, only: fjson_truncate_string

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      character(len=:), allocatable :: message

      character(len=*), parameter :: caller = "fjson_finish_object"

      if (.not.fjson_is_handle_init(fj_h)) return

      if (fj_h%scope(fj_h%n_entries_scope) /= JSON_ARRAY) then
         call fjson_raise_error(fj_h, ARRAY_ERROR)
      end if

      if (fj_h%first_entry) then
         ! The array is empty, which is allowed by the standard
         fj_h%first_entry = .false.
      else
         call fjson_write(fj_h, "")
      end if

      call fjson_truncate_string(fj_h%prefix, 2)
      if (allocated(fj_h%prefix)) then
         message = fj_h%prefix // ']'
      else
         message = ']'
      end if
      call fjson_write(fj_h, message, should_advance=.false.)

      fj_h%scope(fj_h%n_entries_scope) = JSON_NULL
      fj_h%n_entries_scope = fj_h%n_entries_scope - 1
      fj_h%first_entry = .false.

      ! Only one JSON value is allowed as the root entry of the text
      if (fj_h%n_entries_scope == 1) fj_h%finished = .true.

   end subroutine

   !> If an error has occured, write an error message to the file and sets the
   !! relevent code
   subroutine fjson_raise_error(fj_h, error_code)

      use FJSON_CONSTANTS, only: NO_ERROR, FILE_OPEN_ERROR, HANDLE_UNINIT_ERROR

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> Error code being raised
      integer(kind=i4), intent(in) :: error_code

      character(len=:), allocatable :: message

      character(len=*), parameter :: caller = "fjson_raise_error"

      ! Do nothing if error flag if the current operation was successful
      if (error_code == NO_ERROR) return
      fj_h%last_error = error_code

      ! Output error code to file, as long as the nature of the error does not
      ! preclude this
      if (error_code /= FILE_OPEN_ERROR .and. &
          error_code /= HANDLE_UNINIT_ERROR) then
         message = fj_h%base_prefix // " " // fjson_error_message(error_code)
         call fjson_write(fj_h, message)
      end if

   end subroutine

   !> Return an error message based on the provided error code
   pure function fjson_error_message(error_code) result(error_message)

      use FJSON_CONSTANTS, only: NO_ERROR, FILE_OPEN_ERROR, FINISHED_ERROR, &
                                 NAME_VALUE_ERROR, VALUE_ERROR, OBJECT_ERROR, &
                                 ARRAY_ERROR, HANDLE_UNINIT_ERROR, &
                                 FILE_WRITE_ERROR

      implicit none

      !> Error code being raised
      integer(kind=i4), intent(in) :: error_code

      !> Error message for provided error code
      character(len=:), allocatable :: error_message

      character(len=*), parameter :: caller = "fjson_error_message"

      select case (error_code)
         case (NO_ERROR)
            error_message = "FortJSON Error:  No error"
         case (FILE_OPEN_ERROR)
            error_message = "FortJSON Error:  Errors encountered when opening &
                            &files."
         case (FINISHED_ERROR)
            error_message = "FortJSON Error:  Operations were requested on a &
                            &JSON text which was already complete."
         case (NAME_VALUE_ERROR)
            error_message = "FortJSON Error:  Errors encountered when writing &
                            &a name/value pair."
         case (VALUE_ERROR)
            error_message = "FortJSON Error:  Errors encountered when writing &
                            &a value."
         case (OBJECT_ERROR)
            error_message = "FortJSON Error:  Errors encountered when writing &
                            &an object."
         case (ARRAY_ERROR)
            error_message = "FortJSON Error:  Errors encountered when writing &
                            &an array."
         case (HANDLE_UNINIT_ERROR)
            error_message = "FortJSON Error:  Handle was accessed while &
                            &uninitialized."
         case (FILE_WRITE_ERROR)
            error_message = "FortJSON Error:  Errors encountered when writing &
                            &to files."
         case default
            error_message = "FortJSON Error:  Unknown errors encountered."
      end select

   end function

   !> Return the last error code.  The meanings of error codes are defined in
   !! the FJSON_CONSTANTS module, and the associated error message is returned
   !! by fjson_error_message().
   pure function fjson_get_last_error(fj_h) result(last_error)

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(in) :: fj_h

      !> The last error code encountered by the code.
      integer(kind=i4) :: last_error

      character(len=*), parameter :: caller = "fjson_get_last_error"

      last_error = fj_h%last_error

   end function

   !> Checks the initialization status of the handle and returns the result.  If
   !! the handle is not initialized, will also raise an HANDLE_UNINIT_ERROR
   !! error in the handle.
   function fjson_is_handle_init(fj_h) result(handle_init)

      use FJSON_CONSTANTS, only: HANDLE_UNINIT_ERROR

      implicit none

      !> FortJSON handle
      type(fjson_handle), intent(inout) :: fj_h

      !> Is the handle initialized?
      logical :: handle_init

      character(len=*), parameter :: caller = "fjson_is_handle_init"

      handle_init = fj_h%handle_init
      if (.not.handle_init) call fjson_raise_error(fj_h, HANDLE_UNINIT_ERROR)

   end function

   !> Checks if the unit and JSON file associated with the current handle
   !! are open
   subroutine fjson_check_file_io(fj_h)

      use FJSON_CONSTANTS, only: FILE_WRITE_ERROR

      implicit none

      !> FortJSON handle associated with the file/unit
      type(fjson_handle), intent(inout) :: fj_h

      character(len=*), parameter :: caller = "fjson_check_file_io"

      logical :: is_opened

      if (.not.allocated(fj_h%file_name)) then
         call fjson_raise_error(fj_h, FILE_WRITE_ERROR)
         return
      end if

      inquire(file=fj_h%file_name, opened=is_opened)
      if (.not.is_opened) then
         call fjson_raise_error(fj_h, FILE_WRITE_ERROR)
         return
      end if

      inquire(unit=fj_h%print_unit, opened=is_opened)
      if (.not.is_opened) then
         call fjson_raise_error(fj_h, FILE_WRITE_ERROR)
         return
      end if

   end subroutine

end module FJSON_RW
