! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the root directory.

!> Constants used by FortJSON
module FJSON_CONSTANTS

   use FJSON_PRECISION, only: i4

   implicit none

   !> Uninitialized value for integer(kind=i4) values
   integer(kind=i4), parameter :: UNSET = -910910_i4

   !> Uninitialized value for string values
   character(len=*), parameter :: UNSET_STR = "N/A"

   !> Used by JSON scoping stack to denote the element is not set
   integer(kind=i4), parameter :: JSON_NULL = 0

   !> Used by JSON scoping stack to denote the element is the root/top-level of
   !! the JSON file
   integer(kind=i4), parameter :: JSON_ROOT = 1

   !> Used by JSON scoping stack to denote the element is an object
   integer(kind=i4), parameter :: JSON_OBJECT = 2

   !> Used by JSON scoping stack to denote the element is an array
   integer(kind=i4), parameter :: JSON_ARRAY = 3

   !> Used by error handling to denote no error has been encountered
   integer(kind=i4), parameter :: NO_ERROR = 0

   !> Used by error handling to denote that an error occured when opening a file
   integer(kind=i4), parameter :: FILE_OPEN_ERROR = 1

   !> Used by error handling to denote that writing was attempted to a JSON text
   !! which was finished, i.e. already had a completely top-level entry
   integer(kind=i4), parameter :: FINISHED_ERROR = 2

   !> Used by error handling to denote that a name/value pair was
   !! inappropriately written
   integer(kind=i4), parameter :: NAME_VALUE_ERROR = 3

   !> Used by error handling to denote that a value was inappropriately written
   integer(kind=i4), parameter :: VALUE_ERROR = 4

   !> Used by error handling to denote that something went wrong when writing an
   !! object
   integer(kind=i4), parameter :: OBJECT_ERROR = 5

   !> Used by error handling to denote that something went wrong when writing an
   !! array
   integer(kind=i4), parameter :: ARRAY_ERROR = 6

   !> Used by error handling to denote that the handle was accessed while it was
   !! uninitialized.  Note that initializing the handle will clear any error
   !! code, including this one.
   integer(kind=i4), parameter :: HANDLE_UNINIT_ERROR = 7

   !> Used by error handling to denote that an error occurred when writing to a
   !! file
   integer(kind=i4), parameter :: FILE_WRITE_ERROR = 8

   !> Number of characters in RFC3339 datetime string with time zone offset.
   integer(kind=i4), parameter :: DATETIME_LEN = 29

end module FJSON_CONSTANTS
