! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the root directory.

!> Definitions for data types used by FortJSON
module FJSON_DATATYPE

   use FJSON_PRECISION, only: i4

   implicit none

   private

   !> Main data type of FortJSON.  Encapsulates details about the current state
   !! of the JSON file.
   type, public :: fjson_handle

      !> Has handle been initialized?
      logical :: handle_init = .false.

      !> Filename of JSON file associated with this handle (when relevant)
      character(len=:), allocatable :: file_name

      !> Fortran unit associated with JSON file
      integer(kind=i4) :: print_unit

      !> String to append to beginning of each line.  Will be modified to add
      !! indentation of the JSON output to make it "pretty".
      character(len=:), allocatable :: prefix

      !> The base prefix that this handle was initialized with, i.e. with no
      !! indentation to indicate scope.  Used for outputting error messages.
      character(len=:), allocatable :: base_prefix

      !> A simple stack to control the scoping for the JSON file.  Used to
      !! error check the user's modifications of the JSON file (only name/value
      !! pairs can reside in objects, only values can reside in arrays, etc.)
      integer(kind=i4), allocatable :: scope(:)

      !> The number of valid entries on the stack, i.e. index of the top element
      !! of the stack.
      integer(kind=i4) :: n_entries_scope

      !> Whether we've finished writing the unique root entry to the JSON file.
      !! The root of a JSON file should contain a single object or array, but
      !! some parsers accept comma-separated lists of objects or array.
      !! This is an unnecessary extension of the standard (as arrays are already
      !! supported) which will break more careful parsers.
      !! When we've reached the root level of the JSON file, we mark the file as
      !! finished and do not allow any additional writing.
      logical :: finished

      !> Controls placement of commas in objects and arrays.  The first
      !! name/value pair in an object or value in an array will not have a comma
      !! before it; all others will.  There is no need to save the state of the
      !! scope; once we've entered into a new scope, the fact that we've entered
      !! into a new scope means that the previous scope has at least one
      !! name/value pair or value in it.
      logical :: first_entry

      !> The last error encountered.  This is a dirty bit and, once an error is
      !! encountered, it will remain set until a new error is encountered.
      integer(kind=i4) :: last_error

   end type

end module FJSON_DATATYPE
