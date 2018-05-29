! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the root directory.

!> Handles initializing/deinitializing the FortJSON handle and opening/closing
!! files
module FJSON_SETUP

   use FJSON_DATATYPE,  only: fjson_handle
   use FJSON_PRECISION, only: i4

   implicit none

   private

   public :: fjson_init_io
   public :: fjson_reset_fj_handle
   ! We count units as a data structure here
   ! TODO:  Change "file" to "unit"
   public :: fjson_open_file
   public :: fjson_close_file

contains

   !> Initializes a FortJSON handle for reading and writing to files.
   subroutine fjson_init_io(fj_h, print_unit, file_name, prefix)

      use FJSON_CONSTANTS, only: NO_ERROR, JSON_ROOT

      implicit none

      !> FortJSON handle to initialize
      type(fjson_handle), intent(out) :: fj_h

      !> Fortran unit to modify
      integer(kind=i4), intent(in) :: print_unit

      !> Filename for file to modify (if relevant)
      character(len=*), intent(in) :: file_name

      !> Prefix to append to beginning of every line
      character(len=*), intent(in) :: prefix

      character(len=*), parameter :: caller = "fjson_init_io"

      ! For safety
      call fjson_reset_fj_handle(fj_h)

      fj_h%handle_init = .true.
      fj_h%print_unit  = print_unit
      fj_h%file_name   = file_name
      fj_h%prefix      = prefix
      fj_h%base_prefix = prefix
      fj_h%first_entry = .true.
      fj_h%finished    = .false.
      fj_h%last_error  = NO_ERROR

      fj_h%n_entries_scope = 1
      if (allocated(fj_h%scope)) deallocate(fj_h%scope)
      allocate(fj_h%scope(1))
      fj_h%scope(1) = JSON_ROOT

   end subroutine

   !> Resets a FortJSON handle.
   subroutine fjson_reset_fj_handle(fj_h)

      use FJSON_CONSTANTS, only: NO_ERROR, UNSET

      implicit none

      !> FortJSON handle to reset
      type(fjson_handle), intent(inout) :: fj_h

      character(len=*), parameter :: caller = "fjson_reset_fj_handle"

      fj_h%handle_init = .false.
      if (allocated(fj_h%file_name)) deallocate(fj_h%file_name)
      fj_h%print_unit  = UNSET
      fj_h%first_entry = .true.
      fj_h%finished    = .true.
      fj_h%last_error  = NO_ERROR
      if (allocated(fj_h%prefix)) deallocate(fj_h%prefix)
      if (allocated(fj_h%base_prefix)) deallocate(fj_h%base_prefix)

      fj_h%n_entries_scope = 0
      if (allocated(fj_h%scope)) deallocate(fj_h%scope)

   end subroutine

   !> Initialize a FortJSON handle for the JSON file and opens the file.
   subroutine fjson_open_file(fj_h, print_unit, file_name)

      use FJSON_CONSTANTS, only: FILE_OPEN_ERROR
      use FJSON_RW,        only: fjson_raise_error ! TODO: Bad dependency!

      implicit none

      !> FortJSON handle to be initialized
      type(fjson_handle), intent(out) :: fj_h

      !> Fortran unit to modify
      integer(kind=i4), intent(in) :: print_unit

      !> Filename for file to modify (if relevant)
      character(len=*), intent(in) :: file_name

      integer(kind=i4) :: ios

      character(len=*), parameter :: caller = "fjson_open_file"

      call fjson_init_io(fj_h, print_unit, file_name, "")

      open(unit=fj_h%print_unit, file=fj_h%file_name, iostat=ios)

      if (ios /= 0) call fjson_raise_error(fj_h, FILE_OPEN_ERROR)

   end subroutine

   !> Closes the JSON file and resets the FortJSON handle.
   subroutine fjson_close_file(fj_h)

      implicit none

      !> FortJSON handle associated with the file to be closed
      type(fjson_handle), intent(inout) :: fj_h

      character(len=*), parameter :: caller = "fjson_close_file"

      close(fj_h%print_unit)

      call fjson_reset_fj_handle(fj_h)

   end subroutine

end module FJSON_SETUP
