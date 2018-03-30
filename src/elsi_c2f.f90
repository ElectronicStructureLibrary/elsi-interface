! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains various utilities to aid in conversion between C and
!! Fortran.
!!
module ELSI_C2F

   use, intrinsic :: ISO_C_BINDING

   implicit none

   private

   public :: c_int_to_f_logical
   public :: c_string_to_f_string

contains

!>
!! This routine encodes the standard convention that 0 is .false., any other
!! integer is .true..
!!
function c_int_to_f_logical(int_c) result(logical_f)

   implicit none

   integer(kind=c_int), intent(in) :: int_c
   logical                         :: logical_f

   if(int_c == 0) then
     logical_f = .false.
   else
     logical_f = .true.
   endif

end function

!>
!! This routine converts a C string into a Fortran string. A Fortran string is
!! NOT just a character array without a NULL character. A Fortran string (i.e.
!! char(*)) is a separate data type from a character array (i.e. char,
!! dimension(*)) and they are NOT interoperable in interfaces.
!!
function c_string_to_f_string(string_c) result(string_f)

   implicit none

   character(kind=c_char,len=1), intent(in) :: string_c(*)
   character(len=:), allocatable            :: string_f

   integer(kind=c_int) :: string_f_len

   string_f_len = 0

   do
     if(string_c(string_f_len+1) == C_NULL_CHAR) exit
     string_f_len = string_f_len+1
   enddo

   allocate(character(len=string_f_len) :: string_f)

   string_f = transfer(string_c(1:string_f_len),string_f)

end function

end module ELSI_C2F
