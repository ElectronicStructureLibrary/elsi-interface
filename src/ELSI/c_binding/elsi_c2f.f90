! Copyright (c) 2015-2018, the ELSI team. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!  * Neither the name of the "ELectronic Structure Infrastructure" project nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

   character(kind=c_char,len=1), dimension(*), intent(in) :: string_c
   character(len=:), allocatable                          :: string_f

   integer :: string_f_len

   string_f_len = 0

   do
     if(string_c(string_f_len+1) == C_NULL_CHAR) exit
     string_f_len = string_f_len+1
   enddo

   allocate(character(len=string_f_len) :: string_f)
   string_f = transfer(string_c(1:string_f_len),string_f)

end function

end module ELSI_C2F
