! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide stub routines which are only compiled when the actual ELPA-AEO is not
!! available.
!!
module ELPA

   use ELSI_PRECISION, only: i4

   implicit none

   private

   type, public :: elpa_t
      integer(kind=i4) :: dummy

   contains

   procedure, nopass :: elpa_dummy

   end type

   type, public :: elpa_autotune_t
      integer(kind=i4) :: dummy

   contains

   procedure, nopass :: elpa_autotune_dummy

   end type

contains

subroutine elpa_dummy()

   implicit none

   write(*,"(A)") "**Error! An ELPA stub routine was called"
   stop

end subroutine

subroutine elpa_autotune_dummy()

   implicit none

   write(*,"(A)") "**Error! An ELPA stub routine was called"
   stop

end subroutine

end module ELPA
