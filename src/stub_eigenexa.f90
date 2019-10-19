! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide stub routines which are only compiled when the actual EigenExa is not
!! available.
!!
module ELSI_EIGENEXA

   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_PRECISION, only: r8

   implicit none

   private

   public :: elsi_init_eigenexa
   public :: elsi_cleanup_eigenexa
   public :: elsi_solve_eigenexa

   interface elsi_solve_eigenexa
      module procedure elsi_solve_eigenexa_real
   end interface

contains

subroutine elsi_init_eigenexa(ph,bh)

   implicit none

   type(elsi_param_t) :: ph
   type(elsi_basic_t) :: bh

end subroutine

subroutine elsi_solve_eigenexa_real(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t) :: ph
   type(elsi_basic_t) :: bh
   real(kind=r8) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8) :: eval(ph%n_basis)
   real(kind=r8) :: evec(bh%n_lrow,bh%n_lcol)

   write(*,"(A)") "**Error! An EigenExa stub routine was called"
   stop

end subroutine

subroutine elsi_cleanup_eigenexa(ph)

   implicit none

   type(elsi_param_t) :: ph

end subroutine

end module ELSI_EIGENEXA
