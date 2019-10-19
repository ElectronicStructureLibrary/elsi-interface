! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide stub routines which are only compiled when the actual MAGMA is not
!! available.
!!
module ELSI_MAGMA

   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_PRECISION, only: r8

   implicit none

   private

   public :: elsi_init_magma
   public :: elsi_cleanup_magma
   public :: elsi_solve_magma

   interface elsi_solve_magma
      module procedure elsi_solve_magma_real
      module procedure elsi_solve_magma_cmplx
   end interface

contains

subroutine elsi_init_magma(ph)

   implicit none

   type(elsi_param_t) :: ph

end subroutine

subroutine elsi_solve_magma_real(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t) :: ph
   type(elsi_basic_t) :: bh
   real(kind=r8) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8) :: eval(ph%n_basis)
   real(kind=r8) :: evec(bh%n_lrow,bh%n_lcol)

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

subroutine elsi_solve_magma_cmplx(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t) :: ph
   type(elsi_basic_t) :: bh
   complex(kind=r8) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8) :: eval(ph%n_basis)
   complex(kind=r8) :: evec(bh%n_lrow,bh%n_lcol)

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

subroutine elsi_cleanup_magma(ph)

   implicit none

   type(elsi_param_t) :: ph

end subroutine

end module ELSI_MAGMA
