! Copyright (c) 2015-2021, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide stub routines which are only compiled when the actual BSEPACK is not
!! available.
!!
module ELSI_BSEPACK

   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_PRECISION, only: r8

   implicit none

   private

   public :: elsi_solve_bsepack

   interface elsi_solve_bsepack
      module procedure elsi_solve_bsepack_real
      module procedure elsi_solve_bsepack_cmplx
   end interface

contains

subroutine elsi_solve_bsepack_real(ph,bh,mat_a,mat_b,eval,evec)

   implicit none

   type(elsi_param_t) :: ph
   type(elsi_basic_t) :: bh
   real(kind=r8) :: mat_a(bh%n_lrow,bh%n_lcol)
   real(kind=r8) :: mat_b(bh%n_lrow,bh%n_lcol)
   real(kind=r8) :: eval(ph%n_basis)
   real(kind=r8) :: evec(ph%bse_n_lrow,ph%bse_n_lcol)

   write(*,"(A)") "**Error! A BSEPACK stub routine was called"
   stop

end subroutine

subroutine elsi_solve_bsepack_cmplx(ph,bh,mat_a,mat_b,eval,evec)

   implicit none

   type(elsi_param_t) :: ph
   type(elsi_basic_t) :: bh
   complex(kind=r8) :: mat_a(bh%n_lrow,bh%n_lcol)
   complex(kind=r8) :: mat_b(bh%n_lrow,bh%n_lcol)
   real(kind=r8) :: eval(ph%n_basis)
   complex(kind=r8) :: evec(ph%bse_n_lrow,ph%bse_n_lcol)

   write(*,"(A)") "**Error! A BSEPACK stub routine was called"
   stop

end subroutine

end module ELSI_BSEPACK
