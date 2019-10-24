! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide stub routines which are only compiled when the actual SLEPc-SIPs is
!! not available.
!!
module ELSI_SIPS

   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_PRECISION, only: r8,i4

   implicit none

   private

   public :: elsi_init_sips
   public :: elsi_cleanup_sips
   public :: elsi_solve_sips
   public :: elsi_build_dm_sips
   public :: elsi_build_edm_sips

   interface elsi_solve_sips
      module procedure elsi_solve_sips_real
   end interface

   interface elsi_build_dm_sips
      module procedure elsi_build_dm_sips_real
   end interface

   interface elsi_build_edm_sips
      module procedure elsi_build_edm_sips_real
   end interface

contains

subroutine elsi_init_sips(ph,bh)

   implicit none

   type(elsi_param_t) :: ph
   type(elsi_basic_t) :: bh

end subroutine

subroutine elsi_solve_sips_real(ph,bh,row_ind,col_ptr,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t) :: ph
   type(elsi_basic_t) :: bh
   integer(kind=i4) :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4) :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8) :: ham(bh%nnz_l_sp1)
   real(kind=r8) :: ovlp(bh%nnz_l_sp1)
   real(kind=r8) :: eval(ph%n_states)
   real(kind=r8) :: evec(bh%n_lcol_sp1,ph%n_states)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called"
   stop

end subroutine

subroutine elsi_build_dm_sips_real(ph,bh,row_ind,col_ptr,occ,dm)

   implicit none

   type(elsi_param_t) :: ph
   type(elsi_basic_t) :: bh
   integer(kind=i4) :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4) :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8) :: occ(ph%n_states)
   real(kind=r8) :: dm(bh%nnz_l_sp1)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called"
   stop

end subroutine

subroutine elsi_build_edm_sips_real(ph,bh,row_ind,col_ptr,occ,edm)

   implicit none

   type(elsi_param_t) :: ph
   type(elsi_basic_t) :: bh
   integer(kind=i4) :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4) :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8) :: occ(ph%n_states)
   real(kind=r8) :: edm(bh%nnz_l_sp1)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called"
   stop

end subroutine

subroutine elsi_cleanup_sips(ph)

   implicit none

   type(elsi_param_t) :: ph

end subroutine

end module ELSI_SIPS
