! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains wrappers of internal helper routines in order to make
!! these routines easily usable from outside ELSI.
!!
module ELSI_TOOLS

   use ELSI_CONSTANTS, only: ELPA_SOLVER
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_ELPA, only: elsi_update_dm_elpa
   use ELSI_NTPOLY, only: elsi_update_dm_ntpoly,CopyMatrix
   use ELSI_OCC, only: elsi_mu_and_occ,elsi_entropy
   use ELSI_PRECISION, only: i4,r8
   use ELSI_REDIST, only: elsi_blacs_to_ntpoly_hs,elsi_ntpoly_to_blacs_dm
   use ELSI_UTILS, only: elsi_check_init,elsi_gram_schmidt

   implicit none

   private

   public :: elsi_orthonormalize_ev_real
   public :: elsi_orthonormalize_ev_complex
   public :: elsi_extrapolate_dm_real
   public :: elsi_extrapolate_dm_complex
   public :: elsi_compute_mu_and_occ
   public :: elsi_compute_entropy

contains

!>
!! This routine orthonormalizes eigenvectors with respect to an overlap matrix.
!!
subroutine elsi_orthonormalize_ev_real(eh,ovlp,evec)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   real(kind=r8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   character(len=*), parameter :: caller = "elsi_orthonormalize_ev_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_gram_schmidt(eh%ph,eh%bh,ovlp,evec)

end subroutine

!>
!! This routine orthonormalizes eigenvectors with respect to an overlap matrix.
!!
subroutine elsi_orthonormalize_ev_complex(eh,ovlp,evec)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   complex(kind=r8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   character(len=*), parameter :: caller = "elsi_orthonormalize_ev_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_gram_schmidt(eh%ph,eh%bh,ovlp,evec)

end subroutine

!>
!! This routine extrapolates density matrix for a new overlap.
!!
subroutine elsi_extrapolate_dm_real(eh,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< New overlap
   real(kind=r8), intent(inout) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix

   character(len=*), parameter :: caller = "elsi_extrapolate_dm_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_update_dm_elpa(eh%ph,eh%bh,eh%ovlp_real_copy,ovlp,dm)

      eh%ovlp_real_copy = ovlp
   case default
      eh%ph%unit_ovlp = .true.

      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,dm,ovlp,eh%ph%nt_ham,eh%ph%nt_dm)

      eh%ph%first_blacs_to_ntpoly = .true.

      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dm,eh%ph%nt_ovlp,&
           eh%ph%nt_dm)

      eh%ph%unit_ovlp = .false.

      call elsi_update_dm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ovlp_copy,eh%ph%nt_ovlp,&
           eh%ph%nt_ham,eh%ph%nt_dm)
      call elsi_ntpoly_to_blacs_dm(eh%ph,eh%bh,eh%ph%nt_dm,dm)
      call CopyMatrix(eh%ph%nt_ovlp,eh%ph%nt_ovlp_copy)
   end select

end subroutine

!>
!! This routine extrapolates density matrix for a new overlap.
!!
subroutine elsi_extrapolate_dm_complex(eh,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< New overlap
   complex(kind=r8), intent(inout) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix

   character(len=*), parameter :: caller = "elsi_extrapolate_dm_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_update_dm_elpa(eh%ph,eh%bh,eh%ovlp_cmplx_copy,ovlp,dm)

      eh%ovlp_cmplx_copy = ovlp
   case default
      eh%ph%unit_ovlp = .true.

      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,dm,ovlp,eh%ph%nt_ham,eh%ph%nt_dm)

      eh%ph%first_blacs_to_ntpoly = .true.

      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dm,eh%ph%nt_ovlp,&
           eh%ph%nt_dm)

      eh%ph%unit_ovlp = .false.

      call elsi_update_dm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ovlp_copy,eh%ph%nt_ovlp,&
           eh%ph%nt_ham,eh%ph%nt_dm)
      call elsi_ntpoly_to_blacs_dm(eh%ph,eh%bh,eh%ph%nt_dm,dm)
      call CopyMatrix(eh%ph%nt_ovlp,eh%ph%nt_ovlp_copy)
   end select

end subroutine

!>
!! This routine computes the chemical potential and occupation numbers.
!!
subroutine elsi_compute_mu_and_occ(eh,n_electron,n_state,n_spin,n_kpt,k_wt,&
   eval,occ,mu)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   real(kind=r8), intent(in) :: n_electron !< Number of electrons
   integer(kind=i4), intent(in) :: n_state !< Number of states
   integer(kind=i4), intent(in) :: n_spin !< Number of spins
   integer(kind=i4), intent(in) :: n_kpt !< Number of k-points
   real(kind=r8), intent(in) :: k_wt(n_kpt) !< K-points weights
   real(kind=r8), intent(in) :: eval(n_state,n_spin,n_kpt) !< Eigenvalues
   real(kind=r8), intent(out) :: occ(n_state,n_spin,n_kpt) !< Occupation members
   real(kind=r8), intent(out) :: mu !< Chemical potential

   character(len=*), parameter :: caller = "elsi_compute_mu_and_occ"

   call elsi_mu_and_occ(eh%ph,eh%bh,n_electron,n_state,n_spin,n_kpt,k_wt,eval,&
        occ,mu)

end subroutine

!>
!! This routine computes the entropy.
!!
subroutine elsi_compute_entropy(eh,n_state,n_spin,n_kpt,k_wt,eval,occ,mu,ts)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_state !< Number of states
   integer(kind=i4), intent(in) :: n_spin !< Number of spins
   integer(kind=i4), intent(in) :: n_kpt !< Number of k-points
   real(kind=r8), intent(in) :: k_wt(n_kpt) !< K-points weights
   real(kind=r8), intent(in) :: eval(n_state,n_spin,n_kpt) !< Eigenvalues
   real(kind=r8), intent(in) :: occ(n_state,n_spin,n_kpt) !< Occupation numbers
   real(kind=r8), intent(in) :: mu !< Input chemical potential
   real(kind=r8), intent(out) :: ts !< Entropy

   character(len=*), parameter :: caller = "elsi_compute_entropy"

   call elsi_entropy(eh%ph,n_state,n_spin,n_kpt,k_wt,eval,occ,mu,ts)

end subroutine

end module
