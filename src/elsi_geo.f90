! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Perform wavefunction and density matrix extrapolation between geometry steps.
!!
module ELSI_GEO

   use ELSI_CONSTANTS, only: ELPA_SOLVER
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_ELPA, only: elsi_update_dm_elpa
   use ELSI_NTPOLY, only: elsi_update_dm_ntpoly,CopyMatrix
   use ELSI_PRECISION, only: i4,r8
   use ELSI_REDIST, only: elsi_blacs_to_ntpoly_hs,elsi_ntpoly_to_blacs_dm
   use ELSI_UTIL, only: elsi_check_init,elsi_gram_schmidt

   implicit none

   private

   public :: elsi_orthonormalize_ev_real
   public :: elsi_orthonormalize_ev_complex
   public :: elsi_extrapolate_dm_real
   public :: elsi_extrapolate_dm_complex

contains

!>
!! Orthonormalize eigenvectors with respect to an overlap matrix.
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
!! Orthonormalize eigenvectors with respect to an overlap matrix.
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
!! Extrapolate density matrix for a new overlap.
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
!! Extrapolate density matrix for a new overlap.
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

end module ELSI_GEO
