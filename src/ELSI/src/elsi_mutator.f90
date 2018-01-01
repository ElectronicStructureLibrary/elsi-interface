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
!! This module provides routines for modifying ELSI and solver parameters.
!!
module ELSI_MUTATOR

   use ELSI_CONSTANTS, only: ELPA_SOLVER,OMM_SOLVER,PEXSI_SOLVER,CHESS_SOLVER,&
                             SIPS_SOLVER
   use ELSI_DATATYPE
   use ELSI_ELPA,      only: elsi_compute_edm_elpa_real,&
                             elsi_compute_edm_elpa_cmplx
   use ELSI_IO,        only: elsi_say
   use ELSI_MALLOC
   use ELSI_MPI,       only: elsi_stop
   use ELSI_MATCONV,   only: elsi_pexsi_to_blacs_dm_real,&
                             elsi_pexsi_to_blacs_dm_cmplx,&
                             elsi_blacs_to_sips_dm_real,&
                             elsi_blacs_to_sips_dm_cmplx
   use ELSI_OMM,       only: elsi_compute_edm_omm_real,&
                             elsi_compute_edm_omm_cmplx
   use ELSI_PEXSI,     only: elsi_compute_edm_pexsi_real,&
                             elsi_compute_edm_pexsi_cmplx
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS

   implicit none

   private

   ! Mutator
   public :: elsi_set_output
   public :: elsi_set_write_unit
   public :: elsi_set_unit_ovlp
   public :: elsi_set_zero_def
   public :: elsi_set_sing_check
   public :: elsi_set_sing_tol
   public :: elsi_set_sing_stop
   public :: elsi_set_uplo
   public :: elsi_set_elpa_solver
   public :: elsi_set_elpa_n_single
   public :: elsi_set_omm_flavor
   public :: elsi_set_omm_n_elpa
   public :: elsi_set_omm_tol
   public :: elsi_set_omm_ev_shift
   public :: elsi_set_omm_psp
   public :: elsi_set_pexsi_n_mu
   public :: elsi_set_pexsi_n_pole
   public :: elsi_set_pexsi_np_per_pole
   public :: elsi_set_pexsi_np_symbo
   public :: elsi_set_pexsi_ordering
   public :: elsi_set_pexsi_temp
   public :: elsi_set_pexsi_gap
   public :: elsi_set_pexsi_delta_e
   public :: elsi_set_pexsi_mu_min
   public :: elsi_set_pexsi_mu_max
   public :: elsi_set_pexsi_inertia_tol
   public :: elsi_set_chess_erf_decay
   public :: elsi_set_chess_erf_decay_min
   public :: elsi_set_chess_erf_decay_max
   public :: elsi_set_chess_ev_ham_min
   public :: elsi_set_chess_ev_ham_max
   public :: elsi_set_chess_ev_ovlp_min
   public :: elsi_set_chess_ev_ovlp_max
   public :: elsi_set_sips_n_elpa
   public :: elsi_set_sips_slice_type
   public :: elsi_set_sips_n_slice
   public :: elsi_set_sips_inertia
   public :: elsi_set_sips_left_bound
   public :: elsi_set_sips_slice_buf
   public :: elsi_set_sips_ev_min
   public :: elsi_set_sips_ev_max
   public :: elsi_set_dmp_method
   public :: elsi_set_dmp_max_step
   public :: elsi_set_dmp_tol
   public :: elsi_set_mu_broaden_scheme
   public :: elsi_set_mu_broaden_width
   public :: elsi_set_mu_tol
   public :: elsi_set_mu_spin_degen
   public :: elsi_set_output_solver_timings
   public :: elsi_set_solver_timings_unit
   public :: elsi_set_solver_timings_file
   public :: elsi_set_solver_timing_tag
   public :: elsi_get_pexsi_mu_min
   public :: elsi_get_pexsi_mu_max
   public :: elsi_get_ovlp_sing
   public :: elsi_get_n_sing
   public :: elsi_get_mu
   public :: elsi_get_edm_real
   public :: elsi_get_edm_complex
   public :: elsi_get_edm_real_sparse
   public :: elsi_get_edm_complex_sparse

contains

!>
!! This routine sets ELSI output level.
!!
subroutine elsi_set_output(e_h,out_level)

   implicit none

   type(elsi_handle), intent(inout) :: e_h       !< Handle
   integer(kind=i4),  intent(in)    :: out_level !< Output level

   character*40, parameter :: caller = "elsi_set_output"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(out_level <= 0) then
      e_h%stdio%print_info        = .false.
      e_h%print_mem               = .false.
      e_h%omm_output              = .false.
      e_h%pexsi_options%verbosity = 1
      e_h%elpa_output             = .false.
   elseif(out_level == 1) then
      e_h%stdio%print_info        = .true.
      e_h%print_mem               = .false.
      e_h%omm_output              = .false.
      e_h%pexsi_options%verbosity = 1
      e_h%elpa_output             = .false.
   elseif(out_level == 2) then
      e_h%stdio%print_info        = .true.
      e_h%print_mem               = .false.
      e_h%omm_output              = .true.
      e_h%pexsi_options%verbosity = 2
      e_h%elpa_output             = .true.
   else
      e_h%stdio%print_info        = .true.
      e_h%print_mem               = .true.
      e_h%omm_output              = .true.
      e_h%pexsi_options%verbosity = 2
      e_h%elpa_output             = .true.
   endif

end subroutine

!>
!! This routine sets the unit to be used by ELSI output.
!!
subroutine elsi_set_write_unit(e_h,write_unit)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: write_unit !< Unit

   character*40, parameter :: caller = "elsi_set_write_unit"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%stdio%print_unit = write_unit

end subroutine

!>
!! This routine sets the overlap matrix to be identity.
!!
subroutine elsi_set_unit_ovlp(e_h,unit_ovlp)

   implicit none

   type(elsi_handle), intent(inout) :: e_h       !< Handle
   integer(kind=i4),  intent(in)    :: unit_ovlp !< Overlap is identity?

   character*40, parameter :: caller = "elsi_set_unit_ovlp"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(unit_ovlp == 0) then
      e_h%ovlp_is_unit = .false.
   else
      e_h%ovlp_is_unit = .true.
   endif

end subroutine

!>
!! This routine sets the threshold to define "zero".
!!
subroutine elsi_set_zero_def(e_h,zero_def)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   real(kind=r8),     intent(in)    :: zero_def !< Zero tolerance

   character*40, parameter :: caller = "elsi_set_zero_def"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%zero_def = zero_def

end subroutine

!>
!! This routine switches on/off the singularity check of the overlap matrix.
!!
subroutine elsi_set_sing_check(e_h,sing_check)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: sing_check !< Check singularity?

   character*40, parameter :: caller = "elsi_set_sing_check"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(sing_check == 0) then
      e_h%check_sing   = .false.
      e_h%ovlp_is_sing = .false.
      e_h%n_nonsing    = e_h%n_basis
   else
      e_h%check_sing = .true.
   endif

end subroutine

!>
!! This routine sets the tolerance of the singularity check.
!!
subroutine elsi_set_sing_tol(e_h,sing_tol)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   real(kind=r8),     intent(in)    :: sing_tol !< Singularity tolerance

   character*40, parameter :: caller = "elsi_set_sing_tol"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%sing_tol = sing_tol

end subroutine

!>
!! This routine sets whether to stop in case of singular overlap matrix.
!!
subroutine elsi_set_sing_stop(e_h,sing_stop)

   implicit none

   type(elsi_handle), intent(inout) :: e_h       !< Handle
   integer(kind=i4),  intent(in)    :: sing_stop !< Stop if overlap is singular

   character*40, parameter :: caller = "elsi_set_sing_stop"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(sing_stop == 0) then
      e_h%stop_sing = .false.
   else
      e_h%stop_sing = .true.
   endif

end subroutine

!>
!! This routine sets the input matrices to be full, upper triangular, or lower
!! triangular.
!!
subroutine elsi_set_uplo(e_h,uplo)

   implicit none

   type(elsi_handle), intent(inout) :: e_h  !< Handle
   integer(kind=i4),  intent(in)    :: uplo !< Input matrix triangular?

   character*40, parameter :: caller = "elsi_set_uplo"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%uplo = uplo

end subroutine

!>
!! This routine sets the ELPA solver.
!!
subroutine elsi_set_elpa_solver(e_h,elpa_solver)

   implicit none

   type(elsi_handle), intent(inout) :: e_h         !< Handle
   integer(kind=i4),  intent(in)    :: elpa_solver !< ELPA solver

   character*40, parameter :: caller = "elsi_set_elpa_solver"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(elpa_solver < 1 .or. elpa_solver > 2) then
      call elsi_stop(" Unsupported elpa_solver.",e_h,caller)
   endif

   e_h%elpa_solver = elpa_solver

end subroutine

!>
!! This routine sets the number of steps using single precision ELPA.
!!
subroutine elsi_set_elpa_n_single(e_h,n_single)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   integer(kind=i4),  intent(in)    :: n_single !< Single precision steps

   character*40, parameter :: caller = "elsi_set_elpa_n_single"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%n_single_steps = n_single

end subroutine

!>
!! This routine sets the flavor of libOMM.
!!
subroutine elsi_set_omm_flavor(e_h,omm_flavor)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: omm_flavor !< libOMM flavor

   character*40, parameter :: caller = "elsi_set_omm_flavor"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(omm_flavor /= 0 .and. omm_flavor /= 2) then
      call elsi_stop(" Unsupported omm_flavor.",e_h,caller)
   endif

   e_h%omm_flavor = omm_flavor

end subroutine

!>
!! This routine sets the number of ELPA steps when using libOMM.
!!
subroutine elsi_set_omm_n_elpa(e_h,n_elpa)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   integer(kind=i4),  intent(in)    :: n_elpa !< ELPA steps

   character*40, parameter :: caller = "elsi_set_omm_n_elpa"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%omm_n_elpa = n_elpa

end subroutine

!>
!! This routine sets the tolerance of OMM minimization.
!!
subroutine elsi_set_omm_tol(e_h,min_tol)

   implicit none

   type(elsi_handle), intent(inout) :: e_h     !< Handle
   real(kind=r8),     intent(in)    :: min_tol !< Tolerance of OMM minimization

   character*40, parameter :: caller = "elsi_set_omm_tol"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%min_tol = min_tol

end subroutine

!>
!! This routine sets the shift of the eigenspectrum.
!!
subroutine elsi_set_omm_ev_shift(e_h,ev_shift)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   real(kind=r8),     intent(in)    :: ev_shift !< Shift of the eigenspectrum

   character*40, parameter :: caller = "elsi_set_omm_ev_shift"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%eta = ev_shift

end subroutine

!>
!! This routine switches on/off the matrix multiplications using PSP.
!!
subroutine elsi_set_omm_psp(e_h,use_psp)

   implicit none

   type(elsi_handle), intent(inout) :: e_h     !< Handle
   integer(kind=i4),  intent(in)    :: use_psp !< Use pspBLAS?

   character*40, parameter :: caller = "elsi_set_omm_psp"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(use_psp == 0) then
      e_h%use_psp = .false.
   else
      e_h%use_psp = .false.
   endif

end subroutine

!>
!! This routine sets the number of mu points when using PEXSI driver 2.
!!
subroutine elsi_set_pexsi_n_mu(e_h,n_mu)

   implicit none

   type(elsi_handle), intent(inout) :: e_h  !< Handle
   integer(kind=i4),  intent(in)    :: n_mu !< Number of mu points

   character*40, parameter :: caller = "elsi_set_pexsi_n_mu"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(n_mu < 1) then
      call elsi_stop(" Number of mu points should be at least 1.",e_h,caller)
   endif

   e_h%pexsi_options%nPoints = n_mu

end subroutine

!>
!! This routine sets the number of poles in the pole expansion.
!!
subroutine elsi_set_pexsi_n_pole(e_h,n_pole)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   integer(kind=i4),  intent(in)    :: n_pole !< Number of poles

   character*40, parameter :: caller = "elsi_set_pexsi_n_pole"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(n_pole < 1) then
      call elsi_stop(" Number of poles should be at least 1.",e_h,caller)
   endif

   e_h%pexsi_options%numPole = n_pole

end subroutine

!>
!! This routine sets the number of MPI tasks assigned for one pole.
!!
subroutine elsi_set_pexsi_np_per_pole(e_h,np_per_pole)

   implicit none

   type(elsi_handle), intent(inout) :: e_h         !< Handle
   integer(kind=i4),  intent(in)    :: np_per_pole !< Number of tasks per pole

   character*40, parameter :: caller = "elsi_set_pexsi_np_per_pole"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%np_per_pole = np_per_pole

end subroutine

!>
!! This routine sets the number of MPI tasks for the symbolic factorization.
!!
subroutine elsi_set_pexsi_np_symbo(e_h,np_symbo)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   integer(kind=i4),  intent(in)    :: np_symbo !< Number of tasks for symbolic

   character*40, parameter :: caller = "elsi_set_pexsi_np_symbo"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(np_symbo < 1) then
      call elsi_stop(" Number of MPI tasks for symbolic factorization should"//&
              " be at least 1.",e_h,caller)
   endif

   e_h%pexsi_options%npSymbFact = np_symbo

end subroutine

!>
!! This routine sets the matrix reordering method.
!!
subroutine elsi_set_pexsi_ordering(e_h,ordering)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   integer(kind=i4),  intent(in)    :: ordering !< Matrix reordering method

   character*40, parameter :: caller = "elsi_set_pexsi_ordering"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%pexsi_options%ordering = ordering

end subroutine

!>
!! This routine sets the temperature parameter in PEXSI.
!!
subroutine elsi_set_pexsi_temp(e_h,temp)

   implicit none

   type(elsi_handle), intent(inout) :: e_h  !< Handle
   real(kind=r8),     intent(in)    :: temp !< Temperature

   character*40, parameter :: caller = "elsi_set_pexsi_temp"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%pexsi_options%temperature = temp

end subroutine

!>
!! This routine sets the spectral gap in PEXSI.
!!
subroutine elsi_set_pexsi_gap(e_h,gap)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle
   real(kind=r8),     intent(in)    :: gap !< Gap

   character*40, parameter :: caller = "elsi_set_pexsi_gap"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(gap < 0.0_r8) then
      call elsi_stop(" Gap cannot be negative.",e_h,caller)
   endif

   e_h%pexsi_options%gap = gap

end subroutine

!>
!! This routine sets the spectrum width in PEXSI.
!!
subroutine elsi_set_pexsi_delta_e(e_h,delta_e)

   implicit none

   type(elsi_handle), intent(inout) :: e_h     !< Handle
   real(kind=r8),     intent(in)    :: delta_e !< Spectrum width

   character*40, parameter :: caller = "elsi_set_pexsi_delta_e"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(delta_e < 0.0_r8) then
      call elsi_stop(" Spectrum width cannot be negative.",e_h,caller)
   endif

   e_h%pexsi_options%deltaE = delta_e

end subroutine

!>
!! This routine sets the lower bound of the chemical potential in PEXSI.
!!
subroutine elsi_set_pexsi_mu_min(e_h,mu_min)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(in)    :: mu_min !< Lower bound of mu

   character*40, parameter :: caller = "elsi_set_pexsi_mu_min"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%pexsi_options%muMin0 = mu_min

end subroutine

!>
!! This routine sets the upper bound of the chemical potential in PEXSI.
!!
subroutine elsi_set_pexsi_mu_max(e_h,mu_max)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(in)    :: mu_max !< Upper bound of mu

   character*40, parameter :: caller = "elsi_set_pexsi_mu_max"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%pexsi_options%muMax0 = mu_max

end subroutine

!>
!! This routine sets the tolerance of the estimation of the chemical potential
!! in the inertia counting procedure.
!!
subroutine elsi_set_pexsi_inertia_tol(e_h,inertia_tol)

   implicit none

   type(elsi_handle), intent(inout) :: e_h         !< Handle
   real(kind=r8),     intent(in)    :: inertia_tol !< Tolerance of inertia counting

   character*40, parameter :: caller = "elsi_set_pexsi_inertia_tol"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(inertia_tol < 0.0_r8) then
      call elsi_stop(" Inertia counting tolerance cannot be negative.",e_h,&
              caller)
   endif

   e_h%pexsi_options%muInertiaTolerance = inertia_tol

end subroutine

!>
!! This routine sets the initial guess of the error function decay length in
!! CheSS.
!!
subroutine elsi_set_chess_erf_decay(e_h,decay)

   implicit none

   type(elsi_handle), intent(inout) :: e_h   !< Handle
   real(kind=r8),     intent(in)    :: decay !< Decay length

   character*40, parameter :: caller = "elsi_set_chess_erf_decay"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%erf_decay = decay

end subroutine

!>
!! This routine sets the lower bound of the decay length in CheSS.
!!
subroutine elsi_set_chess_erf_decay_min(e_h,decay_min)

   implicit none

   type(elsi_handle), intent(inout) :: e_h       !< Handle
   real(kind=r8),     intent(in)    :: decay_min !< Minimum decay length

   character*40, parameter :: caller = "elsi_set_chess_erf_decay_min"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%erf_decay_min = decay_min

end subroutine

!>
!! This routine sets the upper bound of the decay length in CheSS.
!!
subroutine elsi_set_chess_erf_decay_max(e_h,decay_max)

   implicit none

   type(elsi_handle), intent(inout) :: e_h       !< Handle
   real(kind=r8),     intent(in)    :: decay_max !< Maximum decay length

   character*40, parameter :: caller = "elsi_set_chess_erf_decay_max"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%erf_decay_max = decay_max

end subroutine

!>
!! This routine sets the lower bound of the eigenvalues of the Hamiltonian.
!!
subroutine elsi_set_chess_ev_ham_min(e_h,ev_min)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(in)    :: ev_min !< Minimum eigenvalue

   character*40, parameter :: caller = "elsi_set_chess_ev_ham_min"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%ev_ham_min = ev_min

end subroutine

!>
!! This routine sets the upper bound of the eigenvalues of the Hamiltonian.
!!
subroutine elsi_set_chess_ev_ham_max(e_h,ev_max)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(in)    :: ev_max !< Maximum eigenvalue

   character*40, parameter :: caller = "elsi_set_chess_ev_ham_max"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%ev_ham_max = ev_max

end subroutine

!>
!! This routine sets the lower bound of the eigenvalues of the overlap matrix.
!!
subroutine elsi_set_chess_ev_ovlp_min(e_h,ev_min)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(in)    :: ev_min !< Minimum eigenvalue

   character*40, parameter :: caller = "elsi_set_chess_ev_ovlp_min"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%ev_ovlp_min = ev_min

end subroutine

!>
!! This routine sets the upper bound of the eigenvalues of the overlap matrix.
!!
subroutine elsi_set_chess_ev_ovlp_max(e_h,ev_max)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(in)    :: ev_max !< Maximum eigenvalue

   character*40, parameter :: caller = "elsi_set_chess_ev_ovlp_max"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%ev_ovlp_max = ev_max

end subroutine

!>
!! This routine sets the number of ELPA steps when using SIPs.
!!
subroutine elsi_set_sips_n_elpa(e_h,n_elpa)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   integer(kind=i4),  intent(in)    :: n_elpa !< ELPA steps

   character*40, parameter :: caller = "elsi_set_sips_n_elpa"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%sips_n_elpa = n_elpa

end subroutine

!>
!! This routine switches on and off inertia counting in SIPs.
!!
subroutine elsi_set_sips_inertia(e_h,do_inertia)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: do_inertia !< Inertia counting option

   character*40, parameter :: caller = "elsi_set_sips_inertia"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(do_inertia == 0) then
      e_h%inertia_option = 0
   else
      e_h%inertia_option = 1
   endif

end subroutine

!>
!! This routine sets the slicing method when using SIPs.
!!
subroutine elsi_set_sips_slice_type(e_h,slice_type)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: slice_type !< Method of slicing

   character*40, parameter :: caller = "elsi_set_sips_slice_type"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(slice_type < 2 .or. slice_type > 3) then
      call elsi_stop(" Unsupported slice_type.",e_h,caller)
   endif

   e_h%slicing_method = slice_type

end subroutine

!>
!! This routine sets the number of slices in SIPs.
!!
subroutine elsi_set_sips_n_slice(e_h,n_slice)

   implicit none

   type(elsi_handle), intent(inout) :: e_h     !< Handle
   integer(kind=i4),  intent(in)    :: n_slice !< Number of slices

   character*40, parameter :: caller = "elsi_set_sips_n_slice"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(mod(e_h%n_procs,n_slice) == 0) then
      e_h%n_slices     = n_slice
      e_h%np_per_slice = e_h%n_procs/n_slice
   else
      call elsi_stop(" The total number of MPI tasks must be a multiple of"//&
              " the number of slices.",e_h,caller)
   endif

end subroutine

!>
!! This routine sets the method to bound the left side of the eigenvalue
!! interval in SIPs.
!!
subroutine elsi_set_sips_left_bound(e_h,left_bound)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: left_bound !< How to bound left side?

   character*40, parameter :: caller = "elsi_set_sips_left_bound"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%unbound = left_bound

end subroutine

!>
!! This routine sets a small buffer to expand the eigenvalue interval in SIPs.
!!
subroutine elsi_set_sips_slice_buf(e_h,slice_buffer)

   implicit none

   type(elsi_handle), intent(inout) :: e_h          !< Handle
   real(kind=r8),     intent(in)    :: slice_buffer !< Buffer to expand interval

   character*40, parameter :: caller = "elsi_set_sips_slice_buf"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%slice_buffer = slice_buffer

end subroutine

!>
!! This routine sets the lower bound of eigenvalues to be solved by SIPs.
!!
subroutine elsi_set_sips_ev_min(e_h,ev_min)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(in)    :: ev_min !< Minimum eigenvalue

   character*40, parameter :: caller = "elsi_set_sips_ev_min"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%ev_min = ev_min

end subroutine

!>
!! This routine sets the upper bound of eigenvalues to be solved by SIPs.
!!
subroutine elsi_set_sips_ev_max(e_h,ev_max)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(in)    :: ev_max !< Maximum eigenvalue

   character*40, parameter :: caller = "elsi_set_sips_ev_max"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%ev_max = ev_max

end subroutine

!>
!! This routine sets the density matrix purification method.
!!
subroutine elsi_set_dmp_method(e_h,dmp_method)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: dmp_method !< Purification method

   character*40, parameter :: caller = "elsi_set_dmp_method"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%dmp_method = dmp_method

end subroutine

!>
!! This routine sets the maximum number of density matrix purification steps.
!!
subroutine elsi_set_dmp_max_step(e_h,max_step)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   integer(kind=i4),  intent(in)    :: max_step !< Maximum number of steps

   character*40, parameter :: caller = "elsi_set_dmp_max_step"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%max_dmp_iter = max_step

end subroutine

!>
!! This routine sets the tolerance of the density matrix purification.
!!
subroutine elsi_set_dmp_tol(e_h,dmp_tol)

   implicit none

   type(elsi_handle), intent(inout) :: e_h     !< Handle
   real(kind=r8),     intent(in)    :: dmp_tol !< Tolerance

   character*40, parameter :: caller = "elsi_set_dmp_tol"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%dmp_tol = dmp_tol

end subroutine

!>
!! This routine sets the broadening scheme to determine the chemical potential
!! and the occupation numbers.
!!
subroutine elsi_set_mu_broaden_scheme(e_h,broaden_scheme)

   implicit none

   type(elsi_handle), intent(inout) :: e_h            !< Handle
   integer(kind=i4),  intent(in)    :: broaden_scheme !< Broadening method

   character*40, parameter :: caller = "elsi_set_mu_broaden_method"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%broaden_scheme = broaden_scheme

end subroutine

!>
!! This routine sets the broadening width to determine the chemical potential
!! and the occupation numbers.
!!
subroutine elsi_set_mu_broaden_width(e_h,broaden_width)

   implicit none

   type(elsi_handle), intent(inout) :: e_h           !< Handle
   real(kind=r8),     intent(in)    :: broaden_width !< Broadening width

   character*40, parameter :: caller = "elsi_set_mu_broaden_width"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(broaden_width < 0.0_r8) then
      call elsi_stop(" Broadening width cannot be negative.",e_h,caller)
   endif

   e_h%broaden_width = broaden_width

end subroutine

!>
!! This routine sets the desired accuracy of the determination of the chemical
!! potential and the occupation numbers.
!!
subroutine elsi_set_mu_tol(e_h,mu_tol)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(in)    :: mu_tol !< Accuracy of mu

   character*40, parameter :: caller = "elsi_set_mu_tol"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(mu_tol < 0.0_r8) then
      call elsi_stop(" Occupation number accuracy cannot be negative.",e_h,&
              caller)
   endif

   e_h%occ_tolerance = mu_tol

end subroutine

!>
!! This routine sets the spin degeneracy in the determination of the chemical
!! potential and the occupation numbers.
!!
subroutine elsi_set_mu_spin_degen(e_h,spin_degen)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   real(kind=r8),     intent(in)    :: spin_degen !< Spin degeneracy

   character*40, parameter :: caller = "elsi_set_mu_spin_degen"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%spin_degen  = spin_degen
   e_h%spin_is_set = .true.

end subroutine

!>
!! This routine sets whether the detailed solver timings file should be output.
!!
subroutine elsi_set_output_solver_timings(e_h,output_solver_timings)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                   !< Handle
   integer(kind=i4),  intent(in)    :: output_solver_timings !< Output timings?

   character*40, parameter :: caller = "elsi_set_output_solver_timings"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   if(output_solver_timings == 0) then
      e_h%output_solver_timings = .false.
   else
      e_h%output_solver_timings = .true.
   endif

end subroutine

!>
!! This routine sets the unit to which detailed solver timings are output.
!!
subroutine elsi_set_solver_timings_unit(e_h,solver_timings_unit)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                 !< Handle
   integer(kind=i4),  intent(in)    :: solver_timings_unit !< Unit

   character*40, parameter :: caller = "elsi_set_solver_timings_unit"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%solver_timings_file%print_unit = solver_timings_unit

end subroutine

!>
!! This routine sets the file to which detailed solver timings are output.
!!
subroutine elsi_set_solver_timings_file(e_h,solver_timings_file)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                 !< Handle
   character(len=*),  intent(in)    :: solver_timings_file !< File

   character*40, parameter :: caller = "elsi_set_solver_timings_file"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%solver_timings_file%file_name = solver_timings_file

end subroutine

!>
!! This routine sets the next user_tag for the solver timings.
!!
subroutine elsi_set_solver_timing_tag(e_h,user_tag)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   character(len=*),  intent(in)    :: user_tag !< Tag

   character*40, parameter :: caller = "elsi_set_solver_timing_tag"

   e_h%solver_timings%next_user_tag = user_tag

end subroutine

!>
!! This routine gets the lower bound of the chemical potential returned by the
!! inertia counting in PEXSI.
!!
subroutine elsi_get_pexsi_mu_min(e_h,mu_min)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(out)   :: mu_min !< Lower bound of mu

   character*40, parameter :: caller = "elsi_get_pexsi_mu_min"

   call elsi_check_handle(e_h,caller)

   mu_min = e_h%pexsi_options%muMin0

end subroutine

!>
!! This routine gets the upper bound of the chemical potential returned by the
!! inertia counting in PEXSI.
!!
subroutine elsi_get_pexsi_mu_max(e_h,mu_max)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(out)   :: mu_max !< Upper bound of mu

   character*40, parameter :: caller = "elsi_get_pexsi_mu_max"

   call elsi_check_handle(e_h,caller)

   mu_max = e_h%pexsi_options%muMax0

end subroutine

!>
!! This routine gets the result of the singularity check of the overlap matrix.
!!
subroutine elsi_get_ovlp_sing(e_h,ovlp_sing)

   implicit none

   type(elsi_handle), intent(inout) :: e_h       !< Handle
   integer(kind=i4),  intent(out)   :: ovlp_sing !< Is overlap singular?

   character*40, parameter :: caller = "elsi_get_ovlp_sing"

   call elsi_check_handle(e_h,caller)

   if(e_h%ovlp_is_sing) then
      ovlp_sing = 1
   else
      ovlp_sing = 0
   endif

end subroutine

!>
!! This routine gets the number of basis functions that are removed due to
!! overlap singularity.
!!
subroutine elsi_get_n_sing(e_h,n_sing)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   integer(kind=i4),  intent(out)   :: n_sing !< Number of singular basis

   character*40, parameter :: caller = "elsi_get_n_sing"

   call elsi_check_handle(e_h,caller)

   n_sing = e_h%n_basis-e_h%n_nonsing

end subroutine

!>
!! This routine gets the chemical potential.
!!
subroutine elsi_get_mu(e_h,mu)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle
   real(kind=r8),     intent(out)   :: mu  !< Chemical potential

   character*40, parameter :: caller = "elsi_get_mu"

   call elsi_check_handle(e_h,caller)

   mu = e_h%mu

   if(.not. e_h%mu_ready) then
      call elsi_say(e_h,"  ATTENTION! The return value of mu may be 0, since"//&
              " mu has not been computed.")
   endif

   e_h%mu_ready = .false.

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_real(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                        !< Handle
   real(kind=r8),     intent(out)   :: edm(e_h%n_lrow,e_h%n_lcol) !< Energy density matrix

   real(kind=r8), allocatable :: tmp_real(:,:)

   character*40, parameter :: caller = "elsi_get_edm_real"

   call elsi_check_handle(e_h,caller)

   if(e_h%edm_ready_real) then
      select case(e_h%solver)
      case(ELPA_SOLVER)
         call elsi_allocate(e_h,tmp_real,e_h%n_lrow,e_h%n_lcol,"tmp_real",&
                 caller)
         call elsi_compute_edm_elpa_real(e_h,e_h%eval_elpa,e_h%evec_real_elpa,&
                 edm,tmp_real)
         call elsi_deallocate(e_h,tmp_real,"tmp_real")
      case(OMM_SOLVER)
         call elsi_compute_edm_omm_real(e_h,edm)
      case(PEXSI_SOLVER)
         call elsi_compute_edm_pexsi_real(e_h,e_h%dm_real_pexsi)
         call elsi_pexsi_to_blacs_dm_real(e_h,edm)
      case(CHESS_SOLVER)
         call elsi_stop(" CHESS not yet implemented.",e_h,caller)
      case(SIPS_SOLVER)
         call elsi_stop(" SIPS not yet implemented.",e_h,caller)
      case(DMP_SOLVER)
         call elsi_stop(" DMP not yet implemented.",e_h,caller)
      case default
         call elsi_stop(" Unsupported solver.",e_h,caller)
      end select

      e_h%edm_ready_real = .false.
   else
      call elsi_stop(" Energy-weighted density matrix not computed.",e_h,caller)
   endif

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_real_sparse(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h               !< Handle
   real(kind=r8),     intent(out)   :: edm(e_h%nnz_l_sp) !< Energy density matrix

   real(kind=r8), allocatable :: tmp_real(:,:)

   character*40, parameter :: caller = "elsi_get_edm_real_sparse"

   call elsi_check_handle(e_h,caller)

   if(e_h%edm_ready_real) then
      select case(e_h%solver)
      case(ELPA_SOLVER)
         call elsi_allocate(e_h,tmp_real,e_h%n_lrow,e_h%n_lcol,"tmp_real",&
                 caller)
         call elsi_compute_edm_elpa_real(e_h,e_h%eval_elpa,e_h%evec_real_elpa,&
                 e_h%dm_real_elpa,tmp_real)
         call elsi_deallocate(e_h,tmp_real,"tmp_real")
         call elsi_blacs_to_sips_dm_real(e_h,edm)
      case(OMM_SOLVER)
         call elsi_compute_edm_omm_real(e_h,e_h%dm_real_elpa)
         call elsi_blacs_to_sips_dm_real(e_h,edm)
      case(PEXSI_SOLVER)
         call elsi_compute_edm_pexsi_real(e_h,edm)
      case(CHESS_SOLVER)
         call elsi_stop(" CHESS not yet implemented.",e_h,caller)
      case(SIPS_SOLVER)
         call elsi_stop(" SIPS not yet implemented.",e_h,caller)
      case(DMP_SOLVER)
         call elsi_stop(" DMP not yet implemented.",e_h,caller)
      case default
         call elsi_stop(" Unsupported solver.",e_h,caller)
      end select

      e_h%edm_ready_real = .false.
   else
      call elsi_stop(" Energy-weighted density matrix not computed.",e_h,caller)
   endif

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_complex(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                        !< Handle
   complex(kind=r8),  intent(out)   :: edm(e_h%n_lrow,e_h%n_lcol) !< Energy density matrix
   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character*40, parameter :: caller = "elsi_get_edm_complex"

   call elsi_check_handle(e_h,caller)

   if(e_h%edm_ready_cmplx) then
      select case(e_h%solver)
      case(ELPA_SOLVER)
         call elsi_allocate(e_h,tmp_cmplx,e_h%n_lrow,e_h%n_lcol,"tmp_cmplx",&
                 caller)
         call elsi_compute_edm_elpa_cmplx(e_h,e_h%eval_elpa,&
                 e_h%evec_cmplx_elpa,edm,tmp_cmplx)
         call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")
      case(OMM_SOLVER)
         call elsi_compute_edm_omm_cmplx(e_h,edm)
      case(PEXSI_SOLVER)
         call elsi_compute_edm_pexsi_cmplx(e_h,e_h%dm_cmplx_pexsi)
         call elsi_pexsi_to_blacs_dm_cmplx(e_h,edm)
      case(CHESS_SOLVER)
         call elsi_stop(" CHESS not yet implemented.",e_h,caller)
      case(SIPS_SOLVER)
         call elsi_stop(" SIPS not yet implemented.",e_h,caller)
      case(DMP_SOLVER)
         call elsi_stop(" DMP not yet implemented.",e_h,caller)
      case default
         call elsi_stop(" Unsupported solver.",e_h,caller)
      end select

      e_h%edm_ready_cmplx = .false.
   else
      call elsi_stop(" Energy-weighted density matrix not computed.",e_h,caller)
   endif

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_complex_sparse(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h               !< Handle
   complex(kind=r8),  intent(out)   :: edm(e_h%nnz_l_sp) !< Energy density matrix

   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character*40, parameter :: caller = "elsi_get_edm_complex_sparse"

   call elsi_check_handle(e_h,caller)

   if(e_h%edm_ready_cmplx) then
      select case(e_h%solver)
      case(ELPA_SOLVER)
         call elsi_allocate(e_h,tmp_cmplx,e_h%n_lrow,e_h%n_lcol,"tmp_cmplx",&
                 caller)
         call elsi_compute_edm_elpa_cmplx(e_h,e_h%eval_elpa,&
                 e_h%evec_cmplx_elpa,e_h%dm_cmplx_elpa,tmp_cmplx)
         call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")
         call elsi_blacs_to_sips_dm_cmplx(e_h,edm)
      case(OMM_SOLVER)
         call elsi_compute_edm_omm_cmplx(e_h,e_h%dm_cmplx_elpa)
         call elsi_blacs_to_sips_dm_cmplx(e_h,edm)
      case(PEXSI_SOLVER)
         call elsi_compute_edm_pexsi_cmplx(e_h,edm)
      case(CHESS_SOLVER)
         call elsi_stop(" CHESS not yet implemented.",e_h,caller)
      case(SIPS_SOLVER)
         call elsi_stop(" SIPS not yet implemented.",e_h,caller)
      case(DMP_SOLVER)
         call elsi_stop(" DMP not yet implemented.",e_h,caller)
      case default
         call elsi_stop(" Unsupported solver.",e_h,caller)
      end select

      e_h%edm_ready_cmplx = .false.
   else
      call elsi_stop(" Energy-weighted density matrix not computed.",e_h,caller)
   endif

end subroutine

end module ELSI_MUTATOR
