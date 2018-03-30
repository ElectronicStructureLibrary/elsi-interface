! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides routines for modifying ELSI and solver parameters.
!!
module ELSI_MUTATOR

   use ELSI_CONSTANTS,  only: ELPA_SOLVER,OMM_SOLVER,PEXSI_SOLVER,SIPS_SOLVER,&
                              DMP_SOLVER,PEXSI_CSC,SIESTA_CSC
   use ELSI_DATATYPE,   only: elsi_handle
   use ELSI_ELPA,       only: elsi_compute_edm_elpa_real,&
                              elsi_compute_edm_elpa_cmplx
   use ELSI_IO,         only: elsi_say
   use ELSI_MALLOC,     only: elsi_allocate,elsi_deallocate
   use ELSI_MAT_REDIST, only: elsi_blacs_to_siesta_dm_cmplx,&
                              elsi_blacs_to_siesta_dm_real,&
                              elsi_blacs_to_sips_dm_cmplx,&
                              elsi_blacs_to_sips_dm_real,&
                              elsi_pexsi_to_blacs_dm_cmplx,&
                              elsi_pexsi_to_blacs_dm_real,&
                              elsi_pexsi_to_siesta_dm_cmplx,&
                              elsi_pexsi_to_siesta_dm_real,&
                              elsi_sips_to_blacs_dm_real,&
                              elsi_sips_to_siesta_dm_real
   use ELSI_MPI,        only: elsi_stop
   use ELSI_OMM,        only: elsi_compute_edm_omm_real,&
                              elsi_compute_edm_omm_cmplx
   use ELSI_PEXSI,      only: elsi_compute_edm_pexsi_real,&
                              elsi_compute_edm_pexsi_cmplx
   use ELSI_PRECISION,  only: r8,i4
   use ELSI_SIPS,       only: elsi_compute_edm_sips_real
   use ELSI_UTILS,      only: elsi_check_handle

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
   public :: elsi_set_csc_blk
   public :: elsi_set_elpa_solver
   public :: elsi_set_omm_flavor
   public :: elsi_set_omm_n_elpa
   public :: elsi_set_omm_tol
   public :: elsi_set_omm_ev_shift
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
   public :: elsi_set_sips_n_elpa
   public :: elsi_set_sips_n_slice
   public :: elsi_set_sips_slice_type
   public :: elsi_set_sips_buffer
   public :: elsi_set_sips_inertia_tol
   public :: elsi_set_sips_interval
   public :: elsi_set_sips_first_ev
   public :: elsi_set_dmp_method
   public :: elsi_set_dmp_max_step
   public :: elsi_set_dmp_tol
   public :: elsi_set_mu_broaden_scheme
   public :: elsi_set_mu_broaden_width
   public :: elsi_set_mu_tol
   public :: elsi_set_mu_spin_degen
   public :: elsi_set_mu_mp_order
   public :: elsi_set_output_log
   public :: elsi_set_log_unit
   public :: elsi_set_log_file
   public :: elsi_set_log_tag
   public :: elsi_set_uuid
   public :: elsi_set_calling_code
   public :: elsi_get_pexsi_mu_min
   public :: elsi_get_pexsi_mu_max
   public :: elsi_get_n_sing
   public :: elsi_get_mu
   public :: elsi_get_entropy
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

   character(len=40), parameter :: caller = "elsi_set_output"

   call elsi_check_handle(e_h,caller)

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

   character(len=40), parameter :: caller = "elsi_set_write_unit"

   call elsi_check_handle(e_h,caller)

   e_h%stdio%print_unit = write_unit

end subroutine

!>
!! This routine sets the overlap matrix to be identity.
!!
subroutine elsi_set_unit_ovlp(e_h,unit_ovlp)

   implicit none

   type(elsi_handle), intent(inout) :: e_h       !< Handle
   integer(kind=i4),  intent(in)    :: unit_ovlp !< Overlap is identity?

   character(len=40), parameter :: caller = "elsi_set_unit_ovlp"

   call elsi_check_handle(e_h,caller)

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

   character(len=40), parameter :: caller = "elsi_set_zero_def"

   call elsi_check_handle(e_h,caller)

   e_h%zero_def = zero_def

end subroutine

!>
!! This routine switches on/off the singularity check of the overlap matrix.
!!
subroutine elsi_set_sing_check(e_h,sing_check)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: sing_check !< Check singularity?

   character(len=40), parameter :: caller = "elsi_set_sing_check"

   call elsi_check_handle(e_h,caller)

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

   character(len=40), parameter :: caller = "elsi_set_sing_tol"

   call elsi_check_handle(e_h,caller)

   e_h%sing_tol = sing_tol

end subroutine

!>
!! This routine sets whether to stop in case of singular overlap matrix.
!!
subroutine elsi_set_sing_stop(e_h,sing_stop)

   implicit none

   type(elsi_handle), intent(inout) :: e_h       !< Handle
   integer(kind=i4),  intent(in)    :: sing_stop !< Stop if overlap is singular

   character(len=40), parameter :: caller = "elsi_set_sing_stop"

   call elsi_check_handle(e_h,caller)

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

   character(len=40), parameter :: caller = "elsi_set_uplo"

   call elsi_check_handle(e_h,caller)

   e_h%uplo = uplo

end subroutine

!>
!! This routine sets the block size in 1D block-cyclic distributed CSC format.
!!
subroutine elsi_set_csc_blk(e_h,blk)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle
   integer(kind=i4),  intent(in)    :: blk !< Block size

   character(len=40), parameter :: caller = "elsi_set_csc_blk"

   call elsi_check_handle(e_h,caller)

   e_h%blk_sp2 = blk

end subroutine

!>
!! This routine sets the ELPA solver.
!!
subroutine elsi_set_elpa_solver(e_h,elpa_solver)

   implicit none

   type(elsi_handle), intent(inout) :: e_h         !< Handle
   integer(kind=i4),  intent(in)    :: elpa_solver !< ELPA solver

   character(len=40), parameter :: caller = "elsi_set_elpa_solver"

   call elsi_check_handle(e_h,caller)

   if(elpa_solver < 1 .or. elpa_solver > 2) then
      call elsi_stop(e_h,"Unsupported elpa_solver.",caller)
   endif

   e_h%elpa_solver = elpa_solver

end subroutine

!>
!! This routine sets the flavor of libOMM.
!!
subroutine elsi_set_omm_flavor(e_h,omm_flavor)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: omm_flavor !< libOMM flavor

   character(len=40), parameter :: caller = "elsi_set_omm_flavor"

   call elsi_check_handle(e_h,caller)

   if(omm_flavor /= 0 .and. omm_flavor /= 2) then
      call elsi_stop(e_h,"Unsupported omm_flavor.",caller)
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

   character(len=40), parameter :: caller = "elsi_set_omm_n_elpa"

   call elsi_check_handle(e_h,caller)

   e_h%omm_n_elpa = n_elpa

end subroutine

!>
!! This routine sets the tolerance of OMM minimization.
!!
subroutine elsi_set_omm_tol(e_h,min_tol)

   implicit none

   type(elsi_handle), intent(inout) :: e_h     !< Handle
   real(kind=r8),     intent(in)    :: min_tol !< Tolerance of OMM minimization

   character(len=40), parameter :: caller = "elsi_set_omm_tol"

   call elsi_check_handle(e_h,caller)

   e_h%omm_tol = min_tol

end subroutine

!>
!! This routine sets the shift of the eigenspectrum.
!!
subroutine elsi_set_omm_ev_shift(e_h,ev_shift)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   real(kind=r8),     intent(in)    :: ev_shift !< Shift of the eigenspectrum

   character(len=40), parameter :: caller = "elsi_set_omm_ev_shift"

   call elsi_check_handle(e_h,caller)

   e_h%omm_ev_shift = ev_shift

end subroutine

!>
!! This routine sets the number of mu points when using PEXSI driver 2.
!!
subroutine elsi_set_pexsi_n_mu(e_h,n_mu)

   implicit none

   type(elsi_handle), intent(inout) :: e_h  !< Handle
   integer(kind=i4),  intent(in)    :: n_mu !< Number of mu points

   character(len=40), parameter :: caller = "elsi_set_pexsi_n_mu"

   call elsi_check_handle(e_h,caller)

   if(n_mu < 1) then
      call elsi_stop(e_h,"Number of mu points should be at least 1.",caller)
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

   character(len=40), parameter :: caller = "elsi_set_pexsi_n_pole"

   call elsi_check_handle(e_h,caller)

   if(n_pole < 1) then
      call elsi_stop(e_h,"Number of poles should be at least 1.",caller)
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

   character(len=40), parameter :: caller = "elsi_set_pexsi_np_per_pole"

   call elsi_check_handle(e_h,caller)

   e_h%pexsi_np_per_pole = np_per_pole

end subroutine

!>
!! This routine sets the number of MPI tasks for the symbolic factorization.
!!
subroutine elsi_set_pexsi_np_symbo(e_h,np_symbo)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   integer(kind=i4),  intent(in)    :: np_symbo !< Number of tasks for symbolic

   character(len=40), parameter :: caller = "elsi_set_pexsi_np_symbo"

   call elsi_check_handle(e_h,caller)

   if(np_symbo < 1) then
      call elsi_stop(e_h,"Number of MPI tasks for symbolic factorization"//&
              " should be at least 1.",caller)
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

   character(len=40), parameter :: caller = "elsi_set_pexsi_ordering"

   call elsi_check_handle(e_h,caller)

   e_h%pexsi_options%ordering = ordering

end subroutine

!>
!! This routine sets the temperature parameter in PEXSI.
!!
subroutine elsi_set_pexsi_temp(e_h,temp)

   implicit none

   type(elsi_handle), intent(inout) :: e_h  !< Handle
   real(kind=r8),     intent(in)    :: temp !< Temperature

   character(len=40), parameter :: caller = "elsi_set_pexsi_temp"

   call elsi_check_handle(e_h,caller)

   e_h%pexsi_options%temperature = temp

end subroutine

!>
!! This routine sets the spectral gap in PEXSI.
!!
subroutine elsi_set_pexsi_gap(e_h,gap)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle
   real(kind=r8),     intent(in)    :: gap !< Gap

   character(len=40), parameter :: caller = "elsi_set_pexsi_gap"

   call elsi_check_handle(e_h,caller)

   if(gap < 0.0_r8) then
      call elsi_stop(e_h,"Gap cannot be negative.",caller)
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

   character(len=40), parameter :: caller = "elsi_set_pexsi_delta_e"

   call elsi_check_handle(e_h,caller)

   if(delta_e < 0.0_r8) then
      call elsi_stop(e_h,"Spectrum width cannot be negative.",caller)
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

   character(len=40), parameter :: caller = "elsi_set_pexsi_mu_min"

   call elsi_check_handle(e_h,caller)

   e_h%pexsi_options%muMin0 = mu_min

end subroutine

!>
!! This routine sets the upper bound of the chemical potential in PEXSI.
!!
subroutine elsi_set_pexsi_mu_max(e_h,mu_max)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(in)    :: mu_max !< Upper bound of mu

   character(len=40), parameter :: caller = "elsi_set_pexsi_mu_max"

   call elsi_check_handle(e_h,caller)

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

   character(len=40), parameter :: caller = "elsi_set_pexsi_inertia_tol"

   call elsi_check_handle(e_h,caller)

   if(inertia_tol < 0.0_r8) then
      call elsi_stop(e_h,"Inertia counting tolerance cannot be negative.",&
              caller)
   endif

   e_h%pexsi_options%muInertiaTolerance = inertia_tol

end subroutine

!>
!! This routine sets the number of ELPA steps when using SIPS.
!!
subroutine elsi_set_sips_n_elpa(e_h,n_elpa)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   integer(kind=i4),  intent(in)    :: n_elpa !< ELPA steps

   character(len=40), parameter :: caller = "elsi_set_sips_n_elpa"

   call elsi_check_handle(e_h,caller)

   e_h%sips_n_elpa = n_elpa

end subroutine

!>
!! This routine sets the number of slices in SIPS.
!!
subroutine elsi_set_sips_n_slice(e_h,n_slice)

   implicit none

   type(elsi_handle), intent(inout) :: e_h     !< Handle
   integer(kind=i4),  intent(in)    :: n_slice !< Number of slices

   character(len=40), parameter :: caller = "elsi_set_sips_n_slice"

   call elsi_check_handle(e_h,caller)

   if(mod(e_h%n_procs,n_slice) == 0) then
      e_h%sips_n_slices     = n_slice
      e_h%sips_np_per_slice = e_h%n_procs/n_slice
   else
      call elsi_stop(e_h," The total number of MPI tasks must be a multiple"//&
              " of the number of slices.",caller)
   endif

end subroutine

!>
!! This routine sets the type of slices to be used in SIPS.
!!
subroutine elsi_set_sips_slice_type(e_h,slice_type)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: slice_type !< Slice type

   character(len=40), parameter :: caller = "elsi_set_sips_slice_type"

   call elsi_check_handle(e_h,caller)

   if(slice_type /= 0 .and. slice_type /= 2 .and. slice_type /= 4) then
      call elsi_stop(e_h,"Unsupported slice type.",caller)
   endif

   e_h%sips_slice_type = slice_type

end subroutine

!>
!! This routine sets a small buffer to expand the eigenvalue interval in SIPS.
!!
subroutine elsi_set_sips_buffer(e_h,buffer)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(in)    :: buffer !< Buffer to expand interval

   character(len=40), parameter :: caller = "elsi_set_sips_buffer"

   call elsi_check_handle(e_h,caller)

   e_h%sips_buffer = buffer

end subroutine

!>
!! This routine sets the tolerance to stop inertia counting in SIPS.
!!
subroutine elsi_set_sips_inertia_tol(e_h,inertia_tol)

   implicit none

   type(elsi_handle), intent(inout) :: e_h         !< Handle
   real(kind=r8),     intent(in)    :: inertia_tol !< Stopping criterion

   character(len=40), parameter :: caller = "elsi_set_sips_inertia_tol"

   call elsi_check_handle(e_h,caller)

   e_h%sips_inertia_tol = inertia_tol

end subroutine

!>
!! This routine sets the global interval to be solved by SIPS.
!!
subroutine elsi_set_sips_interval(e_h,lower,upper)

   implicit none

   type(elsi_handle), intent(inout) :: e_h   !< Handle
   real(kind=r8),     intent(in)    :: lower !< Lower bound
   real(kind=r8),     intent(in)    :: upper !< Upper bound

   character(len=40), parameter :: caller = "elsi_set_sips_interval"

   call elsi_check_handle(e_h,caller)

   e_h%sips_interval(1) = lower
   e_h%sips_interval(2) = upper

end subroutine

!>
!! This routine sets the index of the first eigensolution to be solved.
!!
subroutine elsi_set_sips_first_ev(e_h,first_ev)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   integer(kind=i4),  intent(in)    :: first_ev !< Index of first eigensolution

   character(len=40), parameter :: caller = "elsi_set_sips_first_ev"

   call elsi_check_handle(e_h,caller)

   if(first_ev < 1) then
      e_h%sips_first_ev = 1
   elseif(first_ev > e_h%n_basis-e_h%n_states+1) then
      e_h%sips_first_ev = e_h%n_basis-e_h%n_states+1
   else
      e_h%sips_first_ev = first_ev
   endif

end subroutine

!>
!! This routine sets the density matrix purification method.
!!
subroutine elsi_set_dmp_method(e_h,dmp_method)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: dmp_method !< Purification method

   character(len=40), parameter :: caller = "elsi_set_dmp_method"

   call elsi_check_handle(e_h,caller)

   e_h%dmp_method = dmp_method

end subroutine

!>
!! This routine sets the maximum number of density matrix purification steps.
!!
subroutine elsi_set_dmp_max_step(e_h,max_step)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   integer(kind=i4),  intent(in)    :: max_step !< Maximum number of steps

   character(len=40), parameter :: caller = "elsi_set_dmp_max_step"

   call elsi_check_handle(e_h,caller)

   e_h%dmp_max_iter = max_step

end subroutine

!>
!! This routine sets the tolerance of the density matrix purification.
!!
subroutine elsi_set_dmp_tol(e_h,dmp_tol)

   implicit none

   type(elsi_handle), intent(inout) :: e_h     !< Handle
   real(kind=r8),     intent(in)    :: dmp_tol !< Tolerance

   character(len=40), parameter :: caller = "elsi_set_dmp_tol"

   call elsi_check_handle(e_h,caller)

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

   character(len=40), parameter :: caller = "elsi_set_mu_broaden_method"

   call elsi_check_handle(e_h,caller)

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

   character(len=40), parameter :: caller = "elsi_set_mu_broaden_width"

   call elsi_check_handle(e_h,caller)

   if(broaden_width < 0.0_r8) then
      call elsi_stop(e_h,"Broadening width cannot be negative.",caller)
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

   character(len=40), parameter :: caller = "elsi_set_mu_tol"

   call elsi_check_handle(e_h,caller)

   if(mu_tol < 0.0_r8) then
      call elsi_stop(e_h,"Occupation number accuracy cannot be negative.",&
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

   character(len=40), parameter :: caller = "elsi_set_mu_spin_degen"

   call elsi_check_handle(e_h,caller)

   e_h%spin_degen  = spin_degen
   e_h%spin_is_set = .true.

end subroutine

!>
!! This routine sets the order of the Methfessel-Paxton broadening scheme.
!!
subroutine elsi_set_mu_mp_order(e_h,mp_order)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   integer(kind=i4),  intent(in)    :: mp_order !< Order

   character(len=40), parameter :: caller = "elsi_set_mu_mp_order"

   call elsi_check_handle(e_h,caller)

   e_h%mp_order = mp_order

end subroutine

!>
!! This routine sets whether a log file should be output.
!!
subroutine elsi_set_output_log(e_h,output_log)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: output_log !< Output log

   character(len=40), parameter :: caller = "elsi_set_output_log"

   call elsi_check_handle(e_h,caller)

   if(output_log == 0) then
      e_h%log_file%print_info = .false.
   else
      e_h%log_file%print_info = .true.
   endif

end subroutine

!>
!! This routine sets the unit to which the log is output.
!!
subroutine elsi_set_log_unit(e_h,log_unit)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   integer(kind=i4),  intent(in)    :: log_unit !< Unit

   character(len=40), parameter :: caller = "elsi_set_log_unit"

   call elsi_check_handle(e_h,caller)

   e_h%log_file%print_unit = log_unit

end subroutine

!>
!! This routine sets the name of the log file.
!!
subroutine elsi_set_log_file(e_h,file_name)

   implicit none

   type(elsi_handle), intent(inout) :: e_h       !< Handle
   character(len=*),  intent(in)    :: file_name !< File name

   character(len=40), parameter :: caller = "elsi_set_log_file"

   e_h%log_file%file_name = file_name

end subroutine

!>
!! This routine sets the user_tag for the log.
!!
subroutine elsi_set_log_tag(e_h,user_tag)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   character(len=*),  intent(in)    :: user_tag !< Tag

   character(len=40), parameter :: caller = "elsi_set_log_tag"

   e_h%log_file%user_tag = user_tag

end subroutine

!>
!! This routine sets the UUID.
!!
subroutine elsi_set_uuid(e_h,uuid)

   implicit none

   type(elsi_handle), intent(inout) :: e_h  !< Handle
   character(len=*),  intent(in)    :: uuid !< UUID

   character(len=40), parameter :: caller = "elsi_set_uuid"

   e_h%uuid_exists = .true.
   e_h%uuid        = uuid

end subroutine

!>
!! This routine sets the name of the calling code.
!!
subroutine elsi_set_calling_code(e_h,code)

   implicit none

   type(elsi_handle), intent(inout) :: e_h  !< Handle
   character(len=*),  intent(in)    :: code !< Name of calling code

   character(len=40), parameter :: caller = "elsi_set_calling_code"

   e_h%caller = code

end subroutine

!>
!! This routine gets the lower bound of the chemical potential returned by the
!! inertia counting in PEXSI.
!!
subroutine elsi_get_pexsi_mu_min(e_h,mu_min)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(out)   :: mu_min !< Lower bound of mu

   character(len=40), parameter :: caller = "elsi_get_pexsi_mu_min"

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

   character(len=40), parameter :: caller = "elsi_get_pexsi_mu_max"

   call elsi_check_handle(e_h,caller)

   mu_max = e_h%pexsi_options%muMax0

end subroutine

!>
!! This routine gets the number of basis functions that are removed due to
!! overlap singularity.
!!
subroutine elsi_get_n_sing(e_h,n_sing)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   integer(kind=i4),  intent(out)   :: n_sing !< Number of singular basis

   character(len=40), parameter :: caller = "elsi_get_n_sing"

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

   character(len=40), parameter :: caller = "elsi_get_mu"

   call elsi_check_handle(e_h,caller)

   mu = e_h%mu

   if(.not. e_h%mu_ready) then
      call elsi_say(e_h,"  ATTENTION! The return value of mu may be 0, since"//&
              " mu has not been computed.",e_h%stdio)
   endif

   e_h%mu_ready = .false.

end subroutine

!>
!! This routine gets the entropy.
!!
subroutine elsi_get_entropy(e_h,ts)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle
   real(kind=r8),     intent(out)   :: ts  !< Entropy

   character(len=40), parameter :: caller = "elsi_get_entropy"

   call elsi_check_handle(e_h,caller)

   ts = e_h%ts

   if(.not. e_h%ts_ready) then
      call elsi_say(e_h,"  ATTENTION! The return value of ts may be 0, since"//&
              " ts has not been computed.",e_h%stdio)
   endif

   e_h%ts_ready = .false.

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_real(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                        !< Handle
   real(kind=r8),     intent(out)   :: edm(e_h%n_lrow,e_h%n_lcol) !< Energy density matrix

   integer(kind=i4) :: solver_save

   real(kind=r8), allocatable :: tmp_real(:,:)

   character(len=40), parameter :: caller = "elsi_get_edm_real"

   call elsi_check_handle(e_h,caller)

   solver_save = e_h%solver

   if(e_h%solver == OMM_SOLVER .and. e_h%n_elsi_calls <= e_h%omm_n_elpa) then
      solver_save = OMM_SOLVER
      e_h%solver  = ELPA_SOLVER
   endif

   if(e_h%solver == SIPS_SOLVER .and. e_h%n_elsi_calls <= e_h%sips_n_elpa) then
      solver_save = SIPS_SOLVER
      e_h%solver  = ELPA_SOLVER
   endif

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
      case(SIPS_SOLVER)
         call elsi_compute_edm_sips_real(e_h,e_h%dm_real_pexsi)
         call elsi_sips_to_blacs_dm_real(e_h,edm)
      case default
         call elsi_stop(e_h,"Unsupported density matrix solver.",caller)
      end select

      e_h%edm_ready_real = .false.
   else
      call elsi_stop(e_h,"Energy-weighted density matrix not computed.",caller)
   endif

   e_h%solver = solver_save

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_real_sparse(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h               !< Handle
   real(kind=r8),     intent(out)   :: edm(e_h%nnz_l_sp) !< Energy density matrix

   integer(kind=i4) :: solver_save

   real(kind=r8), allocatable :: tmp_real(:,:)

   character(len=40), parameter :: caller = "elsi_get_edm_real_sparse"

   call elsi_check_handle(e_h,caller)

   solver_save = e_h%solver

   if(e_h%solver == OMM_SOLVER .and. e_h%n_elsi_calls <= e_h%omm_n_elpa) then
      solver_save = OMM_SOLVER
      e_h%solver  = ELPA_SOLVER
   endif

   if(e_h%solver == SIPS_SOLVER .and. e_h%n_elsi_calls <= e_h%sips_n_elpa) then
      solver_save = SIPS_SOLVER
      e_h%solver  = ELPA_SOLVER
   endif

   if(e_h%edm_ready_real) then
      select case(e_h%solver)
      case(ELPA_SOLVER)
         call elsi_allocate(e_h,tmp_real,e_h%n_lrow,e_h%n_lcol,"tmp_real",&
                 caller)
         call elsi_compute_edm_elpa_real(e_h,e_h%eval_elpa,e_h%evec_real_elpa,&
                 e_h%dm_real_elpa,tmp_real)
         call elsi_deallocate(e_h,tmp_real,"tmp_real")

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_real(e_h,edm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_real(e_h,edm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select
      case(OMM_SOLVER)
         call elsi_compute_edm_omm_real(e_h,e_h%dm_real_elpa)

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_real(e_h,edm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_real(e_h,edm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select
      case(PEXSI_SOLVER)
         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_compute_edm_pexsi_real(e_h,edm)
         case(SIESTA_CSC)
            call elsi_compute_edm_pexsi_real(e_h,e_h%dm_real_pexsi)
            call elsi_pexsi_to_siesta_dm_real(e_h,edm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select
      case(SIPS_SOLVER)
         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_compute_edm_sips_real(e_h,edm)
         case(SIESTA_CSC)
            call elsi_compute_edm_sips_real(e_h,e_h%dm_real_pexsi)
            call elsi_sips_to_siesta_dm_real(e_h,edm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select
      case default
         call elsi_stop(e_h,"Unsupported density matrix solver.",caller)
      end select

      e_h%edm_ready_real = .false.
   else
      call elsi_stop(e_h,"Energy-weighted density matrix not computed.",caller)
   endif

   e_h%solver = solver_save

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_complex(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                        !< Handle
   complex(kind=r8),  intent(out)   :: edm(e_h%n_lrow,e_h%n_lcol) !< Energy density matrix

   integer(kind=i4) :: solver_save

   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character(len=40), parameter :: caller = "elsi_get_edm_complex"

   call elsi_check_handle(e_h,caller)

   solver_save = e_h%solver

   if(e_h%solver == OMM_SOLVER .and. e_h%n_elsi_calls <= e_h%omm_n_elpa) then
      solver_save = OMM_SOLVER
      e_h%solver  = ELPA_SOLVER
   endif

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
      case default
         call elsi_stop(e_h,"Unsupported density matrix solver.",caller)
      end select

      e_h%edm_ready_cmplx = .false.
   else
      call elsi_stop(e_h,"Energy-weighted density matrix not computed.",caller)
   endif

   e_h%solver = solver_save

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_complex_sparse(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h               !< Handle
   complex(kind=r8),  intent(out)   :: edm(e_h%nnz_l_sp) !< Energy density matrix

   integer(kind=i4) :: solver_save

   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character(len=40), parameter :: caller = "elsi_get_edm_complex_sparse"

   call elsi_check_handle(e_h,caller)

   solver_save = e_h%solver

   if(e_h%solver == OMM_SOLVER .and. e_h%n_elsi_calls <= e_h%omm_n_elpa) then
      solver_save = OMM_SOLVER
      e_h%solver  = ELPA_SOLVER
   endif

   if(e_h%edm_ready_cmplx) then
      select case(e_h%solver)
      case(ELPA_SOLVER)
         call elsi_allocate(e_h,tmp_cmplx,e_h%n_lrow,e_h%n_lcol,"tmp_cmplx",&
                 caller)
         call elsi_compute_edm_elpa_cmplx(e_h,e_h%eval_elpa,&
                 e_h%evec_cmplx_elpa,e_h%dm_cmplx_elpa,tmp_cmplx)
         call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_cmplx(e_h,edm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_cmplx(e_h,edm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select
      case(OMM_SOLVER)
         call elsi_compute_edm_omm_cmplx(e_h,e_h%dm_cmplx_elpa)

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_cmplx(e_h,edm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_cmplx(e_h,edm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select
      case(PEXSI_SOLVER)
         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_compute_edm_pexsi_cmplx(e_h,edm)
         case(SIESTA_CSC)
            call elsi_compute_edm_pexsi_cmplx(e_h,e_h%dm_cmplx_pexsi)
            call elsi_pexsi_to_siesta_dm_cmplx(e_h,edm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select
      case default
         call elsi_stop(e_h,"Unsupported density matrix solver.",caller)
      end select

      e_h%edm_ready_cmplx = .false.
   else
      call elsi_stop(e_h,"Energy-weighted density matrix not computed.",caller)
   endif

   e_h%solver = solver_save

end subroutine

end module ELSI_MUTATOR
