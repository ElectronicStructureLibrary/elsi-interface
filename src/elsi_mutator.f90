! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides routines for modifying ELSI and solver parameters.
!!
module ELSI_MUTATOR

   use ELSI_CONSTANTS, only: ELPA_SOLVER,OMM_SOLVER,PEXSI_SOLVER,SIPS_SOLVER,&
       NTPOLY_SOLVER,PEXSI_CSC,SIESTA_CSC
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_ELPA, only: elsi_compute_edm_elpa
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_REDIST, only: elsi_blacs_to_siesta_dm,elsi_blacs_to_sips_dm,&
       elsi_ntpoly_to_blacs_dm,elsi_ntpoly_to_siesta_dm,elsi_ntpoly_to_sips_dm,&
       elsi_pexsi_to_blacs_dm,elsi_pexsi_to_siesta_dm,elsi_sips_to_blacs_dm,&
       elsi_sips_to_siesta_dm
   use ELSI_MPI, only: elsi_stop
   use ELSI_NTPOLY, only: elsi_compute_edm_ntpoly
   use ELSI_OMM, only: elsi_compute_edm_omm
   use ELSI_PEXSI, only: elsi_compute_edm_pexsi
   use ELSI_PRECISION, only: r8,i4
   use ELSI_SIPS, only: elsi_compute_edm_sips
   use ELSI_UTILS, only: elsi_check_init

   implicit none

   private

   public :: elsi_set_output
   public :: elsi_set_write_unit
   public :: elsi_set_unit_ovlp
   public :: elsi_set_zero_def
   public :: elsi_set_sing_check
   public :: elsi_set_sing_tol
   public :: elsi_set_sing_stop
   public :: elsi_set_csc_blk
   public :: elsi_set_elpa_solver
   public :: elsi_set_elpa_n_single
   public :: elsi_set_elpa_gpu
   public :: elsi_set_elpa_gpu_kernels
   public :: elsi_set_elpa_autotune
   public :: elsi_set_omm_flavor
   public :: elsi_set_omm_n_elpa
   public :: elsi_set_omm_tol
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
   public :: elsi_set_ntpoly_method
   public :: elsi_set_ntpoly_isr
   public :: elsi_set_ntpoly_tol
   public :: elsi_set_ntpoly_filter
   public :: elsi_set_ntpoly_max_iter
   public :: elsi_set_mu_broaden_scheme
   public :: elsi_set_mu_broaden_width
   public :: elsi_set_mu_tol
   public :: elsi_set_mu_spin_degen
   public :: elsi_set_mu_mp_order
   public :: elsi_set_output_log
   public :: elsi_set_log_tag
   public :: elsi_set_uuid
   public :: elsi_get_pexsi_mu_min
   public :: elsi_get_pexsi_mu_max
   public :: elsi_get_initialized
   public :: elsi_get_version
   public :: elsi_get_datestamp
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
subroutine elsi_set_output(eh,out_level)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: out_level !< Output level

   character(len=*), parameter :: caller = "elsi_set_output"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(out_level <= 0) then
      eh%bh%print_info = 0
      eh%ph%omm_output = .false.
      eh%ph%pexsi_options%verbosity = 1
      eh%ph%elpa_output = .false.
      eh%ph%nt_output = .false.
   else if(out_level == 1) then
      eh%bh%print_info = 1
      eh%ph%omm_output = .false.
      eh%ph%pexsi_options%verbosity = 1
      eh%ph%elpa_output = .false.
      eh%ph%nt_output = .false.
   else if(out_level == 2) then
      eh%bh%print_info = 2
      eh%ph%omm_output = .true.
      eh%ph%pexsi_options%verbosity = 2
      eh%ph%elpa_output = .true.
      eh%ph%nt_output = .true.
   else
      eh%bh%print_info = 3
      eh%ph%omm_output = .true.
      eh%ph%pexsi_options%verbosity = 2
      eh%ph%elpa_output = .true.
      eh%ph%nt_output = .true.
   end if

end subroutine

!>
!! This routine sets the unit to be used by ELSI output.
!!
subroutine elsi_set_write_unit(eh,write_unit)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: write_unit !< Unit

   character(len=*), parameter :: caller = "elsi_set_write_unit"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%bh%print_unit = write_unit

end subroutine

!>
!! This routine sets the overlap matrix to be identity.
!!
subroutine elsi_set_unit_ovlp(eh,unit_ovlp)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: unit_ovlp !< Overlap is identity?

   character(len=*), parameter :: caller = "elsi_set_unit_ovlp"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(unit_ovlp == 0) then
      eh%ph%ovlp_is_unit = .false.
   else
      eh%ph%ovlp_is_unit = .true.
   end if

end subroutine

!>
!! This routine sets the threshold to define "zero".
!!
subroutine elsi_set_zero_def(eh,zero_def)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: zero_def !< Zero tolerance

   character(len=*), parameter :: caller = "elsi_set_zero_def"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%bh%def0 = zero_def

end subroutine

!>
!! This routine switches on/off the singularity check of the overlap matrix.
!!
subroutine elsi_set_sing_check(eh,sing_check)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: sing_check !< Check singularity?

   character(len=*), parameter :: caller = "elsi_set_sing_check"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(sing_check == 0) then
      eh%ph%check_sing = .false.
      eh%ph%ovlp_is_sing = .false.
      eh%ph%n_good = eh%ph%n_basis
   else
      eh%ph%check_sing = .true.
   end if

end subroutine

!>
!! This routine sets the tolerance of the singularity check.
!!
subroutine elsi_set_sing_tol(eh,sing_tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: sing_tol !< Singularity tolerance

   character(len=*), parameter :: caller = "elsi_set_sing_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%sing_tol = sing_tol

end subroutine

!>
!! This routine sets whether to stop in case of singular overlap matrix.
!!
subroutine elsi_set_sing_stop(eh,sing_stop)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: sing_stop !< Stop if overlap is singular

   character(len=*), parameter :: caller = "elsi_set_sing_stop"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(sing_stop == 0) then
      eh%ph%stop_sing = .false.
   else
      eh%ph%stop_sing = .true.
   end if

end subroutine

!>
!! This routine sets the block size in 1D block-cyclic distributed CSC format.
!!
subroutine elsi_set_csc_blk(eh,blk)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: blk !< Block size

   character(len=*), parameter :: caller = "elsi_set_csc_blk"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%bh%blk_sp2 = blk

end subroutine

!>
!! This routine sets the ELPA solver.
!!
subroutine elsi_set_elpa_solver(eh,elpa_solver)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: elpa_solver !< ELPA solver

   character(len=*), parameter :: caller = "elsi_set_elpa_solver"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(elpa_solver < 1 .or. elpa_solver > 2) then
      call elsi_stop(eh%bh,"Invalid choice.",caller)
   end if

   eh%ph%elpa_solver = elpa_solver

end subroutine

!>
!! This routine sets the number of single precision steps with ELPA.
!!
subroutine elsi_set_elpa_n_single(eh,n_single)

   implicit none

   type(elsi_handle), intent(inout) :: eh
   integer(kind=i4), intent(in) :: n_single

   character(len=*), parameter :: caller = "elsi_set_elpa_n_single"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(n_single < 0) then
      eh%ph%elpa_n_single = 0
   else
      eh%ph%elpa_n_single = n_single
   end if

end subroutine

!>
!! This routine sets whether GPU acceleration (not including GPU kernels for
!! back-transforming eigenvectors) should be enabled in ELPA. No effect if no
!! GPU acceleration available.
!!
subroutine elsi_set_elpa_gpu(eh,use_gpu)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: use_gpu !< Use GPU acceleration?

   character(len=*), parameter :: caller = "elsi_set_elpa_gpu"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(use_gpu == 0) then
      eh%ph%elpa_gpu = .false.
   else
      eh%ph%elpa_gpu = .true.
   end if

end subroutine

!>
!! This routine sets whether GPU acceleration (including GPU kernels for back-
!! transforming eigenvectors) should be enabled in ELPA. No effect if no GPU
!! acceleration available.
!!
subroutine elsi_set_elpa_gpu_kernels(eh,use_gpu_kernels)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: use_gpu_kernels !< Use GPU kernels?

   character(len=*), parameter :: caller = "elsi_set_elpa_gpu_kernels"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(use_gpu_kernels == 0) then
      eh%ph%elpa_gpu_kernels = .false.
   else
      call elsi_set_elpa_gpu(eh,1)

      eh%ph%elpa_gpu_kernels = .true.
   end if

end subroutine

!>
!! This routine sets whether auto-tuning should be enabled in ELPA. No effect if
!! no auto-tuning available.
!!
subroutine elsi_set_elpa_autotune(eh,use_autotune)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: use_autotune !< Use auto-tuning?

   character(len=*), parameter :: caller = "elsi_set_elpa_autotune"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(use_autotune == 0) then
      eh%ph%elpa_autotune = .false.
   else
      eh%ph%elpa_autotune = .true.
   end if

end subroutine

!>
!! This routine sets the flavor of libOMM.
!!
subroutine elsi_set_omm_flavor(eh,omm_flavor)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: omm_flavor !< libOMM flavor

   character(len=*), parameter :: caller = "elsi_set_omm_flavor"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(omm_flavor /= 0 .and. omm_flavor /= 2) then
      call elsi_stop(eh%bh,"Invalid choice.",caller)
   end if

   eh%ph%omm_flavor = omm_flavor

end subroutine

!>
!! This routine sets the number of ELPA steps when using libOMM.
!!
subroutine elsi_set_omm_n_elpa(eh,n_elpa)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_elpa !< ELPA steps

   character(len=*), parameter :: caller = "elsi_set_omm_n_elpa"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(n_elpa < 0) then
      eh%ph%omm_n_elpa = 0
   else
      eh%ph%omm_n_elpa = n_elpa
   end if

end subroutine

!>
!! This routine sets the tolerance of OMM minimization.
!!
subroutine elsi_set_omm_tol(eh,min_tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: min_tol !< Tolerance of OMM minimization

   character(len=*), parameter :: caller = "elsi_set_omm_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(min_tol <= 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should be positive.",caller)
   end if

   eh%ph%omm_tol = min_tol

end subroutine

!>
!! This routine sets the number of mu points when using PEXSI driver 2.
!!
subroutine elsi_set_pexsi_n_mu(eh,n_mu)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_mu !< Number of mu points

   character(len=*), parameter :: caller = "elsi_set_pexsi_n_mu"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(n_mu < 1) then
      call elsi_stop(eh%bh,"Input value should be at least 1.",caller)
   end if

   eh%ph%pexsi_options%nPoints = n_mu

end subroutine

!>
!! This routine sets the number of poles in the pole expansion.
!!
subroutine elsi_set_pexsi_n_pole(eh,n_pole)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_pole !< Number of poles

   character(len=*), parameter :: caller = "elsi_set_pexsi_n_pole"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(n_pole < 1) then
      call elsi_stop(eh%bh,"Input value should be at least 1.",caller)
   end if

   eh%ph%pexsi_options%numPole = n_pole

end subroutine

!>
!! This routine sets the number of MPI tasks assigned for one pole.
!!
subroutine elsi_set_pexsi_np_per_pole(eh,np_per_pole)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: np_per_pole !< Number of tasks per pole

   character(len=*), parameter :: caller = "elsi_set_pexsi_np_per_pole"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(np_per_pole < 1) then
      call elsi_stop(eh%bh,"Input value should be at least 1.",caller)
   end if

   eh%ph%pexsi_np_per_pole = np_per_pole

end subroutine

!>
!! This routine sets the number of MPI tasks for the symbolic factorization.
!!
subroutine elsi_set_pexsi_np_symbo(eh,np_symbo)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: np_symbo !< Number of tasks for symbolic

   character(len=*), parameter :: caller = "elsi_set_pexsi_np_symbo"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(np_symbo < 1) then
      call elsi_stop(eh%bh,"Input value should be at least 1.",caller)
   end if

   eh%ph%pexsi_options%npSymbFact = np_symbo

end subroutine

!>
!! This routine sets the matrix reordering method.
!!
subroutine elsi_set_pexsi_ordering(eh,ordering)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: ordering !< Matrix reordering method

   character(len=*), parameter :: caller = "elsi_set_pexsi_ordering"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%pexsi_options%ordering = ordering

end subroutine

!>
!! This routine sets the temperature parameter in PEXSI.
!!
subroutine elsi_set_pexsi_temp(eh,temp)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: temp !< Temperature

   character(len=*), parameter :: caller = "elsi_set_pexsi_temp"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(temp <= 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should be positive.",caller)
   end if

   eh%ph%pexsi_options%temperature = temp

end subroutine

!>
!! This routine sets the spectral gap in PEXSI.
!!
subroutine elsi_set_pexsi_gap(eh,gap)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: gap !< Gap

   character(len=*), parameter :: caller = "elsi_set_pexsi_gap"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(gap < 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should not be negative.",caller)
   end if

   eh%ph%pexsi_options%gap = gap

end subroutine

!>
!! This routine sets the spectrum width in PEXSI.
!!
subroutine elsi_set_pexsi_delta_e(eh,delta_e)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: delta_e !< Spectrum width

   character(len=*), parameter :: caller = "elsi_set_pexsi_delta_e"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(delta_e < 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should not be negative.",caller)
   end if

   eh%ph%pexsi_options%deltaE = delta_e

end subroutine

!>
!! This routine sets the lower bound of the chemical potential in PEXSI.
!!
subroutine elsi_set_pexsi_mu_min(eh,mu_min)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: mu_min !< Lower bound of mu

   character(len=*), parameter :: caller = "elsi_set_pexsi_mu_min"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%pexsi_options%muMin0 = mu_min

end subroutine

!>
!! This routine sets the upper bound of the chemical potential in PEXSI.
!!
subroutine elsi_set_pexsi_mu_max(eh,mu_max)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: mu_max !< Upper bound of mu

   character(len=*), parameter :: caller = "elsi_set_pexsi_mu_max"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%pexsi_options%muMax0 = mu_max

end subroutine

!>
!! This routine sets the tolerance of the estimation of the chemical potential
!! in the inertia counting procedure.
!!
subroutine elsi_set_pexsi_inertia_tol(eh,inertia_tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: inertia_tol !< Inertia counting tolerance

   character(len=*), parameter :: caller = "elsi_set_pexsi_inertia_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(inertia_tol <= 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should be positive.",caller)
   end if

   eh%ph%pexsi_options%muInertiaTolerance = inertia_tol

end subroutine

!>
!! This routine sets the number of ELPA steps when using SLEPc-SIPs.
!!
subroutine elsi_set_sips_n_elpa(eh,n_elpa)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_elpa !< ELPA steps

   character(len=*), parameter :: caller = "elsi_set_sips_n_elpa"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(n_elpa < 0) then
      eh%ph%sips_n_elpa = 0
   else
      eh%ph%sips_n_elpa = n_elpa
   end if

end subroutine

!>
!! This routine sets the number of slices in SLEPc-SIPs.
!!
subroutine elsi_set_sips_n_slice(eh,n_slice)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_slice !< Number of slices

   character(len=*), parameter :: caller = "elsi_set_sips_n_slice"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(mod(eh%bh%n_procs,n_slice) == 0) then
      eh%ph%sips_n_slices = n_slice
      eh%ph%sips_np_per_slice = eh%bh%n_procs/n_slice
   else
      call elsi_stop(eh%bh,"Input value should be a divisor of total number"//&
              " of MPI tasks.",caller)
   end if

end subroutine

!>
!! This routine sets the type of slices to be used in SLEPc-SIPs.
!!
subroutine elsi_set_sips_slice_type(eh,slice_type)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: slice_type !< Slice type

   character(len=*), parameter :: caller = "elsi_set_sips_slice_type"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(slice_type /= 0 .and. slice_type /= 2 .and. slice_type /= 4) then
      call elsi_stop(eh%bh,"Invalid choice.",caller)
   end if

   eh%ph%sips_slice_type = slice_type

end subroutine

!>
!! This routine sets a small buffer to expand the eigenvalue interval in
!! SLEPc-SIPs.
!!
subroutine elsi_set_sips_buffer(eh,buffer)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: buffer !< Buffer to expand interval

   character(len=*), parameter :: caller = "elsi_set_sips_buffer"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(buffer <= 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should be positive.",caller)
   end if

   eh%ph%sips_buffer = buffer

end subroutine

!>
!! This routine sets the tolerance to stop inertia counting in SLEPc-SIPs.
!!
subroutine elsi_set_sips_inertia_tol(eh,inertia_tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: inertia_tol !< Stopping criterion

   character(len=*), parameter :: caller = "elsi_set_sips_inertia_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(inertia_tol <= 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should be positive.",caller)
   end if

   eh%ph%sips_inertia_tol = inertia_tol

end subroutine

!>
!! This routine sets the global interval to be solved by SLEPc-SIPs.
!!
subroutine elsi_set_sips_interval(eh,lower,upper)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: lower !< Lower bound
   real(kind=r8), intent(in) :: upper !< Upper bound

   character(len=*), parameter :: caller = "elsi_set_sips_interval"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(lower > upper) then
      call elsi_stop(eh%bh,"Lower bound should not be larger than upper"//&
              " bound.",caller)
   end if

   eh%ph%sips_interval(1) = lower
   eh%ph%sips_interval(2) = upper

end subroutine

!>
!! This routine sets the index of the first eigensolution to be solved.
!!
subroutine elsi_set_sips_first_ev(eh,first_ev)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: first_ev !< Index of first eigensolution

   character(len=*), parameter :: caller = "elsi_set_sips_first_ev"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(first_ev < 1) then
      eh%ph%sips_first_ev = 1
   else if(first_ev > eh%ph%n_basis-eh%ph%n_states+1) then
      eh%ph%sips_first_ev = eh%ph%n_basis-eh%ph%n_states+1
   else
      eh%ph%sips_first_ev = first_ev
   end if

end subroutine

!>
!! This routine sets the density matrix purification method in NTPoly.
!!
subroutine elsi_set_ntpoly_method(eh,nt_method)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: nt_method !< Purification method

   character(len=*), parameter :: caller = "elsi_set_ntpoly_method"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(nt_method < 0 .or. nt_method > 3) then
      call elsi_stop(eh%bh,"Invalid choice.",caller)
   end if

   eh%ph%nt_method = nt_method

end subroutine

!>
!! This routine sets the inverse square root method in NTPoly.
!!
subroutine elsi_set_ntpoly_isr(eh,nt_isr)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: nt_isr !< Inverse square root method

   character(len=*), parameter :: caller = "elsi_set_ntpoly_isr"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(nt_isr /= 2 .and. nt_isr /= 3 .and. nt_isr /= 5) then
      call elsi_stop(eh%bh,"Invalid choice.",caller)
   end if

   eh%ph%nt_isr = nt_isr

end subroutine

!>
!! This routine sets the tolerance of the density matrix purification.
!!
subroutine elsi_set_ntpoly_tol(eh,nt_tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: nt_tol !< Tolerance

   character(len=*), parameter :: caller = "elsi_set_ntpoly_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(nt_tol <= 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should be positive.",caller)
   end if

   eh%ph%nt_tol = nt_tol

end subroutine

!>
!! This routine sets the threshold to filter intermediate results in NTPoly.
!!
subroutine elsi_set_ntpoly_filter(eh,filter)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: filter !< Filter

   character(len=*), parameter :: caller = "elsi_set_ntpoly_filter"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%nt_filter = filter

end subroutine

!>
!! This routine sets the maximum number of density matrix purification steps.
!!
subroutine elsi_set_ntpoly_max_iter(eh,max_iter)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: max_iter !< Maximum number of iterations

   character(len=*), parameter :: caller = "elsi_set_ntpoly_max_iter"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(max_iter < 1) then
      call elsi_stop(eh%bh,"Input value should be at least 1.",caller)
   end if

   eh%ph%nt_max_iter = max_iter

end subroutine

!>
!! This routine sets the broadening scheme to determine the chemical potential
!! and the occupation numbers.
!!
subroutine elsi_set_mu_broaden_scheme(eh,broaden_scheme)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: broaden_scheme !< Broadening method

   character(len=*), parameter :: caller = "elsi_set_mu_broaden_method"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(elpa_solver < 0 .or. elpa_solver > 4) then
      call elsi_stop(eh%bh,"Invalid choice.",caller)
   end if

   eh%ph%mu_scheme = broaden_scheme

end subroutine

!>
!! This routine sets the broadening width to determine the chemical potential
!! and the occupation numbers.
!!
subroutine elsi_set_mu_broaden_width(eh,broaden_width)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: broaden_width !< Broadening width

   character(len=*), parameter :: caller = "elsi_set_mu_broaden_width"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(broaden_width <= 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should be positive.",caller)
   end if

   eh%ph%mu_width = broaden_width

end subroutine

!>
!! This routine sets the desired accuracy of the determination of the chemical
!! potential and the occupation numbers.
!!
subroutine elsi_set_mu_tol(eh,mu_tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: mu_tol !< Accuracy of mu

   character(len=*), parameter :: caller = "elsi_set_mu_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(mu_tol < 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should not be negative.",caller)
   end if

   eh%ph%mu_tol = mu_tol

end subroutine

!>
!! This routine sets the spin degeneracy in the determination of the chemical
!! potential and the occupation numbers.
!!
subroutine elsi_set_mu_spin_degen(eh,spin_degen)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: spin_degen !< Spin degeneracy

   character(len=*), parameter :: caller = "elsi_set_mu_spin_degen"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%spin_degen = spin_degen
   eh%ph%spin_is_set = .true.

end subroutine

!>
!! This routine sets the order of the Methfessel-Paxton broadening scheme.
!!
subroutine elsi_set_mu_mp_order(eh,mp_order)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: mp_order !< Order

   character(len=*), parameter :: caller = "elsi_set_mu_mp_order"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(mp_order < 0) then
      call elsi_stop(eh%bh,"Input value should be at least 0.",caller)
   end if

   eh%ph%mu_mp_order = mp_order

end subroutine

!>
!! This routine sets whether a log file should be output.
!!
subroutine elsi_set_output_log(eh,output_log)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: output_log !< Output log

   character(len=*), parameter :: caller = "elsi_set_output_log"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%bh%print_json = output_log

end subroutine

!>
!! This routine sets the user_tag for the log.
!!
subroutine elsi_set_log_tag(eh,user_tag)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   character(len=*), intent(in) :: user_tag !< Tag

   character(len=*), parameter :: caller = "elsi_set_log_tag"

   eh%bh%user_tag = user_tag

end subroutine

!>
!! This routine sets the UUID.
!!
subroutine elsi_set_uuid(eh,uuid)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   character(len=*), intent(in) :: uuid !< UUID

   character(len=*), parameter :: caller = "elsi_set_uuid"

   eh%bh%uuid_ready = .true.
   eh%bh%uuid = uuid

end subroutine

!>
!! This routine gets the lower bound of the chemical potential returned by the
!! inertia counting in PEXSI.
!!
subroutine elsi_get_pexsi_mu_min(eh,mu_min)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: mu_min !< Lower bound of mu

   character(len=*), parameter :: caller = "elsi_get_pexsi_mu_min"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   mu_min = eh%ph%pexsi_options%muMin0

end subroutine

!>
!! This routine gets the upper bound of the chemical potential returned by the
!! inertia counting in PEXSI.
!!
subroutine elsi_get_pexsi_mu_max(eh,mu_max)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: mu_max !< Upper bound of mu

   character(len=*), parameter :: caller = "elsi_get_pexsi_mu_max"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   mu_max = eh%ph%pexsi_options%muMax0

end subroutine

!>
!! This routine returns 0 if the input handle has not been initialized; returns
!! 1 if it has been initialized.
!!
subroutine elsi_get_initialized(eh,handle_init)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(out) :: handle_init !< Handle initialized?

   character(len=*), parameter :: caller = "elsi_get_initialized"

   if(eh%handle_init) then
      handle_init = 1
   else
      handle_init = 0
   end if

end subroutine

!>
!! This routine returns version number.
!!
subroutine elsi_get_version(major,minor,patch)

   implicit none

   integer(kind=i4), intent(out) :: major !< Major version number
   integer(kind=i4), intent(out) :: minor !< Minor version number
   integer(kind=i4), intent(out) :: patch !< Patch level

   integer(kind=i4) :: i
   integer(kind=i4) :: j
   character(len=8) :: s1
   character(len=8) :: s2
   character(len=8) :: s3
   character(len=40) :: s4
   character(len=20) :: s5

   character(len=*), parameter :: caller = "elsi_get_version"

   call elsi_version_info(s1,s2,s3,s4,s5)

   i = index(s1,".",.false.)
   j = index(s1,".",.true.)

   read(s1(1:i-1),*) major
   read(s1(i+1:j-1),*) minor
   read(s1(j+1:8),*) patch

end subroutine

!>
!! This routine returns date stamp.
!!
subroutine elsi_get_datestamp(datestamp)

   implicit none

   integer(kind=i4), intent(out) :: datestamp !< Date stamp

   character(len=8) :: s1
   character(len=8) :: s2
   character(len=8) :: s3
   character(len=40) :: s4
   character(len=20) :: s5

   character(len=*), parameter :: caller = "elsi_get_datestamp"

   call elsi_version_info(s1,s2,s3,s4,s5)

   read(s2,*) datestamp

end subroutine

!>
!! This routine gets the number of basis functions that are removed due to
!! overlap singularity.
!!
subroutine elsi_get_n_sing(eh,n_sing)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(out) :: n_sing !< Number of singular basis

   character(len=*), parameter :: caller = "elsi_get_n_sing"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   n_sing = eh%ph%n_basis-eh%ph%n_good

end subroutine

!>
!! This routine gets the chemical potential.
!!
subroutine elsi_get_mu(eh,mu)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: mu !< Chemical potential

   character(len=*), parameter :: caller = "elsi_get_mu"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   mu = eh%ph%mu

end subroutine

!>
!! This routine gets the entropy.
!!
subroutine elsi_get_entropy(eh,ts)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: ts !< Entropy

   character(len=*), parameter :: caller = "elsi_get_entropy"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   ts = eh%ph%ts

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_real(eh,edm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: edm(eh%bh%n_lrow,eh%bh%n_lcol) !< Energy density matrix

   integer(kind=i4) :: solver_save

   real(kind=r8), allocatable :: tmp_real(:,:)

   character(len=*), parameter :: caller = "elsi_get_edm_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   solver_save = eh%ph%solver

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver_save = OMM_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(eh%ph%solver == SIPS_SOLVER .and. eh%ph%n_calls <= eh%ph%sips_n_elpa) then
      solver_save = SIPS_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(eh%ph%edm_ready_real) then
      select case(eh%ph%solver)
      case(ELPA_SOLVER)
         call elsi_allocate(eh%bh,tmp_real,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "tmp_real",caller)

         call elsi_compute_edm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,eh%eval,&
                 eh%evec_real,eh%occ,edm,tmp_real)

         call elsi_deallocate(eh%bh,tmp_real,"tmp_real")
      case(OMM_SOLVER)
         call elsi_compute_edm_omm(eh%ph,eh%bh,eh%omm_c_real,edm)
      case(PEXSI_SOLVER)
         call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,eh%dm_real_csc)
         call elsi_pexsi_to_blacs_dm(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%dm_real_csc,edm)
      case(SIPS_SOLVER)
         call elsi_compute_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%occ,eh%dm_real_csc)
         call elsi_sips_to_blacs_dm(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%dm_real_csc,edm)
      case(NTPOLY_SOLVER)
         call elsi_compute_edm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_dm)
         call elsi_ntpoly_to_blacs_dm(eh%bh,eh%ph%nt_dm,edm)
      case default
         call elsi_stop(eh%bh,"Unsupported density matrix solver.",caller)
      end select

      eh%ph%edm_ready_real = .false.
   else
      call elsi_stop(eh%bh,"Energy-weighted density matrix cannot be"//&
              " computed before density matrix.",caller)
   end if

   eh%ph%solver = solver_save

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_real_sparse(eh,edm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: edm(eh%bh%nnz_l_sp) !< Energy density matrix

   integer(kind=i4) :: solver_save

   real(kind=r8), allocatable :: tmp_real(:,:)

   character(len=*), parameter :: caller = "elsi_get_edm_real_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   solver_save = eh%ph%solver

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver_save = OMM_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(eh%ph%solver == SIPS_SOLVER .and. eh%ph%n_calls <= eh%ph%sips_n_elpa) then
      solver_save = SIPS_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(eh%ph%edm_ready_real) then
      select case(eh%ph%solver)
      case(ELPA_SOLVER)
         call elsi_allocate(eh%bh,tmp_real,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "tmp_real",caller)

         call elsi_compute_edm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,eh%eval,&
                 eh%evec_real,eh%occ,eh%dm_real_den,tmp_real)

         call elsi_deallocate(eh%bh,tmp_real,"tmp_real")

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_real_den,edm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                    eh%dm_real_den,edm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select
      case(OMM_SOLVER)
         call elsi_compute_edm_omm(eh%ph,eh%bh,eh%omm_c_real,eh%dm_real_den)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_real_den,edm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                    eh%dm_real_den,edm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select
      case(PEXSI_SOLVER)
         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,edm)
         case(SIESTA_CSC)
            call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,&
                    eh%dm_real_csc)
            call elsi_pexsi_to_siesta_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_real_csc,eh%row_ind_sp2,&
                    eh%col_ptr_sp2,edm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select
      case(SIPS_SOLVER)
         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_compute_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%occ,edm)
         case(SIESTA_CSC)
            call elsi_compute_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%occ,eh%dm_real_csc)
            call elsi_sips_to_siesta_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_real_csc,eh%row_ind_sp2,&
                    eh%col_ptr_sp2,edm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select
      case(NTPOLY_SOLVER)
         call elsi_compute_edm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_dm)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%ph%nt_dm,edm)
         case(SIESTA_CSC)
            call elsi_ntpoly_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                    eh%ph%nt_dm,edm)
         end select
      case default
         call elsi_stop(eh%bh,"Unsupported density matrix solver.",caller)
      end select

      eh%ph%edm_ready_real = .false.
   else
      call elsi_stop(eh%bh,"Energy-weighted density matrix cannot be"//&
              " computed before density matrix.",caller)
   end if

   eh%ph%solver = solver_save

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_complex(eh,edm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(out) :: edm(eh%bh%n_lrow,eh%bh%n_lcol) !< Energy density matrix

   integer(kind=i4) :: solver_save

   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character(len=*), parameter :: caller = "elsi_get_edm_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   solver_save = eh%ph%solver

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver_save = OMM_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(eh%ph%edm_ready_cmplx) then
      select case(eh%ph%solver)
      case(ELPA_SOLVER)
         call elsi_allocate(eh%bh,tmp_cmplx,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "tmp_cmplx",caller)

         call elsi_compute_edm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,eh%eval,&
                 eh%evec_cmplx,eh%occ,edm,tmp_cmplx)

         call elsi_deallocate(eh%bh,tmp_cmplx,"tmp_cmplx")
      case(OMM_SOLVER)
         call elsi_compute_edm_omm(eh%ph,eh%bh,eh%omm_c_cmplx,edm)
      case(PEXSI_SOLVER)
         call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,&
                 eh%dm_cmplx_csc)
         call elsi_pexsi_to_blacs_dm(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%dm_cmplx_csc,edm)
      case(NTPOLY_SOLVER)
         call elsi_compute_edm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_dm)
         call elsi_ntpoly_to_blacs_dm(eh%bh,eh%ph%nt_dm,edm)
      case default
         call elsi_stop(eh%bh,"Unsupported density matrix solver.",caller)
      end select

      eh%ph%edm_ready_cmplx = .false.
   else
      call elsi_stop(eh%bh,"Energy-weighted density matrix cannot be"//&
              " computed before density matrix.",caller)
   end if

   eh%ph%solver = solver_save

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_complex_sparse(eh,edm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(out) :: edm(eh%bh%nnz_l_sp) !< Energy density matrix

   integer(kind=i4) :: solver_save

   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character(len=*), parameter :: caller = "elsi_get_edm_complex_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   solver_save = eh%ph%solver

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver_save = OMM_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(eh%ph%edm_ready_cmplx) then
      select case(eh%ph%solver)
      case(ELPA_SOLVER)
         call elsi_allocate(eh%bh,tmp_cmplx,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "tmp_cmplx",caller)

         call elsi_compute_edm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,eh%eval,&
                 eh%evec_cmplx,eh%occ,eh%dm_cmplx_den,tmp_cmplx)

         call elsi_deallocate(eh%bh,tmp_cmplx,"tmp_cmplx")

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_cmplx_den,edm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                    eh%dm_cmplx_den,edm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select
      case(OMM_SOLVER)
         call elsi_compute_edm_omm(eh%ph,eh%bh,eh%omm_c_cmplx,eh%dm_cmplx_den)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_cmplx_den,edm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                    eh%dm_cmplx_den,edm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select
      case(PEXSI_SOLVER)
         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,edm)
         case(SIESTA_CSC)
            call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,&
                    eh%dm_cmplx_csc)
            call elsi_pexsi_to_siesta_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_cmplx_csc,eh%row_ind_sp2,&
                    eh%col_ptr_sp2,edm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select
      case(NTPOLY_SOLVER)
         call elsi_compute_edm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_dm)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%ph%nt_dm,edm)
         case(SIESTA_CSC)
            call elsi_ntpoly_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                    eh%ph%nt_dm,edm)
         end select
      case default
         call elsi_stop(eh%bh,"Unsupported density matrix solver.",caller)
      end select

      eh%ph%edm_ready_cmplx = .false.
   else
      call elsi_stop(eh%bh,"Energy-weighted density matrix cannot be"//&
              " computed before density matrix.",caller)
   end if

   eh%ph%solver = solver_save

end subroutine

end module ELSI_MUTATOR
