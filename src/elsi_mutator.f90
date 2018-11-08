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
   use ELSI_MPI, only: elsi_stop
   use ELSI_NTPOLY, only: elsi_compute_edm_ntpoly
   use ELSI_OMM, only: elsi_compute_edm_omm
   use ELSI_PEXSI, only: elsi_compute_edm_pexsi
   use ELSI_PRECISION, only: r8,i4
   use ELSI_REDIST, only: elsi_blacs_to_siesta_dm,elsi_blacs_to_sips_dm,&
       elsi_ntpoly_to_blacs_dm,elsi_ntpoly_to_siesta_dm,elsi_ntpoly_to_sips_dm,&
       elsi_pexsi_to_blacs_dm,elsi_pexsi_to_siesta_dm,elsi_sips_to_blacs_dm,&
       elsi_sips_to_siesta_dm
   use ELSI_SIPS, only: elsi_build_edm_sips
   use ELSI_UTILS, only: elsi_check_init,elsi_build_edm

   implicit none

   private

   public :: elsi_set_output
   public :: elsi_set_output_unit
   public :: elsi_set_output_log
   public :: elsi_set_output_tag
   public :: elsi_set_uuid
   public :: elsi_set_unit_ovlp
   public :: elsi_set_zero_def
   public :: elsi_set_illcond_check
   public :: elsi_set_illcond_tol
   public :: elsi_set_illcond_abort
   public :: elsi_set_csc_blk
   public :: elsi_set_elpa_solver
   public :: elsi_set_elpa_cholesky
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
   public :: elsi_set_sips_ev_min
   public :: elsi_set_sips_ev_max
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
   public :: elsi_get_initialized
   public :: elsi_get_version
   public :: elsi_get_datestamp
   public :: elsi_get_n_illcond
   public :: elsi_get_pexsi_mu_min
   public :: elsi_get_pexsi_mu_max
   public :: elsi_get_mu
   public :: elsi_get_entropy
   public :: elsi_get_edm_real
   public :: elsi_get_edm_complex
   public :: elsi_get_edm_real_sparse
   public :: elsi_get_edm_complex_sparse

   ! Deprecated
   public :: elsi_set_write_unit
   public :: elsi_set_sing_check
   public :: elsi_set_sing_tol
   public :: elsi_set_sing_stop
   public :: elsi_set_sips_interval
   public :: elsi_get_n_sing

   interface elsi_set_write_unit
      module procedure elsi_set_output_unit
   end interface

   interface elsi_set_sing_check
      module procedure elsi_set_illcond_check
   end interface

   interface elsi_set_sing_tol
      module procedure elsi_set_illcond_tol
   end interface

   interface elsi_set_sing_stop
      module procedure elsi_set_illcond_abort
   end interface

   interface elsi_get_n_sing
      module procedure elsi_get_n_illcond
   end interface

contains

!>
!! This routine sets ELSI output level.
!!
subroutine elsi_set_output(eh,output)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: output !< Output level

   character(len=*), parameter :: caller = "elsi_set_output"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(output <= 0) then
      eh%bh%print_info = 0
      eh%ph%omm_output = .false.
      eh%ph%pexsi_options%verbosity = 1
      eh%ph%elpa_output = .false.
      eh%ph%nt_output = .false.
   else if(output == 1) then
      eh%bh%print_info = 1
      eh%ph%omm_output = .false.
      eh%ph%pexsi_options%verbosity = 1
      eh%ph%elpa_output = .false.
      eh%ph%nt_output = .false.
   else if(output == 2) then
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
subroutine elsi_set_output_unit(eh,print_unit)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: print_unit !< Unit

   character(len=*), parameter :: caller = "elsi_set_output_unit"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%bh%print_unit = print_unit

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
subroutine elsi_set_output_tag(eh,output_tag)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   character(len=*), intent(in) :: output_tag !< Tag

   character(len=*), parameter :: caller = "elsi_set_output_tag"

   eh%bh%user_tag = output_tag

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
!! This routine sets the overlap matrix to be identity.
!!
subroutine elsi_set_unit_ovlp(eh,unit_ovlp)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: unit_ovlp !< Overlap is identity?

   character(len=*), parameter :: caller = "elsi_set_unit_ovlp"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(unit_ovlp == 0) then
      eh%ph%unit_ovlp = .false.
   else
      eh%ph%unit_ovlp = .true.
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
!! This routine sets whether ill-conditioning should be checked.
!!
subroutine elsi_set_illcond_check(eh,illcond_check)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: illcond_check !< Check singularity?

   character(len=*), parameter :: caller = "elsi_set_illcond_check"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(illcond_check == 0) then
      eh%ph%ill_check = .false.
      eh%ph%ill_ovlp = .false.
      eh%ph%n_good = eh%ph%n_basis
   else
      eh%ph%ill_check = .true.
   end if

end subroutine

!>
!! This routine sets the tolerance of ill-conditioning.
!!
subroutine elsi_set_illcond_tol(eh,illcond_tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: illcond_tol !< Singularity tolerance

   character(len=*), parameter :: caller = "elsi_set_illcond_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%ill_tol = illcond_tol

end subroutine

!>
!! This routine sets whether to abort in case of ill-conditioning.
!!
subroutine elsi_set_illcond_abort(eh,illcond_abort)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: illcond_abort !< Abort if ill-conditioned

   character(len=*), parameter :: caller = "elsi_set_illcond_abort"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(illcond_abort == 0) then
      eh%ph%ill_abort = .false.
   else
      eh%ph%ill_abort = .true.
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
subroutine elsi_set_elpa_solver(eh,solver)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: solver !< ELPA solver

   character(len=*), parameter :: caller = "elsi_set_elpa_solver"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(solver < 1 .or. solver > 2) then
      call elsi_stop(eh%bh,"Invalid choice.",caller)
   end if

   eh%ph%elpa_solver = solver

end subroutine

!>
!! This routine sets whether the Cholesky factorization step in ELPA has been
!! performed externally or should be performed internally.
!!
subroutine elsi_set_elpa_cholesky(eh,cholesky)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: cholesky !< Cholesky factorized?

   character(len=*), parameter :: caller = "elsi_set_elpa_cholesky"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(cholesky == 0) then
      eh%ph%elpa_first = .true.
   else
      eh%ph%elpa_first = .false.
   end if

end subroutine

!>
!! This routine sets the number of single precision steps with ELPA.
!!
subroutine elsi_set_elpa_n_single(eh,n_single)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_single !< Single-precision steps

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
subroutine elsi_set_elpa_gpu(eh,gpu)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: gpu !< Use GPU acceleration?

   character(len=*), parameter :: caller = "elsi_set_elpa_gpu"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(gpu == 0) then
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
subroutine elsi_set_elpa_gpu_kernels(eh,gpu_kernels)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: gpu_kernels !< Use GPU kernels?

   character(len=*), parameter :: caller = "elsi_set_elpa_gpu_kernels"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(gpu_kernels == 0) then
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
subroutine elsi_set_elpa_autotune(eh,autotune)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: autotune !< Use auto-tuning?

   character(len=*), parameter :: caller = "elsi_set_elpa_autotune"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(autotune == 0) then
      eh%ph%elpa_autotune = .false.
   else
      eh%ph%elpa_autotune = .true.
   end if

end subroutine

!>
!! This routine sets the flavor of libOMM.
!!
subroutine elsi_set_omm_flavor(eh,flavor)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: flavor !< libOMM flavor

   character(len=*), parameter :: caller = "elsi_set_omm_flavor"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(flavor /= 0 .and. flavor /= 2) then
      call elsi_stop(eh%bh,"Invalid choice.",caller)
   end if

   eh%ph%omm_flavor = flavor

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
subroutine elsi_set_omm_tol(eh,tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: tol !< Tolerance of OMM minimization

   character(len=*), parameter :: caller = "elsi_set_omm_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(tol <= 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should be positive.",caller)
   end if

   eh%ph%omm_tol = tol

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
!! This routine sets the electronic temperature in PEXSI.
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
!! This routine sets the lower bound of the interval to be solved by SLEPc-SIPs.
!!
subroutine elsi_set_sips_ev_min(eh,ev_min)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ev_min !< Lower bound

   character(len=*), parameter :: caller = "elsi_set_sips_ev_min"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%sips_interval(1) = ev_min

end subroutine

!>
!! This routine sets the upper bound of the interval to be solved by SLEPc-SIPs.
!!
subroutine elsi_set_sips_ev_max(eh,ev_max)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ev_max !< Upper bound

   character(len=*), parameter :: caller = "elsi_set_sips_ev_max"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%sips_interval(2) = ev_max

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
!! This routine sets the density matrix purification method in NTPoly.
!!
subroutine elsi_set_ntpoly_method(eh,method)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: method !< Purification method

   character(len=*), parameter :: caller = "elsi_set_ntpoly_method"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(method < 0 .or. method > 3) then
      call elsi_stop(eh%bh,"Invalid choice.",caller)
   end if

   eh%ph%nt_method = method

end subroutine

!>
!! This routine sets the inverse square root method in NTPoly.
!!
subroutine elsi_set_ntpoly_isr(eh,isr)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: isr !< Inverse square root method

   character(len=*), parameter :: caller = "elsi_set_ntpoly_isr"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(isr /= 2 .and. isr /= 3 .and. isr /= 5) then
      call elsi_stop(eh%bh,"Invalid choice.",caller)
   end if

   eh%ph%nt_isr = isr

end subroutine

!>
!! This routine sets the tolerance of the density matrix purification.
!!
subroutine elsi_set_ntpoly_tol(eh,tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: tol !< Tolerance

   character(len=*), parameter :: caller = "elsi_set_ntpoly_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(tol <= 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should be positive.",caller)
   end if

   eh%ph%nt_tol = tol

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
   integer(kind=i4), intent(in) :: max_iter !< Number of iterations

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

   if(broaden_scheme < 0 .or. broaden_scheme > 4) then
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
subroutine elsi_set_mu_tol(eh,tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: tol !< Accuracy of mu

   character(len=*), parameter :: caller = "elsi_set_mu_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(tol < 0.0_r8) then
      call elsi_stop(eh%bh,"Input value should not be negative.",caller)
   end if

   eh%ph%mu_tol = tol

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
subroutine elsi_get_initialized(eh,initialized)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(out) :: initialized !< Handle initialized?

   character(len=*), parameter :: caller = "elsi_get_initialized"

   if(eh%handle_init) then
      initialized = 1
   else
      initialized = 0
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
!! ill-conditioning.
!!
subroutine elsi_get_n_illcond(eh,n_illcond)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(out) :: n_illcond !< Number of removed functions

   character(len=*), parameter :: caller = "elsi_get_n_illcond"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   n_illcond = eh%ph%n_basis-eh%ph%n_good

end subroutine

!>
!! This routine gets the Fermi level (chemical potential).
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
subroutine elsi_get_entropy(eh,entropy)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: entropy !< Entropy

   character(len=*), parameter :: caller = "elsi_get_entropy"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   entropy = eh%ph%ts

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_real(eh,edm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: edm(eh%bh%n_lrow,eh%bh%n_lcol) !< Energy density matrix

   integer(kind=i4) :: solver_save

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
         call elsi_build_edm(eh%ph,eh%bh,eh%row_map,eh%col_map,&
              eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%eval(1:eh%ph%n_states),&
              eh%evec_real,edm)
      case(OMM_SOLVER)
         call elsi_compute_edm_omm(eh%ph,eh%bh,eh%omm_c_real,edm)
      case(PEXSI_SOLVER)
         call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,eh%dm_real_csc)
         call elsi_pexsi_to_blacs_dm(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%dm_real_csc,edm)
      case(SIPS_SOLVER)
         call elsi_build_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%dm_real_csc)
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
         call elsi_build_edm(eh%ph,eh%bh,eh%row_map,eh%col_map,&
              eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%eval(1:eh%ph%n_states),&
              eh%evec_real,eh%dm_real_den)

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
                 eh%col_ptr_sp1,eh%dm_real_csc,eh%row_ind_sp2,eh%col_ptr_sp2,&
                 edm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select
      case(SIPS_SOLVER)
         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_build_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),edm)
         case(SIESTA_CSC)
            call elsi_build_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%dm_real_csc)
            call elsi_sips_to_siesta_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                 eh%col_ptr_sp1,eh%dm_real_csc,eh%row_ind_sp2,eh%col_ptr_sp2,&
                 edm)
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
         call elsi_build_edm(eh%ph,eh%bh,eh%row_map,eh%col_map,&
              eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%eval(1:eh%ph%n_states),&
              eh%evec_cmplx,edm)
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
         call elsi_build_edm(eh%ph,eh%bh,eh%row_map,eh%col_map,&
              eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%eval(1:eh%ph%n_states),&
              eh%evec_cmplx,eh%dm_cmplx_den)

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
                 eh%col_ptr_sp1,eh%dm_cmplx_csc,eh%row_ind_sp2,eh%col_ptr_sp2,&
                 edm)
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
