! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide routines to modify parameters of ELSI and the solvers.
!!
module ELSI_SET

   use ELSI_CONSTANT, only: SINGLE_PROC
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_MPI, only: elsi_stop
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTIL, only: elsi_check_init

   implicit none

   private

   public :: elsi_set_output
   public :: elsi_set_output_unit
   public :: elsi_set_output_log
   public :: elsi_set_save_ovlp
   public :: elsi_set_unit_ovlp
   public :: elsi_set_zero_def
   public :: elsi_set_illcond_check
   public :: elsi_set_illcond_tol
   public :: elsi_set_spin_degeneracy
   public :: elsi_set_energy_gap
   public :: elsi_set_spectrum_width
   public :: elsi_set_dimensionality
   public :: elsi_set_extrapolation
   public :: elsi_set_elpa_solver
   public :: elsi_set_elpa_n_single
   public :: elsi_set_elpa_gpu
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
   public :: elsi_set_pexsi_mu_min
   public :: elsi_set_pexsi_mu_max
   public :: elsi_set_pexsi_inertia_tol
   public :: elsi_set_eigenexa_method
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
   public :: elsi_set_ntpoly_n_layer
   public :: elsi_set_magma_solver
   public :: elsi_set_mu_broaden_scheme
   public :: elsi_set_mu_broaden_width
   public :: elsi_set_mu_tol
   public :: elsi_set_mu_mp_order

   ! Deprecated
   public :: elsi_set_write_unit
   public :: elsi_set_sing_check
   public :: elsi_set_sing_tol
   public :: elsi_set_elpa_gpu_kernels
   public :: elsi_set_pexsi_gap
   public :: elsi_set_pexsi_delta_e
   public :: elsi_set_sips_interval
   public :: elsi_set_mu_spin_degen
   public :: elsi_set_uuid

   interface elsi_set_write_unit
      module procedure elsi_set_output_unit
   end interface

   interface elsi_set_sing_check
      module procedure elsi_set_illcond_check
   end interface

   interface elsi_set_sing_tol
      module procedure elsi_set_illcond_tol
   end interface

   interface elsi_set_mu_spin_degen
      module procedure elsi_set_spin_degeneracy
   end interface

   interface elsi_set_pexsi_gap
      module procedure elsi_set_energy_gap
   end interface

   interface elsi_set_pexsi_delta_e
      module procedure elsi_set_spectrum_width
   end interface

contains

!>
!! Set the output level of ELSI.
!!
subroutine elsi_set_output(eh,output)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: output !< Output level

   character(len=*), parameter :: caller = "elsi_set_output"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   select case(output)
   case(1)
      eh%bh%print_info = 1
      eh%ph%omm_output = .false.
      eh%ph%pexsi_options%verbosity = 1
      eh%ph%elpa_output = .false.
      eh%ph%nt_output = .false.
   case(2)
      eh%bh%print_info = 2
      eh%ph%omm_output = .true.
      eh%ph%pexsi_options%verbosity = 2
      eh%ph%elpa_output = .true.
      eh%ph%nt_output = .true.
   case(3)
      eh%bh%print_info = 3
      eh%ph%omm_output = .true.
      eh%ph%pexsi_options%verbosity = 2
      eh%ph%elpa_output = .true.
      eh%ph%nt_output = .true.
   case default
      eh%bh%print_info = 0
      eh%ph%omm_output = .false.
      eh%ph%pexsi_options%verbosity = 1
      eh%ph%elpa_output = .false.
      eh%ph%nt_output = .false.
   end select

end subroutine

!>
!! Set the unit to be used by ELSI output.
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
!! Set whether a log file should be output.
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
!! Set whether to save the overlap matrix or its Cholesky factor for density
!! matrix extrapolation.
!!
subroutine elsi_set_save_ovlp(eh,save_ovlp)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: save_ovlp !< Save overlap?

   character(len=*), parameter :: caller = "elsi_set_save_ovlp"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(save_ovlp == 0) then
      eh%ph%save_ovlp = .false.
   else
      eh%ph%save_ovlp = .true.
   end if

end subroutine

!>
!! Set the overlap matrix to be identity.
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
!! Set the threshold to define "zero".
!!
subroutine elsi_set_zero_def(eh,zero_def)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: zero_def !< Zero threshold

   character(len=*), parameter :: caller = "elsi_set_zero_def"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%bh%def0 = zero_def

end subroutine

!>
!! Set whether ill-conditioning should be checked.
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
!! Set the tolerance of ill-conditioning.
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
!! Set the spin degeneracy that controls the maximum number of electrons on a
!! state.
!!
subroutine elsi_set_spin_degeneracy(eh,spin_degeneracy)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: spin_degeneracy !< Spin degeneracy

   character(len=*), parameter :: caller = "elsi_set_spin_degeneracy"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%ph%spin_degen = spin_degeneracy
   eh%ph%spin_is_set = .true.

end subroutine

!>
!! Set an estimate of the band (or HOMO/LUMO) gap.
!!
subroutine elsi_set_energy_gap(eh,gap)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: gap !< Energy gap

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_energy_gap"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(gap < 0.0_r8) then
      write(msg,"(A)") "Input value cannot be negative"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%energy_gap = gap
   eh%ph%pexsi_options%gap = gap

end subroutine

!>
!! Set an estimate of the eigenspectrum width.
!!
subroutine elsi_set_spectrum_width(eh,width)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: width !< Spectrum width

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_spectrum_width"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(width < 0.0_r8) then
      write(msg,"(A)") "Input value cannot be negative"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%spectrum_width = width
   eh%ph%pexsi_options%deltaE = width

end subroutine

!>
!! Set the dimensionality of the simulating system.
!!
subroutine elsi_set_dimensionality(eh,dimensionality)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: dimensionality !< Dimensionality

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_dimensionality"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(dimensionality < 1 .or. dimensionality > 3) then
      write(msg,"(A)") "Input value should be 1, 2, or 3"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%dimensionality = dimensionality

end subroutine

!>
!! Set the method for density matrix extrapolation.
!!
subroutine elsi_set_extrapolation(eh,extrapolation)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: extrapolation !< extrapolation

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_extrapolation"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(extrapolation < 0 .or. extrapolation > 1) then
      write(msg,"(A)") "Input value should be 0 or 1"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%extrapolation = extrapolation

end subroutine

!>
!! Set the ELPA solver.
!!
subroutine elsi_set_elpa_solver(eh,solver)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: solver !< ELPA solver

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_elpa_solver"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(solver < 1 .or. solver > 2) then
      write(msg,"(A)") "Input value should be 1 or 2"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%elpa_solver = solver

end subroutine

!>
!! Set the number of single precision steps with ELPA.
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
!! Set whether GPU acceleration (not including GPU kernels for back-transforming
!! eigenvectors) should be enabled in ELPA. No effect if no GPU acceleration
!! available.
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
!! Set whether auto-tuning should be enabled in ELPA. No effect if no
!! auto-tuning available.
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
!! Set the flavor of libOMM.
!!
subroutine elsi_set_omm_flavor(eh,flavor)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: flavor !< libOMM flavor

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_omm_flavor"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(flavor /= 0 .and. flavor /= 2) then
      write(msg,"(A)") "Input value should be 0 or 2"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%omm_flavor = flavor

end subroutine

!>
!! Set the number of ELPA steps when using libOMM.
!!
subroutine elsi_set_omm_n_elpa(eh,n_elpa)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_elpa !< ELPA steps

   character(len=*), parameter :: caller = "elsi_set_omm_n_elpa"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(n_elpa < 1) then
      eh%ph%omm_n_elpa = 1
   else
      eh%ph%omm_n_elpa = n_elpa
   end if

end subroutine

!>
!! Set the tolerance of OMM minimization.
!!
subroutine elsi_set_omm_tol(eh,tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: tol !< Tolerance of OMM minimization

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_omm_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(tol <= 0.0_r8) then
      write(msg,"(A)") "Input value should be positive"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%omm_tol = tol

end subroutine

!>
!! Set the number of mu points when using PEXSI driver 2.
!!
subroutine elsi_set_pexsi_n_mu(eh,n_mu)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_mu !< Number of mu points

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_pexsi_n_mu"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(n_mu < 1) then
      write(msg,"(A)") "Input value cannot be smaller than 1"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%pexsi_options%nPoints = n_mu

end subroutine

!>
!! Set the number of poles in the pole expansion.
!!
subroutine elsi_set_pexsi_n_pole(eh,n_pole)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_pole !< Number of poles

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_pexsi_n_pole"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(n_pole < 1) then
      write(msg,"(A)") "Input value cannot be smaller than 1"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%pexsi_options%numPole = n_pole

end subroutine

!>
!! Set the number of MPI tasks assigned for one pole.
!!
subroutine elsi_set_pexsi_np_per_pole(eh,np_per_pole)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: np_per_pole !< Number of tasks per pole

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_pexsi_np_per_pole"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(np_per_pole < 1) then
      write(msg,"(A)") "Input value cannot be smaller than 1"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%pexsi_np_per_pole = np_per_pole

end subroutine

!>
!! Set the number of MPI tasks for the symbolic factorization.
!!
subroutine elsi_set_pexsi_np_symbo(eh,np_symbo)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: np_symbo !< Number of tasks for symbolic

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_pexsi_np_symbo"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(np_symbo < 1) then
      write(msg,"(A)") "Input value cannot be smaller than 1"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%pexsi_options%npSymbFact = np_symbo

end subroutine

!>
!! Set the matrix reordering method.
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
!! Set the electronic temperature in PEXSI.
!!
subroutine elsi_set_pexsi_temp(eh,temp)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: temp !< Temperature

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_pexsi_temp"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(temp <= 0.0_r8) then
      write(msg,"(A)") "Input value should be positive"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%pexsi_options%temperature = temp

end subroutine

!>
!! Set the lower bound of the chemical potential in PEXSI.
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
!! Set the upper bound of the chemical potential in PEXSI.
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
!! Set the tolerance of the estimation of the chemical potential in the inertia
!! counting procedure.
!!
subroutine elsi_set_pexsi_inertia_tol(eh,inertia_tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: inertia_tol !< Inertia counting tolerance

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_pexsi_inertia_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(inertia_tol <= 0.0_r8) then
      write(msg,"(A)") "Input value should be positive"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%pexsi_options%muInertiaTolerance = inertia_tol

end subroutine

!>
!! Set the EigenExa algorithm.
!!
subroutine elsi_set_eigenexa_method(eh,method)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: method !< Method

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_eigenexa_method"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(method < 1 .or. method > 2) then
      write(msg,"(A)") "Input value should be 1 or 2"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%exa_method = method

end subroutine

!>
!! Set the number of ELPA steps when using SLEPc-SIPs.
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
!! Set the number of slices in SLEPc-SIPs.
!!
subroutine elsi_set_sips_n_slice(eh,n_slice)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_slice !< Number of slices

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_sips_n_slice"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(mod(eh%bh%n_procs,n_slice) == 0) then
      eh%ph%sips_n_slices = n_slice
   else
      write(msg,"(A)") "Input value should be a divisor of total number of"//&
         " MPI tasks"
      call elsi_stop(eh%bh,msg,caller)
   end if

end subroutine

!>
!! Set the type of slices to be used in SLEPc-SIPs.
!!
subroutine elsi_set_sips_slice_type(eh,slice_type)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: slice_type !< Slice type

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_sips_slice_type"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(slice_type /= 0 .and. slice_type /= 2 .and. slice_type /= 4) then
      write(msg,"(A)") "Input value should be 0, 2, or 4"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%sips_slice_type = slice_type

end subroutine

!>
!! Set a small buffer to expand the eigenvalue interval in SLEPc-SIPs.
!!
subroutine elsi_set_sips_buffer(eh,buffer)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: buffer !< Buffer to expand interval

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_sips_buffer"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(buffer <= 0.0_r8) then
      write(msg,"(A)") "Input value should be positive"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%sips_buffer = buffer

end subroutine

!>
!! Set the tolerance to stop inertia counting in SLEPc-SIPs.
!!
subroutine elsi_set_sips_inertia_tol(eh,inertia_tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: inertia_tol !< Stopping criterion

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_sips_inertia_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(inertia_tol <= 0.0_r8) then
      write(msg,"(A)") "Input value should be positive"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%sips_inertia_tol = inertia_tol

end subroutine

!>
!! Set the lower bound of the interval to be solved by SLEPc-SIPs.
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
!! Set the upper bound of the interval to be solved by SLEPc-SIPs.
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
!! Set the global interval to be solved by SLEPc-SIPs.
!!
subroutine elsi_set_sips_interval(eh,lower,upper)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: lower !< Lower bound
   real(kind=r8), intent(in) :: upper !< Upper bound

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_sips_interval"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(lower >= upper) then
      write(msg,"(A)") "Lower bound must be smaller than upper bound"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%sips_interval(1) = lower
   eh%ph%sips_interval(2) = upper

end subroutine

!>
!! Set the density matrix purification method in NTPoly.
!!
subroutine elsi_set_ntpoly_method(eh,method)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: method !< Purification method

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_ntpoly_method"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(method < 0 .or. method > 3) then
      write(msg,"(A)") "Input value should be 0, 1, 2, or 3"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%nt_method = method

end subroutine

!>
!! Set the inverse square root method in NTPoly.
!!
subroutine elsi_set_ntpoly_isr(eh,isr)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: isr !< Inverse square root method

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_ntpoly_isr"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(isr /= 2 .and. isr /= 3 .and. isr /= 5) then
      write(msg,"(A)") "Input value should be 2, 3, or 5"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%nt_isr = isr

end subroutine

!>
!! Set the tolerance of the density matrix purification.
!!
subroutine elsi_set_ntpoly_tol(eh,tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: tol !< Tolerance

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_ntpoly_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(tol <= 0.0_r8) then
      write(msg,"(A)") "Input value should be positive"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%nt_tol = tol

end subroutine

!>
!! Set the threshold to filter intermediate results in NTPoly.
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
!! Set the maximum number of density matrix purification steps.
!!
subroutine elsi_set_ntpoly_max_iter(eh,max_iter)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: max_iter !< Number of iterations

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_ntpoly_max_iter"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(max_iter < 1) then
      write(msg,"(A)") "Input value cannot be smaller than 1"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%nt_max_iter = max_iter

end subroutine

!>
!! Set the number of layers in the 3D process grid of NTPoly.
!!
subroutine elsi_set_ntpoly_n_layer(eh,n_layer)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_layer !< Number of process layers

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_ntpoly_n_layer"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(mod(eh%bh%n_procs,n_layer) == 0) then
      eh%ph%nt_n_layers = n_layer
   else
      write(msg,"(A)") "Input value should be a divisor of total number of"//&
         " MPI tasks"
      call elsi_stop(eh%bh,msg,caller)
   end if

end subroutine

!>
!! Set the solver of MAGMA.
!!
subroutine elsi_set_magma_solver(eh,solver)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: solver !< MAGMA solver

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_magma_solver"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(solver < 1 .or. solver > 2) then
      write(msg,"(A)") "Input value should be 1 or 2"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%magma_solver = solver

end subroutine

!>
!! Set the broadening scheme to determine the chemical potential and the
!! occupation numbers.
!!
subroutine elsi_set_mu_broaden_scheme(eh,broaden_scheme)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: broaden_scheme !< Broadening method

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_mu_broaden_method"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(broaden_scheme < 0 .or. broaden_scheme > 4) then
      write(msg,"(A)") "Input value should be 0, 1, 2, 3, or 4"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%mu_scheme = broaden_scheme

end subroutine

!>
!! Set the broadening width to determine the chemical potential and the
!! occupation numbers.
!!
subroutine elsi_set_mu_broaden_width(eh,broaden_width)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: broaden_width !< Broadening width

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_mu_broaden_width"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(broaden_width <= 0.0_r8) then
      write(msg,"(A)") "Input value should be positive"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%mu_width = broaden_width

end subroutine

!>
!! Set the desired accuracy of the determination of the chemical potential and
!! the occupation numbers.
!!
subroutine elsi_set_mu_tol(eh,tol)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: tol !< Accuracy of mu

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_mu_tol"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(tol < 0.0_r8) then
      write(msg,"(A)") "Input value cannot be negative"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%mu_tol = tol

end subroutine

!>
!! Set the order of the Methfessel-Paxton broadening scheme.
!!
subroutine elsi_set_mu_mp_order(eh,mp_order)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: mp_order !< Order

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_set_mu_mp_order"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(mp_order < 0) then
      write(msg,"(A)") "Input value cannot be negative"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eh%ph%mu_mp_order = mp_order

end subroutine

!>
!! Set a UUID.
!! (Deprecated)
!!
subroutine elsi_set_uuid(eh,uuid)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   character(len=*), intent(in) :: uuid !< UUID

   character(len=36) :: tmp

   character(len=*), parameter :: caller = "elsi_set_uuid"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   tmp = uuid

end subroutine

!>
!! Set whether GPU kernels for back-transforming eigenvectors should be enabled
!! in ELPA.
!! (Deprecated)
!!
subroutine elsi_set_elpa_gpu_kernels(eh,gpu_kernels)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: gpu_kernels !< Use GPU kernels?

   integer(kind=i4) :: tmp

   character(len=*), parameter :: caller = "elsi_set_elpa_gpu_kernels"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   tmp = gpu_kernels

end subroutine

end module ELSI_SET
