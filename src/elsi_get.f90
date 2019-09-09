! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide routines to retrieve information from ELSI and the solvers.
!!
module ELSI_GET

   use ELSI_CONSTANT, only: PEXSI_CSC,SIESTA_CSC,GENERIC_COO,ELPA_SOLVER,&
       OMM_SOLVER,PEXSI_SOLVER,EIGENEXA_SOLVER,SIPS_SOLVER,NTPOLY_SOLVER
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_MPI, only: elsi_stop
   use ELSI_NTPOLY, only: elsi_compute_edm_ntpoly
   use ELSI_OMM, only: elsi_compute_edm_omm
   use ELSI_PEXSI, only: elsi_compute_edm_pexsi
   use ELSI_PRECISION, only: r8,i4
   use ELSI_REDIST, only: elsi_blacs_to_generic_dm,elsi_blacs_to_siesta_dm,&
       elsi_blacs_to_sips_dm,elsi_ntpoly_to_blacs_dm,elsi_ntpoly_to_generic_dm,&
       elsi_ntpoly_to_siesta_dm,elsi_ntpoly_to_sips_dm,elsi_pexsi_to_blacs_dm,&
       elsi_pexsi_to_generic_dm,elsi_pexsi_to_siesta_dm,elsi_sips_to_blacs_dm,&
       elsi_sips_to_generic_dm,elsi_sips_to_siesta_dm
   use ELSI_SIPS, only: elsi_build_edm_sips
   use ELSI_UTIL, only: elsi_check_init,elsi_build_edm

   implicit none

   private

   public :: elsi_get_initialized
   public :: elsi_get_version
   public :: elsi_get_datestamp
   public :: elsi_get_n_illcond
   public :: elsi_get_ovlp_ev_min
   public :: elsi_get_ovlp_ev_max
   public :: elsi_get_pexsi_mu_min
   public :: elsi_get_pexsi_mu_max
   public :: elsi_get_mu
   public :: elsi_get_entropy
   public :: elsi_get_edm_real
   public :: elsi_get_edm_complex
   public :: elsi_get_edm_real_sparse
   public :: elsi_get_edm_complex_sparse
   public :: elsi_get_eval
   public :: elsi_get_evec_real
   public :: elsi_get_evec_complex

   ! Deprecated
   public :: elsi_get_n_sing

   interface elsi_get_n_sing
      module procedure elsi_get_n_illcond
   end interface

contains

!>
!! Return 0 if the input handle has not been initialized; return 1 otherwise.
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
!! Return the version number of ELSI.
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
!! Return the date stamp of ELSI source code.
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
!! Get the number of basis functions that are removed due to ill-conditioning.
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
!! Get the lowest eigenvalue of the overlap matrix.
!!
subroutine elsi_get_ovlp_ev_min(eh,ev_min)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: ev_min !< Lowest eigenvalue

   character(len=*), parameter :: caller = "elsi_get_ovlp_ev_min"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   ev_min = eh%ph%ovlp_ev_min

end subroutine

!>
!! Get the highest eigenvalue of the overlap matrix.
!!
subroutine elsi_get_ovlp_ev_max(eh,ev_max)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: ev_max !< Highest eigenvalue

   character(len=*), parameter :: caller = "elsi_get_ovlp_ev_max"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   ev_max = eh%ph%ovlp_ev_max

end subroutine

!>
!! Get the lower bound of the chemical potential returned by the inertia
!! counting in PEXSI.
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
!! Get the upper bound of the chemical potential returned by the inertia
!! counting in PEXSI.
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
!! Get the Fermi level (chemical potential).
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
!! Get the electronic entropy.
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
!! Get the energy-weighted density matrix.
!!
subroutine elsi_get_edm_real(eh,edm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: edm(eh%bh%n_lrow,eh%bh%n_lcol) !< Energy density matrix

   integer(kind=i4) :: solver_save
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_get_edm_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   solver_save = eh%ph%solver

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver_save = OMM_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(eh%ph%solver == EIGENEXA_SOLVER) then
      solver_save = EIGENEXA_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(eh%ph%solver == SIPS_SOLVER .and. eh%ph%n_calls <= eh%ph%sips_n_elpa) then
      solver_save = SIPS_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(.not. eh%ph%edm_ready) then
      write(msg,"(A)") "Energy-weighted density matrix cannot be computed"//&
         " before density matrix"
      call elsi_stop(eh%bh,msg,caller)
   end if

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_build_edm(eh%ph,eh%bh,eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),&
           eh%eval(1:eh%ph%n_states),eh%evec_real,edm)
   case(OMM_SOLVER)
      call elsi_compute_edm_omm(eh%ph,eh%bh,eh%omm_c_real,edm)
   case(PEXSI_SOLVER)
      call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,eh%dm_real_sp)
      call elsi_pexsi_to_blacs_dm(eh%ph,eh%bh,eh%dm_real_sp,eh%row_ind_sp1,&
           eh%col_ptr_sp1,edm)
   case(SIPS_SOLVER)
      call elsi_build_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
           eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%dm_real_sp)
      call elsi_sips_to_blacs_dm(eh%ph,eh%bh,eh%dm_real_sp,eh%row_ind_sp1,&
           eh%col_ptr_sp1,edm)
   case(NTPOLY_SOLVER)
      call elsi_compute_edm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_dm)
      call elsi_ntpoly_to_blacs_dm(eh%ph,eh%bh,eh%ph%nt_dm,edm)
   case default
      write(msg,"(A)") "Unsupported density matrix solver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%edm_ready = .false.
   eh%ph%solver = solver_save

end subroutine

!>
!! Get the energy-weighted density matrix.
!!
subroutine elsi_get_edm_real_sparse(eh,edm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: edm(eh%bh%nnz_l_sp) !< Energy density matrix

   integer(kind=i4) :: solver_save
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_get_edm_real_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   solver_save = eh%ph%solver

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver_save = OMM_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(eh%ph%solver == EIGENEXA_SOLVER) then
      solver_save = EIGENEXA_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(eh%ph%solver == SIPS_SOLVER .and. eh%ph%n_calls <= eh%ph%sips_n_elpa) then
      solver_save = SIPS_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(.not. eh%ph%edm_ready) then
      write(msg,"(A)") "Energy-weighted density matrix cannot be computed"//&
         " before density matrix"
      call elsi_stop(eh%bh,msg,caller)
   end if

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_build_edm(eh%ph,eh%bh,eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),&
           eh%eval(1:eh%ph%n_states),eh%evec_real,eh%dm_real_den)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%dm_real_den,edm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%dm_real_den,edm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_blacs_to_generic_dm(eh%ph,eh%bh,eh%dm_real_den,eh%map_den,&
              edm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select
   case(OMM_SOLVER)
      call elsi_compute_edm_omm(eh%ph,eh%bh,eh%omm_c_real,eh%dm_real_den)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%dm_real_den,edm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%dm_real_den,edm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_blacs_to_generic_dm(eh%ph,eh%bh,eh%dm_real_den,eh%map_den,&
              edm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select
   case(PEXSI_SOLVER)
      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,edm)
      case(SIESTA_CSC)
         call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,eh%dm_real_sp)
         call elsi_pexsi_to_siesta_dm(eh%ph,eh%bh,eh%dm_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,edm,eh%row_ind_sp2,eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,eh%dm_real_sp)
         call elsi_pexsi_to_generic_dm(eh%ph,eh%bh,eh%dm_real_sp,&
              eh%row_ind_sp1,eh%col_ptr_sp1,eh%map_sp1,edm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select
   case(SIPS_SOLVER)
      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_build_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),edm)
      case(SIESTA_CSC)
         call elsi_build_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%dm_real_sp)
         call elsi_sips_to_siesta_dm(eh%ph,eh%bh,eh%dm_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,edm,eh%row_ind_sp2,eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_build_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%dm_real_sp)
         call elsi_sips_to_generic_dm(eh%ph,eh%bh,eh%dm_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%map_sp1,edm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select
   case(NTPOLY_SOLVER)
      call elsi_compute_edm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_dm)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%ph%nt_dm,edm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_ntpoly_to_siesta_dm(eh%ph,eh%bh,eh%ph%nt_dm,edm,&
              eh%row_ind_sp2,eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_ntpoly_to_generic_dm(eh%ph,eh%bh,eh%ph%nt_dm,eh%ph%nt_map,&
              edm,eh%perm_sp3)
      end select
   case default
      write(msg,"(A)") "Unsupported density matrix solver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%edm_ready = .false.
   eh%ph%solver = solver_save

end subroutine

!>
!! Get the energy-weighted density matrix.
!!
subroutine elsi_get_edm_complex(eh,edm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(out) :: edm(eh%bh%n_lrow,eh%bh%n_lcol) !< Energy density matrix

   integer(kind=i4) :: solver_save
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_get_edm_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   solver_save = eh%ph%solver

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver_save = OMM_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(.not. eh%ph%edm_ready) then
      write(msg,"(A)") "Energy-weighted density matrix cannot be computed"//&
         " before density matrix"
      call elsi_stop(eh%bh,msg,caller)
   end if

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_build_edm(eh%ph,eh%bh,eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),&
           eh%eval(1:eh%ph%n_states),eh%evec_cmplx,edm)
   case(OMM_SOLVER)
      call elsi_compute_edm_omm(eh%ph,eh%bh,eh%omm_c_cmplx,edm)
   case(PEXSI_SOLVER)
      call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,eh%dm_cmplx_sp)
      call elsi_pexsi_to_blacs_dm(eh%ph,eh%bh,eh%dm_cmplx_sp,eh%row_ind_sp1,&
           eh%col_ptr_sp1,edm)
   case(NTPOLY_SOLVER)
      call elsi_compute_edm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_dm)
      call elsi_ntpoly_to_blacs_dm(eh%ph,eh%bh,eh%ph%nt_dm,edm)
   case default
      write(msg,"(A)") "Unsupported density matrix solver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%edm_ready = .false.
   eh%ph%solver = solver_save

end subroutine

!>
!! Get the energy-weighted density matrix.
!!
subroutine elsi_get_edm_complex_sparse(eh,edm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(out) :: edm(eh%bh%nnz_l_sp) !< Energy density matrix

   integer(kind=i4) :: solver_save
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_get_edm_complex_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   solver_save = eh%ph%solver

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver_save = OMM_SOLVER
      eh%ph%solver = ELPA_SOLVER
   end if

   if(.not. eh%ph%edm_ready) then
      write(msg,"(A)") "Energy-weighted density matrix cannot be computed"//&
         " before density matrix"
      call elsi_stop(eh%bh,msg,caller)
   end if

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_build_edm(eh%ph,eh%bh,eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),&
           eh%eval(1:eh%ph%n_states),eh%evec_cmplx,eh%dm_cmplx_den)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%dm_cmplx_den,edm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%dm_cmplx_den,edm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_blacs_to_generic_dm(eh%ph,eh%bh,eh%dm_cmplx_den,eh%map_den,&
              edm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select
   case(OMM_SOLVER)
      call elsi_compute_edm_omm(eh%ph,eh%bh,eh%omm_c_cmplx,eh%dm_cmplx_den)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%dm_cmplx_den,edm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%dm_cmplx_den,edm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_blacs_to_generic_dm(eh%ph,eh%bh,eh%dm_cmplx_den,eh%map_den,&
              edm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select
   case(PEXSI_SOLVER)
      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,edm)
      case(SIESTA_CSC)
         call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,eh%dm_cmplx_sp)
         call elsi_pexsi_to_siesta_dm(eh%ph,eh%bh,eh%dm_cmplx_sp,&
              eh%row_ind_sp1,eh%col_ptr_sp1,edm,eh%row_ind_sp2,eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_compute_edm_pexsi(eh%ph,eh%bh,eh%pexsi_ne_vec,eh%dm_cmplx_sp)
         call elsi_pexsi_to_generic_dm(eh%ph,eh%bh,eh%dm_cmplx_sp,&
              eh%row_ind_sp1,eh%col_ptr_sp1,eh%map_sp1,edm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select
   case(NTPOLY_SOLVER)
      call elsi_compute_edm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_dm)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%ph%nt_dm,edm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_ntpoly_to_siesta_dm(eh%ph,eh%bh,eh%ph%nt_dm,edm,&
              eh%row_ind_sp2,eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_ntpoly_to_generic_dm(eh%ph,eh%bh,eh%ph%nt_dm,eh%ph%nt_map,&
              edm,eh%perm_sp3)
      end select
   case default
      write(msg,"(A)") "Unsupported density matrix solver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%edm_ready = .false.
   eh%ph%solver = solver_save

end subroutine

!>
!! Get eigenvalues when elsi_dm_{real|complex}_{sparse} has been called with
!! either ELPA or SLEPc-SIPs.
!!
subroutine elsi_get_eval(eh,eval)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: eval(eh%ph%n_basis) !< Eigenvalues

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_get_eval"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(.not. eh%ph%eval_ready) then
      write(msg,"(A)") "Eigenvalues not available"
      call elsi_stop(eh%bh,msg,caller)
   end if

   eval = eh%eval
   eh%ph%eval_ready = .false.

end subroutine

!>
!! Get eigenvectors when elsi_dm has been called with ELPA or SLEPc-SIPs.
!!
subroutine elsi_get_evec_real(eh,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(out) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_get_evec_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(.not. eh%ph%evec_ready) then
      write(msg,"(A)") "Eigenvectors not available"
      call elsi_stop(eh%bh,msg,caller)
   end if

   evec = eh%evec_real
   eh%ph%evec_ready = .false.

end subroutine

!>
!! Get eigenvectors when elsi_dm has been called with ELPA or SLEPc-SIPs.
!!
subroutine elsi_get_evec_complex(eh,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(out) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_get_evec_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(.not. eh%ph%evec_ready) then
      write(msg,"(A)") "Eigenvectors not available"
      call elsi_stop(eh%bh,msg,caller)
   end if

   evec = eh%evec_cmplx
   eh%ph%evec_ready = .false.

end subroutine

end module ELSI_GET
