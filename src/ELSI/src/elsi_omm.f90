! Copyright (c) 2015-2017, the ELSI team. All rights reserved.
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
!! This module provides interfaces to libOMM.
!!
module ELSI_OMM

   use ISO_C_BINDING
   use ELSI_CONSTANTS, only: REAL_VALUES,COMPLEX_VALUES
   use ELSI_DIMENSIONS, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4
   use ELSI_TIMERS
   use ELSI_UTILS
   use ELPA1
   use MATRIXSWITCH, only: m_add,m_deallocate,ms_scalapack_setup,matrix

   implicit none
   private

   public :: elsi_solve_evp_omm
   public :: elsi_compute_edm_omm
   public :: elsi_set_omm_default
   public :: elsi_print_omm_options
   public :: ms_scalapack_setup_no_opt
   public :: tomato_TB_get_dims
   public :: tomato_TB_real

contains

!>
!! This routine interfaces to libOMM.
!!
subroutine elsi_solve_evp_omm(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   logical          :: success
   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_solve_evp_omm"

   if(elsi_h%overlap_is_singular) then
      call elsi_stop(" libOMM cannot treat singular overlap matrix yet."//&
                     " Exiting...",elsi_h,caller)
   endif

   ! Move to elsi_check?
   if((elsi_h%omm_flavor /= 0) .and. (elsi_h%omm_flavor /= 2)) then
      call elsi_stop(" libOMM supports flavor = 0 or 2. Exiting...",elsi_h,caller)
   endif

   call elsi_start_density_matrix_time(elsi_h)

   if(.not. elsi_h%overlap_is_unit) then
      if(elsi_h%omm_flavor == 2) then
         if(elsi_h%n_elsi_calls == 1) then
            call elsi_start_cholesky_time(elsi_h)

            ! Cholesky factorization
            select case(elsi_h%matrix_data_type)
            case(COMPLEX_VALUES)
               ! Compute S = (U^T)U, U -> S
               success = elpa_cholesky_complex_double(elsi_h%n_g_size,elsi_h%ovlp_omm%zval,&
                            elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                            elsi_h%mpi_comm_col,.false.)

               success = elpa_invert_trm_complex_double(elsi_h%n_g_size,elsi_h%ovlp_omm%zval,&
                            elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                            elsi_h%mpi_comm_col,.false.)

            case(REAL_VALUES)
               ! Compute S = (U^T)U, U -> S
               success = elpa_cholesky_real_double(elsi_h%n_g_size,elsi_h%ovlp_omm%dval,&
                            elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                            elsi_h%mpi_comm_col,.false.)

               success = elpa_invert_trm_real_double(elsi_h%n_g_size,elsi_h%ovlp_omm%dval,&
                            elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                            elsi_h%mpi_comm_col,.false.)
            end select

            call elsi_stop_cholesky_time(elsi_h)
         endif

         if(elsi_h%n_elsi_calls > elsi_h%n_elpa_steps+1) then
            ! Invert one more time
            select case(elsi_h%matrix_data_type)
            case(COMPLEX_VALUES)
               success = elpa_invert_trm_complex_double(elsi_h%n_g_size,elsi_h%ovlp_omm%zval,&
                            elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                            elsi_h%mpi_comm_col,.false.)
            case(REAL_VALUES)
               success = elpa_invert_trm_real_double(elsi_h%n_g_size,elsi_h%ovlp_omm%dval,&
                            elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                            elsi_h%mpi_comm_col,.false.)
            end select
         endif
      endif ! omm_flavor == 2
   endif ! overlap_is_unit

   if(elsi_h%n_elsi_calls == 1) then
      elsi_h%coeff_initialized = .false.
   else
      elsi_h%coeff_initialized = .true.
   endif

   if(elsi_h%n_elsi_calls == elsi_h%n_elpa_steps+1) then
      elsi_h%new_overlap = .true.
   else
      elsi_h%new_overlap = .false.
   endif

   ! Shift eigenvalue spectrum
   if(elsi_h%eta .ne. 0.0_r8) then
      call m_add(elsi_h%ovlp_omm,'N',elsi_h%ham_omm,-elsi_h%eta,1.0_r8,"lap")
   endif

   ! Solve the eigenvalue problem
   call elsi_statement_print("  Starting OMM density matrix solver",elsi_h)

   select case(elsi_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call omm(elsi_h%n_g_size,elsi_h%n_states,elsi_h%ham_omm,elsi_h%ovlp_omm,&
               elsi_h%new_overlap,elsi_h%total_energy,elsi_h%den_mat_omm,elsi_h%calc_ed,&
               elsi_h%eta,elsi_h%coeff_omm,elsi_h%coeff_initialized,elsi_h%t_den_mat_omm,&
               elsi_h%scale_kinetic,elsi_h%omm_flavor,elsi_h%nk_times_nspin,elsi_h%i_k_spin,&
               elsi_h%min_tol,elsi_h%omm_output,elsi_h%do_dealloc,"pzdbc","lap")

   case(REAL_VALUES)
      if(elsi_h%use_psp) then
         call omm(elsi_h%n_g_size,elsi_h%n_states,elsi_h%ham_omm,elsi_h%ovlp_omm,&
                  elsi_h%new_overlap,elsi_h%total_energy,elsi_h%den_mat_omm,elsi_h%calc_ed,&
                  elsi_h%eta,elsi_h%coeff_omm,elsi_h%coeff_initialized,elsi_h%t_den_mat_omm,&
                  elsi_h%scale_kinetic,elsi_h%omm_flavor,elsi_h%nk_times_nspin,elsi_h%i_k_spin,&
                  elsi_h%min_tol,elsi_h%omm_output,elsi_h%do_dealloc,"pddbc","psp")
      else
         call omm(elsi_h%n_g_size,elsi_h%n_states,elsi_h%ham_omm,elsi_h%ovlp_omm,&
                  elsi_h%new_overlap,elsi_h%total_energy,elsi_h%den_mat_omm,elsi_h%calc_ed,&
                  elsi_h%eta,elsi_h%coeff_omm,elsi_h%coeff_initialized,elsi_h%t_den_mat_omm,&
                  elsi_h%scale_kinetic,elsi_h%omm_flavor,elsi_h%nk_times_nspin,elsi_h%i_k_spin,&
                  elsi_h%min_tol,elsi_h%omm_output,elsi_h%do_dealloc,"pddbc","lap")
      endif
   end select

   call MPI_Barrier(elsi_h%mpi_comm,mpierr)
   call elsi_stop_density_matrix_time(elsi_h)

end subroutine

!> 
!! This routine computes the energy-weighted density matrix.
!! 
subroutine elsi_compute_edm_omm(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_compute_edm_omm"

   elsi_h%calc_ed = .true.

   select case(elsi_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call omm(elsi_h%n_g_size,elsi_h%n_states,elsi_h%ham_omm,elsi_h%ovlp_omm,&
               elsi_h%new_overlap,elsi_h%total_energy,elsi_h%den_mat_omm,elsi_h%calc_ed,&
               elsi_h%eta,elsi_h%coeff_omm,elsi_h%coeff_initialized,elsi_h%t_den_mat_omm,&
               elsi_h%scale_kinetic,elsi_h%omm_flavor,elsi_h%nk_times_nspin,elsi_h%i_k_spin,&
               elsi_h%min_tol,elsi_h%omm_output,elsi_h%do_dealloc,"pzdbc","lap")

   case(REAL_VALUES)
      call omm(elsi_h%n_g_size,elsi_h%n_states,elsi_h%ham_omm,elsi_h%ovlp_omm,&
               elsi_h%new_overlap,elsi_h%total_energy,elsi_h%den_mat_omm,elsi_h%calc_ed,&
               elsi_h%eta,elsi_h%coeff_omm,elsi_h%coeff_initialized,elsi_h%t_den_mat_omm,&
               elsi_h%scale_kinetic,elsi_h%omm_flavor,elsi_h%nk_times_nspin,elsi_h%i_k_spin,&
               elsi_h%min_tol,elsi_h%omm_output,elsi_h%do_dealloc,"pddbc","lap")
   end select

   elsi_h%calc_ed = .false.

end subroutine

!> 
!! This routine sets default libOMM parameters.
!! 
subroutine elsi_set_omm_default(elsi_h)
   
   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_set_omm_default"

   ! How many steps of ELPA to run before OMM
   elsi_h%n_elpa_steps = 6

   ! How do we perform the calculation
   ! 0 = Basic
   ! 1 = Cholesky factorisation of S requested
   ! 2 = Cholesky already performed, U is provided in S
   ! 3 = Use preconditioning based on the energy density
   elsi_h%omm_flavor = 0

   ! How to scale the kinetic energy matrix
   elsi_h%scale_kinetic = 5.0_r8

   ! Calculate the energy density matrix
   elsi_h%calc_ed = .false.

   ! Eigenspectrum shift parameter
   elsi_h%eta = 0.0_r8

   ! Tolerance for minimization
   elsi_h%min_tol = 1.0e-10_r8

   ! n_k_points * n_spin
   elsi_h%nk_times_nspin = 1

   ! Combined k_point spin index
   elsi_h%i_k_spin = 1

   ! Output level?
   elsi_h%omm_output = .false.

   ! Deallocate temporary arrays?
   elsi_h%do_dealloc = .false.

   ! Use pspBLAS sparse linear algebra?
   elsi_h%use_psp = .false.

end subroutine

!>
!! This routine prints libOMM settings.
!!          
subroutine elsi_print_omm_options(elsi_h)

   implicit none

   type(elsi_handle), intent(in) :: elsi_h !< Handle

   character*200 :: info_str

   character*40, parameter :: caller = "elsi_print_omm_options"

   write(info_str,"(A)") "  libOMM settings (in the same unit of Hamiltonian):"
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | ELPA steps before using libOMM ',I2)") elsi_h%n_elpa_steps
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | OMM flavor ',I2)") elsi_h%omm_flavor
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Eigenspectrum shift parameter ',F10.4)") elsi_h%eta
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Tolerance of OMM minimization ',E10.1)") elsi_h%min_tol
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Use pspBLAS for sparse linear algebra? ',L1)") elsi_h%use_psp
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! The following is a wrapper around ms_scalapack_setup without optional variables.
!! This is done by using *_present variables to mimic "present" built-in.
!! This is needed for language interoperability.
!!
subroutine ms_scalapack_setup_no_opt(mpi_comm,nprow,order,bs_def,bs_list_present,&
                                     bs_list,icontxt_present,icontxt)

   implicit none

   character(1),     intent(in) :: order
   integer(kind=i4), intent(in) :: mpi_comm
   integer(kind=i4), intent(in) :: nprow
   integer(kind=i4), intent(in) :: bs_def
   logical,          intent(in) :: bs_list_present
   integer(kind=i4), intent(in) :: bs_list(:)
   logical,          intent(in) :: icontxt_present
   integer(kind=i4), intent(in) :: icontxt

   character*40, parameter :: caller = "ms_scalapack_setup_no_opt"

   if(bs_list_present) then
      if(icontxt_present) then
         call ms_scalapack_setup(mpi_comm,nprow,order,bs_def,&
                                 bs_list=bs_list,icontxt=icontxt)
      else
         call ms_scalapack_setup(mpi_comm,nprow,order,bs_def,bs_list=bs_list)
      endif
   else
      if(icontxt_present) then
         call ms_scalapack_setup(mpi_comm,nprow,order,bs_def,icontxt=icontxt)
      else
         call ms_scalapack_setup(mpi_comm,nprow,order,bs_def)
      endif
   endif

end subroutine

!>
!! Wrapper around tomato_TB to return array dimension, but no matrix data.
!! This is needed for interoperability between languages.
!!
subroutine tomato_TB_get_dims(template_basedir,system_label,&
                              switch1,frac_occ,num_orbs_per_atom,&
                              switch2,num_orbs,num_cells_dir,&
                              switch3,sparsity,orb_r_cut,&
                              num_occ_states,&
                              gamma_point,k_point,&
                              defect,defect_perturbation,&
                              n_rows_H,n_cols_H,n_rows_S,n_cols_S,&
                              m_storage,build_matrix)

   implicit none

   character(5),     intent(in)    :: m_storage
   character(*),     intent(in)    :: template_basedir
   character(*),     intent(in)    :: system_label
   logical,          intent(in)    :: switch1
   logical,          intent(in)    :: switch2
   logical,          intent(in)    :: switch3
   logical,          intent(in)    :: gamma_point
   logical,          intent(in)    :: build_matrix
   logical,          intent(in)    :: defect
   real(kind=r8),    intent(in)    :: defect_perturbation
   integer(kind=i4), intent(out)   :: n_rows_H
   integer(kind=i4), intent(out)   :: n_cols_H
   integer(kind=i4), intent(out)   :: n_rows_S
   integer(kind=i4), intent(out)   :: n_cols_S
   integer(kind=i4), intent(out)   :: num_occ_states
   integer(kind=i4), intent(inout) :: num_orbs_per_atom
   integer(kind=i4), intent(inout) :: num_orbs
   integer(kind=i4), intent(inout) :: num_cells_dir(3)
   real(kind=r8),    intent(inout) :: frac_occ
   real(kind=r8),    intent(inout) :: sparsity
   real(kind=r8),    intent(inout) :: orb_r_cut
   real(kind=r8),    intent(inout) :: k_point(3)

   type(matrix) :: H
   type(matrix) :: S

   character*40, parameter :: caller = "tomato_TB_get_dims"

   call tomato_TB(template_basedir,system_label,&
                  switch1,frac_occ,num_orbs_per_atom,&
                  switch2,num_orbs,num_cells_dir,&
                  switch3,sparsity,orb_r_cut,&
                  num_occ_states,&
                  gamma_point,k_point,&
                  defect,defect_perturbation,&
                  H,S,m_storage,&
                  build_matrix)

   n_rows_H = H%dim1
   n_cols_H = H%dim2
   n_rows_S = S%dim1
   n_cols_S = S%dim2

   call m_deallocate(S)
   call m_deallocate(H)

end subroutine

!>
!! Wrapper around tomato_TB for the real case to eliminate the usage of the
!! type(matrix) data type. This is needed for interoperability between
!! languages. There is currently a matrix copy in this subroutine that can
!! likely be eliminated.
!!
subroutine tomato_TB_real(template_basedir,system_label,&
                          switch1,frac_occ,num_orbs_per_atom,&
                          switch2,num_orbs,num_cells_dir,&
                          switch3,sparsity,orb_r_cut,&
                          num_occ_states,&
                          gamma_point,k_point,&
                          defect,defect_perturbation,&
                          n_rows_H,n_cols_H,H_dval,&
                          n_rows_S,n_cols_S,S_dval,&
                          m_storage,build_matrix)

   implicit none

   character(5),     intent(in)    :: m_storage
   character(*),     intent(in)    :: template_basedir
   character(*),     intent(in)    :: system_label
   logical,          intent(in)    :: switch1
   logical,          intent(in)    :: switch2
   logical,          intent(in)    :: switch3
   logical,          intent(in)    :: gamma_point
   logical,          intent(in)    :: build_matrix
   logical,          intent(in)    :: defect
   real(kind=r8),    intent(in)    :: defect_perturbation
   integer(kind=i4), intent(in)    :: n_rows_H 
   integer(kind=i4), intent(in)    :: n_cols_H
   integer(kind=i4), intent(in)    :: n_rows_S
   integer(kind=i4), intent(in)    :: n_cols_S
   integer(kind=i4), intent(out)   :: num_occ_states
   real(kind=r8),    intent(out)   :: H_dval(n_rows_H,n_cols_H)
   real(kind=r8),    intent(out)   :: S_dval(n_rows_S,n_cols_S)
   integer(kind=i4), intent(inout) :: num_orbs_per_atom
   integer(kind=i4), intent(inout) :: num_orbs
   integer(kind=i4), intent(inout) :: num_cells_dir(3)
   real(kind=r8),    intent(inout) :: frac_occ
   real(kind=r8),    intent(inout) :: sparsity
   real(kind=r8),    intent(inout) :: orb_r_cut
   real(kind=r8),    intent(inout) :: k_point(3)

   type(matrix) :: H
   type(matrix) :: S

   character*40, parameter :: caller = "tomato_TB_real"

   call tomato_TB(template_basedir,system_label,&
                  switch1,frac_occ,num_orbs_per_atom,&
                  switch2,num_orbs,num_cells_dir,&
                  switch3,sparsity,orb_r_cut,&
                  num_occ_states,&
                  gamma_point,k_point,&
                  defect,defect_perturbation,&
                  H,S,m_storage,&
                  build_matrix)

   H_dval = H%dval
   S_dval = S%dval

   call m_deallocate(S)
   call m_deallocate(H)

end subroutine

end module ELSI_OMM
