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

   use iso_c_binding
   use ELSI_PRECISION, only : dp
   use ELSI_CONSTANTS, only : REAL_VALUES, COMPLEX_VALUES
   use ELSI_DIMENSIONS
   use ELSI_TIMERS
   use ELSI_UTILS
   use ELPA1
   use MatrixSwitch

   implicit none
   private

   public :: elsi_solve_evp_omm
   public :: elsi_set_omm_default_options
   public :: elsi_print_omm_options
   public :: ms_scalapack_setup_no_opt
   public :: tomato_TB_get_dims
   public :: tomato_TB_real

contains

!=======================
! ELSI routines for OMM
!=======================

!>
!! This routine interfaces to libOMM.
!!
subroutine elsi_solve_evp_omm()

   implicit none
   include "mpif.h"

   logical :: success
   character*40, parameter :: caller = "elsi_solve_evp_omm"

   if(overlap_is_singular) then
      call elsi_stop(" libOMM cannot treat singular overlap matrix yet."//&
                     " Exiting...",caller)
   endif

   if((omm_flavor /= 0) .and. (omm_flavor /= 2)) then
      call elsi_stop(" libOMM supports flavor = 0 or 2. Exiting...",caller)
   endif

   call elsi_start_density_matrix_time()

   if(.not.overlap_is_unit) then
      if(omm_flavor == 2) then
         if(n_elsi_calls == 1) then
            call elsi_start_cholesky_time()

            ! Cholesky factorization
            select case (mode)
               case (COMPLEX_VALUES)
                  ! Compute S = (U^T)U, U -> S
                  success = elpa_cholesky_complex_double(n_g_size,ovlp_omm%zval,n_l_rows,&
                               n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)

                  success = elpa_invert_trm_complex_double(n_g_size,ovlp_omm%zval,n_l_rows,&
                               n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)

               case (REAL_VALUES)
                  ! Compute S = (U^T)U, U -> S
                  success = elpa_cholesky_real_double(n_g_size,ovlp_omm%dval,n_l_rows,&
                               n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)

                  success = elpa_invert_trm_real_double(n_g_size,ovlp_omm%dval,n_l_rows,&
                               n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)

            end select

            call elsi_stop_cholesky_time()
         endif

         if(n_elsi_calls > n_elpa_steps+1) then
            ! Invert one more time
            select case (mode)
               case (COMPLEX_VALUES)
                  success = elpa_invert_trm_complex_double(n_g_size,ovlp_omm%zval,n_l_rows,&
                               n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
               case (REAL_VALUES)
                  success = elpa_invert_trm_real_double(n_g_size,ovlp_omm%dval,n_l_rows,&
                               n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
            end select
         endif
      endif ! omm_flavor == 2
   endif ! .not.overlap_is_unit

   if(n_elsi_calls == 1) then
      coeff_initialized = .false.
   else
      coeff_initialized = .true.
   endif

   if(n_elsi_calls == n_elpa_steps+1) then
      new_overlap = .true.
   else
      new_overlap = .false.
   endif

   ! Shift eigenvalue spectrum
   if(eta .ne. 0.0_dp) then
      call m_add(ovlp_omm,'N',ham_omm,-eta,1.0_dp,"lap")
   endif

   ! Solve the eigenvalue problem
   call elsi_statement_print("  Starting OMM density matrix solver")

   select case (mode)
      case (COMPLEX_VALUES)
         call omm(n_g_size,n_states,ham_omm,ovlp_omm,new_overlap,total_energy,&
                  den_mat_omm,calc_ed,eta,coeff_omm,coeff_initialized,t_den_mat_omm,&
                  scale_kinetic,omm_flavor,nk_times_nspin,i_k_spin,min_tol,omm_verbose,&
                  do_dealloc,"pzdbc","lap")

      case (REAL_VALUES)
         if(use_psp) then
            call omm(n_g_size,n_states,ham_omm,ovlp_omm,new_overlap,total_energy,&
                     den_mat_omm,calc_ed,eta,coeff_omm,coeff_initialized,t_den_mat_omm,&
                     scale_kinetic,omm_flavor,nk_times_nspin,i_k_spin,min_tol,&
                     omm_verbose,do_dealloc,"pddbc","psp")

         else
            call omm(n_g_size,n_states,ham_omm,ovlp_omm,new_overlap,total_energy,&
                     den_mat_omm,calc_ed,eta,coeff_omm,coeff_initialized,t_den_mat_omm,&
                     scale_kinetic,omm_flavor,nk_times_nspin,i_k_spin,min_tol,&
                     omm_verbose,do_dealloc,"pddbc","lap")
         endif
   end select

   call MPI_Barrier(mpi_comm_global,mpierr)
   call elsi_stop_density_matrix_time()

end subroutine

!> 
!! Set OMM variables to ELSI default.
!! 
subroutine elsi_set_omm_default_options()
   
   implicit none
   
   !< How many steps of ELPA to run before OMM
   n_elpa_steps = 3

   !< How do we perform the calculation
   !! 0 = Basic
   !! 1 = Cholesky factorisation of S requested
   !! 2 = Cholesky already performed, U is provided in S
   !! 3 = Use preconditioning based on the energy density
   omm_flavor = 2

   !< How to scale the kinetic energy matrix
   scale_kinetic = 5.0_dp

   !< Calculate the energy density matrix
   calc_ed = .false.

   !< Eigenspectrum shift parameter
   eta = 0.0_dp

   !< Tolerance for minimization
   min_tol = 1.0e-9_dp

   !< n_k_points * n_spin
   nk_times_nspin = 1

   !< Combined k_point spin index
   i_k_spin = 1

   !< Output level?
   omm_verbose = .true.

   !< Deallocate temporary arrays?
   do_dealloc = .false.

   !< Use pspBLAS sparse linear algebra?
   use_psp = .false.
      
end subroutine

!>
!! Print OMM settings.
!!          
subroutine elsi_print_omm_options()

   implicit none

   character*200 :: info_str

   write(info_str,"(A)") "  libOMM settings (in the same unit of Hamiltonian):"
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | ELPA steps before using libOMM ',I2)") n_elpa_steps
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | OMM flavor ',I2)") omm_flavor
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Eigenspectrum shift parameter ',F10.4)") eta
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Tolerance of OMM minimization ',E10.1)") min_tol
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Use pspBLAS for sparse linear algebra? ',L1)") use_psp
   call elsi_statement_print(info_str)

end subroutine

!> The following is a wrapper around ms_scalapack_setup without optional variables
!! This is done by using *_present variables to mimic "present" built-in
!! This is needed for language interoperability
subroutine ms_scalapack_setup_no_opt(mpi_comm,nprow,order,bs_def,bs_list_present,&
                                     bs_list,icontxt_present,icontxt)

  implicit none

  character(1), intent(in) :: order 
  integer, intent(in)      :: mpi_comm
  integer, intent(in)      :: nprow 
  integer, intent(in)      :: bs_def 
  logical, intent(in)      :: bs_list_present
  integer, intent(in)      :: bs_list(:)
  logical, intent(in)      :: icontxt_present
  integer, intent(in)      :: icontxt

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

!> Wrapper around tomato_TB to return array dimension, but no matrix data.
!! This is needed for interoperability between languages.
subroutine tomato_TB_get_dims(template_basedir,system_label,&
                              switch1,frac_occ,num_orbs_per_atom,&
                              switch2,num_orbs,num_cells_dir,&
                              switch3,sparsity,orb_r_cut,&
                              num_occ_states,&
                              gamma_point,k_point,&
                              defect,defect_perturbation,&
                              n_rows_H, n_cols_H, n_rows_S, n_cols_S,&
                              m_storage, build_matrix)

  implicit none

  character(5),  intent(in)    :: m_storage
  character(*),  intent(in)    :: template_basedir
  character(*),  intent(in)    :: system_label
  logical,       intent(in)    :: switch1
  logical,       intent(in)    :: switch2
  logical,       intent(in)    :: switch3
  logical,       intent(in)    :: gamma_point
  logical,       intent(in)    :: build_matrix
  logical,       intent(in)    :: defect
  real(kind=dp), intent(in)    :: defect_perturbation
  integer,       intent(out)   :: n_rows_H
  integer,       intent(out)   :: n_cols_H
  integer,       intent(out)   :: n_rows_S
  integer,       intent(out)   :: n_cols_S
  integer,       intent(out)   :: num_occ_states
  integer,       intent(inout) :: num_orbs_per_atom
  integer,       intent(inout) :: num_orbs
  integer,       intent(inout) :: num_cells_dir(3)
  real(kind=dp), intent(inout) :: frac_occ
  real(kind=dp), intent(inout) :: sparsity
  real(kind=dp), intent(inout) :: orb_r_cut
  real(kind=dp), intent(inout) :: k_point(3)

  type(matrix) :: H
  type(matrix) :: S

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

!> Wrapper around tomato_TB for the real case to eliminate the
!! usage of the type(matrix) data type.
!! This is needed for interoperability between languages.
!! There is currently a matrix copy in this subroutine that can
!! likely be eliminated.
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

  character(5),  intent(in)    :: m_storage
  character(*),  intent(in)    :: template_basedir
  character(*),  intent(in)    :: system_label
  logical,       intent(in)    :: switch1
  logical,       intent(in)    :: switch2
  logical,       intent(in)    :: switch3
  logical,       intent(in)    :: gamma_point
  logical,       intent(in)    :: build_matrix
  logical,       intent(in)    :: defect
  real(kind=dp), intent(in)    :: defect_perturbation
  integer,       intent(in)    :: n_rows_H 
  integer,       intent(in)    :: n_cols_H
  integer,       intent(in)    :: n_rows_S
  integer,       intent(in)    :: n_cols_S
  integer,       intent(out)   :: num_occ_states
  real(kind=dp), intent(out)   :: H_dval(n_rows_H,n_cols_H)
  real(kind=dp), intent(out)   :: S_dval(n_rows_S,n_cols_S)
  integer,       intent(inout) :: num_orbs_per_atom
  integer,       intent(inout) :: num_orbs
  integer,       intent(inout) :: num_cells_dir(3)
  real(kind=dp), intent(inout) :: frac_occ
  real(kind=dp), intent(inout) :: sparsity
  real(kind=dp), intent(inout) :: orb_r_cut
  real(kind=dp), intent(inout) :: k_point(3)

  type(matrix) :: H
  type(matrix) :: S

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
