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
                     " Exiting...", caller)
   endif

   call elsi_start_solve_evp_time()

   if(n_elsi_calls == 1) then
      coeff_initialized = .false.
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
   else
      coeff_initialized = .true.
   endif

   if(n_elsi_calls == n_elpa_steps+1) then
      new_overlap = .true.
   else
      new_overlap = .false.
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

   ! Shift eigenvalue spectrum
   if(eta .ne. 0.0d0) then
      call m_add(ovlp_omm,'N',ham_omm,-eta,1.0d0,"lap")
   endif

   ! Solve the eigenvalue problem
   call elsi_statement_print("  Starting OMM density matrix solver")

   select case (mode)
      case (COMPLEX_VALUES)
         call omm(n_g_size,n_states,ham_omm,ovlp_omm,new_overlap,total_energy,&
                  den_mat_omm,calc_ED,eta,coeff_omm,coeff_initialized,t_den_mat_omm,&
                  scale_kinetic,omm_flavour,nk_times_nspin,i_k_spin,min_tol,omm_verbose,&
                  do_dealloc,"pzdbc","lap")

      case (REAL_VALUES)
         if(use_psp) then
            call omm(n_g_size,n_states,ham_omm,ovlp_omm,new_overlap,total_energy,&
                     den_mat_omm,calc_ED,eta,coeff_omm,coeff_initialized,t_den_mat_omm,&
                     scale_kinetic,omm_flavour,nk_times_nspin,i_k_spin,min_tol,&
                     omm_verbose,do_dealloc,"pddbc","psp")

         else
            call omm(n_g_size,n_states,ham_omm,ovlp_omm,new_overlap,total_energy,&
                     den_mat_omm,calc_ED,eta,coeff_omm,coeff_initialized,t_den_mat_omm,&
                     scale_kinetic,omm_flavour,nk_times_nspin,i_k_spin,min_tol,&
                     omm_verbose,do_dealloc,"pddbc","lap")
         endif
   end select

   call MPI_Barrier(mpi_comm_global,mpierr)
   call elsi_stop_solve_evp_time()

end subroutine ! elsi_solve_evp_omm

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
   omm_flavour = 2

   !< Scaling of the kinetic energy matrix
   scale_kinetic = 5.0d0

   !< Calculate the energy weighted density matrix
   calc_ed = .false.

   !< Eigenspectrum shift parameter
   eta = 0.0d0

   !< Tolerance for minimization
   min_tol = 1.0d-9

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

   write(info_str,"(1X,' | Eigenspectrum shift parameter ',F10.4)") eta
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Tolerance of OMM minimization ',E10.1)") min_tol
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Use pspBLAS for sparse linear algebra? ',L1)") use_psp
   call elsi_statement_print(info_str)

end subroutine

end module ELSI_OMM
