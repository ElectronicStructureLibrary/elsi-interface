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

   use ELSI_CONSTANTS, only: REAL_VALUES,COMPLEX_VALUES
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS
   use ELPA1
   use MATRIXSWITCH, only: m_add

   implicit none

   private

   public :: elsi_solve_evp_omm
   public :: elsi_compute_edm_omm
   public :: elsi_set_omm_default
   public :: elsi_print_omm_options

contains

!>
!! This routine interfaces to libOMM.
!!
subroutine elsi_solve_evp_omm(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   logical          :: success
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   integer(kind=i4) :: mpierr
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_solve_evp_omm"

   if(elsi_h%ovlp_is_sing) then
      call elsi_stop(" libOMM cannot treat singular overlap matrix yet."//&
              " Exiting...",elsi_h,caller)
   endif

   if(.not. elsi_h%ovlp_is_unit) then
      if(elsi_h%omm_flavor == 2) then
         if(elsi_h%n_elsi_calls == 1) then
            call elsi_get_time(elsi_h,t0)

            ! Cholesky factorization
            select case(elsi_h%matrix_data_type)
            case(COMPLEX_VALUES)
               ! Compute S = (U^T)U, U -> S
               success = elpa_cholesky_complex_double(elsi_h%n_basis,&
                            elsi_h%ovlp_omm%zval,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                            elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                            .false.)

               success = elpa_invert_trm_complex_double(elsi_h%n_basis,&
                            elsi_h%ovlp_omm%zval,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                            elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                            .false.)
            case(REAL_VALUES)
               ! Compute S = (U^T)U, U -> S
               success = elpa_cholesky_real_double(elsi_h%n_basis,&
                            elsi_h%ovlp_omm%dval,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                            elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                            .false.)

               success = elpa_invert_trm_real_double(elsi_h%n_basis,&
                            elsi_h%ovlp_omm%dval,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                            elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                            .false.)
            end select

            call elsi_get_time(elsi_h,t1)

            write(info_str,"('  Finished Cholesky decomposition')")
            call elsi_statement_print(info_str,elsi_h)
            write(info_str,"('  | Time :',F10.3,' s')") t1-t0
            call elsi_statement_print(info_str,elsi_h)
         endif

         if(elsi_h%n_elsi_calls > elsi_h%n_elpa_steps+1) then
            ! Invert one more time
            select case(elsi_h%matrix_data_type)
            case(COMPLEX_VALUES)
               success = elpa_invert_trm_complex_double(elsi_h%n_basis,&
                            elsi_h%ovlp_omm%zval,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                            elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                            .false.)
            case(REAL_VALUES)
               success = elpa_invert_trm_real_double(elsi_h%n_basis,&
                            elsi_h%ovlp_omm%dval,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                            elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                            .false.)
            end select
         endif
      endif ! omm_flavor == 2
   endif ! ovlp_is_unit

   if(elsi_h%n_elsi_calls == 1) then
      elsi_h%coeff_ready = .false.
   else
      elsi_h%coeff_ready = .true.
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

   call elsi_get_time(elsi_h,t0)

   call elsi_statement_print("  Starting OMM density matrix solver",elsi_h)

   select case(elsi_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call omm(elsi_h%n_basis,elsi_h%n_states_omm,elsi_h%ham_omm,elsi_h%ovlp_omm,&
              elsi_h%new_overlap,elsi_h%energy_hdm,elsi_h%dm_omm,elsi_h%calc_ed,&
              elsi_h%eta,elsi_h%coeff_omm,elsi_h%coeff_ready,elsi_h%tdm_omm,&
              elsi_h%scale_kinetic,elsi_h%omm_flavor,1,1,elsi_h%min_tol,&
              elsi_h%omm_output,elsi_h%do_dealloc,"pzdbc","lap")

   case(REAL_VALUES)
      if(elsi_h%use_psp) then
         call omm(elsi_h%n_basis,elsi_h%n_states_omm,elsi_h%ham_omm,elsi_h%ovlp_omm,&
                 elsi_h%new_overlap,elsi_h%energy_hdm,elsi_h%dm_omm,elsi_h%calc_ed,&
                 elsi_h%eta,elsi_h%coeff_omm,elsi_h%coeff_ready,elsi_h%tdm_omm,&
                 elsi_h%scale_kinetic,elsi_h%omm_flavor,1,1,elsi_h%min_tol,&
                 elsi_h%omm_output,elsi_h%do_dealloc,"pddbc","psp")
      else
         call omm(elsi_h%n_basis,elsi_h%n_states_omm,elsi_h%ham_omm,elsi_h%ovlp_omm,&
                 elsi_h%new_overlap,elsi_h%energy_hdm,elsi_h%dm_omm,elsi_h%calc_ed,&
                 elsi_h%eta,elsi_h%coeff_omm,elsi_h%coeff_ready,elsi_h%tdm_omm,&
                 elsi_h%scale_kinetic,elsi_h%omm_flavor,1,1,elsi_h%min_tol,&
                 elsi_h%omm_output,elsi_h%do_dealloc,"pddbc","lap")
      endif
   end select

   call MPI_Barrier(elsi_h%mpi_comm,mpierr)

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished density matrix calculation')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine computes the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_omm(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character*200 :: info_str

   character*40, parameter :: caller = "elsi_compute_edm_omm"

   call elsi_get_time(elsi_h,t0)

   elsi_h%calc_ed = .true.

   select case(elsi_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call omm(elsi_h%n_basis,elsi_h%n_states_omm,elsi_h%ham_omm,elsi_h%ovlp_omm,&
              elsi_h%new_overlap,elsi_h%energy_hdm,elsi_h%dm_omm,elsi_h%calc_ed,&
              elsi_h%eta,elsi_h%coeff_omm,elsi_h%coeff_ready,elsi_h%tdm_omm,&
              elsi_h%scale_kinetic,elsi_h%omm_flavor,1,1,elsi_h%min_tol,&
              elsi_h%omm_output,elsi_h%do_dealloc,"pzdbc","lap")
   case(REAL_VALUES)
      call omm(elsi_h%n_basis,elsi_h%n_states_omm,elsi_h%ham_omm,elsi_h%ovlp_omm,&
              elsi_h%new_overlap,elsi_h%energy_hdm,elsi_h%dm_omm,elsi_h%calc_ed,&
              elsi_h%eta,elsi_h%coeff_omm,elsi_h%coeff_ready,elsi_h%tdm_omm,&
              elsi_h%scale_kinetic,elsi_h%omm_flavor,1,1,elsi_h%min_tol,&
              elsi_h%omm_output,elsi_h%do_dealloc,"pddbc","lap")
   end select

   elsi_h%calc_ed = .false.

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished energy density matrix calculation')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

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
   elsi_h%min_tol = 1.0e-12_r8

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

   write(info_str,"(1X,' | ELPA steps before OMM         ',I10)") elsi_h%n_elpa_steps
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | libOMM flavor                 ',I10)") elsi_h%omm_flavor
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Tolerance of OMM minimization ',E10.2)") elsi_h%min_tol
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Use PSP sparse linear algebra ',L10)") elsi_h%use_psp
   call elsi_statement_print(info_str,elsi_h)

end subroutine

end module ELSI_OMM
