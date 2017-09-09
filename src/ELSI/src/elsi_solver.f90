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
!! This module provides solver interfaces to ELPA, libOMM, PEXSI, CheSS, SIPs.
!!
module ELSI_SOLVER

   use ELSI_CHESS
   use ELSI_CONSTANTS, only: ELPA,LIBOMM,PEXSI,CHESS,SIPS,REAL_VALUES,&
                             COMPLEX_VALUES,SINGLE_PROC,UNSET
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_ELPA
   use ELSI_MATCONV
   use ELSI_OMM
   use ELSI_PEXSI
   use ELSI_PRECISION, only: r8,i4
   use ELSI_SIPS
   use ELSI_UTILS
   use MATRIXSWITCH, only: m_allocate

   implicit none

   private

   ! Solver
   public :: elsi_ev_real
   public :: elsi_ev_complex
   public :: elsi_ev_real_sparse
   public :: elsi_dm_real
   public :: elsi_dm_complex
   public :: elsi_dm_real_sparse

contains

!>
!! This routine gets the energy.
!!
subroutine elsi_get_energy(e_h,energy)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(out)   :: energy !< Energy of the system

   real(kind=r8)    :: tmp_real
   integer(kind=i4) :: i_state
   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_get_energy"

   select case(e_h%solver)
   case(ELPA)
      energy = 0.0_r8

      do i_state = 1,e_h%n_states_solve
         energy = energy+e_h%i_weight*e_h%eval(i_state)*&
                     e_h%occ_num(i_state,e_h%i_spin,e_h%i_kpt)
      enddo
   case(LIBOMM)
      energy = 2.0_r8*e_h%energy_hdm*e_h%i_weight
   case(PEXSI)
      energy = e_h%energy_hdm*e_h%i_weight
   case(CHESS)
      energy = e_h%energy_hdm*e_h%i_weight
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",e_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen. Please choose"//&
               " ELPA, LIBOMM, or PEXSI solver. Exiting...",e_h,caller)
   end select

   if(e_h%n_spins*e_h%n_kpts > 1) then
      if(e_h%myid /= 0) then
         energy = 0.0_r8
      endif

      call MPI_Allreduce(energy,tmp_real,1,mpi_real8,mpi_sum,e_h%mpi_comm_all,&
              mpierr)

      energy = tmp_real
   endif

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors.
!!
subroutine elsi_ev_real(e_h,h_in,s_in,e_val_out,e_vec_out)

   implicit none

   type(elsi_handle) :: e_h                                  !< Handle
   real(kind=r8)     :: h_in(e_h%n_l_rows,e_h%n_l_cols)      !< Hamiltonian
   real(kind=r8)     :: s_in(e_h%n_l_rows,e_h%n_l_cols)      !< Overlap
   real(kind=r8)     :: e_val_out(e_h%n_basis)               !< Eigenvalues
   real(kind=r8)     :: e_vec_out(e_h%n_l_rows,e_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_real"

   call elsi_check_handle(e_h,caller)

   ! Update counter
   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   ! REAL case
   e_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(e_h,caller)

   select case(e_h%solver)
   case(ELPA)
      ! Set matrices
      call elsi_set_ham(e_h,h_in)
      call elsi_set_ovlp(e_h,s_in)
      call elsi_set_evec(e_h,e_vec_out)
      call elsi_set_eval(e_h,e_val_out)

      ! Solve
      if(e_h%parallel_mode == SINGLE_PROC) then
         call elsi_solve_evp_elpa_sp(e_h)
      else ! MULTI_PROC
         call elsi_solve_evp_elpa(e_h)
      endif
   case(LIBOMM)
      call elsi_stop(" LIBOMM is not an eigensolver. Choose ELPA or SIPS if"//&
              " needed. Exiting...",e_h,caller)
   case(PEXSI)
      call elsi_stop(" PEXSI is not an eigensolver. Choose ELPA or SIPS if"//&
              " needed. Exiting...",e_h,caller)
   case(CHESS)
      call elsi_stop(" CHESS is not an eigensolver. Choose ELPA or SIPS if"//&
              " needed. Exiting...",e_h,caller)
   case(SIPS)
      call elsi_print_sips_options(e_h)

      ! Initialize SIPs
      call elsi_init_sips(e_h)

      ! Convert 2D dense to 1D CSC
      if(.not. e_h%ovlp_is_unit .and. e_h%n_elsi_calls == 1) then
         call elsi_set_full_mat(e_h,s_in)
      endif
      call elsi_set_full_mat(e_h,h_in)
      call elsi_blacs_to_sips_hs(e_h,h_in,s_in)

      ! Set matrices
      call elsi_set_sparse_ham(e_h,e_h%ham_real_sips)
      call elsi_set_sparse_ovlp(e_h,e_h%ovlp_real_sips)
      call elsi_set_row_ind(e_h,e_h%row_ind_sips)
      call elsi_set_col_ptr(e_h,e_h%col_ptr_sips)
      call elsi_set_eval(e_h,e_val_out)
      call elsi_set_evec(e_h,e_vec_out)

      ! Solve
      call elsi_solve_evp_sips(e_h)

      ! Convert non-distributed dense to 2D dense
      call elsi_sips_to_blacs_ev(e_h)
   case default
      call elsi_stop(" No supported solver has been chosen. Exiting...",e_h,&
              caller)
   end select

   e_h%matrix_data_type = UNSET

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors.
!!
subroutine elsi_ev_complex(e_h,h_in,s_in,e_val_out,e_vec_out)

   implicit none

   type(elsi_handle) :: e_h                                  !< Handle
   complex(kind=r8)  :: h_in(e_h%n_l_rows,e_h%n_l_cols)      !< Hamiltonian
   complex(kind=r8)  :: s_in(e_h%n_l_rows,e_h%n_l_cols)      !< Overlap
   real(kind=r8)     :: e_val_out(e_h%n_basis)               !< Eigenvalues
   complex(kind=r8)  :: e_vec_out(e_h%n_l_rows,e_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_complex"

   call elsi_check_handle(e_h,caller)

   ! Update counter
   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   ! COMPLEX case
   e_h%matrix_data_type = COMPLEX_VALUES

   ! Safety check
   call elsi_check(e_h,caller)

   select case(e_h%solver)
   case(ELPA)
      ! Set matrices
      call elsi_set_ham(e_h,h_in)
      call elsi_set_ovlp(e_h,s_in)
      call elsi_set_evec(e_h,e_vec_out)
      call elsi_set_eval(e_h,e_val_out)

      ! Solve
      if(e_h%parallel_mode == SINGLE_PROC) then
         call elsi_solve_evp_elpa_sp(e_h)
      else ! MULTI_PROC
         call elsi_solve_evp_elpa(e_h)
      endif
   case(LIBOMM)
      call elsi_stop(" LIBOMM is not an eigensolver. Choose ELPA if needed."//&
              " Exiting...",e_h,caller)
   case(PEXSI)
      call elsi_stop(" PEXSI is not an eigensolver. Choose ELPA if needed."//&
              " Exiting...",e_h,caller)
   case(CHESS)
      call elsi_stop(" CHESS is not an eigensolver. Choose ELPA if needed."//&
              " Exiting...",e_h,caller)
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",e_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen. Exiting...",e_h,&
              caller)
   end select

   e_h%matrix_data_type = UNSET

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors.
!!
subroutine elsi_ev_real_sparse(e_h,h_in,s_in,e_val_out,e_vec_out)

   implicit none

   type(elsi_handle) :: e_h                                  !< Handle
   real(kind=r8)     :: h_in(e_h%nnz_l_sp)                   !< Hamiltonian
   real(kind=r8)     :: s_in(e_h%nnz_l_sp)                   !< Overlap
   real(kind=r8)     :: e_val_out(e_h%n_basis)               !< Eigenvalues
   real(kind=r8)     :: e_vec_out(e_h%n_l_rows,e_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_real_sparse"

   call elsi_check_handle(e_h,caller)

   ! Update counter
   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   ! REAL case
   e_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(e_h,caller)

   select case(e_h%solver)
   case(ELPA)
      ! Convert 1D CSC to 2D dense
      call elsi_sips_to_blacs_hs(e_h,h_in,s_in)

      ! Set matrices
      call elsi_set_ham(e_h,e_h%ham_real_elpa)
      call elsi_set_ovlp(e_h,e_h%ovlp_real_elpa)
      call elsi_set_evec(e_h,e_vec_out)
      call elsi_set_eval(e_h,e_val_out)

      ! Solve
      call elsi_solve_evp_elpa(e_h)
   case(LIBOMM)
      call elsi_stop(" LIBOMM is not an eigensolver. Choose ELPA if needed."//&
              " Exiting...",e_h,caller)
   case(PEXSI)
      call elsi_stop(" PEXSI is not an eigensolver. Choose ELPA if needed."//&
              " Exiting...",e_h,caller)
   case(CHESS)
      call elsi_stop(" CHESS is not an eigensolver. Choose ELPA if needed."//&
              " Exiting...",e_h,caller)
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",e_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen. Exiting...",e_h,&
              caller)
   end select

   e_h%matrix_data_type = UNSET

end subroutine

!>
!! This routine computes the density matrix.
!!
subroutine elsi_dm_real(e_h,h_in,s_in,d_out,energy_out)

   implicit none

   type(elsi_handle) :: e_h                              !< Handle
   real(kind=r8)     :: h_in(e_h%n_l_rows,e_h%n_l_cols)  !< Hamiltonian
   real(kind=r8)     :: s_in(e_h%n_l_rows,e_h%n_l_cols)  !< Overlap
   real(kind=r8)     :: d_out(e_h%n_l_rows,e_h%n_l_cols) !< Density matrix
   real(kind=r8)     :: energy_out                       !< Energy

   character*40, parameter :: caller = "elsi_dm_real"

   call elsi_check_handle(e_h,caller)

   ! Update counter
   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   ! REAL case
   e_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(e_h,caller)

   select case(e_h%solver)
   case(ELPA)
      ! Allocate
      if(.not. allocated(e_h%eval_elpa)) then
         call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(e_h%evec_real_elpa)) then
         call elsi_allocate(e_h,e_h%evec_real_elpa,e_h%n_l_rows,e_h%n_l_cols,&
                 "evec_real_elpa",caller)
      endif

      ! Set matrices
      call elsi_set_ham(e_h,h_in)
      call elsi_set_ovlp(e_h,s_in)
      call elsi_set_evec(e_h,e_h%evec_real_elpa)
      call elsi_set_eval(e_h,e_h%eval_elpa)
      call elsi_set_dm(e_h,d_out)

      ! Solve
      call elsi_solve_evp_elpa(e_h)

      ! Compute density matrix
      call elsi_compute_occ_elpa(e_h)
      call elsi_compute_dm_elpa(e_h)
      call elsi_get_energy(e_h,energy_out)

      e_h%mu_ready = .true.
   case(LIBOMM)
      call elsi_print_omm_options(e_h)

      if(e_h%n_elsi_calls <= e_h%n_elpa_steps) then
         if(e_h%n_elsi_calls == 1 .and. e_h%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(e_h,e_h%ovlp_real_omm,e_h%n_l_rows,e_h%n_l_cols,&
                    "ovlp_real_omm",caller)
            e_h%ovlp_real_omm = s_in
         endif

         ! Compute libOMM initial guess by ELPA
         e_h%solver = ELPA

         ! Allocate
         if(.not. allocated(e_h%eval_elpa)) then
            call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
         endif
         if(.not. allocated(e_h%evec_real_elpa)) then
            call elsi_allocate(e_h,e_h%evec_real_elpa,e_h%n_l_rows,&
                    e_h%n_l_cols,"evec_real_elpa",caller)
         endif

         ! Set matrices
         call elsi_set_ham(e_h,h_in)
         call elsi_set_ovlp(e_h,s_in)
         call elsi_set_evec(e_h,e_h%evec_real_elpa)
         call elsi_set_eval(e_h,e_h%eval_elpa)
         call elsi_set_dm(e_h,d_out)

         ! Solve
         call elsi_solve_evp_elpa(e_h)

         ! Compute density matrix
         call elsi_compute_occ_elpa(e_h)
         call elsi_compute_dm_elpa(e_h)
         call elsi_get_energy(e_h,energy_out)

         ! Switch back to libOMM
         e_h%solver = LIBOMM
      else ! ELPA is done
         if(allocated(e_h%ovlp_real_omm)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            s_in = e_h%ovlp_real_omm
            call elsi_deallocate(e_h,e_h%ovlp_real_omm,"ovlp_real_omm")
         endif

         ! Allocate
         if(.not. e_h%coeff%is_initialized) then
            call m_allocate(e_h%coeff,e_h%n_states_omm,e_h%n_basis,"pddbc")
         endif

         ! Set matrices
         call elsi_set_ham(e_h,h_in)
         call elsi_set_ovlp(e_h,s_in)
         call elsi_set_dm(e_h,d_out)

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(e_h%n_elpa_steps > 0 .and. &
            e_h%n_elsi_calls == e_h%n_elpa_steps+1) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,e_h%evec_real,1,1,&
                    e_h%sc_desc,0.0_r8,d_out,1,1,e_h%sc_desc)

            e_h%coeff%dval(1:e_h%coeff%iaux2(1),1:e_h%coeff%iaux2(2)) = &
               d_out(1:e_h%coeff%iaux2(1),1:e_h%coeff%iaux2(2))

            ! ELPA matrices are no longer needed
            if(associated(e_h%ham_real)) then
               nullify(e_h%ham_real)
            endif
            if(associated(e_h%ovlp_real)) then
               nullify(e_h%ovlp_real)
            endif
            if(associated(e_h%evec_real)) then
               nullify(e_h%evec_real)
            endif
            if(associated(e_h%eval)) then
               nullify(e_h%eval)
            endif
            if(associated(e_h%dm_real)) then
               nullify(e_h%dm_real)
            endif
            if(allocated(e_h%evec_real_elpa)) then
               call elsi_deallocate(e_h,e_h%evec_real_elpa,"evec_real_elpa")
            endif
            if(allocated(e_h%eval_elpa)) then
               call elsi_deallocate(e_h,e_h%eval_elpa,"eval_elpa")
            endif
            if(allocated(e_h%eval_all)) then
               call elsi_deallocate(e_h,e_h%eval_all,"eval_all")
            endif
            if(allocated(e_h%occ_num)) then
               call elsi_deallocate(e_h,e_h%occ_num,"occ_num")
            endif
            if(allocated(e_h%k_weight)) then
               call elsi_deallocate(e_h,e_h%k_weight,"k_weight")
            endif
         endif

         ! Solve
         call elsi_solve_evp_omm(e_h)

         e_h%dm_omm%dval = 2.0_r8*e_h%dm_omm%dval
         call elsi_get_energy(e_h,energy_out)
      endif
   case(PEXSI)
      call elsi_print_pexsi_options(e_h)

      ! Initialize PEXSI
      call elsi_init_pexsi(e_h)

      ! Convert 2D dense to 1D CSC
      if(.not. e_h%ovlp_is_unit .and. e_h%n_elsi_calls == 1) then
         call elsi_set_full_mat(e_h,s_in)
      endif
      call elsi_set_full_mat(e_h,h_in)
      call elsi_blacs_to_pexsi_hs(e_h,h_in,s_in)

      ! Allocate
      if(.not. allocated(e_h%dm_real_pexsi)) then
         call elsi_allocate(e_h,e_h%dm_real_pexsi,e_h%nnz_l_sp,"dm_real_pexsi",&
                 caller)
      endif
      e_h%dm_real_pexsi = 0.0_r8

      ! Set matrices
      call elsi_set_sparse_ham(e_h,e_h%ham_real_pexsi)
      call elsi_set_sparse_ovlp(e_h,e_h%ovlp_real_pexsi)
      call elsi_set_row_ind(e_h,e_h%row_ind_pexsi)
      call elsi_set_col_ptr(e_h,e_h%col_ptr_pexsi)
      call elsi_set_sparse_dm(e_h,e_h%dm_real_pexsi)

      ! Solve
      call elsi_solve_evp_pexsi(e_h)

      call elsi_pexsi_to_blacs_dm(e_h,d_out)
      call elsi_get_energy(e_h,energy_out)

      e_h%mu_ready = .true.
   case(CHESS)
      call elsi_print_chess_options(e_h)

      ! Convert 2D dense to non-distributed CSC
      if(.not. e_h%ovlp_is_unit .and. e_h%n_elsi_calls == 1) then
         call elsi_set_full_mat(e_h,s_in)
      endif
      call elsi_set_full_mat(e_h,h_in)
      call elsi_blacs_to_chess_hs(e_h,h_in,s_in)

      ! Initialize CheSS
      call elsi_init_chess(e_h)

      ! Set matrices
      call elsi_set_sparse_ham(e_h,e_h%ham_real_chess)
      call elsi_set_sparse_ovlp(e_h,e_h%ovlp_real_chess)
      call elsi_set_row_ind(e_h,e_h%row_ind_chess)
      call elsi_set_col_ptr(e_h,e_h%col_ptr_chess)

      ! Solve
      call elsi_solve_evp_chess(e_h)

      call elsi_chess_to_blacs_dm(e_h,d_out)
      call elsi_get_energy(e_h,energy_out)

      e_h%mu_ready = .true.
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",e_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen. Exiting...",e_h,&
              caller)
   end select

   e_h%matrix_data_type = UNSET
   e_h%edm_ready_real = .true.

end subroutine

!>
!! This routine computes the density matrix.
!!
subroutine elsi_dm_complex(e_h,h_in,s_in,d_out,energy_out)

   implicit none

   type(elsi_handle) :: e_h                              !< Handle
   complex(kind=r8)  :: h_in(e_h%n_l_rows,e_h%n_l_cols)  !< Hamiltonian
   complex(kind=r8)  :: s_in(e_h%n_l_rows,e_h%n_l_cols)  !< Overlap
   complex(kind=r8)  :: d_out(e_h%n_l_rows,e_h%n_l_cols) !< Density matrix
   real(kind=r8)     :: energy_out                       !< Energy

   character*40, parameter :: caller = "elsi_dm_complex"

   call elsi_check_handle(e_h,caller)

   ! Update counter
   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   ! COMPLEX case
   e_h%matrix_data_type = COMPLEX_VALUES

   ! Safety check
   call elsi_check(e_h,caller)

   select case(e_h%solver)
   case(ELPA)
      ! Allocate
      if(.not. allocated(e_h%eval_elpa)) then
         call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(e_h%evec_cmplx_elpa)) then
         call elsi_allocate(e_h,e_h%evec_cmplx_elpa,e_h%n_l_rows,e_h%n_l_cols,&
                 "evec_cmplx_elpa",caller)
      endif

      ! Set matrices
      call elsi_set_ham(e_h,h_in)
      call elsi_set_ovlp(e_h,s_in)
      call elsi_set_evec(e_h,e_h%evec_cmplx_elpa)
      call elsi_set_eval(e_h,e_h%eval_elpa)
      call elsi_set_dm(e_h,d_out)

      ! Solve
      call elsi_solve_evp_elpa(e_h)

      ! Compute density matrix
      call elsi_compute_occ_elpa(e_h)
      call elsi_compute_dm_elpa(e_h)
      call elsi_get_energy(e_h,energy_out)

      e_h%mu_ready = .true.
   case(LIBOMM)
      call elsi_print_omm_options(e_h)

      if(e_h%n_elsi_calls <= e_h%n_elpa_steps) then
         if(e_h%n_elsi_calls == 1 .and. e_h%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(e_h,e_h%ovlp_cmplx_omm,e_h%n_l_rows,&
                    e_h%n_l_cols,"ovlp_cmplx_omm",caller)
            e_h%ovlp_cmplx_omm = s_in
         endif

         ! Compute libOMM initial guess by ELPA
         e_h%solver = ELPA

         ! Allocate
         if(.not. allocated(e_h%eval_elpa)) then
            call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
         endif
         if(.not. allocated(e_h%evec_cmplx_elpa)) then
            call elsi_allocate(e_h,e_h%evec_cmplx_elpa,e_h%n_l_rows,&
                    e_h%n_l_cols,"evec_cmplx_elpa",caller)
         endif

         ! Set matrices
         call elsi_set_ham(e_h,h_in)
         call elsi_set_ovlp(e_h,s_in)
         call elsi_set_evec(e_h,e_h%evec_cmplx_elpa)
         call elsi_set_eval(e_h,e_h%eval_elpa)
         call elsi_set_dm(e_h,d_out)

         ! Solve
         call elsi_solve_evp_elpa(e_h)

         ! Compute density matrix
         call elsi_compute_occ_elpa(e_h)
         call elsi_compute_dm_elpa(e_h)
         call elsi_get_energy(e_h,energy_out)

         ! Switch back to libOMM
         e_h%solver = LIBOMM
      else ! ELPA is done
         if(allocated(e_h%ovlp_cmplx_omm)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            s_in = e_h%ovlp_cmplx_omm
            call elsi_deallocate(e_h,e_h%ovlp_cmplx_omm,"ovlp_cmplx_omm")
         endif

         ! Allocate
         if(.not. e_h%coeff%is_initialized) then
            call m_allocate(e_h%coeff,e_h%n_states_omm,e_h%n_basis,"pzdbc")
         endif

         ! Set matrices
         call elsi_set_ham(e_h,h_in)
         call elsi_set_ovlp(e_h,s_in)
         call elsi_set_dm(e_h,d_out)

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(e_h%n_elpa_steps > 0 .and. &
            e_h%n_elsi_calls == e_h%n_elpa_steps+1) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),&
                    e_h%evec_cmplx,1,1,e_h%sc_desc,(0.0_r8,0.0_r8),d_out,1,1,&
                    e_h%sc_desc)

            e_h%coeff%zval(1:e_h%coeff%iaux2(1),1:e_h%coeff%iaux2(2)) = &
               d_out(1:e_h%coeff%iaux2(1),1:e_h%coeff%iaux2(2))

            ! ELPA matrices are no longer needed
            if(associated(e_h%ham_cmplx)) then
               nullify(e_h%ham_cmplx)
            endif
            if(associated(e_h%ovlp_cmplx)) then
               nullify(e_h%ovlp_cmplx)
            endif
            if(associated(e_h%evec_cmplx)) then
               nullify(e_h%evec_cmplx)
            endif
            if(associated(e_h%eval)) then
               nullify(e_h%eval)
            endif
            if(associated(e_h%dm_cmplx)) then
               nullify(e_h%dm_cmplx)
            endif
            if(allocated(e_h%evec_cmplx_elpa)) then
               call elsi_deallocate(e_h,e_h%evec_cmplx_elpa,"evec_cmplx_elpa")
            endif
            if(allocated(e_h%eval_elpa)) then
               call elsi_deallocate(e_h,e_h%eval_elpa,"eval_elpa")
            endif
            if(allocated(e_h%eval_all)) then
               call elsi_deallocate(e_h,e_h%eval_all,"eval_all")
            endif
            if(allocated(e_h%occ_num)) then
               call elsi_deallocate(e_h,e_h%occ_num,"occ_num")
            endif
            if(allocated(e_h%k_weight)) then
               call elsi_deallocate(e_h,e_h%k_weight,"k_weight")
            endif
         endif

         ! Solve
         call elsi_solve_evp_omm(e_h)

         e_h%dm_omm%zval = 2.0_r8*e_h%dm_omm%zval
         call elsi_get_energy(e_h,energy_out)
      endif
   case(PEXSI)
      call elsi_print_pexsi_options(e_h)

      ! Initialize PEXSI
      call elsi_init_pexsi(e_h)

      ! Convert 2D dense to 1D CSC
      if(.not. e_h%ovlp_is_unit .and. e_h%n_elsi_calls == 1) then
         call elsi_set_full_mat(e_h,s_in)
      endif
      call elsi_set_full_mat(e_h,h_in)
      call elsi_blacs_to_pexsi_hs(e_h,h_in,s_in)

      ! Allocate
      if(.not. allocated(e_h%dm_cmplx_pexsi)) then
         call elsi_allocate(e_h,e_h%dm_cmplx_pexsi,e_h%nnz_l_sp,&
                 "dm_cmplx_pexsi",caller)
      endif
      e_h%dm_cmplx_pexsi = (0.0_r8,0.0_r8)

      ! Set matrices
      call elsi_set_sparse_ham(e_h,e_h%ham_cmplx_pexsi)
      call elsi_set_sparse_ovlp(e_h,e_h%ovlp_cmplx_pexsi)
      call elsi_set_row_ind(e_h,e_h%row_ind_pexsi)
      call elsi_set_col_ptr(e_h,e_h%col_ptr_pexsi)
      call elsi_set_sparse_dm(e_h,e_h%dm_cmplx_pexsi)

      ! Solve
      call elsi_solve_evp_pexsi(e_h)

      call elsi_pexsi_to_blacs_dm(e_h,d_out)
      call elsi_get_energy(e_h,energy_out)

      e_h%mu_ready = .true.
   case(CHESS)
      call elsi_stop(" CHESS not yet implemented. Exiting...",e_h,caller)
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",e_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen. Exiting...",e_h,&
              caller)
   end select

   e_h%matrix_data_type = UNSET
   e_h%edm_ready_cmplx = .true.

end subroutine

!>
!! This routine computes the density matrix.
!!
subroutine elsi_dm_real_sparse(e_h,h_in,s_in,d_out,energy_out)

   implicit none

   type(elsi_handle) :: e_h                 !< Handle
   real(kind=r8)     :: h_in(e_h%nnz_l_sp)  !< Hamiltonian
   real(kind=r8)     :: s_in(e_h%nnz_l_sp)  !< Overlap
   real(kind=r8)     :: d_out(e_h%nnz_l_sp) !< Density matrix
   real(kind=r8)     :: energy_out          !< Energy

   character*40, parameter :: caller = "elsi_dm_real_sparse"

   call elsi_check_handle(e_h,caller)

   ! Update counter
   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   ! REAL case
   e_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(e_h,caller)

   select case(e_h%solver)
   case(ELPA)
      ! Convert 1D CSC to 2D dense
      call elsi_sips_to_blacs_hs(e_h,h_in,s_in)

      ! Allocate
      if(.not. allocated(e_h%eval_elpa)) then
         call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(e_h%evec_real_elpa)) then
         call elsi_allocate(e_h,e_h%evec_real_elpa,e_h%n_l_rows,e_h%n_l_cols,&
                 "evec_real_elpa",caller)
      endif
      if(.not. allocated(e_h%dm_real_elpa)) then
         call elsi_allocate(e_h,e_h%dm_real_elpa,e_h%n_l_rows,e_h%n_l_cols,&
                 "dm_real_elpa",caller)
      endif

      ! Set matrices
      call elsi_set_ham(e_h,e_h%ham_real_elpa)
      call elsi_set_ovlp(e_h,e_h%ovlp_real_elpa)
      call elsi_set_evec(e_h,e_h%evec_real_elpa)
      call elsi_set_eval(e_h,e_h%eval_elpa)
      call elsi_set_dm(e_h,e_h%dm_real_elpa)

      ! Solve eigenvalue problem
      call elsi_solve_evp_elpa(e_h)

      ! Compute density matrix
      call elsi_compute_occ_elpa(e_h)
      call elsi_compute_dm_elpa(e_h)
      call elsi_blacs_to_sips_dm(e_h,d_out)
      call elsi_get_energy(e_h,energy_out)

      e_h%mu_ready = .true.
   case(LIBOMM)
      call elsi_print_omm_options(e_h)

      ! Convert 1D CSC to 2D dense
      call elsi_sips_to_blacs_hs(e_h,h_in,s_in)

      if(e_h%n_elsi_calls <= e_h%n_elpa_steps) then
         if(e_h%n_elsi_calls == 1 .and. e_h%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(e_h,e_h%ovlp_real_omm,e_h%n_l_rows,&
                    e_h%n_l_cols,"ovlp_real_omm",caller)
            e_h%ovlp_real_omm = e_h%ovlp_real_elpa
         endif

         ! Compute libOMM initial guess by ELPA
         e_h%solver = ELPA

         ! Allocate
         if(.not. allocated(e_h%eval_elpa)) then
            call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
         endif
         if(.not. allocated(e_h%evec_real_elpa)) then
            call elsi_allocate(e_h,e_h%evec_real_elpa,e_h%n_l_rows,&
                    e_h%n_l_cols,"evec_real_elpa",caller)
         endif
         if(.not. allocated(e_h%dm_real_elpa)) then
            call elsi_allocate(e_h,e_h%dm_real_elpa,e_h%n_l_rows,e_h%n_l_cols,&
                    "dm_real_elpa",caller)
         endif

         ! Set matrices
         call elsi_set_ham(e_h,e_h%ham_real_elpa)
         call elsi_set_ovlp(e_h,e_h%ovlp_real_elpa)
         call elsi_set_evec(e_h,e_h%evec_real_elpa)
         call elsi_set_eval(e_h,e_h%eval_elpa)
         call elsi_set_dm(e_h,e_h%dm_real_elpa)

         ! Solve eigenvalue problem
         call elsi_solve_evp_elpa(e_h)

         ! Compute density matrix
         call elsi_compute_occ_elpa(e_h)
         call elsi_compute_dm_elpa(e_h)
         call elsi_blacs_to_sips_dm(e_h,d_out)
         call elsi_get_energy(e_h,energy_out)

         ! Switch back to libOMM
         e_h%solver = LIBOMM
      else ! ELPA is done
         if(allocated(e_h%ovlp_real_omm)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            e_h%ovlp_real_elpa = e_h%ovlp_real_omm
            call elsi_deallocate(e_h,e_h%ovlp_real_omm,"ovlp_real_omm")
         endif

         ! Allocate
         if(.not. e_h%coeff%is_initialized) then
            call m_allocate(e_h%coeff,e_h%n_states_omm,e_h%n_basis,"pddbc")
         endif
         if(.not. allocated(e_h%dm_real_elpa)) then
            call elsi_allocate(e_h,e_h%dm_real_elpa,e_h%n_l_rows,e_h%n_l_cols,&
                    "dm_real_elpa",caller)
         endif

         ! Set matrices
         call elsi_set_ham(e_h,e_h%ham_real_elpa)
         call elsi_set_ovlp(e_h,e_h%ovlp_real_elpa)
         call elsi_set_dm(e_h,e_h%dm_real_elpa)

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(e_h%n_elpa_steps > 0 .and. &
            e_h%n_elsi_calls == e_h%n_elpa_steps+1) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,e_h%evec_real,1,1,&
                    e_h%sc_desc,0.0_r8,e_h%dm_real_elpa,1,1,e_h%sc_desc)

            e_h%coeff%dval(1:e_h%coeff%iaux2(1),1:e_h%coeff%iaux2(2)) = &
               e_h%dm_real_elpa(1:e_h%coeff%iaux2(1),1:e_h%coeff%iaux2(2))

            ! ELPA matrices are no longer needed
            if(associated(e_h%ham_real)) then
               nullify(e_h%ham_real)
            endif
            if(associated(e_h%ovlp_real)) then
               nullify(e_h%ovlp_real)
            endif
            if(associated(e_h%evec_real)) then
               nullify(e_h%evec_real)
            endif
            if(associated(e_h%dm_real)) then
               nullify(e_h%dm_real)
            endif
            if(associated(e_h%eval)) then
               nullify(e_h%eval)
            endif
            if(allocated(e_h%evec_real_elpa)) then
               call elsi_deallocate(e_h,e_h%evec_real_elpa,"evec_real_elpa")
            endif
            if(allocated(e_h%eval_elpa)) then
               call elsi_deallocate(e_h,e_h%eval_elpa,"eval_elpa")
            endif
            if(allocated(e_h%occ_num)) then
               call elsi_deallocate(e_h,e_h%occ_num,"occ_num")
            endif
         endif

         ! Solve
         call elsi_solve_evp_omm(e_h)

         e_h%dm_omm%dval = 2.0_r8*e_h%dm_omm%dval
         call elsi_blacs_to_sips_dm(e_h,d_out)
         call elsi_get_energy(e_h,energy_out)
      endif
   case(PEXSI)
      call elsi_print_pexsi_options(e_h)

      ! Set matrices
      call elsi_set_sparse_ham(e_h,h_in)
      call elsi_set_sparse_ovlp(e_h,s_in)
      call elsi_set_sparse_dm(e_h,d_out)

      ! Initialize PEXSI
      call elsi_init_pexsi(e_h)

      ! Solve
      call elsi_solve_evp_pexsi(e_h)

      call elsi_get_energy(e_h,energy_out)

      e_h%mu_ready = .true.
   case(CHESS)
      call elsi_stop(" CHESS not yet implemented. Exiting...",e_h,caller)
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",e_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen. Exiting...",e_h,&
              caller)
   end select

   e_h%matrix_data_type = UNSET
   e_h%edm_ready_real = .true.

end subroutine

end module ELSI_SOLVER
