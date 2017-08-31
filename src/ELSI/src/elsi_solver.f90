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
subroutine elsi_get_energy(elsi_h,energy)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(out)   :: energy !< Energy of the system

   real(kind=r8)    :: tmp_real
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: mpierr
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_get_energy"

   select case(elsi_h%solver)
   case(ELPA)
      energy = 0.0_r8

      do i_state = 1,elsi_h%n_states_solve
         energy = energy+elsi_h%i_weight*elsi_h%eval(i_state)*&
                     elsi_h%occ_num(i_state,elsi_h%i_spin,elsi_h%i_kpt)
      enddo
   case(LIBOMM)
      energy = 2.0_r8*elsi_h%energy_hdm*elsi_h%i_weight
   case(PEXSI)
      energy = elsi_h%energy_hdm*elsi_h%i_weight
   case(CHESS)
      energy = elsi_h%energy_hdm*elsi_h%i_weight
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen."//&
               " Please choose ELPA, LIBOMM, or PEXSI solver."//&
               " Exiting...",elsi_h,caller)
   end select

   if(elsi_h%n_spins*elsi_h%n_kpts > 1) then
      if(elsi_h%myid /= 0) then
         energy = 0.0_r8
      endif

      call MPI_Allreduce(energy,tmp_real,1,mpi_real8,mpi_sum,&
              elsi_h%mpi_comm_all,mpierr)

      energy = tmp_real
   endif

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors.
!!
subroutine elsi_ev_real(elsi_h,h_in,s_in,e_val_out,e_vec_out)

   implicit none

   type(elsi_handle) :: elsi_h                                     !< Handle
   real(kind=r8)     :: h_in(elsi_h%n_l_rows,elsi_h%n_l_cols)      !< Hamiltonian
   real(kind=r8)     :: s_in(elsi_h%n_l_rows,elsi_h%n_l_cols)      !< Overlap
   real(kind=r8)     :: e_val_out(elsi_h%n_basis)                  !< Eigenvalues
   real(kind=r8)     :: e_vec_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_real"

   call elsi_check_handle(elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! REAL case
   elsi_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   select case(elsi_h%solver)
   case(ELPA)
      ! Set matrices
      call elsi_set_ham(elsi_h,h_in)
      call elsi_set_ovlp(elsi_h,s_in)
      call elsi_set_evec(elsi_h,e_vec_out)
      call elsi_set_eval(elsi_h,e_val_out)

      ! Solve
      if(elsi_h%parallel_mode == SINGLE_PROC) then
         call elsi_solve_evp_elpa_sp(elsi_h)
      else ! MULTI_PROC
         call elsi_solve_evp_elpa(elsi_h)
      endif
   case(LIBOMM)
      call elsi_stop(" LIBOMM is not an eigensolver. Choose ELPA or"//&
              " SIPS if needed. Exiting...",elsi_h,caller)
   case(PEXSI)
      call elsi_stop(" PEXSI is not an eigensolver. Choose ELPA or"//&
              " SIPS if needed. Exiting...",elsi_h,caller)
   case(CHESS)
      call elsi_stop(" CHESS is not an eigensolver. Choose ELPA or"//&
              " SIPS if needed. Exiting...",elsi_h,caller)
   case(SIPS)
      call elsi_print_sips_options(elsi_h)

      ! Initialize SIPs
      call elsi_init_sips(elsi_h)

      ! Convert BLACS H and S to SIPs format
      if((.not. elsi_h%ovlp_is_unit) .and. (elsi_h%n_elsi_calls == 1)) then
         call elsi_set_full_mat(elsi_h,s_in)
      endif
      call elsi_set_full_mat(elsi_h,h_in)
      call elsi_blacs_to_sips_hs(elsi_h,h_in,s_in)

      ! Set matrices
      call elsi_set_sparse_ham(elsi_h,elsi_h%ham_real_sips)
      call elsi_set_sparse_ovlp(elsi_h,elsi_h%ovlp_real_sips)
      call elsi_set_row_ind(elsi_h,elsi_h%row_ind_sips)
      call elsi_set_col_ptr(elsi_h,elsi_h%col_ptr_sips)
      call elsi_set_eval(elsi_h,e_val_out)
      call elsi_set_evec(elsi_h,e_vec_out)

      ! Solve
      call elsi_solve_evp_sips(elsi_h)

      ! Convert SIPs eigenvectors to BLACS format
      call elsi_sips_to_blacs_ev(elsi_h)
   case default
      call elsi_stop(" No supported solver has been chosen."//&
              " Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors.
!!
subroutine elsi_ev_complex(elsi_h,h_in,s_in,e_val_out,e_vec_out)

   implicit none

   type(elsi_handle) :: elsi_h                                     !< Handle
   complex(kind=r8)  :: h_in(elsi_h%n_l_rows,elsi_h%n_l_cols)      !< Hamiltonian
   complex(kind=r8)  :: s_in(elsi_h%n_l_rows,elsi_h%n_l_cols)      !< Overlap
   real(kind=r8)     :: e_val_out(elsi_h%n_basis)                  !< Eigenvalues
   complex(kind=r8)  :: e_vec_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_complex"

   call elsi_check_handle(elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! COMPLEX case
   elsi_h%matrix_data_type = COMPLEX_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   select case(elsi_h%solver)
   case(ELPA)
      ! Set matrices
      call elsi_set_ham(elsi_h,h_in)
      call elsi_set_ovlp(elsi_h,s_in)
      call elsi_set_evec(elsi_h,e_vec_out)
      call elsi_set_eval(elsi_h,e_val_out)

      ! Solve
      if(elsi_h%parallel_mode == SINGLE_PROC) then
         call elsi_solve_evp_elpa_sp(elsi_h)
      else ! MULTI_PROC
         call elsi_solve_evp_elpa(elsi_h)
      endif
   case(LIBOMM)
      call elsi_stop(" LIBOMM is not an eigensolver. Choose ELPA if"//&
              " needed. Exiting...",elsi_h,caller)
   case(PEXSI)
      call elsi_stop(" PEXSI is not an eigensolver. Choose ELPA if"//&
              " needed. Exiting...",elsi_h,caller)
   case(CHESS)
      call elsi_stop(" CHESS is not an eigensolver. Choose ELPA if"//&
              " needed. Exiting...",elsi_h,caller)
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen."//&
              " Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors.
!!
subroutine elsi_ev_real_sparse(elsi_h,h_in,s_in,e_val_out,e_vec_out)

   implicit none

   type(elsi_handle) :: elsi_h                                     !< Handle
   real(kind=r8)     :: h_in(elsi_h%nnz_l_sp)                      !< Hamiltonian
   real(kind=r8)     :: s_in(elsi_h%nnz_l_sp)                      !< Overlap
   real(kind=r8)     :: e_val_out(elsi_h%n_basis)                  !< Eigenvalues
   real(kind=r8)     :: e_vec_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_real_sparse"

   call elsi_check_handle(elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! REAL case
   elsi_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   select case(elsi_h%solver)
   case(ELPA)
      ! Convert PEXSI H and S to BLACS format
      call elsi_pexsi_to_blacs_hs(elsi_h,h_in,s_in)

      ! Set matrices
      call elsi_set_ham(elsi_h,elsi_h%ham_real_elpa)
      call elsi_set_ovlp(elsi_h,elsi_h%ovlp_real_elpa)
      call elsi_set_evec(elsi_h,e_vec_out)
      call elsi_set_eval(elsi_h,e_val_out)

      ! Solve
      call elsi_solve_evp_elpa(elsi_h)
   case(LIBOMM)
      call elsi_stop(" LIBOMM is not an eigensolver. Choose ELPA if"//&
              " needed. Exiting...",elsi_h,caller)
   case(PEXSI)
      call elsi_stop(" PEXSI is not an eigensolver. Choose ELPA if"//&
              " needed. Exiting...",elsi_h,caller)
   case(CHESS)
      call elsi_stop(" CHESS is not an eigensolver. Choose ELPA if"//&
              " needed. Exiting...",elsi_h,caller)
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen."//&
              " Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET

end subroutine

!>
!! This routine computes the density matrix.
!!
subroutine elsi_dm_real(elsi_h,h_in,s_in,d_out,energy_out)

   implicit none

   type(elsi_handle) :: elsi_h                                 !< Handle
   real(kind=r8)     :: h_in(elsi_h%n_l_rows,elsi_h%n_l_cols)  !< Hamiltonian
   real(kind=r8)     :: s_in(elsi_h%n_l_rows,elsi_h%n_l_cols)  !< Overlap
   real(kind=r8)     :: d_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix
   real(kind=r8)     :: energy_out                             !< Energy

   character*40, parameter :: caller = "elsi_dm_real"

   call elsi_check_handle(elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! REAL case
   elsi_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   select case(elsi_h%solver)
   case(ELPA)
      ! Allocate
      if(.not. allocated(elsi_h%eval_elpa)) then
         call elsi_allocate(elsi_h,elsi_h%eval_elpa,elsi_h%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(elsi_h%evec_real_elpa)) then
         call elsi_allocate(elsi_h,elsi_h%evec_real_elpa,elsi_h%n_l_rows,&
                 elsi_h%n_l_cols,"evec_real_elpa",caller)
      endif

      ! Set matrices
      call elsi_set_ham(elsi_h,h_in)
      call elsi_set_ovlp(elsi_h,s_in)
      call elsi_set_evec(elsi_h,elsi_h%evec_real_elpa)
      call elsi_set_eval(elsi_h,elsi_h%eval_elpa)
      call elsi_set_dm(elsi_h,d_out)

      ! Solve
      call elsi_solve_evp_elpa(elsi_h)

      ! Compute density matrix
      call elsi_compute_occ_elpa(elsi_h)
      call elsi_compute_dm_elpa(elsi_h)
      call elsi_get_energy(elsi_h,energy_out)

      elsi_h%mu_ready = .true.
   case(LIBOMM)
      call elsi_print_omm_options(elsi_h)

      if(elsi_h%n_elsi_calls <= elsi_h%n_elpa_steps) then
         if((elsi_h%n_elsi_calls == 1) .and. (elsi_h%omm_flavor == 0)) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(elsi_h,elsi_h%ovlp_real_omm,elsi_h%n_l_rows,&
                    elsi_h%n_l_cols,"ovlp_real_omm",caller)
            elsi_h%ovlp_real_omm = s_in
         endif

         ! Compute libOMM initial guess by ELPA
         elsi_h%solver = ELPA

         ! Allocate
         if(.not. allocated(elsi_h%eval_elpa)) then
            call elsi_allocate(elsi_h,elsi_h%eval_elpa,elsi_h%n_basis,&
                    "eval_elpa",caller)
         endif
         if(.not. allocated(elsi_h%evec_real_elpa)) then
            call elsi_allocate(elsi_h,elsi_h%evec_real_elpa,elsi_h%n_l_rows,&
                    elsi_h%n_l_cols,"evec_real_elpa",caller)
         endif

         ! Set matrices
         call elsi_set_ham(elsi_h,h_in)
         call elsi_set_ovlp(elsi_h,s_in)
         call elsi_set_evec(elsi_h,elsi_h%evec_real_elpa)
         call elsi_set_eval(elsi_h,elsi_h%eval_elpa)
         call elsi_set_dm(elsi_h,d_out)

         ! Solve
         call elsi_solve_evp_elpa(elsi_h)

         ! Compute density matrix
         call elsi_compute_occ_elpa(elsi_h)
         call elsi_compute_dm_elpa(elsi_h)
         call elsi_get_energy(elsi_h,energy_out)

         ! Switch back to libOMM here to guarantee elsi_customize_omm
         elsi_h%solver = LIBOMM
      else ! ELPA is done
         if(allocated(elsi_h%ovlp_real_omm)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            s_in = elsi_h%ovlp_real_omm
            call elsi_deallocate(elsi_h,elsi_h%ovlp_real_omm,"ovlp_real_omm")

            call elsi_set_full_mat(elsi_h,s_in)
         endif

         ! Allocate
         if(.not. elsi_h%coeff_omm%is_initialized) then
            call m_allocate(elsi_h%coeff_omm,elsi_h%n_states_omm,elsi_h%n_basis,"pddbc")
         endif

         ! Set matrices
         call elsi_set_ham(elsi_h,h_in)
         call elsi_set_ovlp(elsi_h,s_in)
         call elsi_set_dm(elsi_h,d_out)

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if((elsi_h%n_elpa_steps > 0) .and. &
            (elsi_h%n_elsi_calls == elsi_h%n_elpa_steps+1)) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pdtran(elsi_h%n_basis,elsi_h%n_basis,1.0_r8,elsi_h%evec_real,&
                    1,1,elsi_h%sc_desc,0.0_r8,d_out,1,1,elsi_h%sc_desc)

            elsi_h%coeff_omm%dval(1:elsi_h%coeff_omm%iaux2(1),1:elsi_h%coeff_omm%iaux2(2)) = &
               d_out(1:elsi_h%coeff_omm%iaux2(1),1:elsi_h%coeff_omm%iaux2(2))

            ! ELPA matrices are no longer needed
            if(associated(elsi_h%ham_real)) then
               nullify(elsi_h%ham_real)
            endif
            if(associated(elsi_h%ovlp_real)) then
               nullify(elsi_h%ovlp_real)
            endif
            if(associated(elsi_h%evec_real)) then
               nullify(elsi_h%evec_real)
            endif
            if(associated(elsi_h%eval)) then
               nullify(elsi_h%eval)
            endif
            if(associated(elsi_h%dm_real)) then
               nullify(elsi_h%dm_real)
            endif
            if(allocated(elsi_h%evec_real_elpa)) then
               call elsi_deallocate(elsi_h,elsi_h%evec_real_elpa,"evec_real_elpa")
            endif
            if(allocated(elsi_h%eval_elpa)) then
               call elsi_deallocate(elsi_h,elsi_h%eval_elpa,"eval_elpa")
            endif
            if(allocated(elsi_h%eval_all)) then
               call elsi_deallocate(elsi_h,elsi_h%eval_all,"eval_all")
            endif
            if(allocated(elsi_h%occ_num)) then
               call elsi_deallocate(elsi_h,elsi_h%occ_num,"occ_num")
            endif
            if(allocated(elsi_h%k_weight)) then
               call elsi_deallocate(elsi_h,elsi_h%k_weight,"k_weight")
            endif
         endif

         ! Solve
         call elsi_solve_evp_omm(elsi_h)

         elsi_h%dm_omm%dval = 2.0_r8*elsi_h%dm_omm%dval
         call elsi_get_energy(elsi_h,energy_out)
      endif
   case(PEXSI)
      call elsi_print_pexsi_options(elsi_h)

      ! Initialize PEXSI
      call elsi_init_pexsi(elsi_h)

      ! Convert BLACS H and S to PEXSI format
      if((.not. elsi_h%ovlp_is_unit) .and. (elsi_h%n_elsi_calls == 1)) then
         call elsi_set_full_mat(elsi_h,s_in)
      endif
      call elsi_set_full_mat(elsi_h,h_in)
      call elsi_blacs_to_pexsi_hs(elsi_h,h_in,s_in)

      ! Allocate
      if(.not. allocated(elsi_h%dm_real_pexsi)) then
         call elsi_allocate(elsi_h,elsi_h%dm_real_pexsi,elsi_h%nnz_l_sp,&
                 "dm_real_pexsi",caller)
      endif
      elsi_h%dm_real_pexsi = 0.0_r8

      ! Set matrices
      call elsi_set_sparse_ham(elsi_h,elsi_h%ham_real_pexsi)
      call elsi_set_sparse_ovlp(elsi_h,elsi_h%ovlp_real_pexsi)
      call elsi_set_row_ind(elsi_h,elsi_h%row_ind_pexsi)
      call elsi_set_col_ptr(elsi_h,elsi_h%col_ptr_pexsi)
      call elsi_set_sparse_dm(elsi_h,elsi_h%dm_real_pexsi)

      ! Solve
      call elsi_solve_evp_pexsi(elsi_h)

      ! Convert PEXSI density matrix to BLACS format
      call elsi_pexsi_to_blacs_dm(elsi_h,d_out)
      call elsi_get_energy(elsi_h,energy_out)

      elsi_h%mu_ready = .true.
   case(CHESS)
      call elsi_print_chess_options(elsi_h)

      ! Convert BLACS H and S to CheSS format
      if((.not. elsi_h%ovlp_is_unit) .and. (elsi_h%n_elsi_calls == 1)) then
         call elsi_set_full_mat(elsi_h,s_in)
      endif
      call elsi_set_full_mat(elsi_h,h_in)
      call elsi_blacs_to_chess_hs(elsi_h,h_in,s_in)

      ! Initialize CheSS
      call elsi_init_chess(elsi_h)

      ! Set matrices
      call elsi_set_sparse_ham(elsi_h,elsi_h%ham_real_chess)
      call elsi_set_sparse_ovlp(elsi_h,elsi_h%ovlp_real_chess)
      call elsi_set_row_ind(elsi_h,elsi_h%row_ind_chess)
      call elsi_set_col_ptr(elsi_h,elsi_h%col_ptr_chess)

      ! Solve
      call elsi_solve_evp_chess(elsi_h)

      ! Convert CheSS density matrix to BLACS format
      call elsi_chess_to_blacs_dm(elsi_h,d_out)
      call elsi_get_energy(elsi_h,energy_out)

      elsi_h%mu_ready = .true.
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen."//&
              " Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET
   elsi_h%edm_ready_real = .true.

end subroutine

!>
!! This routine computes the density matrix.
!!
subroutine elsi_dm_complex(elsi_h,h_in,s_in,d_out,energy_out)

   implicit none

   type(elsi_handle) :: elsi_h                                 !< Handle
   complex(kind=r8)  :: h_in(elsi_h%n_l_rows,elsi_h%n_l_cols)  !< Hamiltonian
   complex(kind=r8)  :: s_in(elsi_h%n_l_rows,elsi_h%n_l_cols)  !< Overlap
   complex(kind=r8)  :: d_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix
   real(kind=r8)     :: energy_out                             !< Energy

   character*40, parameter :: caller = "elsi_dm_complex"

   call elsi_check_handle(elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! COMPLEX case
   elsi_h%matrix_data_type = COMPLEX_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   select case(elsi_h%solver)
   case(ELPA)
      ! Allocate
      if(.not. allocated(elsi_h%eval_elpa)) then
         call elsi_allocate(elsi_h,elsi_h%eval_elpa,elsi_h%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(elsi_h%evec_complex_elpa)) then
         call elsi_allocate(elsi_h,elsi_h%evec_complex_elpa,elsi_h%n_l_rows,&
                 elsi_h%n_l_cols,"evec_complex_elpa",caller)
      endif

      ! Set matrices
      call elsi_set_ham(elsi_h,h_in)
      call elsi_set_ovlp(elsi_h,s_in)
      call elsi_set_evec(elsi_h,elsi_h%evec_complex_elpa)
      call elsi_set_eval(elsi_h,elsi_h%eval_elpa)
      call elsi_set_dm(elsi_h,d_out)

      ! Solve
      call elsi_solve_evp_elpa(elsi_h)

      ! Compute density matrix
      call elsi_compute_occ_elpa(elsi_h)
      call elsi_compute_dm_elpa(elsi_h)
      call elsi_get_energy(elsi_h,energy_out)

      elsi_h%mu_ready = .true.
   case(LIBOMM)
      call elsi_print_omm_options(elsi_h)

      if(elsi_h%n_elsi_calls <= elsi_h%n_elpa_steps) then
         if((elsi_h%n_elsi_calls == 1) .and. (elsi_h%omm_flavor == 0)) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(elsi_h,elsi_h%ovlp_complex_omm,elsi_h%n_l_rows,&
                    elsi_h%n_l_cols,"ovlp_complex_omm",caller)
            elsi_h%ovlp_complex_omm = s_in
         endif

         ! Compute libOMM initial guess by ELPA
         elsi_h%solver = ELPA

         ! Allocate
         if(.not. allocated(elsi_h%eval_elpa)) then
            call elsi_allocate(elsi_h,elsi_h%eval_elpa,elsi_h%n_basis,&
                    "eval_elpa",caller)
         endif
         if(.not. allocated(elsi_h%evec_complex_elpa)) then
            call elsi_allocate(elsi_h,elsi_h%evec_complex_elpa,elsi_h%n_l_rows,&
                    elsi_h%n_l_cols,"evec_complex_elpa",caller)
         endif

         ! Set matrices
         call elsi_set_ham(elsi_h,h_in)
         call elsi_set_ovlp(elsi_h,s_in)
         call elsi_set_evec(elsi_h,elsi_h%evec_complex_elpa)
         call elsi_set_eval(elsi_h,elsi_h%eval_elpa)
         call elsi_set_dm(elsi_h,d_out)

         ! Solve
         call elsi_solve_evp_elpa(elsi_h)

         ! Compute density matrix
         call elsi_compute_occ_elpa(elsi_h)
         call elsi_compute_dm_elpa(elsi_h)
         call elsi_get_energy(elsi_h,energy_out)

         ! Switch back to libOMM here to guarantee elsi_customize_omm
         elsi_h%solver = LIBOMM
      else ! ELPA is done
         if(allocated(elsi_h%ovlp_complex_omm)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            s_in = elsi_h%ovlp_complex_omm
            call elsi_deallocate(elsi_h,elsi_h%ovlp_complex_omm,"ovlp_complex_omm")

            call elsi_set_full_mat(elsi_h,s_in)
         endif

         ! Allocate
         if(.not. elsi_h%coeff_omm%is_initialized) then
            call m_allocate(elsi_h%coeff_omm,elsi_h%n_states_omm,elsi_h%n_basis,"pzdbc")
         endif

         ! Set matrices
         call elsi_set_ham(elsi_h,h_in)
         call elsi_set_ovlp(elsi_h,s_in)
         call elsi_set_dm(elsi_h,d_out)

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if((elsi_h%n_elpa_steps > 0) .and. &
            (elsi_h%n_elsi_calls == elsi_h%n_elpa_steps+1)) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pztranc(elsi_h%n_basis,elsi_h%n_basis,(1.0_r8,0.0_r8),&
                    elsi_h%evec_complex,1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),&
                    d_out,1,1,elsi_h%sc_desc)

            elsi_h%coeff_omm%zval(1:elsi_h%coeff_omm%iaux2(1),1:elsi_h%coeff_omm%iaux2(2)) = &
               d_out(1:elsi_h%coeff_omm%iaux2(1),1:elsi_h%coeff_omm%iaux2(2))

            ! ELPA matrices are no longer needed
            if(associated(elsi_h%ham_complex)) then
               nullify(elsi_h%ham_complex)
            endif
            if(associated(elsi_h%ovlp_complex)) then
               nullify(elsi_h%ovlp_complex)
            endif
            if(associated(elsi_h%evec_complex)) then
               nullify(elsi_h%evec_complex)
            endif
            if(associated(elsi_h%eval)) then
               nullify(elsi_h%eval)
            endif
            if(associated(elsi_h%dm_complex)) then
               nullify(elsi_h%dm_complex)
            endif
            if(allocated(elsi_h%evec_complex_elpa)) then
               call elsi_deallocate(elsi_h,elsi_h%evec_complex_elpa,"evec_complex_elpa")
            endif
            if(allocated(elsi_h%eval_elpa)) then
               call elsi_deallocate(elsi_h,elsi_h%eval_elpa,"eval_elpa")
            endif
            if(allocated(elsi_h%eval_all)) then
               call elsi_deallocate(elsi_h,elsi_h%eval_all,"eval_all")
            endif
            if(allocated(elsi_h%occ_num)) then
               call elsi_deallocate(elsi_h,elsi_h%occ_num,"occ_num")
            endif
            if(allocated(elsi_h%k_weight)) then
               call elsi_deallocate(elsi_h,elsi_h%k_weight,"k_weight")
            endif
         endif

         ! Solve
         call elsi_solve_evp_omm(elsi_h)

         elsi_h%dm_omm%zval = 2.0_r8*elsi_h%dm_omm%zval
         call elsi_get_energy(elsi_h,energy_out)
      endif
   case(PEXSI)
      call elsi_print_pexsi_options(elsi_h)

      ! Initialize PEXSI
      call elsi_init_pexsi(elsi_h)

      ! Convert BLACS H and S to PEXSI format
      if((.not. elsi_h%ovlp_is_unit) .and. (elsi_h%n_elsi_calls == 1)) then
         call elsi_set_full_mat(elsi_h,s_in)
      endif
      call elsi_set_full_mat(elsi_h,h_in)
      call elsi_blacs_to_pexsi_hs(elsi_h,h_in,s_in)

      ! Allocate
      if(.not. allocated(elsi_h%dm_complex_pexsi)) then
         call elsi_allocate(elsi_h,elsi_h%dm_complex_pexsi,elsi_h%nnz_l_sp,&
                 "dm_complex_pexsi",caller)
      endif
      elsi_h%dm_complex_pexsi = (0.0_r8,0.0_r8)

      ! Set matrices
      call elsi_set_sparse_ham(elsi_h,elsi_h%ham_complex_pexsi)
      call elsi_set_sparse_ovlp(elsi_h,elsi_h%ovlp_complex_pexsi)
      call elsi_set_row_ind(elsi_h,elsi_h%row_ind_pexsi)
      call elsi_set_col_ptr(elsi_h,elsi_h%col_ptr_pexsi)
      call elsi_set_sparse_dm(elsi_h,elsi_h%dm_complex_pexsi)

      ! Solve
      call elsi_solve_evp_pexsi(elsi_h)

      ! Convert PEXSI density matrix to BLACS format
      call elsi_pexsi_to_blacs_dm(elsi_h,d_out)
      call elsi_get_energy(elsi_h,energy_out)

      elsi_h%mu_ready = .true.
   case(CHESS)
      call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen."//&
              " Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET
   elsi_h%edm_ready_complex = .true.

end subroutine

!>
!! This routine computes the density matrix.
!!
subroutine elsi_dm_real_sparse(elsi_h,h_in,s_in,d_out,energy_out)

   implicit none

   type(elsi_handle) :: elsi_h                 !< Handle
   real(kind=r8)     :: h_in(elsi_h%nnz_l_sp)  !< Hamiltonian
   real(kind=r8)     :: s_in(elsi_h%nnz_l_sp)  !< Overlap
   real(kind=r8)     :: d_out(elsi_h%nnz_l_sp) !< Density matrix
   real(kind=r8)     :: energy_out             !< Energy

   character*40, parameter :: caller = "elsi_dm_real_sparse"

   call elsi_check_handle(elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! REAL case
   elsi_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   select case(elsi_h%solver)
   case(ELPA)
      ! Convert PEXSI H and S to BLACS format
      call elsi_pexsi_to_blacs_hs(elsi_h,h_in,s_in)

      ! Allocate
      if(.not. allocated(elsi_h%eval_elpa)) then
         call elsi_allocate(elsi_h,elsi_h%eval_elpa,elsi_h%n_basis,&
                 "eval_elpa",caller)
      endif
      if(.not. allocated(elsi_h%evec_real_elpa)) then
         call elsi_allocate(elsi_h,elsi_h%evec_real_elpa,elsi_h%n_l_rows,&
                 elsi_h%n_l_cols,"evec_real_elpa",caller)
      endif
      if(.not. allocated(elsi_h%dm_real_elpa)) then
         call elsi_allocate(elsi_h,elsi_h%dm_real_elpa,elsi_h%n_l_rows,&
                 elsi_h%n_l_cols,"dm_real_elpa",caller)
      endif

      ! Set matrices
      call elsi_set_ham(elsi_h,elsi_h%ham_real_elpa)
      call elsi_set_ovlp(elsi_h,elsi_h%ovlp_real_elpa)
      call elsi_set_evec(elsi_h,elsi_h%evec_real_elpa)
      call elsi_set_eval(elsi_h,elsi_h%eval_elpa)
      call elsi_set_dm(elsi_h,elsi_h%dm_real_elpa)

      ! Solve eigenvalue problem
      call elsi_solve_evp_elpa(elsi_h)

      ! Compute density matrix
      call elsi_compute_occ_elpa(elsi_h)
      call elsi_compute_dm_elpa(elsi_h)
      call elsi_blacs_to_pexsi_dm(elsi_h,d_out)
      call elsi_get_energy(elsi_h,energy_out)

      elsi_h%mu_ready = .true.
   case(LIBOMM)
      call elsi_print_omm_options(elsi_h)

      ! Convert PEXSI H and S to BLACS format
      call elsi_pexsi_to_blacs_hs(elsi_h,h_in,s_in)

      if(elsi_h%n_elsi_calls <= elsi_h%n_elpa_steps) then
         if((elsi_h%n_elsi_calls == 1) .and. (elsi_h%omm_flavor == 0)) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(elsi_h,elsi_h%ovlp_real_omm,elsi_h%n_l_rows,&
                    elsi_h%n_l_cols,"ovlp_real_omm",caller)
            elsi_h%ovlp_real_omm = elsi_h%ovlp_real_elpa
         endif

         ! Compute libOMM initial guess by ELPA
         elsi_h%solver = ELPA

         ! Allocate
         if(.not. allocated(elsi_h%eval_elpa)) then
            call elsi_allocate(elsi_h,elsi_h%eval_elpa,elsi_h%n_basis,&
                    "eval_elpa",caller)
         endif
         if(.not. allocated(elsi_h%evec_real_elpa)) then
            call elsi_allocate(elsi_h,elsi_h%evec_real_elpa,elsi_h%n_l_rows,&
                    elsi_h%n_l_cols,"evec_real_elpa",caller)
         endif
         if(.not. allocated(elsi_h%dm_real_elpa)) then
            call elsi_allocate(elsi_h,elsi_h%dm_real_elpa,elsi_h%n_l_rows,&
                    elsi_h%n_l_cols,"dm_real_elpa",caller)
         endif

         ! Set matrices
         call elsi_set_ham(elsi_h,elsi_h%ham_real_elpa)
         call elsi_set_ovlp(elsi_h,elsi_h%ovlp_real_elpa)
         call elsi_set_evec(elsi_h,elsi_h%evec_real_elpa)
         call elsi_set_eval(elsi_h,elsi_h%eval_elpa)
         call elsi_set_dm(elsi_h,elsi_h%dm_real_elpa)

         ! Solve eigenvalue problem
         call elsi_solve_evp_elpa(elsi_h)

         ! Compute density matrix
         call elsi_compute_occ_elpa(elsi_h)
         call elsi_compute_dm_elpa(elsi_h)
         call elsi_blacs_to_pexsi_dm(elsi_h,d_out)
         call elsi_get_energy(elsi_h,energy_out)

         ! Switch back to libOMM here to guarantee elsi_customize_omm
         elsi_h%solver = LIBOMM
      else ! ELPA is done
         if(allocated(elsi_h%ovlp_real_omm)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            elsi_h%ovlp_real_elpa = elsi_h%ovlp_real_omm
            call elsi_deallocate(elsi_h,elsi_h%ovlp_real_omm,"ovlp_real_omm")

            call elsi_set_full_mat(elsi_h,elsi_h%ovlp_real_elpa)
         endif

         ! Allocate
         if(.not. elsi_h%coeff_omm%is_initialized) then
            call m_allocate(elsi_h%coeff_omm,elsi_h%n_states_omm,elsi_h%n_basis,"pddbc")
         endif

         ! Set matrices
         call elsi_set_ham(elsi_h,elsi_h%ham_real_elpa)
         call elsi_set_ovlp(elsi_h,elsi_h%ovlp_real_elpa)
         call elsi_set_dm(elsi_h,elsi_h%dm_real_elpa)

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if((elsi_h%n_elpa_steps > 0) .and. &
            (elsi_h%n_elsi_calls == elsi_h%n_elpa_steps+1)) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pdtran(elsi_h%n_basis,elsi_h%n_basis,1.0_r8,elsi_h%evec_real,1,1,&
                    elsi_h%sc_desc,0.0_r8,elsi_h%dm_real_elpa,1,1,elsi_h%sc_desc)

            elsi_h%coeff_omm%dval(1:elsi_h%coeff_omm%iaux2(1),1:elsi_h%coeff_omm%iaux2(2)) = &
               elsi_h%dm_real_elpa(1:elsi_h%coeff_omm%iaux2(1),1:elsi_h%coeff_omm%iaux2(2))

            ! ELPA matrices are no longer needed
            if(associated(elsi_h%ham_real)) then
               nullify(elsi_h%ham_real)
            endif
            if(associated(elsi_h%ovlp_real)) then
               nullify(elsi_h%ovlp_real)
            endif
            if(associated(elsi_h%evec_real)) then
               nullify(elsi_h%evec_real)
            endif
            if(associated(elsi_h%dm_real)) then
               nullify(elsi_h%dm_real)
            endif
            if(associated(elsi_h%eval)) then
               nullify(elsi_h%eval)
            endif
            if(allocated(elsi_h%evec_real_elpa)) then
               call elsi_deallocate(elsi_h,elsi_h%evec_real_elpa,"evec_real_elpa")
            endif
            if(allocated(elsi_h%eval_elpa)) then
               call elsi_deallocate(elsi_h,elsi_h%eval_elpa,"eval_elpa")
            endif
            if(allocated(elsi_h%occ_num)) then
               call elsi_deallocate(elsi_h,elsi_h%occ_num,"occ_num")
            endif
         endif

         ! Solve
         call elsi_solve_evp_omm(elsi_h)

         elsi_h%dm_omm%dval = 2.0_r8*elsi_h%dm_omm%dval
         call elsi_blacs_to_pexsi_dm(elsi_h,d_out)
         call elsi_get_energy(elsi_h,energy_out)
      endif
   case(PEXSI)
      call elsi_print_pexsi_options(elsi_h)

      ! Set matrices
      call elsi_set_sparse_ham(elsi_h,h_in)
      call elsi_set_sparse_ovlp(elsi_h,s_in)
      call elsi_set_sparse_dm(elsi_h,d_out)

      ! Initialize PEXSI
      call elsi_init_pexsi(elsi_h)

      ! Solve
      call elsi_solve_evp_pexsi(elsi_h)

      call elsi_get_energy(elsi_h,energy_out)

      elsi_h%mu_ready = .true.
   case(CHESS)
      call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
   case(SIPS)
      call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
   case default
      call elsi_stop(" No supported solver has been chosen."//&
              " Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET
   elsi_h%edm_ready_real = .true.

end subroutine

end module ELSI_SOLVER
