! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains subroutines to solve an eigenproblem or to compute the
!! density matrix, using one of the solvers ELPA, libOMM, PEXSI, SIPS, DMP.
!!
module ELSI_SOLVER

   use ELSI_CONSTANTS,  only: ELPA_SOLVER,OMM_SOLVER,PEXSI_SOLVER,SIPS_SOLVER,&
                              DMP_SOLVER,REAL_DATA,CMPLX_DATA,MULTI_PROC,&
                              SINGLE_PROC,PEXSI_CSC,SIESTA_CSC,UNSET,OUTPUT_EV,&
                              OUTPUT_DM,TIME_LEN
   use ELSI_DATATYPE,   only: elsi_handle
   use ELSI_DMP,        only: elsi_solve_evp_dmp_real
   use ELSI_ELPA,       only: elsi_compute_occ_elpa,elsi_compute_dm_elpa_real,&
                              elsi_solve_evp_elpa_real,&
                              elsi_compute_dm_elpa_cmplx,&
                              elsi_solve_evp_elpa_cmplx
   use ELSI_LAPACK,     only: elsi_solve_evp_lapack_real,&
                              elsi_solve_evp_lapack_cmplx
   use ELSI_MALLOC,     only: elsi_allocate,elsi_deallocate
   use ELSI_MAT_REDIST, only: elsi_blacs_to_pexsi_hs_cmplx,&
                              elsi_blacs_to_pexsi_hs_real,&
                              elsi_blacs_to_siesta_dm_cmplx,&
                              elsi_blacs_to_siesta_dm_real,&
                              elsi_blacs_to_sips_dm_cmplx,&
                              elsi_blacs_to_sips_dm_real,&
                              elsi_blacs_to_sips_hs_cmplx,&
                              elsi_blacs_to_sips_hs_real,&
                              elsi_pexsi_to_blacs_dm_cmplx,&
                              elsi_pexsi_to_blacs_dm_real,&
                              elsi_pexsi_to_siesta_dm_cmplx,&
                              elsi_pexsi_to_siesta_dm_real,&
                              elsi_siesta_to_blacs_hs_cmplx,&
                              elsi_siesta_to_blacs_hs_real,&
                              elsi_siesta_to_pexsi_hs_cmplx,&
                              elsi_siesta_to_pexsi_hs_real,&
                              elsi_siesta_to_sips_hs_real,&
                              elsi_sips_to_blacs_dm_real,&
                              elsi_sips_to_blacs_ev_real,&
                              elsi_sips_to_blacs_hs_cmplx,&
                              elsi_sips_to_blacs_hs_real,&
                              elsi_sips_to_siesta_dm_real
   use ELSI_MPI,        only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8
   use ELSI_OMM,        only: elsi_solve_evp_omm_real,elsi_solve_evp_omm_cmplx
   use ELSI_PEXSI,      only: elsi_init_pexsi,elsi_solve_evp_pexsi_real,&
                              elsi_solve_evp_pexsi_cmplx
   use ELSI_PRECISION,  only: r8,i4
   use ELSI_SETUP,      only: elsi_set_blacs
   use ELSI_SIPS,       only: elsi_init_sips,elsi_solve_evp_sips_real,&
                              elsi_compute_dm_sips_real
   use ELSI_TIMINGS,    only: elsi_get_time,elsi_process_timing,&
                              elsi_get_datetime_rfc3339
   use ELSI_UTILS,      only: elsi_check,elsi_check_handle,elsi_ready_handle
   use MATRIXSWITCH,    only: m_allocate

   implicit none

   private

   ! Solver
   public :: elsi_ev_real
   public :: elsi_ev_complex
   public :: elsi_ev_real_sparse
   public :: elsi_ev_complex_sparse
   public :: elsi_dm_real
   public :: elsi_dm_complex
   public :: elsi_dm_real_sparse
   public :: elsi_dm_complex_sparse

contains

!>
!! This routine gets the energy.
!!
subroutine elsi_get_energy(e_h,energy,solver)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   real(kind=r8),     intent(out)   :: energy !< Energy of the system
   integer(kind=i4),  intent(in)    :: solver !< Solver in use

   real(kind=r8)    :: tmp_real
   integer(kind=i4) :: i_state
   integer(kind=i4) :: ierr

   character(len=40), parameter :: caller = "elsi_get_energy"

   select case(solver)
   case(ELPA_SOLVER)
      energy = 0.0_r8

      do i_state = 1,e_h%n_states_solve
         energy = energy+e_h%i_weight*e_h%eval_elpa(i_state)*&
                     e_h%occ_num(i_state,e_h%i_spin,e_h%i_kpt)
      enddo
   case(OMM_SOLVER)
      energy = e_h%spin_degen*e_h%energy_hdm*e_h%i_weight
   case(PEXSI_SOLVER)
      energy = e_h%energy_hdm*e_h%i_weight
   case(SIPS_SOLVER)
      energy = 0.0_r8

      do i_state = 1,e_h%n_states_solve
         energy = energy+e_h%i_weight*e_h%eval_elpa(i_state)*&
                     e_h%occ_num(i_state,e_h%i_spin,e_h%i_kpt)
      enddo
   case(DMP_SOLVER)
      energy = e_h%spin_degen*e_h%energy_hdm*e_h%i_weight
   case default
      call elsi_stop(e_h,"Unsupported solver.",caller)
   end select

   if(e_h%n_spins*e_h%n_kpts > 1) then
      if(e_h%myid /= 0) then
         energy = 0.0_r8
      endif

      call MPI_Allreduce(energy,tmp_real,1,mpi_real8,mpi_sum,e_h%mpi_comm_all,&
              ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

      energy = tmp_real
   endif

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors. Note the
!! intent(inout) - it is because everything has the potential to be reused in
!! the next call.
!!
subroutine elsi_ev_real(e_h,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                         !< Handle
   real(kind=r8),     intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)  !< Hamiltonian
   real(kind=r8),     intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol) !< Overlap
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)           !< Eigenvalues
   real(kind=r8),     intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol) !< Eigenvectors

   real(kind=r8)           :: t0
   integer(kind=i4)        :: sol = UNSET
   character(len=TIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_ev_real"

   call elsi_check_handle(e_h,caller)
   call elsi_ready_handle(e_h,caller)
   call elsi_check(e_h,caller)
   call elsi_get_time(t0)
   call elsi_get_datetime_rfc3339(dt0)

   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   select case(e_h%solver)
   case(ELPA_SOLVER)
      if(e_h%parallel_mode == SINGLE_PROC) then
         call elsi_solve_evp_lapack_real(e_h,ham,ovlp,eval,evec)
      else
         call elsi_solve_evp_elpa_real(e_h,ham,ovlp,eval,evec)
      endif

      sol = ELPA_SOLVER
   case(SIPS_SOLVER)
      if(e_h%n_elsi_calls <= e_h%sips_n_elpa) then
         if(e_h%n_elsi_calls == 1) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(e_h,e_h%ovlp_real_copy,e_h%n_lrow,e_h%n_lcol,&
                    "ovlp_real_copy",caller)
            e_h%ovlp_real_copy = ovlp
         endif

         call elsi_solve_evp_elpa_real(e_h,ham,ovlp,eval,evec)

         sol = ELPA_SOLVER
      else ! ELPA is done
         if(allocated(e_h%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = e_h%ovlp_real_copy
            call elsi_deallocate(e_h,e_h%ovlp_real_copy,"ovlp_real_copy")
         endif

         call elsi_init_sips(e_h)

         call elsi_blacs_to_sips_hs_real(e_h,ham,ovlp)
         call elsi_solve_evp_sips_real(e_h,e_h%ham_real_pexsi,&
                 e_h%ovlp_real_pexsi,eval)
         call elsi_sips_to_blacs_ev_real(e_h,evec)

         sol = SIPS_SOLVER
      endif
   case default
      call elsi_stop(e_h,"Unsupported eigensolver.",caller)
   end select

   call elsi_process_timing(e_h,OUTPUT_EV,REAL_DATA,sol,dt0,t0)

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors. Note the
!! intent(inout) - it is because everything has the potential to be reused in
!! the next call.
!!
subroutine elsi_ev_complex(e_h,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                         !< Handle
   complex(kind=r8),  intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)  !< Hamiltonian
   complex(kind=r8),  intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol) !< Overlap
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)           !< Eigenvalues
   complex(kind=r8),  intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol) !< Eigenvectors

   real(kind=r8)           :: t0
   integer(kind=i4)        :: sol = UNSET
   character(len=TIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_ev_complex"

   call elsi_check_handle(e_h,caller)
   call elsi_ready_handle(e_h,caller)
   call elsi_check(e_h,caller)
   call elsi_get_time(t0)
   call elsi_get_datetime_rfc3339(dt0)

   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   select case(e_h%solver)
   case(ELPA_SOLVER)
      if(e_h%parallel_mode == SINGLE_PROC) then
         call elsi_solve_evp_lapack_cmplx(e_h,ham,ovlp,eval,evec)
      else
         call elsi_solve_evp_elpa_cmplx(e_h,ham,ovlp,eval,evec)
      endif

      sol = ELPA_SOLVER
   case default
      call elsi_stop(e_h,"Unsupported eigensolver.",caller)
   end select

   call elsi_process_timing(e_h,OUTPUT_EV,CMPLX_DATA,sol,dt0,t0)

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors. Note the
!! intent(inout) - it is because everything has the potential to be reused in
!! the next call.
!!
subroutine elsi_ev_real_sparse(e_h,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                         !< Handle
   real(kind=r8),     intent(inout) :: ham(e_h%nnz_l_sp)           !< Hamiltonian
   real(kind=r8),     intent(inout) :: ovlp(e_h%nnz_l_sp)          !< Overlap
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)           !< Eigenvalues
   real(kind=r8),     intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol) !< Eigenvectors

   real(kind=r8)           :: t0
   integer(kind=i4)        :: sol = UNSET
   character(len=TIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_ev_real_sparse"

   call elsi_check_handle(e_h,caller)
   call elsi_ready_handle(e_h,caller)
   call elsi_check(e_h,caller)
   call elsi_get_time(t0)
   call elsi_get_datetime_rfc3339(dt0)

   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   select case(e_h%solver)
   case(ELPA_SOLVER)
      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_real(e_h,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_real(e_h,ham,ovlp)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      call elsi_solve_evp_elpa_real(e_h,e_h%ham_real_elpa,e_h%ovlp_real_elpa,&
              eval,evec)

      sol = ELPA_SOLVER
   case(SIPS_SOLVER)
      if(e_h%n_elsi_calls <= e_h%sips_n_elpa) then
         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_sips_to_blacs_hs_real(e_h,ham,ovlp)
         case(SIESTA_CSC)
            call elsi_siesta_to_blacs_hs_real(e_h,ham,ovlp)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select

         call elsi_solve_evp_elpa_real(e_h,e_h%ham_real_elpa,&
                 e_h%ovlp_real_elpa,eval,evec)

         sol = ELPA_SOLVER
      else ! ELPA is done
         if(allocated(e_h%ham_real_elpa)) then
            call elsi_deallocate(e_h,e_h%ham_real_elpa,"ham_real_elpa")
         endif
         if(allocated(e_h%ovlp_real_elpa)) then
            call elsi_deallocate(e_h,e_h%ovlp_real_elpa,"ovlp_real_elpa")
         endif

         call elsi_init_sips(e_h)

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_solve_evp_sips_real(e_h,ham,ovlp,eval)
         case(SIESTA_CSC)
            call elsi_siesta_to_sips_hs_real(e_h,ham,ovlp)
            call elsi_solve_evp_sips_real(e_h,e_h%ham_real_pexsi,&
                    e_h%ovlp_real_pexsi,eval)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select

         call elsi_sips_to_blacs_ev_real(e_h,evec)

         sol = SIPS_SOLVER
      endif
   case default
      call elsi_stop(e_h,"Unsupported eigensolver.",caller)
   end select

   call elsi_process_timing(e_h,OUTPUT_EV,REAL_DATA,sol,dt0,t0)

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors. Note the
!! intent(inout) - it is because everything has the potential to be reused in
!! the next call.
!!
subroutine elsi_ev_complex_sparse(e_h,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                         !< Handle
   complex(kind=r8),  intent(inout) :: ham(e_h%nnz_l_sp)           !< Hamiltonian
   complex(kind=r8),  intent(inout) :: ovlp(e_h%nnz_l_sp)          !< Overlap
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)           !< Eigenvalues
   complex(kind=r8),  intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol) !< Eigenvectors

   real(kind=r8)           :: t0
   integer(kind=i4)        :: sol = UNSET
   character(len=TIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_ev_complex_sparse"

   call elsi_check_handle(e_h,caller)
   call elsi_ready_handle(e_h,caller)
   call elsi_check(e_h,caller)
   call elsi_get_time(t0)
   call elsi_get_datetime_rfc3339(dt0)

   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   select case(e_h%solver)
   case(ELPA_SOLVER)
      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_cmplx(e_h,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_cmplx(e_h,ham,ovlp)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      call elsi_solve_evp_elpa_cmplx(e_h,e_h%ham_cmplx_elpa,&
              e_h%ovlp_cmplx_elpa,eval,evec)

      sol = ELPA_SOLVER
   case default
      call elsi_stop(e_h,"Unsupported eigensolver.",caller)
   end select

   call elsi_process_timing(e_h,OUTPUT_EV,CMPLX_DATA,sol,dt0,t0)

end subroutine

!>
!! This routine computes the density matrix. Note the intent(inout) - it is
!! because everything has the potential to be reused in the next call.
!!
subroutine elsi_dm_real(e_h,ham,ovlp,dm,energy)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                         !< Handle
   real(kind=r8),     intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)  !< Hamiltonian
   real(kind=r8),     intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol) !< Overlap
   real(kind=r8),     intent(inout) :: dm(e_h%n_lrow,e_h%n_lcol)   !< Density matrix
   real(kind=r8),     intent(inout) :: energy                      !< Energy

   real(kind=r8)           :: t0
   integer(kind=i4)        :: sol = UNSET
   character(len=TIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_dm_real"

   call elsi_check_handle(e_h,caller)
   call elsi_ready_handle(e_h,caller)
   call elsi_check(e_h,caller)
   call elsi_get_time(t0)
   call elsi_get_datetime_rfc3339(dt0)

   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   select case(e_h%solver)
   case(ELPA_SOLVER)
      if(.not. allocated(e_h%eval_elpa)) then
         call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(e_h%evec_real_elpa)) then
         call elsi_allocate(e_h,e_h%evec_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                 "evec_real_elpa",caller)
      endif

      call elsi_solve_evp_elpa_real(e_h,ham,ovlp,e_h%eval_elpa,&
              e_h%evec_real_elpa)
      call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
      call elsi_compute_dm_elpa_real(e_h,e_h%evec_real_elpa,dm,ham)
      call elsi_get_energy(e_h,energy,ELPA_SOLVER)

      sol = ELPA_SOLVER

      e_h%mu_ready = .true.
      e_h%ts_ready = .true.
   case(OMM_SOLVER)
      if(e_h%n_elsi_calls <= e_h%omm_n_elpa) then
         if(e_h%n_elsi_calls == 1 .and. e_h%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(e_h,e_h%ovlp_real_copy,e_h%n_lrow,e_h%n_lcol,&
                    "ovlp_real_copy",caller)
            e_h%ovlp_real_copy = ovlp
         endif

         if(.not. allocated(e_h%eval_elpa)) then
            call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
         endif
         if(.not. allocated(e_h%evec_real_elpa)) then
            call elsi_allocate(e_h,e_h%evec_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                    "evec_real_elpa",caller)
         endif

         call elsi_solve_evp_elpa_real(e_h,ham,ovlp,e_h%eval_elpa,&
                 e_h%evec_real_elpa)
         call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
         call elsi_compute_dm_elpa_real(e_h,e_h%evec_real_elpa,dm,ham)
         call elsi_get_energy(e_h,energy,ELPA_SOLVER)

         sol = ELPA_SOLVER
      else ! ELPA is done
         if(allocated(e_h%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = e_h%ovlp_real_copy
            call elsi_deallocate(e_h,e_h%ovlp_real_copy,"ovlp_real_copy")
         endif

         if(.not. e_h%c_omm%is_initialized) then
            call m_allocate(e_h%c_omm,e_h%omm_n_states,e_h%n_basis,"pddbc")
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(e_h%omm_n_elpa > 0 .and. e_h%n_elsi_calls == e_h%omm_n_elpa+1) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,e_h%evec_real_elpa,1,1,&
                    e_h%sc_desc,0.0_r8,dm,1,1,e_h%sc_desc)

            e_h%c_omm%dval(1:e_h%c_omm%iaux2(1),1:e_h%c_omm%iaux2(2)) = &
               dm(1:e_h%c_omm%iaux2(1),1:e_h%c_omm%iaux2(2))

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

         call elsi_solve_evp_omm_real(e_h,ham,ovlp,dm)
         call elsi_get_energy(e_h,energy,OMM_SOLVER)

         sol = OMM_SOLVER
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(e_h)
      call elsi_blacs_to_pexsi_hs_real(e_h,ham,ovlp)

      if(.not. allocated(e_h%dm_real_pexsi)) then
         call elsi_allocate(e_h,e_h%dm_real_pexsi,e_h%nnz_l_sp,"dm_real_pexsi",&
                 caller)
      endif
      e_h%dm_real_pexsi = 0.0_r8

      call elsi_solve_evp_pexsi_real(e_h,e_h%ham_real_pexsi,&
              e_h%ovlp_real_pexsi,e_h%dm_real_pexsi)
      call elsi_pexsi_to_blacs_dm_real(e_h,dm)
      call elsi_get_energy(e_h,energy,PEXSI_SOLVER)

      sol = PEXSI_SOLVER

      e_h%mu_ready = .true.
   case(SIPS_SOLVER)
      if(.not. allocated(e_h%eval_elpa)) then
         call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(e_h%evec_real_elpa)) then
         call elsi_allocate(e_h,e_h%evec_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                 "evec_real_elpa",caller)
      endif

      if(e_h%n_elsi_calls <= e_h%sips_n_elpa) then
         if(e_h%n_elsi_calls == 1) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(e_h,e_h%ovlp_real_copy,e_h%n_lrow,e_h%n_lcol,&
                    "ovlp_real_copy",caller)
            e_h%ovlp_real_copy = ovlp
         endif

         call elsi_solve_evp_elpa_real(e_h,ham,ovlp,e_h%eval_elpa,&
                 e_h%evec_real_elpa)

         call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
         call elsi_compute_dm_elpa_real(e_h,e_h%evec_real_elpa,dm,ham)
         call elsi_get_energy(e_h,energy,ELPA_SOLVER)

         sol = ELPA_SOLVER
      else ! ELPA is done
         if(allocated(e_h%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = e_h%ovlp_real_copy
            call elsi_deallocate(e_h,e_h%ovlp_real_copy,"ovlp_real_copy")
            call elsi_deallocate(e_h,e_h%evec_real_elpa,"evec_real_elpa")
         endif

         call elsi_init_sips(e_h)
         call elsi_blacs_to_sips_hs_real(e_h,ham,ovlp)

         if(.not. allocated(e_h%dm_real_pexsi)) then
            call elsi_allocate(e_h,e_h%dm_real_pexsi,e_h%nnz_l_sp,&
                    "dm_real_pexsi",caller)
         endif
         e_h%dm_real_pexsi = 0.0_r8

         call elsi_solve_evp_sips_real(e_h,e_h%ham_real_pexsi,&
                 e_h%ovlp_real_pexsi,e_h%eval_elpa)
         call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
         call elsi_compute_dm_sips_real(e_h,e_h%dm_real_pexsi)
         call elsi_sips_to_blacs_dm_real(e_h,dm)
         call elsi_get_energy(e_h,energy,SIPS_SOLVER)

         sol = SIPS_SOLVER
      endif

      e_h%mu_ready = .true.
      e_h%ts_ready = .true.
   case(DMP_SOLVER)
      ! Save Hamiltonian and overlap
      if(e_h%n_elsi_calls==1) then
         call elsi_allocate(e_h,e_h%ham_real_copy,e_h%n_lrow,e_h%n_lcol,&
                 "ham_real_copy",caller)
         call elsi_allocate(e_h,e_h%ovlp_real_copy,e_h%n_lrow,e_h%n_lcol,&
                 "ovlp_real_copy",caller)
         e_h%ovlp_real_copy = ovlp
      endif
      e_h%ham_real_copy = ham

      call elsi_solve_evp_dmp_real(e_h,ham,ovlp,dm)
      call elsi_get_energy(e_h,energy,DMP_SOLVER)

      sol = DMP_SOLVER
   case default
      call elsi_stop(e_h,"Unsupported density matrix solver.",caller)
   end select

   call elsi_process_timing(e_h,OUTPUT_DM,REAL_DATA,sol,dt0,t0)

   e_h%edm_ready_real = .true.

end subroutine

!>
!! This routine computes the density matrix. Note the intent(inout) - it is
!! because everything has the potential to be reused in the next call.
!!
subroutine elsi_dm_complex(e_h,ham,ovlp,dm,energy)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                         !< Handle
   complex(kind=r8),  intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)  !< Hamiltonian
   complex(kind=r8),  intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol) !< Overlap
   complex(kind=r8),  intent(inout) :: dm(e_h%n_lrow,e_h%n_lcol)   !< Density matrix
   real(kind=r8),     intent(inout) :: energy                      !< Energy

   real(kind=r8)           :: t0
   integer(kind=i4)        :: sol = UNSET
   character(len=TIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_dm_complex"

   call elsi_check_handle(e_h,caller)
   call elsi_ready_handle(e_h,caller)
   call elsi_check(e_h,caller)
   call elsi_get_time(t0)
   call elsi_get_datetime_rfc3339(dt0)

   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   select case(e_h%solver)
   case(ELPA_SOLVER)
      if(.not. allocated(e_h%eval_elpa)) then
         call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(e_h%evec_cmplx_elpa)) then
         call elsi_allocate(e_h,e_h%evec_cmplx_elpa,e_h%n_lrow,e_h%n_lcol,&
                 "evec_cmplx_elpa",caller)
      endif

      call elsi_solve_evp_elpa_cmplx(e_h,ham,ovlp,e_h%eval_elpa,&
              e_h%evec_cmplx_elpa)
      call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
      call elsi_compute_dm_elpa_cmplx(e_h,e_h%evec_cmplx_elpa,dm,ham)
      call elsi_get_energy(e_h,energy,ELPA_SOLVER)

      sol = ELPA_SOLVER

      e_h%mu_ready = .true.
      e_h%ts_ready = .true.
   case(OMM_SOLVER)
      if(e_h%n_elsi_calls <= e_h%omm_n_elpa) then
         if(e_h%n_elsi_calls == 1 .and. e_h%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(e_h,e_h%ovlp_cmplx_copy,e_h%n_lrow,e_h%n_lcol,&
                    "ovlp_cmplx_copy",caller)
            e_h%ovlp_cmplx_copy = ovlp
         endif

         if(.not. allocated(e_h%eval_elpa)) then
            call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
         endif
         if(.not. allocated(e_h%evec_cmplx_elpa)) then
            call elsi_allocate(e_h,e_h%evec_cmplx_elpa,e_h%n_lrow,e_h%n_lcol,&
                    "evec_cmplx_elpa",caller)
         endif

         call elsi_solve_evp_elpa_cmplx(e_h,ham,ovlp,e_h%eval_elpa,&
                 e_h%evec_cmplx_elpa)
         call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
         call elsi_compute_dm_elpa_cmplx(e_h,e_h%evec_cmplx_elpa,dm,ham)
         call elsi_get_energy(e_h,energy,ELPA_SOLVER)

         sol = ELPA_SOLVER
      else ! ELPA is done
         if(allocated(e_h%ovlp_cmplx_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = e_h%ovlp_cmplx_copy
            call elsi_deallocate(e_h,e_h%ovlp_cmplx_copy,"ovlp_cmplx_copy")
         endif

         if(.not. e_h%c_omm%is_initialized) then
            call m_allocate(e_h%c_omm,e_h%omm_n_states,e_h%n_basis,"pzdbc")
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(e_h%omm_n_elpa > 0 .and. e_h%n_elsi_calls == e_h%omm_n_elpa+1) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),&
                    e_h%evec_cmplx_elpa,1,1,e_h%sc_desc,(0.0_r8,0.0_r8),dm,1,1,&
                    e_h%sc_desc)

            e_h%c_omm%zval(1:e_h%c_omm%iaux2(1),1:e_h%c_omm%iaux2(2)) = &
               dm(1:e_h%c_omm%iaux2(1),1:e_h%c_omm%iaux2(2))

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

         call elsi_solve_evp_omm_cmplx(e_h,ham,ovlp,dm)
         call elsi_get_energy(e_h,energy,OMM_SOLVER)

         sol = OMM_SOLVER
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(e_h)
      call elsi_blacs_to_pexsi_hs_cmplx(e_h,ham,ovlp)

      if(.not. allocated(e_h%dm_cmplx_pexsi)) then
         call elsi_allocate(e_h,e_h%dm_cmplx_pexsi,e_h%nnz_l_sp,&
                 "dm_cmplx_pexsi",caller)
      endif
      e_h%dm_cmplx_pexsi = (0.0_r8,0.0_r8)

      call elsi_solve_evp_pexsi_cmplx(e_h,e_h%ham_cmplx_pexsi,&
              e_h%ovlp_cmplx_pexsi,e_h%dm_cmplx_pexsi)
      call elsi_pexsi_to_blacs_dm_cmplx(e_h,dm)
      call elsi_get_energy(e_h,energy,PEXSI_SOLVER)

      sol = PEXSI_SOLVER

      e_h%mu_ready = .true.
   case default
      call elsi_stop(e_h,"Unsupported density matrix solver.",caller)
   end select

   call elsi_process_timing(e_h,OUTPUT_DM,CMPLX_DATA,sol,dt0,t0)

   e_h%edm_ready_cmplx = .true.

end subroutine

!>
!! This routine computes the density matrix. Note the intent(inout) - it is
!! because everything has the potential to be reused in the next call.
!!
subroutine elsi_dm_real_sparse(e_h,ham,ovlp,dm,energy)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                !< Handle
   real(kind=r8),     intent(inout) :: ham(e_h%nnz_l_sp)  !< Hamiltonian
   real(kind=r8),     intent(inout) :: ovlp(e_h%nnz_l_sp) !< Overlap
   real(kind=r8),     intent(inout) :: dm(e_h%nnz_l_sp)   !< Density matrix
   real(kind=r8),     intent(inout) :: energy             !< Energy

   real(kind=r8)           :: t0
   integer(kind=i4)        :: sol = UNSET
   character(len=TIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_dm_real_sparse"

   call elsi_check_handle(e_h,caller)
   call elsi_ready_handle(e_h,caller)
   call elsi_check(e_h,caller)
   call elsi_get_time(t0)
   call elsi_get_datetime_rfc3339(dt0)

   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   select case(e_h%solver)
   case(ELPA_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. e_h%blacs_ready) then
         call elsi_init_blacs(e_h)
      endif

      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_real(e_h,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_real(e_h,ham,ovlp)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      if(.not. allocated(e_h%eval_elpa)) then
         call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(e_h%evec_real_elpa)) then
         call elsi_allocate(e_h,e_h%evec_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                 "evec_real_elpa",caller)
      endif
      if(.not. allocated(e_h%dm_real_elpa)) then
         call elsi_allocate(e_h,e_h%dm_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                 "dm_real_elpa",caller)
      endif

      call elsi_solve_evp_elpa_real(e_h,e_h%ham_real_elpa,e_h%ovlp_real_elpa,&
              e_h%eval_elpa,e_h%evec_real_elpa)
      call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
      call elsi_compute_dm_elpa_real(e_h,e_h%evec_real_elpa,e_h%dm_real_elpa,&
              e_h%ham_real_elpa)

      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm_real(e_h,dm)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm_real(e_h,dm)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(e_h,energy,ELPA_SOLVER)

      sol = ELPA_SOLVER

      e_h%mu_ready = .true.
      e_h%ts_ready = .true.
   case(OMM_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. e_h%blacs_ready) then
         call elsi_init_blacs(e_h)
      endif

      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_real(e_h,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_real(e_h,ham,ovlp)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      if(e_h%n_elsi_calls <= e_h%omm_n_elpa) then
         if(e_h%n_elsi_calls == 1 .and. e_h%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(e_h,e_h%ovlp_real_copy,e_h%n_lrow,e_h%n_lcol,&
                    "ovlp_real_copy",caller)
            e_h%ovlp_real_copy = e_h%ovlp_real_elpa
         endif

         if(.not. allocated(e_h%eval_elpa)) then
            call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
         endif
         if(.not. allocated(e_h%evec_real_elpa)) then
            call elsi_allocate(e_h,e_h%evec_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                    "evec_real_elpa",caller)
         endif
         if(.not. allocated(e_h%dm_real_elpa)) then
            call elsi_allocate(e_h,e_h%dm_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                    "dm_real_elpa",caller)
         endif

         call elsi_solve_evp_elpa_real(e_h,e_h%ham_real_elpa,&
                 e_h%ovlp_real_elpa,e_h%eval_elpa,e_h%evec_real_elpa)
         call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
         call elsi_compute_dm_elpa_real(e_h,e_h%evec_real_elpa,&
                 e_h%dm_real_elpa,e_h%ham_real_elpa)

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_real(e_h,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_real(e_h,dm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(e_h,energy,ELPA_SOLVER)

         sol = ELPA_SOLVER
      else ! ELPA is done
         if(allocated(e_h%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            e_h%ovlp_real_elpa = e_h%ovlp_real_copy
            call elsi_deallocate(e_h,e_h%ovlp_real_copy,"ovlp_real_copy")
         endif

         if(.not. e_h%c_omm%is_initialized) then
            call m_allocate(e_h%c_omm,e_h%omm_n_states,e_h%n_basis,"pddbc")
         endif
         if(.not. allocated(e_h%dm_real_elpa)) then
            call elsi_allocate(e_h,e_h%dm_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                    "dm_real_elpa",caller)
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(e_h%omm_n_elpa > 0 .and. e_h%n_elsi_calls == e_h%omm_n_elpa+1) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,e_h%evec_real_elpa,1,1,&
                    e_h%sc_desc,0.0_r8,e_h%dm_real_elpa,1,1,e_h%sc_desc)

            e_h%c_omm%dval(1:e_h%c_omm%iaux2(1),1:e_h%c_omm%iaux2(2)) = &
               e_h%dm_real_elpa(1:e_h%c_omm%iaux2(1),1:e_h%c_omm%iaux2(2))

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

         call elsi_solve_evp_omm_real(e_h,e_h%ham_real_elpa,e_h%ovlp_real_elpa,&
                 e_h%dm_real_elpa)

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_real(e_h,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_real(e_h,dm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(e_h,energy,OMM_SOLVER)

         sol = OMM_SOLVER
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(e_h)

      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_solve_evp_pexsi_real(e_h,ham,ovlp,dm)
      case(SIESTA_CSC)
         call elsi_siesta_to_pexsi_hs_real(e_h,ham,ovlp)

         if(.not. allocated(e_h%dm_real_pexsi)) then
            call elsi_allocate(e_h,e_h%dm_real_pexsi,e_h%nnz_l_sp1,&
                    "dm_real_pexsi",caller)
         endif
         e_h%dm_real_pexsi = 0.0_r8

         call elsi_solve_evp_pexsi_real(e_h,e_h%ham_real_pexsi,&
                 e_h%ovlp_real_pexsi,e_h%dm_real_pexsi)
         call elsi_pexsi_to_siesta_dm_real(e_h,dm)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(e_h,energy,PEXSI_SOLVER)

      sol = PEXSI_SOLVER

      e_h%mu_ready = .true.
   case(SIPS_SOLVER)
      if(e_h%n_elsi_calls <= e_h%sips_n_elpa) then
         ! Set up BLACS if not done by user
         if(.not. e_h%blacs_ready) then
            call elsi_init_blacs(e_h)
         endif

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_sips_to_blacs_hs_real(e_h,ham,ovlp)
         case(SIESTA_CSC)
            call elsi_siesta_to_blacs_hs_real(e_h,ham,ovlp)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select

         if(.not. allocated(e_h%eval_elpa)) then
            call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
         endif
         if(.not. allocated(e_h%evec_real_elpa)) then
            call elsi_allocate(e_h,e_h%evec_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                    "evec_real_elpa",caller)
         endif
         if(.not. allocated(e_h%dm_real_elpa)) then
            call elsi_allocate(e_h,e_h%dm_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                    "dm_real_elpa",caller)
         endif

         call elsi_solve_evp_elpa_real(e_h,e_h%ham_real_elpa,e_h%ovlp_real_elpa,&
                 e_h%eval_elpa,e_h%evec_real_elpa)
         call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
         call elsi_compute_dm_elpa_real(e_h,e_h%evec_real_elpa,e_h%dm_real_elpa,&
                 e_h%ham_real_elpa)

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_real(e_h,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_real(e_h,dm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(e_h,energy,ELPA_SOLVER)

         sol = ELPA_SOLVER
      else ! ELPA is done
         if(allocated(e_h%ham_real_elpa)) then
            call elsi_deallocate(e_h,e_h%ham_real_elpa,"ham_real_elpa")
         endif
         if(allocated(e_h%ovlp_real_elpa)) then
            call elsi_deallocate(e_h,e_h%ovlp_real_elpa,"ovlp_real_elpa")
         endif
         if(allocated(e_h%dm_real_elpa)) then
            call elsi_deallocate(e_h,e_h%dm_real_elpa,"dm_real_elpa")
         endif
         if(.not. allocated(e_h%eval_elpa)) then
            call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
         endif

         call elsi_init_sips(e_h)

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_solve_evp_sips_real(e_h,ham,ovlp,e_h%eval_elpa)
            call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
            call elsi_compute_dm_sips_real(e_h,dm)
            call elsi_get_energy(e_h,energy,SIPS_SOLVER)
         case(SIESTA_CSC)
            call elsi_siesta_to_sips_hs_real(e_h,ham,ovlp)

            if(.not. allocated(e_h%dm_real_pexsi)) then
               call elsi_allocate(e_h,e_h%dm_real_pexsi,e_h%nnz_l_sp1,&
                       "dm_real_pexsi",caller)
            endif
            e_h%dm_real_pexsi = 0.0_r8

            call elsi_solve_evp_sips_real(e_h,e_h%ham_real_pexsi,&
                    e_h%ovlp_real_pexsi,e_h%eval_elpa)
            call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
            call elsi_compute_dm_sips_real(e_h,e_h%dm_real_pexsi)
            call elsi_sips_to_siesta_dm_real(e_h,dm)
            call elsi_get_energy(e_h,energy,SIPS_SOLVER)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select

         sol = SIPS_SOLVER
      endif

      e_h%mu_ready = .true.
      e_h%ts_ready = .true.
   case(DMP_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. e_h%blacs_ready) then
         call elsi_init_blacs(e_h)
      endif

      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_real(e_h,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_real(e_h,ham,ovlp)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      if(.not. allocated(e_h%dm_real_elpa)) then
         call elsi_allocate(e_h,e_h%dm_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                 "dm_real_elpa",caller)
      endif

      ! Save Hamiltonian and overlap
      if(e_h%n_elsi_calls==1) then
         call elsi_allocate(e_h,e_h%ham_real_copy,e_h%n_lrow,e_h%n_lcol,&
                 "ham_real_copy",caller)
         call elsi_allocate(e_h,e_h%ovlp_real_copy,e_h%n_lrow,e_h%n_lcol,&
                 "ovlp_real_copy",caller)
         e_h%ham_real_copy  = e_h%ham_real_elpa
         e_h%ovlp_real_copy = e_h%ovlp_real_elpa
      endif

      call elsi_solve_evp_dmp_real(e_h,e_h%ham_real_elpa,e_h%ovlp_real_elpa,&
              e_h%dm_real_elpa)

      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm_real(e_h,dm)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm_real(e_h,dm)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(e_h,energy,DMP_SOLVER)

      sol = DMP_SOLVER
   case default
      call elsi_stop(e_h,"Unsupported density matrix solver.",caller)
   end select

   call elsi_process_timing(e_h,OUTPUT_DM,REAL_DATA,sol,dt0,t0)

   e_h%edm_ready_real = .true.

end subroutine

!>
!! This routine computes the density matrix. Note the intent(inout) - it is
!! because everything has the potential to be reused in the next call.
!!
subroutine elsi_dm_complex_sparse(e_h,ham,ovlp,dm,energy)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                !< Handle
   complex(kind=r8),  intent(inout) :: ham(e_h%nnz_l_sp)  !< Hamiltonian
   complex(kind=r8),  intent(inout) :: ovlp(e_h%nnz_l_sp) !< Overlap
   complex(kind=r8),  intent(inout) :: dm(e_h%nnz_l_sp)   !< Density matrix
   real(kind=r8),     intent(inout) :: energy             !< Energy

   real(kind=r8)           :: t0
   integer(kind=i4)        :: sol = UNSET
   character(len=TIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_dm_complex_sparse"

   call elsi_check_handle(e_h,caller)
   call elsi_ready_handle(e_h,caller)
   call elsi_check(e_h,caller)
   call elsi_get_time(t0)
   call elsi_get_datetime_rfc3339(dt0)

   e_h%n_elsi_calls = e_h%n_elsi_calls+1

   select case(e_h%solver)
   case(ELPA_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. e_h%blacs_ready) then
         call elsi_init_blacs(e_h)
      endif

      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_cmplx(e_h,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_cmplx(e_h,ham,ovlp)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      if(.not. allocated(e_h%eval_elpa)) then
         call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(e_h%evec_cmplx_elpa)) then
         call elsi_allocate(e_h,e_h%evec_cmplx_elpa,e_h%n_lrow,e_h%n_lcol,&
                 "evec_cmplx_elpa",caller)
      endif
      if(.not. allocated(e_h%dm_cmplx_elpa)) then
         call elsi_allocate(e_h,e_h%dm_cmplx_elpa,e_h%n_lrow,e_h%n_lcol,&
                 "dm_cmplx_elpa",caller)
      endif

      call elsi_solve_evp_elpa_cmplx(e_h,e_h%ham_cmplx_elpa,&
              e_h%ovlp_cmplx_elpa,e_h%eval_elpa,e_h%evec_cmplx_elpa)
      call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
      call elsi_compute_dm_elpa_cmplx(e_h,e_h%evec_cmplx_elpa,&
              e_h%dm_cmplx_elpa,e_h%ham_cmplx_elpa)

      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm_cmplx(e_h,dm)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm_cmplx(e_h,dm)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(e_h,energy,ELPA_SOLVER)

      sol = ELPA_SOLVER

      e_h%mu_ready = .true.
      e_h%ts_ready = .true.
   case(OMM_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. e_h%blacs_ready) then
         call elsi_init_blacs(e_h)
      endif

      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_cmplx(e_h,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_cmplx(e_h,ham,ovlp)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      if(e_h%n_elsi_calls <= e_h%omm_n_elpa) then
         if(e_h%n_elsi_calls == 1 .and. e_h%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(e_h,e_h%ovlp_cmplx_copy,e_h%n_lrow,e_h%n_lcol,&
                    "ovlp_cmplx_copy",caller)
            e_h%ovlp_cmplx_copy = e_h%ovlp_cmplx_elpa
         endif

         if(.not. allocated(e_h%eval_elpa)) then
            call elsi_allocate(e_h,e_h%eval_elpa,e_h%n_basis,"eval_elpa",caller)
         endif
         if(.not. allocated(e_h%evec_cmplx_elpa)) then
            call elsi_allocate(e_h,e_h%evec_cmplx_elpa,e_h%n_lrow,e_h%n_lcol,&
                    "evec_cmplx_elpa",caller)
         endif
         if(.not. allocated(e_h%dm_cmplx_elpa)) then
            call elsi_allocate(e_h,e_h%dm_cmplx_elpa,e_h%n_lrow,e_h%n_lcol,&
                    "dm_cmplx_elpa",caller)
         endif

         call elsi_solve_evp_elpa_cmplx(e_h,e_h%ham_cmplx_elpa,&
                 e_h%ovlp_cmplx_elpa,e_h%eval_elpa,e_h%evec_cmplx_elpa)
         call elsi_compute_occ_elpa(e_h,e_h%eval_elpa)
         call elsi_compute_dm_elpa_cmplx(e_h,e_h%evec_cmplx_elpa,&
                 e_h%dm_cmplx_elpa,e_h%ham_cmplx_elpa)

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_cmplx(e_h,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_cmplx(e_h,dm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(e_h,energy,ELPA_SOLVER)

         sol = ELPA_SOLVER
      else ! ELPA is done
         if(allocated(e_h%ovlp_cmplx_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            e_h%ovlp_cmplx_elpa = e_h%ovlp_cmplx_copy
            call elsi_deallocate(e_h,e_h%ovlp_cmplx_copy,"ovlp_cmplx_copy")
         endif

         if(.not. e_h%c_omm%is_initialized) then
            call m_allocate(e_h%c_omm,e_h%omm_n_states,e_h%n_basis,"pzdbc")
         endif
         if(.not. allocated(e_h%dm_cmplx_elpa)) then
            call elsi_allocate(e_h,e_h%dm_cmplx_elpa,e_h%n_lrow,e_h%n_lcol,&
                    "dm_cmplx_elpa",caller)
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(e_h%omm_n_elpa > 0 .and. e_h%n_elsi_calls == e_h%omm_n_elpa+1) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),&
                    e_h%evec_cmplx_elpa,1,1,e_h%sc_desc,(0.0_r8,0.0_r8),&
                    e_h%dm_cmplx_elpa,1,1,e_h%sc_desc)

            e_h%c_omm%zval(1:e_h%c_omm%iaux2(1),1:e_h%c_omm%iaux2(2)) = &
               e_h%dm_cmplx_elpa(1:e_h%c_omm%iaux2(1),1:e_h%c_omm%iaux2(2))

            if(allocated(e_h%evec_cmplx_elpa)) then
               call elsi_deallocate(e_h,e_h%evec_cmplx_elpa,"evec_cmplx_elpa")
            endif
            if(allocated(e_h%eval_elpa)) then
               call elsi_deallocate(e_h,e_h%eval_elpa,"eval_elpa")
            endif
            if(allocated(e_h%occ_num)) then
               call elsi_deallocate(e_h,e_h%occ_num,"occ_num")
            endif
         endif

         call elsi_solve_evp_omm_cmplx(e_h,e_h%ham_cmplx_elpa,&
                 e_h%ovlp_cmplx_elpa,e_h%dm_cmplx_elpa)

         select case(e_h%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_cmplx(e_h,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_cmplx(e_h,dm)
         case default
            call elsi_stop(e_h,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(e_h,energy,OMM_SOLVER)

         sol = OMM_SOLVER
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(e_h)

      select case(e_h%matrix_format)
      case(PEXSI_CSC)
         call elsi_solve_evp_pexsi_cmplx(e_h,ham,ovlp,dm)
      case(SIESTA_CSC)
         call elsi_siesta_to_pexsi_hs_cmplx(e_h,ham,ovlp)

         if(.not. allocated(e_h%dm_cmplx_pexsi)) then
            call elsi_allocate(e_h,e_h%dm_cmplx_pexsi,e_h%nnz_l_sp1,&
                    "dm_cmplx_pexsi",caller)
         endif
         e_h%dm_cmplx_pexsi = (0.0_r8,0.0_r8)

         call elsi_solve_evp_pexsi_cmplx(e_h,e_h%ham_cmplx_pexsi,&
                 e_h%ovlp_cmplx_pexsi,e_h%dm_cmplx_pexsi)
         call elsi_pexsi_to_siesta_dm_cmplx(e_h,dm)
      case default
         call elsi_stop(e_h,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(e_h,energy,PEXSI_SOLVER)

      sol = PEXSI_SOLVER

      e_h%mu_ready = .true.
   case default
      call elsi_stop(e_h,"Unsupported density matrix solver.",caller)
   end select

   call elsi_process_timing(e_h,OUTPUT_DM,CMPLX_DATA,sol,dt0,t0)

   e_h%edm_ready_cmplx = .true.

end subroutine

!>
!! This routine initializes BLACS, in case that the user selects a sparse format
!! and BLACS is still used internally.
!!
subroutine elsi_init_blacs(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: block_size

   character(len=40), parameter :: caller = "elsi_init_blacs"

   if(e_h%parallel_mode == MULTI_PROC .and. .not. e_h%blacs_ready) then
      ! Set square-like process grid
      do nprow = nint(sqrt(real(e_h%n_procs))),2,-1
         if(mod(e_h%n_procs,nprow) == 0) exit
      enddo

      npcol = e_h%n_procs/nprow

      if(max(nprow,npcol) > e_h%n_basis) then
         call elsi_stop(e_h,"Matrix size is too small for this number of MPI"//&
                 "tasks.",caller)
      endif

      ! Initialize BLACS
      blacs_ctxt = e_h%mpi_comm

      call BLACS_Gridinit(blacs_ctxt,'r',nprow,npcol)

      ! Find block size
      block_size = 1

      do while(2*block_size*max(nprow,npcol) <= e_h%n_basis)
         block_size = 2*block_size
      enddo

      ! Maximum allowed value: 256
      block_size = min(256,block_size)

      ! ELPA works better with a small block_size
      if(e_h%solver == ELPA_SOLVER) then
         block_size = min(32,block_size)
      endif

      call elsi_set_blacs(e_h,blacs_ctxt,block_size)
   endif

end subroutine

end module ELSI_SOLVER
