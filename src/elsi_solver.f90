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

   use ELSI_CONSTANTS, only: ELPA_SOLVER,OMM_SOLVER,PEXSI_SOLVER,SIPS_SOLVER,&
                             DMP_SOLVER,MULTI_PROC,SINGLE_PROC,PEXSI_CSC,&
                             SIESTA_CSC,DATETIME_LEN
   use ELSI_DATATYPE,  only: elsi_handle,elsi_param_t,elsi_basic_t
   use ELSI_DMP,       only: elsi_init_dmp,elsi_solve_dmp_real
   use ELSI_ELPA,      only: elsi_init_elpa,elsi_compute_occ_elpa,&
                             elsi_compute_dm_elpa_real,elsi_solve_elpa_real,&
                             elsi_compute_dm_elpa_cmplx,elsi_solve_elpa_cmplx
   use ELSI_IO,        only: elsi_add_log,elsi_get_time,&
                             fjson_get_datetime_rfc3339
   use ELSI_LAPACK,    only: elsi_solve_lapack_real,elsi_solve_lapack_cmplx
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_REDIST,    only: elsi_blacs_to_pexsi_hs_cmplx,&
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
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8
   use ELSI_OMM,       only: elsi_init_omm,elsi_solve_omm_real,&
                             elsi_solve_omm_cmplx,elsi_init_coeff_omm_real,&
                             elsi_init_coeff_omm_cmplx
   use ELSI_PEXSI,     only: elsi_init_pexsi,elsi_solve_pexsi_real,&
                             elsi_solve_pexsi_cmplx
   use ELSI_PRECISION, only: r8,i4
   use ELSI_SETUP,     only: elsi_set_blacs
   use ELSI_SIPS,      only: elsi_init_sips,elsi_solve_sips_real,&
                             elsi_compute_dm_sips_real
   use ELSI_UTILS,     only: elsi_check,elsi_check_init

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
subroutine elsi_get_energy(ph,bh,energy,solver)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   integer(kind=i4),   intent(in)    :: solver
   real(kind=r8),      intent(out)   :: energy

   real(kind=r8)    :: tmp_real
   integer(kind=i4) :: i_state
   integer(kind=i4) :: ierr

   character(len=40), parameter :: caller = "elsi_get_energy"

   energy = ph%ebs*ph%i_weight

   if(solver == OMM_SOLVER .or. solver == DMP_SOLVER) then
      energy = ph%spin_degen*energy
   endif

   if(ph%n_spins*ph%n_kpts > 1) then
      if(bh%myid /= 0) then
         energy = 0.0_r8
      endif

      call MPI_Allreduce(energy,tmp_real,1,mpi_real8,mpi_sum,bh%comm_all,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      energy = tmp_real
   endif

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors. Note the
!! intent(inout) - it is because everything has the potential to be reused in
!! the next call.
!!
subroutine elsi_ev_real(eh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh                              !< Handle
   real(kind=r8),     intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol)  !< Hamiltonian
   real(kind=r8),     intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   real(kind=r8),     intent(inout) :: eval(eh%ph%n_basis)             !< Eigenvalues
   real(kind=r8),     intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   real(kind=r8)               :: t0
   character(len=DATETIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_ev_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      if(eh%ph%parallel_mode == SINGLE_PROC) then
         call elsi_solve_lapack_real(eh%ph,eh%bh,ham,ovlp,eval,evec)
      else
         call elsi_init_elpa(eh%ph,eh%bh)
         call elsi_solve_elpa_real(eh,eh%ph,eh%bh,ham,ovlp,eval,evec)
      endif
   case(SIPS_SOLVER)
      if(eh%ph%n_calls <= eh%ph%sips_n_elpa) then
         if(eh%ph%n_calls == 1) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_real_copy",caller)
            eh%ovlp_real_copy = ovlp
         endif

         call elsi_init_elpa(eh%ph,eh%bh)
         call elsi_solve_elpa_real(eh,eh%ph,eh%bh,ham,ovlp,eval,evec)
      else ! ELPA is done
         if(allocated(eh%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = eh%ovlp_real_copy
            call elsi_deallocate(eh%bh,eh%ovlp_real_copy,"ovlp_real_copy")
         endif

         call elsi_init_sips(eh%ph,eh%bh)
         call elsi_blacs_to_sips_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
         call elsi_solve_sips_real(eh,eh%ph,eh%bh,eh%ham_real_pexsi,&
                 eh%ovlp_real_pexsi,eval)
         call elsi_sips_to_blacs_ev_real(eh,eh%ph,eh%bh,evec)
      endif
   case default
      call elsi_stop(eh%bh,"Unsupported eigensolver.",caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors. Note the
!! intent(inout) - it is because everything has the potential to be reused in
!! the next call.
!!
subroutine elsi_ev_complex(eh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh                              !< Handle
   complex(kind=r8),  intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol)  !< Hamiltonian
   complex(kind=r8),  intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   real(kind=r8),     intent(inout) :: eval(eh%ph%n_basis)             !< Eigenvalues
   complex(kind=r8),  intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   real(kind=r8)               :: t0
   character(len=DATETIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_ev_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      if(eh%ph%parallel_mode == SINGLE_PROC) then
         call elsi_solve_lapack_cmplx(eh%ph,eh%bh,ham,ovlp,eval,evec)
      else
         call elsi_init_elpa(eh%ph,eh%bh)
         call elsi_solve_elpa_cmplx(eh,eh%ph,eh%bh,ham,ovlp,eval,evec)
      endif
   case default
      call elsi_stop(eh%bh,"Unsupported eigensolver.",caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors. Note the
!! intent(inout) - it is because everything has the potential to be reused in
!! the next call.
!!
subroutine elsi_ev_real_sparse(eh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh                              !< Handle
   real(kind=r8),     intent(inout) :: ham(eh%bh%nnz_l_sp)             !< Hamiltonian
   real(kind=r8),     intent(inout) :: ovlp(eh%bh%nnz_l_sp)            !< Overlap
   real(kind=r8),     intent(inout) :: eval(eh%ph%n_basis)             !< Eigenvalues
   real(kind=r8),     intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   real(kind=r8)               :: t0
   character(len=DATETIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_ev_real_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_solve_elpa_real(eh,eh%ph,eh%bh,eh%ham_real_elpa,&
              eh%ovlp_real_elpa,eval,evec)
   case(SIPS_SOLVER)
      if(eh%ph%n_calls <= eh%ph%sips_n_elpa) then
         call elsi_init_elpa(eh%ph,eh%bh)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_sips_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
         case(SIESTA_CSC)
            call elsi_siesta_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_solve_elpa_real(eh,eh%ph,eh%bh,eh%ham_real_elpa,&
                 eh%ovlp_real_elpa,eval,evec)
      else ! ELPA is done
         if(allocated(eh%ham_real_elpa)) then
            call elsi_deallocate(eh%bh,eh%ham_real_elpa,"ham_real_elpa")
         endif
         if(allocated(eh%ovlp_real_elpa)) then
            call elsi_deallocate(eh%bh,eh%ovlp_real_elpa,"ovlp_real_elpa")
         endif

         call elsi_init_sips(eh%ph,eh%bh)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_solve_sips_real(eh,eh%ph,eh%bh,ham,ovlp,eval)
         case(SIESTA_CSC)
            call elsi_siesta_to_sips_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
            call elsi_solve_sips_real(eh,eh%ph,eh%bh,eh%ham_real_pexsi,&
                    eh%ovlp_real_pexsi,eval)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_sips_to_blacs_ev_real(eh,eh%ph,eh%bh,evec)
      endif
   case default
      call elsi_stop(eh%bh,"Unsupported eigensolver.",caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! This routine computes the eigenvalues and eigenvectors. Note the
!! intent(inout) - it is because everything has the potential to be reused in
!! the next call.
!!
subroutine elsi_ev_complex_sparse(eh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh                              !< Handle
   complex(kind=r8),  intent(inout) :: ham(eh%bh%nnz_l_sp)             !< Hamiltonian
   complex(kind=r8),  intent(inout) :: ovlp(eh%bh%nnz_l_sp)            !< Overlap
   real(kind=r8),     intent(inout) :: eval(eh%ph%n_basis)             !< Eigenvalues
   complex(kind=r8),  intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   real(kind=r8)               :: t0
   character(len=DATETIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_ev_complex_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_cmplx(eh,eh%ph,eh%bh,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_cmplx(eh,eh%ph,eh%bh,ham,ovlp)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_solve_elpa_cmplx(eh,eh%ph,eh%bh,eh%ham_cmplx_elpa,&
              eh%ovlp_cmplx_elpa,eval,evec)
   case default
      call elsi_stop(eh%bh,"Unsupported eigensolver.",caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! This routine computes the density matrix. Note the intent(inout) - it is
!! because everything has the potential to be reused in the next call.
!!
subroutine elsi_dm_real(eh,ham,ovlp,dm,energy)

   implicit none

   type(elsi_handle), intent(inout) :: eh                              !< Handle
   real(kind=r8),     intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol)  !< Hamiltonian
   real(kind=r8),     intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   real(kind=r8),     intent(inout) :: dm(eh%bh%n_lrow,eh%bh%n_lcol)   !< Density matrix
   real(kind=r8),     intent(inout) :: energy                          !< Energy

   real(kind=r8)               :: t0
   character(len=DATETIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_dm_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(.not. allocated(eh%eval_elpa)) then
         call elsi_allocate(eh%bh,eh%eval_elpa,eh%ph%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(eh%evec_real_elpa)) then
         call elsi_allocate(eh%bh,eh%evec_real_elpa,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "evec_real_elpa",caller)
      endif

      call elsi_solve_elpa_real(eh,eh%ph,eh%bh,ham,ovlp,eh%eval_elpa,&
              eh%evec_real_elpa)
      call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
      call elsi_compute_dm_elpa_real(eh,eh%ph,eh%bh,eh%evec_real_elpa,dm,ham)
      call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
   case(OMM_SOLVER)
      if(eh%ph%n_calls <= eh%ph%omm_n_elpa) then
         call elsi_init_elpa(eh%ph,eh%bh)

         if(eh%ph%n_calls == 1 .and. eh%ph%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_real_copy",caller)
            eh%ovlp_real_copy = ovlp
         endif

         if(.not. allocated(eh%eval_elpa)) then
            call elsi_allocate(eh%bh,eh%eval_elpa,eh%ph%n_basis,"eval_elpa",&
                    caller)
         endif
         if(.not. allocated(eh%evec_real_elpa)) then
            call elsi_allocate(eh%bh,eh%evec_real_elpa,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"evec_real_elpa",caller)
         endif

         call elsi_solve_elpa_real(eh,eh%ph,eh%bh,ham,ovlp,eh%eval_elpa,&
                 eh%evec_real_elpa)
         call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
         call elsi_compute_dm_elpa_real(eh,eh%ph,eh%bh,eh%evec_real_elpa,dm,ham)
         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         call elsi_init_omm(eh%ph,eh%bh)

         if(allocated(eh%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = eh%ovlp_real_copy
            call elsi_deallocate(eh%bh,eh%ovlp_real_copy,"ovlp_real_copy")
         endif

         if(.not. eh%c_omm%is_initialized) then
            call elsi_init_coeff_omm_real(eh,eh%ph)
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(eh%ph%omm_n_elpa > 0 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pdtran(eh%ph%n_basis,eh%ph%n_basis,1.0_r8,eh%evec_real_elpa,1,&
                    1,eh%bh%desc,0.0_r8,dm,1,1,eh%bh%desc)

            eh%c_omm%dval(1:eh%c_omm%iaux2(1),1:eh%c_omm%iaux2(2)) = &
               dm(1:eh%c_omm%iaux2(1),1:eh%c_omm%iaux2(2))

            if(allocated(eh%evec_real_elpa)) then
               call elsi_deallocate(eh%bh,eh%evec_real_elpa,"evec_real_elpa")
            endif
            if(allocated(eh%eval_elpa)) then
               call elsi_deallocate(eh%bh,eh%eval_elpa,"eval_elpa")
            endif
            if(allocated(eh%eval_all)) then
               call elsi_deallocate(eh%bh,eh%eval_all,"eval_all")
            endif
            if(allocated(eh%occ_num)) then
               call elsi_deallocate(eh%bh,eh%occ_num,"occ_num")
            endif
            if(allocated(eh%k_weight)) then
               call elsi_deallocate(eh%bh,eh%k_weight,"k_weight")
            endif
         endif

         call elsi_solve_omm_real(eh,eh%ph,eh%bh,ham,ovlp,dm)
         call elsi_get_energy(eh%ph,eh%bh,energy,OMM_SOLVER)
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)
      call elsi_blacs_to_pexsi_hs_real(eh,eh%ph,eh%bh,ham,ovlp)

      if(.not. allocated(eh%dm_real_pexsi)) then
         call elsi_allocate(eh%bh,eh%dm_real_pexsi,eh%bh%nnz_l_sp,&
                 "dm_real_pexsi",caller)
      endif
      eh%dm_real_pexsi = 0.0_r8

      call elsi_solve_pexsi_real(eh,eh%ph,eh%bh,eh%ham_real_pexsi,&
              eh%ovlp_real_pexsi,eh%dm_real_pexsi)
      call elsi_pexsi_to_blacs_dm_real(eh,eh%ph,eh%bh,dm)
      call elsi_get_energy(eh%ph,eh%bh,energy,PEXSI_SOLVER)
   case(SIPS_SOLVER)
      if(.not. allocated(eh%eval_elpa)) then
         call elsi_allocate(eh%bh,eh%eval_elpa,eh%ph%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(eh%evec_real_elpa)) then
         call elsi_allocate(eh%bh,eh%evec_real_elpa,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "evec_real_elpa",caller)
      endif

      if(eh%ph%n_calls <= eh%ph%sips_n_elpa) then
         call elsi_init_elpa(eh%ph,eh%bh)

         if(eh%ph%n_calls == 1) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_real_copy",caller)
            eh%ovlp_real_copy = ovlp
         endif

         call elsi_solve_elpa_real(eh,eh%ph,eh%bh,ham,ovlp,eh%eval_elpa,&
                 eh%evec_real_elpa)
         call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
         call elsi_compute_dm_elpa_real(eh,eh%ph,eh%bh,eh%evec_real_elpa,dm,ham)
         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         if(allocated(eh%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = eh%ovlp_real_copy
            call elsi_deallocate(eh%bh,eh%ovlp_real_copy,"ovlp_real_copy")
            call elsi_deallocate(eh%bh,eh%evec_real_elpa,"evec_real_elpa")
         endif

         call elsi_init_sips(eh%ph,eh%bh)
         call elsi_blacs_to_sips_hs_real(eh,eh%ph,eh%bh,ham,ovlp)

         if(.not. allocated(eh%dm_real_pexsi)) then
            call elsi_allocate(eh%bh,eh%dm_real_pexsi,eh%bh%nnz_l_sp,&
                    "dm_real_pexsi",caller)
         endif
         eh%dm_real_pexsi = 0.0_r8

         call elsi_solve_sips_real(eh,eh%ph,eh%bh,eh%ham_real_pexsi,&
                 eh%ovlp_real_pexsi,eh%eval_elpa)
         call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
         call elsi_compute_dm_sips_real(eh,eh%ph,eh%bh,eh%dm_real_pexsi)
         call elsi_sips_to_blacs_dm_real(eh,eh%ph,eh%bh,dm)
         call elsi_get_energy(eh%ph,eh%bh,energy,SIPS_SOLVER)
      endif
   case(DMP_SOLVER)
      ! Save Hamiltonian and overlap
      if(eh%ph%n_calls==1) then
         call elsi_allocate(eh%bh,eh%ham_real_copy,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "ham_real_copy",caller)
         call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "ovlp_real_copy",caller)
         eh%ovlp_real_copy = ovlp
      endif
      eh%ham_real_copy = ham

      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_init_dmp(eh%ph)
      call elsi_solve_dmp_real(eh,eh%ph,eh%bh,ham,ovlp,dm)
      call elsi_get_energy(eh%ph,eh%bh,energy,DMP_SOLVER)
   case default
      call elsi_stop(eh%bh,"Unsupported density matrix solver.",caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

   eh%ph%edm_ready_real = .true.

end subroutine

!>
!! This routine computes the density matrix. Note the intent(inout) - it is
!! because everything has the potential to be reused in the next call.
!!
subroutine elsi_dm_complex(eh,ham,ovlp,dm,energy)

   implicit none

   type(elsi_handle), intent(inout) :: eh                              !< Handle
   complex(kind=r8),  intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol)  !< Hamiltonian
   complex(kind=r8),  intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   complex(kind=r8),  intent(inout) :: dm(eh%bh%n_lrow,eh%bh%n_lcol)   !< Density matrix
   real(kind=r8),     intent(inout) :: energy                          !< Energy

   real(kind=r8)               :: t0
   character(len=DATETIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_dm_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      if(.not. allocated(eh%eval_elpa)) then
         call elsi_allocate(eh%bh,eh%eval_elpa,eh%ph%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(eh%evec_cmplx_elpa)) then
         call elsi_allocate(eh%bh,eh%evec_cmplx_elpa,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "evec_cmplx_elpa",caller)
      endif

      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_solve_elpa_cmplx(eh,eh%ph,eh%bh,ham,ovlp,eh%eval_elpa,&
              eh%evec_cmplx_elpa)
      call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
      call elsi_compute_dm_elpa_cmplx(eh,eh%ph,eh%bh,eh%evec_cmplx_elpa,dm,ham)
      call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
   case(OMM_SOLVER)
      if(eh%ph%n_calls <= eh%ph%omm_n_elpa) then
         call elsi_init_elpa(eh%ph,eh%bh)

         if(eh%ph%n_calls == 1 .and. eh%ph%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_copy,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_cmplx_copy",caller)
            eh%ovlp_cmplx_copy = ovlp
         endif

         if(.not. allocated(eh%eval_elpa)) then
            call elsi_allocate(eh%bh,eh%eval_elpa,eh%ph%n_basis,"eval_elpa",&
                    caller)
         endif
         if(.not. allocated(eh%evec_cmplx_elpa)) then
            call elsi_allocate(eh%bh,eh%evec_cmplx_elpa,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"evec_cmplx_elpa",caller)
         endif

         call elsi_solve_elpa_cmplx(eh,eh%ph,eh%bh,ham,ovlp,eh%eval_elpa,&
                 eh%evec_cmplx_elpa)
         call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
         call elsi_compute_dm_elpa_cmplx(eh,eh%ph,eh%bh,eh%evec_cmplx_elpa,dm,&
                 ham)
         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         call elsi_init_omm(eh%ph,eh%bh)

         if(allocated(eh%ovlp_cmplx_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = eh%ovlp_cmplx_copy
            call elsi_deallocate(eh%bh,eh%ovlp_cmplx_copy,"ovlp_cmplx_copy")
         endif

         if(.not. eh%c_omm%is_initialized) then
            call elsi_init_coeff_omm_cmplx(eh,eh%ph)
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(eh%ph%omm_n_elpa > 0 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pztranc(eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),&
                    eh%evec_cmplx_elpa,1,1,eh%bh%desc,(0.0_r8,0.0_r8),dm,1,1,&
                    eh%bh%desc)

            eh%c_omm%zval(1:eh%c_omm%iaux2(1),1:eh%c_omm%iaux2(2)) = &
               dm(1:eh%c_omm%iaux2(1),1:eh%c_omm%iaux2(2))

            if(allocated(eh%evec_cmplx_elpa)) then
               call elsi_deallocate(eh%bh,eh%evec_cmplx_elpa,"evec_cmplx_elpa")
            endif
            if(allocated(eh%eval_elpa)) then
               call elsi_deallocate(eh%bh,eh%eval_elpa,"eval_elpa")
            endif
            if(allocated(eh%eval_all)) then
               call elsi_deallocate(eh%bh,eh%eval_all,"eval_all")
            endif
            if(allocated(eh%occ_num)) then
               call elsi_deallocate(eh%bh,eh%occ_num,"occ_num")
            endif
            if(allocated(eh%k_weight)) then
               call elsi_deallocate(eh%bh,eh%k_weight,"k_weight")
            endif
         endif

         call elsi_solve_omm_cmplx(eh,eh%ph,eh%bh,ham,ovlp,dm)
         call elsi_get_energy(eh%ph,eh%bh,energy,OMM_SOLVER)
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)
      call elsi_blacs_to_pexsi_hs_cmplx(eh,eh%ph,eh%bh,ham,ovlp)

      if(.not. allocated(eh%dm_cmplx_pexsi)) then
         call elsi_allocate(eh%bh,eh%dm_cmplx_pexsi,eh%bh%nnz_l_sp,&
                 "dm_cmplx_pexsi",caller)
      endif
      eh%dm_cmplx_pexsi = (0.0_r8,0.0_r8)

      call elsi_solve_pexsi_cmplx(eh,eh%ph,eh%bh,eh%ham_cmplx_pexsi,&
              eh%ovlp_cmplx_pexsi,eh%dm_cmplx_pexsi)
      call elsi_pexsi_to_blacs_dm_cmplx(eh,eh%ph,eh%bh,dm)
      call elsi_get_energy(eh%ph,eh%bh,energy,PEXSI_SOLVER)
   case default
      call elsi_stop(eh%bh,"Unsupported density matrix solver.",caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

   eh%ph%edm_ready_cmplx = .true.

end subroutine

!>
!! This routine computes the density matrix. Note the intent(inout) - it is
!! because everything has the potential to be reused in the next call.
!!
subroutine elsi_dm_real_sparse(eh,ham,ovlp,dm,energy)

   implicit none

   type(elsi_handle), intent(inout) :: eh                   !< Handle
   real(kind=r8),     intent(inout) :: ham(eh%bh%nnz_l_sp)  !< Hamiltonian
   real(kind=r8),     intent(inout) :: ovlp(eh%bh%nnz_l_sp) !< Overlap
   real(kind=r8),     intent(inout) :: dm(eh%bh%nnz_l_sp)   !< Density matrix
   real(kind=r8),     intent(inout) :: energy               !< Energy

   real(kind=r8)               :: t0
   character(len=DATETIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_dm_real_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. eh%bh%blacs_ready) then
         call elsi_init_blacs(eh)
      endif

      call elsi_init_elpa(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      if(.not. allocated(eh%eval_elpa)) then
         call elsi_allocate(eh%bh,eh%eval_elpa,eh%ph%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(eh%evec_real_elpa)) then
         call elsi_allocate(eh%bh,eh%evec_real_elpa,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "evec_real_elpa",caller)
      endif
      if(.not. allocated(eh%dm_real_elpa)) then
         call elsi_allocate(eh%bh,eh%dm_real_elpa,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "dm_real_elpa",caller)
      endif

      call elsi_solve_elpa_real(eh,eh%ph,eh%bh,eh%ham_real_elpa,&
              eh%ovlp_real_elpa,eh%eval_elpa,eh%evec_real_elpa)
      call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
      call elsi_compute_dm_elpa_real(eh,eh%ph,eh%bh,eh%evec_real_elpa,&
              eh%dm_real_elpa,eh%ham_real_elpa)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm_real(eh,eh%ph,eh%bh,dm)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm_real(eh,eh%bh,dm)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
   case(OMM_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. eh%bh%blacs_ready) then
         call elsi_init_blacs(eh)
      endif

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      if(eh%ph%n_calls <= eh%ph%omm_n_elpa) then
         call elsi_init_elpa(eh%ph,eh%bh)

         if(eh%ph%n_calls == 1 .and. eh%ph%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_real_copy",caller)
            eh%ovlp_real_copy = eh%ovlp_real_elpa
         endif

         if(.not. allocated(eh%eval_elpa)) then
            call elsi_allocate(eh%bh,eh%eval_elpa,eh%ph%n_basis,"eval_elpa",&
                    caller)
         endif
         if(.not. allocated(eh%evec_real_elpa)) then
            call elsi_allocate(eh%bh,eh%evec_real_elpa,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"evec_real_elpa",caller)
         endif
         if(.not. allocated(eh%dm_real_elpa)) then
            call elsi_allocate(eh%bh,eh%dm_real_elpa,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "dm_real_elpa",caller)
         endif

         call elsi_solve_elpa_real(eh,eh%ph,eh%bh,eh%ham_real_elpa,&
                 eh%ovlp_real_elpa,eh%eval_elpa,eh%evec_real_elpa)
         call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
         call elsi_compute_dm_elpa_real(eh,eh%ph,eh%bh,eh%evec_real_elpa,&
                 eh%dm_real_elpa,eh%ham_real_elpa)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_real(eh,eh%ph,eh%bh,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_real(eh,eh%bh,dm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         call elsi_init_omm(eh%ph,eh%bh)

         if(allocated(eh%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            eh%ovlp_real_elpa = eh%ovlp_real_copy
            call elsi_deallocate(eh%bh,eh%ovlp_real_copy,"ovlp_real_copy")
         endif

         if(.not. eh%c_omm%is_initialized) then
            call elsi_init_coeff_omm_real(eh,eh%ph)
         endif
         if(.not. allocated(eh%dm_real_elpa)) then
            call elsi_allocate(eh%bh,eh%dm_real_elpa,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "dm_real_elpa",caller)
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(eh%ph%omm_n_elpa > 0 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pdtran(eh%ph%n_basis,eh%ph%n_basis,1.0_r8,eh%evec_real_elpa,1,&
                    1,eh%bh%desc,0.0_r8,eh%dm_real_elpa,1,1,eh%bh%desc)

            eh%c_omm%dval(1:eh%c_omm%iaux2(1),1:eh%c_omm%iaux2(2)) = &
               eh%dm_real_elpa(1:eh%c_omm%iaux2(1),1:eh%c_omm%iaux2(2))

            if(allocated(eh%evec_real_elpa)) then
               call elsi_deallocate(eh%bh,eh%evec_real_elpa,"evec_real_elpa")
            endif
            if(allocated(eh%eval_elpa)) then
               call elsi_deallocate(eh%bh,eh%eval_elpa,"eval_elpa")
            endif
            if(allocated(eh%occ_num)) then
               call elsi_deallocate(eh%bh,eh%occ_num,"occ_num")
            endif
         endif

         call elsi_solve_omm_real(eh,eh%ph,eh%bh,eh%ham_real_elpa,&
                 eh%ovlp_real_elpa,eh%dm_real_elpa)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_real(eh,eh%ph,eh%bh,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_real(eh,eh%bh,dm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(eh%ph,eh%bh,energy,OMM_SOLVER)
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_solve_pexsi_real(eh,eh%ph,eh%bh,ham,ovlp,dm)
      case(SIESTA_CSC)
         call elsi_siesta_to_pexsi_hs_real(eh,eh%ph,eh%bh,ham,ovlp)

         if(.not. allocated(eh%dm_real_pexsi)) then
            call elsi_allocate(eh%bh,eh%dm_real_pexsi,eh%bh%nnz_l_sp1,&
                    "dm_real_pexsi",caller)
         endif
         eh%dm_real_pexsi = 0.0_r8

         call elsi_solve_pexsi_real(eh,eh%ph,eh%bh,eh%ham_real_pexsi,&
                 eh%ovlp_real_pexsi,eh%dm_real_pexsi)
         call elsi_pexsi_to_siesta_dm_real(eh,eh%ph,eh%bh,dm)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(eh%ph,eh%bh,energy,PEXSI_SOLVER)
   case(SIPS_SOLVER)
      if(eh%ph%n_calls <= eh%ph%sips_n_elpa) then
         ! Set up BLACS if not done by user
         if(.not. eh%bh%blacs_ready) then
            call elsi_init_blacs(eh)
         endif

         call elsi_init_elpa(eh%ph,eh%bh)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_sips_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
         case(SIESTA_CSC)
            call elsi_siesta_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         if(.not. allocated(eh%eval_elpa)) then
            call elsi_allocate(eh%bh,eh%eval_elpa,eh%ph%n_basis,"eval_elpa",&
                    caller)
         endif
         if(.not. allocated(eh%evec_real_elpa)) then
            call elsi_allocate(eh%bh,eh%evec_real_elpa,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"evec_real_elpa",caller)
         endif
         if(.not. allocated(eh%dm_real_elpa)) then
            call elsi_allocate(eh%bh,eh%dm_real_elpa,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "dm_real_elpa",caller)
         endif

         call elsi_solve_elpa_real(eh,eh%ph,eh%bh,eh%ham_real_elpa,&
                 eh%ovlp_real_elpa,eh%eval_elpa,eh%evec_real_elpa)
         call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
         call elsi_compute_dm_elpa_real(eh,eh%ph,eh%bh,eh%evec_real_elpa,&
                 eh%dm_real_elpa,eh%ham_real_elpa)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_real(eh,eh%ph,eh%bh,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_real(eh,eh%bh,dm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         if(allocated(eh%ham_real_elpa)) then
            call elsi_deallocate(eh%bh,eh%ham_real_elpa,"ham_real_elpa")
         endif
         if(allocated(eh%ovlp_real_elpa)) then
            call elsi_deallocate(eh%bh,eh%ovlp_real_elpa,"ovlp_real_elpa")
         endif
         if(allocated(eh%dm_real_elpa)) then
            call elsi_deallocate(eh%bh,eh%dm_real_elpa,"dm_real_elpa")
         endif
         if(.not. allocated(eh%eval_elpa)) then
            call elsi_allocate(eh%bh,eh%eval_elpa,eh%ph%n_basis,"eval_elpa",&
                    caller)
         endif

         call elsi_init_sips(eh%ph,eh%bh)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_solve_sips_real(eh,eh%ph,eh%bh,ham,ovlp,eh%eval_elpa)
            call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
            call elsi_compute_dm_sips_real(eh,eh%ph,eh%bh,dm)
            call elsi_get_energy(eh%ph,eh%bh,energy,SIPS_SOLVER)
         case(SIESTA_CSC)
            call elsi_siesta_to_sips_hs_real(eh,eh%ph,eh%bh,ham,ovlp)

            if(.not. allocated(eh%dm_real_pexsi)) then
               call elsi_allocate(eh%bh,eh%dm_real_pexsi,eh%bh%nnz_l_sp1,&
                       "dm_real_pexsi",caller)
            endif
            eh%dm_real_pexsi = 0.0_r8

            call elsi_solve_sips_real(eh,eh%ph,eh%bh,eh%ham_real_pexsi,&
                    eh%ovlp_real_pexsi,eh%eval_elpa)
            call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
            call elsi_compute_dm_sips_real(eh,eh%ph,eh%bh,eh%dm_real_pexsi)
            call elsi_sips_to_siesta_dm_real(eh,eh%ph,eh%bh,dm)
            call elsi_get_energy(eh%ph,eh%bh,energy,SIPS_SOLVER)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select
      endif
   case(DMP_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. eh%bh%blacs_ready) then
         call elsi_init_blacs(eh)
      endif

      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_init_dmp(eh%ph)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_real(eh,eh%ph,eh%bh,ham,ovlp)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      if(.not. allocated(eh%dm_real_elpa)) then
         call elsi_allocate(eh%bh,eh%dm_real_elpa,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "dm_real_elpa",caller)
      endif

      ! Save Hamiltonian and overlap
      if(eh%ph%n_calls==1) then
         call elsi_allocate(eh%bh,eh%ham_real_copy,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "ham_real_copy",caller)
         call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "ovlp_real_copy",caller)
         eh%ham_real_copy  = eh%ham_real_elpa
         eh%ovlp_real_copy = eh%ovlp_real_elpa
      endif

      call elsi_solve_dmp_real(eh,eh%ph,eh%bh,eh%ham_real_elpa,&
              eh%ovlp_real_elpa,eh%dm_real_elpa)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm_real(eh,eh%ph,eh%bh,dm)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm_real(eh,eh%bh,dm)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(eh%ph,eh%bh,energy,DMP_SOLVER)
   case default
      call elsi_stop(eh%bh,"Unsupported density matrix solver.",caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

   eh%ph%edm_ready_real = .true.

end subroutine

!>
!! This routine computes the density matrix. Note the intent(inout) - it is
!! because everything has the potential to be reused in the next call.
!!
subroutine elsi_dm_complex_sparse(eh,ham,ovlp,dm,energy)

   implicit none

   type(elsi_handle), intent(inout) :: eh                   !< Handle
   complex(kind=r8),  intent(inout) :: ham(eh%bh%nnz_l_sp)  !< Hamiltonian
   complex(kind=r8),  intent(inout) :: ovlp(eh%bh%nnz_l_sp) !< Overlap
   complex(kind=r8),  intent(inout) :: dm(eh%bh%nnz_l_sp)   !< Density matrix
   real(kind=r8),     intent(inout) :: energy               !< Energy

   real(kind=r8)               :: t0
   character(len=DATETIME_LEN) :: dt0

   character(len=40), parameter :: caller = "elsi_dm_complex_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. eh%bh%blacs_ready) then
         call elsi_init_blacs(eh)
      endif

      call elsi_init_elpa(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_cmplx(eh,eh%ph,eh%bh,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_cmplx(eh,eh%ph,eh%bh,ham,ovlp)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      if(.not. allocated(eh%eval_elpa)) then
         call elsi_allocate(eh%bh,eh%eval_elpa,eh%ph%n_basis,"eval_elpa",caller)
      endif
      if(.not. allocated(eh%evec_cmplx_elpa)) then
         call elsi_allocate(eh%bh,eh%evec_cmplx_elpa,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "evec_cmplx_elpa",caller)
      endif
      if(.not. allocated(eh%dm_cmplx_elpa)) then
         call elsi_allocate(eh%bh,eh%dm_cmplx_elpa,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "dm_cmplx_elpa",caller)
      endif

      call elsi_solve_elpa_cmplx(eh,eh%ph,eh%bh,eh%ham_cmplx_elpa,&
              eh%ovlp_cmplx_elpa,eh%eval_elpa,eh%evec_cmplx_elpa)
      call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
      call elsi_compute_dm_elpa_cmplx(eh,eh%ph,eh%bh,eh%evec_cmplx_elpa,&
              eh%dm_cmplx_elpa,eh%ham_cmplx_elpa)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm_cmplx(eh,eh%ph,eh%bh,dm)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm_cmplx(eh,eh%bh,dm)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
   case(OMM_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. eh%bh%blacs_ready) then
         call elsi_init_blacs(eh)
      endif

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs_cmplx(eh,eh%ph,eh%bh,ham,ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs_cmplx(eh,eh%ph,eh%bh,ham,ovlp)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      if(eh%ph%n_calls <= eh%ph%omm_n_elpa) then
         call elsi_init_elpa(eh%ph,eh%bh)

         if(eh%ph%n_calls == 1 .and. eh%ph%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_copy,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_cmplx_copy",caller)
            eh%ovlp_cmplx_copy = eh%ovlp_cmplx_elpa
         endif

         if(.not. allocated(eh%eval_elpa)) then
            call elsi_allocate(eh%bh,eh%eval_elpa,eh%ph%n_basis,"eval_elpa",&
                    caller)
         endif
         if(.not. allocated(eh%evec_cmplx_elpa)) then
            call elsi_allocate(eh%bh,eh%evec_cmplx_elpa,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"evec_cmplx_elpa",caller)
         endif
         if(.not. allocated(eh%dm_cmplx_elpa)) then
            call elsi_allocate(eh%bh,eh%dm_cmplx_elpa,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"dm_cmplx_elpa",caller)
         endif

         call elsi_solve_elpa_cmplx(eh,eh%ph,eh%bh,eh%ham_cmplx_elpa,&
                 eh%ovlp_cmplx_elpa,eh%eval_elpa,eh%evec_cmplx_elpa)
         call elsi_compute_occ_elpa(eh,eh%ph,eh%bh,eh%eval_elpa)
         call elsi_compute_dm_elpa_cmplx(eh,eh%ph,eh%bh,eh%evec_cmplx_elpa,&
                 eh%dm_cmplx_elpa,eh%ham_cmplx_elpa)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_cmplx(eh,eh%ph,eh%bh,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_cmplx(eh,eh%bh,dm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         call elsi_init_omm(eh%ph,eh%bh)

         if(allocated(eh%ovlp_cmplx_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            eh%ovlp_cmplx_elpa = eh%ovlp_cmplx_copy
            call elsi_deallocate(eh%bh,eh%ovlp_cmplx_copy,"ovlp_cmplx_copy")
         endif

         if(.not. eh%c_omm%is_initialized) then
            call elsi_init_coeff_omm_cmplx(eh,eh%ph)
         endif
         if(.not. allocated(eh%dm_cmplx_elpa)) then
            call elsi_allocate(eh%bh,eh%dm_cmplx_elpa,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"dm_cmplx_elpa",caller)
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(eh%ph%omm_n_elpa > 0 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
            ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
            call pztranc(eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),&
                    eh%evec_cmplx_elpa,1,1,eh%bh%desc,(0.0_r8,0.0_r8),&
                    eh%dm_cmplx_elpa,1,1,eh%bh%desc)

            eh%c_omm%zval(1:eh%c_omm%iaux2(1),1:eh%c_omm%iaux2(2)) = &
               eh%dm_cmplx_elpa(1:eh%c_omm%iaux2(1),1:eh%c_omm%iaux2(2))

            if(allocated(eh%evec_cmplx_elpa)) then
               call elsi_deallocate(eh%bh,eh%evec_cmplx_elpa,"evec_cmplx_elpa")
            endif
            if(allocated(eh%eval_elpa)) then
               call elsi_deallocate(eh%bh,eh%eval_elpa,"eval_elpa")
            endif
            if(allocated(eh%occ_num)) then
               call elsi_deallocate(eh%bh,eh%occ_num,"occ_num")
            endif
         endif

         call elsi_solve_omm_cmplx(eh,eh%ph,eh%bh,eh%ham_cmplx_elpa,&
                 eh%ovlp_cmplx_elpa,eh%dm_cmplx_elpa)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm_cmplx(eh,eh%ph,eh%bh,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm_cmplx(eh,eh%bh,dm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(eh%ph,eh%bh,energy,OMM_SOLVER)
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_solve_pexsi_cmplx(eh,eh%ph,eh%bh,ham,ovlp,dm)
      case(SIESTA_CSC)
         call elsi_siesta_to_pexsi_hs_cmplx(eh,eh%ph,eh%bh,ham,ovlp)

         if(.not. allocated(eh%dm_cmplx_pexsi)) then
            call elsi_allocate(eh%bh,eh%dm_cmplx_pexsi,eh%bh%nnz_l_sp1,&
                    "dm_cmplx_pexsi",caller)
         endif
         eh%dm_cmplx_pexsi = (0.0_r8,0.0_r8)

         call elsi_solve_pexsi_cmplx(eh,eh%ph,eh%bh,eh%ham_cmplx_pexsi,&
                 eh%ovlp_cmplx_pexsi,eh%dm_cmplx_pexsi)
         call elsi_pexsi_to_siesta_dm_cmplx(eh,eh%ph,eh%bh,dm)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(eh%ph,eh%bh,energy,PEXSI_SOLVER)
   case default
      call elsi_stop(eh%bh,"Unsupported density matrix solver.",caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

   eh%ph%edm_ready_cmplx = .true.

end subroutine

!>
!! This routine initializes BLACS, in case that the user selects a sparse format
!! and BLACS is still used internally.
!!
subroutine elsi_init_blacs(eh)

   implicit none

   type(elsi_handle), intent(inout) :: eh

   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: block_size

   character(len=40), parameter :: caller = "elsi_init_blacs"

   if(eh%ph%parallel_mode == MULTI_PROC .and. .not. eh%bh%blacs_ready) then
      ! Set square-like process grid
      do nprow = nint(sqrt(real(eh%bh%n_procs,kind=r8))),2,-1
         if(mod(eh%bh%n_procs,nprow) == 0) then
            exit
         endif
      enddo

      npcol = eh%bh%n_procs/nprow

      if(max(nprow,npcol) > eh%ph%n_basis) then
         call elsi_stop(eh%bh,"Matrix size is too small for this number of"//&
                 " MPI tasks.",caller)
      endif

      ! Initialize BLACS
      blacs_ctxt = eh%bh%comm

      call BLACS_Gridinit(blacs_ctxt,'r',nprow,npcol)

      ! Find block size
      block_size = 1

      do while(2*block_size*max(nprow,npcol) <= eh%ph%n_basis)
         block_size = 2*block_size
      enddo

      ! Maximum allowed value: 256
      block_size = min(256,block_size)

      ! ELPA works better with a small block_size
      if(eh%ph%solver == ELPA_SOLVER) then
         block_size = min(32,block_size)
      endif

      call elsi_set_blacs(eh,blacs_ctxt,block_size)
   endif

end subroutine

end module ELSI_SOLVER
