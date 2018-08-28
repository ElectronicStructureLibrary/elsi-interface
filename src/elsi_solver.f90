! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains subroutines to solve an eigenproblem or to compute the
!! density matrix, using one of the solvers: ELPA, libOMM, PEXSI, SLEPc-SIPs,
!! and NTPoly.
!!
module ELSI_SOLVER

   use ELSI_CONSTANTS, only: ELPA_SOLVER,OMM_SOLVER,PEXSI_SOLVER,SIPS_SOLVER,&
                             NTPOLY_SOLVER,MULTI_PROC,SINGLE_PROC,PEXSI_CSC,&
                             SIESTA_CSC
   use ELSI_DATATYPE,  only: elsi_handle,elsi_param_t,elsi_basic_t
   use ELSI_NTPOLY,    only: elsi_init_ntpoly,elsi_solve_ntpoly
   use ELSI_ELPA,      only: elsi_init_elpa,elsi_compute_occ_elpa,&
                             elsi_compute_dm_elpa,elsi_solve_elpa
   use ELSI_IO,        only: elsi_add_log,elsi_get_time,&
                             fjson_get_datetime_rfc3339
   use ELSI_LAPACK,    only: elsi_solve_lapack
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_REDIST,    only: elsi_blacs_to_ntpoly_hs,&
                             elsi_blacs_to_pexsi_hs_dim,elsi_blacs_to_pexsi_hs,&
                             elsi_blacs_to_siesta_dm,elsi_blacs_to_sips_dm,&
                             elsi_blacs_to_sips_hs_dim,elsi_blacs_to_sips_hs,&
                             elsi_ntpoly_to_blacs_dm,elsi_ntpoly_to_siesta_dm,&
                             elsi_ntpoly_to_sips_dm,elsi_pexsi_to_blacs_dm,&
                             elsi_pexsi_to_siesta_dm,elsi_siesta_to_blacs_hs,&
                             elsi_siesta_to_ntpoly_hs,&
                             elsi_siesta_to_pexsi_hs_dim,&
                             elsi_siesta_to_pexsi_hs,&
                             elsi_siesta_to_sips_hs_dim,elsi_siesta_to_sips_hs,&
                             elsi_sips_to_blacs_dm,elsi_sips_to_blacs_ev,&
                             elsi_sips_to_blacs_hs,elsi_sips_to_ntpoly_hs,&
                             elsi_sips_to_siesta_dm
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8
   use ELSI_OMM,       only: elsi_init_omm,elsi_solve_omm
   use ELSI_PEXSI,     only: elsi_init_pexsi,elsi_solve_pexsi
   use ELSI_PRECISION, only: r8,i4
   use ELSI_SETUP,     only: elsi_set_blacs
   use ELSI_SIPS,      only: elsi_init_sips,elsi_solve_sips,elsi_compute_dm_sips
   use ELSI_UTILS,     only: elsi_check,elsi_check_init

   implicit none

   private

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
   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_get_energy"

   energy = ph%ebs*ph%i_weight

   if(solver == OMM_SOLVER) then
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

   real(kind=r8)     :: t0
   character(len=29) :: dt0

   character(len=*), parameter :: caller = "elsi_ev_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      if(eh%ph%parallel_mode == SINGLE_PROC) then
         call elsi_solve_lapack(eh%ph,eh%bh,ham,ovlp,eval,evec)
      else
         call elsi_init_elpa(eh%ph,eh%bh)
         call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,ham,ovlp,eval,&
                 evec)
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
         call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,ham,ovlp,eval,&
                 evec)
      else ! ELPA is done
         if(allocated(eh%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = eh%ovlp_real_copy
            call elsi_deallocate(eh%bh,eh%ovlp_real_copy,"ovlp_real_copy")
         endif

         call elsi_init_sips(eh%ph,eh%bh)

         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_blacs_to_sips_hs_dim(eh%ph,eh%bh,ham,ovlp)

            if(eh%ph%ovlp_is_unit) then
               call elsi_allocate(eh%bh,eh%ovlp_real_csc,1,"ovlp_real_csc",&
                       caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_real_csc,eh%bh%nnz_l_sp1,&
                       "ovlp_real_csc",caller)
            endif
            call elsi_allocate(eh%bh,eh%ham_real_csc,eh%bh%nnz_l_sp1,&
                    "ham_real_csc",caller)
            call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lcol_sp1,&
                    eh%ph%n_states,"evec_real",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                    "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                    "col_ptr_sp1",caller)
         endif

         call elsi_blacs_to_sips_hs(eh%ph,eh%bh,ham,ovlp,eh%ham_real_csc,&
                 eh%ovlp_real_csc,eh%row_ind_sp1,eh%col_ptr_sp1)
         call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%ham_real_csc,eh%ovlp_real_csc,eval,eh%evec_real)
         call elsi_sips_to_blacs_ev(eh%ph,eh%bh,eh%evec_real,evec)
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

   real(kind=r8)     :: t0
   character(len=29) :: dt0

   character(len=*), parameter :: caller = "elsi_ev_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      if(eh%ph%parallel_mode == SINGLE_PROC) then
         call elsi_solve_lapack(eh%ph,eh%bh,ham,ovlp,eval,evec)
      else
         call elsi_init_elpa(eh%ph,eh%bh)
         call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,ham,ovlp,eval,&
                 evec)
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

   real(kind=r8)     :: t0
   character(len=29) :: dt0

   character(len=*), parameter :: caller = "elsi_ev_real_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(.not. allocated(eh%ham_real_den)) then
         call elsi_allocate(eh%bh,eh%ham_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "ham_real_den",caller)
      endif
      if(.not. allocated(eh%ovlp_real_den)) then
         if(.not. eh%ph%ovlp_is_unit) then
            call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_real_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_den,1,1,"ovlp_real_den",&
                    caller)
         endif
      endif

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 ham,ovlp,eh%ham_real_den,eh%ovlp_real_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp2,&
                 eh%col_ptr_sp2,ham,ovlp,eh%ham_real_den,eh%ovlp_real_den)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,eh%ham_real_den,&
              eh%ovlp_real_den,eval,evec)
   case(SIPS_SOLVER)
      if(eh%ph%n_calls <= eh%ph%sips_n_elpa) then
         call elsi_init_elpa(eh%ph,eh%bh)

         if(.not. allocated(eh%ham_real_den)) then
            call elsi_allocate(eh%bh,eh%ham_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "ham_real_den",caller)
         endif
         if(.not. allocated(eh%ovlp_real_den)) then
            if(.not. eh%ph%ovlp_is_unit) then
               call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,&
                       eh%bh%n_lcol,"ovlp_real_den",caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_real_den,1,1,"ovlp_real_den",&
                       caller)
            endif
         endif

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_sips_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,ham,ovlp,eh%ham_real_den,eh%ovlp_real_den)
         case(SIESTA_CSC)
            call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp2,&
                    eh%col_ptr_sp2,ham,ovlp,eh%ham_real_den,eh%ovlp_real_den)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
                 eh%ham_real_den,eh%ovlp_real_den,eval,evec)
      else ! ELPA is done
         call elsi_init_sips(eh%ph,eh%bh)

         if(allocated(eh%ham_real_den)) then
            call elsi_deallocate(eh%bh,eh%ham_real_den,"ham_real_den")
         endif
         if(allocated(eh%ovlp_real_den)) then
            call elsi_deallocate(eh%bh,eh%ovlp_real_den,"ovlp_real_den")
         endif
         if(.not. allocated(eh%evec_real)) then
            call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lcol_sp1,&
                    eh%ph%n_states,"evec_real",caller)
         endif

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,ham,&
                    ovlp,eval,eh%evec_real)
         case(SIESTA_CSC)
            if(.not. allocated(eh%row_ind_sp1)) then
               call elsi_siesta_to_sips_hs_dim(eh%ph,eh%bh,eh%col_ptr_sp2)

               if(eh%ph%ovlp_is_unit) then
                  call elsi_allocate(eh%bh,eh%ovlp_real_csc,1,"ovlp_real_csc",&
                          caller)
               else
                  call elsi_allocate(eh%bh,eh%ovlp_real_csc,eh%bh%nnz_l_sp1,&
                          "ovlp_real_csc",caller)
               endif
               call elsi_allocate(eh%bh,eh%ham_real_csc,eh%bh%nnz_l_sp1,&
                       "ham_real_csc",caller)
               call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                       "row_ind_sp1",caller)
               call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                       "col_ptr_sp1",caller)
            endif

            call elsi_siesta_to_sips_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
                    eh%col_ptr_sp2,eh%ham_real_csc,eh%ovlp_real_csc,&
                    eh%row_ind_sp1,eh%col_ptr_sp1)
            call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                    eh%ham_real_csc,eh%ovlp_real_csc,eval,eh%evec_real)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_sips_to_blacs_ev(eh%ph,eh%bh,eh%evec_real,evec)
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

   real(kind=r8)     :: t0
   character(len=29) :: dt0

   character(len=*), parameter :: caller = "elsi_ev_complex_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(.not. allocated(eh%ham_cmplx_den)) then
         call elsi_allocate(eh%bh,eh%ham_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "ham_cmplx_den",caller)
      endif
      if(.not. allocated(eh%ovlp_cmplx_den)) then
         if(.not. eh%ph%ovlp_is_unit) then
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_cmplx_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,1,1,"ovlp_cmplx_den",&
                    caller)
         endif
      endif

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 ham,ovlp,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp2,&
                 eh%col_ptr_sp2,ham,ovlp,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,eh%ham_cmplx_den,&
              eh%ovlp_cmplx_den,eval,evec)
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

   real(kind=r8)     :: t0
   character(len=29) :: dt0

   character(len=*), parameter :: caller = "elsi_dm_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      endif
      if(.not. allocated(eh%evec_real)) then
         call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "evec_real",caller)
      endif
      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                 eh%ph%n_kpts,"occ",caller)
      endif

      call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,ham,ovlp,eh%eval,&
              eh%evec_real)
      call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
      call elsi_compute_dm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,eh%evec_real,&
              eh%occ,dm,ham)
      call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
   case(OMM_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(eh%ph%n_calls <= eh%ph%omm_n_elpa) then
         if(eh%ph%n_calls == 1 .and. eh%ph%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_real_copy",caller)
            eh%ovlp_real_copy = ovlp
         endif
         if(.not. allocated(eh%eval)) then
            call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
         endif
         if(.not. allocated(eh%evec_real)) then
            call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "evec_real",caller)
         endif
         if(.not. allocated(eh%occ)) then
            call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                    eh%ph%n_kpts,"occ",caller)
         endif

         call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,ham,ovlp,&
                 eh%eval,eh%evec_real)
         call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
         call elsi_compute_dm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
                 eh%evec_real,eh%occ,dm,ham)
         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         call elsi_init_omm(eh%ph,eh%bh)

         if(allocated(eh%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = eh%ovlp_real_copy
            call elsi_deallocate(eh%bh,eh%ovlp_real_copy,"ovlp_real_copy")
         endif
         if(.not. allocated(eh%omm_c_real)) then
            call elsi_allocate(eh%bh,eh%omm_c_real,eh%ph%omm_n_lrow,&
                    eh%bh%n_lcol,"omm_c_real",caller)
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(eh%ph%omm_n_elpa > 0 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
            call pdtran(eh%ph%n_basis,eh%ph%n_basis,1.0_r8,eh%evec_real,1,1,&
                    eh%bh%desc,0.0_r8,dm,1,1,eh%bh%desc)

            eh%omm_c_real(1:eh%ph%omm_n_lrow,:) = dm(1:eh%ph%omm_n_lrow,:)

            if(allocated(eh%evec_real)) then
               call elsi_deallocate(eh%bh,eh%evec_real,"evec_real")
            endif
            if(allocated(eh%eval)) then
               call elsi_deallocate(eh%bh,eh%eval,"eval")
            endif
            if(allocated(eh%occ)) then
               call elsi_deallocate(eh%bh,eh%occ,"occ")
            endif
         endif

         call elsi_solve_omm(eh%ph,eh%bh,ham,ovlp,eh%omm_c_real,dm)
         call elsi_get_energy(eh%ph,eh%bh,energy,OMM_SOLVER)
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)

      if(.not. allocated(eh%row_ind_sp1)) then
         call elsi_blacs_to_pexsi_hs_dim(eh%ph,eh%bh,ham,ovlp)

         if(eh%ph%ovlp_is_unit) then
            call elsi_allocate(eh%bh,eh%ovlp_real_csc,1,"ovlp_real_csc",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_csc,eh%bh%nnz_l_sp,&
                    "ovlp_real_csc",caller)
         endif
         call elsi_allocate(eh%bh,eh%ham_real_csc,eh%bh%nnz_l_sp,&
                 "ham_real_csc",caller)
         call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp,"row_ind_sp1",&
                 caller)
         call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp+1,&
                 "col_ptr_sp1",caller)
      endif

      call elsi_blacs_to_pexsi_hs(eh%ph,eh%bh,ham,ovlp,eh%ham_real_csc,&
              eh%ovlp_real_csc,eh%row_ind_sp1,eh%col_ptr_sp1)

      if(.not. allocated(eh%pexsi_ne_vec)) then
         call elsi_allocate(eh%bh,eh%pexsi_ne_vec,eh%ph%pexsi_options%nPoints,&
                 "pexsi_ne_vec",caller)
      endif
      if(.not. allocated(eh%dm_real_csc)) then
         call elsi_allocate(eh%bh,eh%dm_real_csc,eh%bh%nnz_l_sp,"dm_real_csc",&
                 caller)
      endif
      eh%dm_real_csc = 0.0_r8

      call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%pexsi_ne_vec,eh%ham_real_csc,eh%ovlp_real_csc,eh%dm_real_csc)
      call elsi_pexsi_to_blacs_dm(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%dm_real_csc,dm)
      call elsi_get_energy(eh%ph,eh%bh,energy,PEXSI_SOLVER)
   case(SIPS_SOLVER)
      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      endif
      if(.not. allocated(eh%evec_real)) then
         call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "evec_real",caller)
      endif

      if(eh%ph%n_calls <= eh%ph%sips_n_elpa) then
         call elsi_init_elpa(eh%ph,eh%bh)

         if(eh%ph%n_calls == 1) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_real_copy",caller)
            eh%ovlp_real_copy = ovlp
         endif
         if(.not. allocated(eh%occ)) then
            call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                    eh%ph%n_kpts,"occ",caller)
         endif

         call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,ham,ovlp,&
                 eh%eval,eh%evec_real)
         call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
         call elsi_compute_dm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
                 eh%evec_real,eh%occ,dm,ham)
         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         if(allocated(eh%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = eh%ovlp_real_copy
            call elsi_deallocate(eh%bh,eh%ovlp_real_copy,"ovlp_real_copy")
            call elsi_deallocate(eh%bh,eh%evec_real,"evec_real")
         endif

         call elsi_init_sips(eh%ph,eh%bh)

         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_blacs_to_sips_hs_dim(eh%ph,eh%bh,ham,ovlp)

            if(eh%ph%ovlp_is_unit) then
               call elsi_allocate(eh%bh,eh%ovlp_real_csc,1,"ovlp_real_csc",&
                       caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_real_csc,eh%bh%nnz_l_sp1,&
                       "ovlp_real_csc",caller)
            endif
            call elsi_allocate(eh%bh,eh%ham_real_csc,eh%bh%nnz_l_sp1,&
                    "ham_real_csc",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                    "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                    "col_ptr_sp1",caller)
         endif

         call elsi_blacs_to_sips_hs(eh%ph,eh%bh,ham,ovlp,eh%ham_real_csc,&
                 eh%ovlp_real_csc,eh%row_ind_sp1,eh%col_ptr_sp1)

         if(.not. allocated(eh%evec_real)) then
            call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lcol_sp1,&
                    eh%ph%n_states,"evec_real",caller)
         endif
         if(.not. allocated(eh%occ)) then
            call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                    eh%ph%n_kpts,"occ",caller)
         endif
         if(.not. allocated(eh%dm_real_csc)) then
            call elsi_allocate(eh%bh,eh%dm_real_csc,eh%bh%nnz_l_sp,&
                    "dm_real_csc",caller)
         endif
         eh%dm_real_csc = 0.0_r8

         call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%ham_real_csc,eh%ovlp_real_csc,eh%eval,eh%evec_real)
         call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
         call elsi_compute_dm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%occ,eh%dm_real_csc)
         call elsi_sips_to_blacs_dm(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%dm_real_csc,dm)
         call elsi_get_energy(eh%ph,eh%bh,energy,SIPS_SOLVER)
      endif
   case(NTPOLY_SOLVER)
      call elsi_init_ntpoly(eh%ph,eh%bh)
      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,ham,ovlp,eh%ph%nt_ham,&
              eh%ph%nt_ovlp)
      call elsi_solve_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_ovlp,eh%ph%nt_dm)
      call elsi_ntpoly_to_blacs_dm(eh%bh,eh%ph%nt_dm,dm)
      call elsi_get_energy(eh%ph,eh%bh,energy,NTPOLY_SOLVER)
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

   real(kind=r8)     :: t0
   character(len=29) :: dt0

   character(len=*), parameter :: caller = "elsi_dm_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      endif
      if(.not. allocated(eh%evec_cmplx)) then
         call elsi_allocate(eh%bh,eh%evec_cmplx,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "evec_cmplx",caller)
      endif
      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                 eh%ph%n_kpts,"occ",caller)
      endif

      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,ham,ovlp,eh%eval,&
              eh%evec_cmplx)
      call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
      call elsi_compute_dm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
              eh%evec_cmplx,eh%occ,dm,ham)
      call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
   case(OMM_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(eh%ph%n_calls <= eh%ph%omm_n_elpa) then
         if(eh%ph%n_calls == 1 .and. eh%ph%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_copy,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_cmplx_copy",caller)
            eh%ovlp_cmplx_copy = ovlp
         endif
         if(.not. allocated(eh%eval)) then
            call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
         endif
         if(.not. allocated(eh%evec_cmplx)) then
            call elsi_allocate(eh%bh,eh%evec_cmplx,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "evec_cmplx",caller)
         endif
         if(.not. allocated(eh%occ)) then
            call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                    eh%ph%n_kpts,"occ",caller)
         endif

         call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,ham,ovlp,&
                 eh%eval,eh%evec_cmplx)
         call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
         call elsi_compute_dm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
                 eh%evec_cmplx,eh%occ,dm,ham)
         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         call elsi_init_omm(eh%ph,eh%bh)

         if(allocated(eh%ovlp_cmplx_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            ovlp = eh%ovlp_cmplx_copy
            call elsi_deallocate(eh%bh,eh%ovlp_cmplx_copy,"ovlp_cmplx_copy")
         endif
         if(.not. allocated(eh%omm_c_cmplx)) then
            call elsi_allocate(eh%bh,eh%omm_c_cmplx,eh%ph%omm_n_lrow,&
                    eh%bh%n_lcol,"omm_c_cmplx",caller)
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(eh%ph%omm_n_elpa > 0 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
            call pztranc(eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),&
                    eh%evec_cmplx,1,1,eh%bh%desc,(0.0_r8,0.0_r8),dm,1,1,&
                    eh%bh%desc)

            eh%omm_c_cmplx(1:eh%ph%omm_n_lrow,:) = dm(1:eh%ph%omm_n_lrow,:)

            if(allocated(eh%evec_cmplx)) then
               call elsi_deallocate(eh%bh,eh%evec_cmplx,"evec_cmplx")
            endif
            if(allocated(eh%eval)) then
               call elsi_deallocate(eh%bh,eh%eval,"eval")
            endif
            if(allocated(eh%occ)) then
               call elsi_deallocate(eh%bh,eh%occ,"occ")
            endif
         endif

         call elsi_solve_omm(eh%ph,eh%bh,ham,ovlp,eh%omm_c_cmplx,dm)
         call elsi_get_energy(eh%ph,eh%bh,energy,OMM_SOLVER)
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)

      if(.not. allocated(eh%row_ind_sp1)) then
         call elsi_blacs_to_pexsi_hs_dim(eh%ph,eh%bh,ham,ovlp)

         if(eh%ph%ovlp_is_unit) then
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_csc,1,"ovlp_cmplx_csc",&
                    caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_csc,eh%bh%nnz_l_sp,&
                    "ovlp_cmplx_csc",caller)
         endif
         call elsi_allocate(eh%bh,eh%ham_cmplx_csc,eh%bh%nnz_l_sp,&
                 "ham_cmplx_csc",caller)
         call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp,"row_ind_sp1",&
                 caller)
         call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp+1,&
                 "col_ptr_sp1",caller)
      endif

      call elsi_blacs_to_pexsi_hs(eh%ph,eh%bh,ham,ovlp,eh%ham_cmplx_csc,&
              eh%ovlp_cmplx_csc,eh%row_ind_sp1,eh%col_ptr_sp1)

      if(.not. allocated(eh%pexsi_ne_vec)) then
         call elsi_allocate(eh%bh,eh%pexsi_ne_vec,eh%ph%pexsi_options%nPoints,&
                 "pexsi_ne_vec",caller)
      endif
      if(.not. allocated(eh%dm_cmplx_csc)) then
         call elsi_allocate(eh%bh,eh%dm_cmplx_csc,eh%bh%nnz_l_sp,&
                 "dm_cmplx_csc",caller)
      endif
      eh%dm_cmplx_csc = (0.0_r8,0.0_r8)

      call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%pexsi_ne_vec,eh%ham_cmplx_csc,eh%ovlp_cmplx_csc,&
              eh%dm_cmplx_csc)
      call elsi_pexsi_to_blacs_dm(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%dm_cmplx_csc,dm)
      call elsi_get_energy(eh%ph,eh%bh,energy,PEXSI_SOLVER)
   case(NTPOLY_SOLVER)
      call elsi_init_ntpoly(eh%ph,eh%bh)
      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,ham,ovlp,eh%ph%nt_ham,&
              eh%ph%nt_ovlp)
      call elsi_solve_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_ovlp,eh%ph%nt_dm)
      call elsi_ntpoly_to_blacs_dm(eh%bh,eh%ph%nt_dm,dm)
      call elsi_get_energy(eh%ph,eh%bh,energy,NTPOLY_SOLVER)
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

   real(kind=r8)     :: t0
   character(len=29) :: dt0

   character(len=*), parameter :: caller = "elsi_dm_real_sparse"

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

      if(.not. allocated(eh%ham_real_den)) then
         call elsi_allocate(eh%bh,eh%ham_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "ham_real_den",caller)
      endif
      if(.not. allocated(eh%ovlp_real_den)) then
         if(.not. eh%ph%ovlp_is_unit) then
            call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_real_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_den,1,1,"ovlp_real_den",&
                    caller)
         endif
      endif

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 ham,ovlp,eh%ham_real_den,eh%ovlp_real_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp2,&
                 eh%col_ptr_sp2,ham,ovlp,eh%ham_real_den,eh%ovlp_real_den)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      endif
      if(.not. allocated(eh%evec_real)) then
         call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "evec_real",caller)
      endif
      if(.not. allocated(eh%dm_real_den)) then
         call elsi_allocate(eh%bh,eh%dm_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "dm_real_den",caller)
      endif
      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                 eh%ph%n_kpts,"occ",caller)
      endif

      call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,eh%ham_real_den,&
              eh%ovlp_real_den,eh%eval,eh%evec_real)
      call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
      call elsi_compute_dm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,eh%evec_real,&
              eh%occ,eh%dm_real_den,eh%ham_real_den)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%dm_real_den,dm)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                 eh%dm_real_den,dm)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
   case(OMM_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. eh%bh%blacs_ready) then
         call elsi_init_blacs(eh)
      endif

      if(.not. allocated(eh%ham_real_den)) then
         call elsi_allocate(eh%bh,eh%ham_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "ham_real_den",caller)
      endif
      if(.not. allocated(eh%ovlp_real_den)) then
         if(.not. eh%ph%ovlp_is_unit) then
            call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_real_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_den,1,1,"ovlp_real_den",&
                    caller)
         endif
      endif

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 ham,ovlp,eh%ham_real_den,eh%ovlp_real_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp2,&
                 eh%col_ptr_sp2,ham,ovlp,eh%ham_real_den,eh%ovlp_real_den)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_init_elpa(eh%ph,eh%bh)

      if(eh%ph%n_calls <= eh%ph%omm_n_elpa) then
         if(eh%ph%n_calls == 1 .and. eh%ph%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_real_copy",caller)
            eh%ovlp_real_copy = eh%ovlp_real_den
         endif
         if(.not. allocated(eh%eval)) then
            call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
         endif
         if(.not. allocated(eh%evec_real)) then
            call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "evec_real",caller)
         endif
         if(.not. allocated(eh%dm_real_den)) then
            call elsi_allocate(eh%bh,eh%dm_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "dm_real_den",caller)
         endif
         if(.not. allocated(eh%occ)) then
            call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                    eh%ph%n_kpts,"occ",caller)
         endif

         call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
                 eh%ham_real_den,eh%ovlp_real_den,eh%eval,eh%evec_real)
         call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
         call elsi_compute_dm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
                 eh%evec_real,eh%occ,eh%dm_real_den,eh%ham_real_den)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_real_den,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                    eh%dm_real_den,dm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         call elsi_init_omm(eh%ph,eh%bh)

         if(allocated(eh%ovlp_real_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            eh%ovlp_real_den = eh%ovlp_real_copy
            call elsi_deallocate(eh%bh,eh%ovlp_real_copy,"ovlp_real_copy")
         endif
         if(.not. allocated(eh%dm_real_den)) then
            call elsi_allocate(eh%bh,eh%dm_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "dm_real_den",caller)
         endif
         if(.not. allocated(eh%omm_c_real)) then
            call elsi_allocate(eh%bh,eh%omm_c_real,eh%ph%omm_n_lrow,&
                    eh%bh%n_lcol,"omm_c_real",caller)
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(eh%ph%omm_n_elpa > 0 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
            call pdtran(eh%ph%n_basis,eh%ph%n_basis,1.0_r8,eh%evec_real,1,1,&
                    eh%bh%desc,0.0_r8,eh%dm_real_den,1,1,eh%bh%desc)

            eh%omm_c_real(1:eh%ph%omm_n_lrow,:) =&
               eh%dm_real_den(1:eh%ph%omm_n_lrow,:)

            if(allocated(eh%evec_real)) then
               call elsi_deallocate(eh%bh,eh%evec_real,"evec_real")
            endif
            if(allocated(eh%eval)) then
               call elsi_deallocate(eh%bh,eh%eval,"eval")
            endif
            if(allocated(eh%occ)) then
               call elsi_deallocate(eh%bh,eh%occ,"occ")
            endif
         endif

         call elsi_solve_omm(eh%ph,eh%bh,eh%ham_real_den,eh%ovlp_real_den,&
                 eh%omm_c_real,eh%dm_real_den)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_real_den,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                    eh%dm_real_den,dm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(eh%ph,eh%bh,energy,OMM_SOLVER)
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)

      if(.not. allocated(eh%pexsi_ne_vec)) then
         call elsi_allocate(eh%bh,eh%pexsi_ne_vec,eh%ph%pexsi_options%nPoints,&
                 "pexsi_ne_vec",caller)
      endif

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%pexsi_ne_vec,ham,ovlp,dm)
      case(SIESTA_CSC)
         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_siesta_to_pexsi_hs_dim(eh%ph,eh%bh,eh%col_ptr_sp2)

            if(eh%ph%ovlp_is_unit) then
               call elsi_allocate(eh%bh,eh%ovlp_real_csc,1,"ovlp_real_csc",&
                       caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_real_csc,eh%bh%nnz_l_sp1,&
                       "ovlp_real_csc",caller)
            endif
            call elsi_allocate(eh%bh,eh%ham_real_csc,eh%bh%nnz_l_sp1,&
                    "ham_real_csc",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                    "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                    "col_ptr_sp1",caller)
         endif

         call elsi_siesta_to_pexsi_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
                 eh%col_ptr_sp2,eh%ham_real_csc,eh%ovlp_real_csc,&
                 eh%row_ind_sp1,eh%col_ptr_sp1)

         if(.not. allocated(eh%dm_real_csc)) then
            call elsi_allocate(eh%bh,eh%dm_real_csc,eh%bh%nnz_l_sp1,&
                    "dm_real_csc",caller)
         endif
         eh%dm_real_csc = 0.0_r8

         call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%pexsi_ne_vec,eh%ham_real_csc,eh%ovlp_real_csc,&
                 eh%dm_real_csc)
         call elsi_pexsi_to_siesta_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                 eh%col_ptr_sp1,eh%dm_real_csc,eh%row_ind_sp2,eh%col_ptr_sp2,dm)
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

         if(.not. allocated(eh%ham_real_den)) then
            call elsi_allocate(eh%bh,eh%ham_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "ham_real_den",caller)
         endif
         if(.not. allocated(eh%ovlp_real_den)) then
            if(.not. eh%ph%ovlp_is_unit) then
               call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,&
                       eh%bh%n_lcol,"ovlp_real_den",caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_real_den,1,1,"ovlp_real_den",&
                       caller)
            endif
         endif

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_sips_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,ham,ovlp,eh%ham_real_den,eh%ovlp_real_den)
         case(SIESTA_CSC)
            call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp2,&
                    eh%col_ptr_sp2,ham,ovlp,eh%ham_real_den,eh%ovlp_real_den)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         if(.not. allocated(eh%eval)) then
            call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
         endif
         if(.not. allocated(eh%evec_real)) then
            call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "evec_real",caller)
         endif
         if(.not. allocated(eh%dm_real_den)) then
            call elsi_allocate(eh%bh,eh%dm_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "dm_real_den",caller)
         endif
         if(.not. allocated(eh%occ)) then
            call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                    eh%ph%n_kpts,"occ",caller)
         endif

         call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
                 eh%ham_real_den,eh%ovlp_real_den,eh%eval,eh%evec_real)
         call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
         call elsi_compute_dm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
                 eh%evec_real,eh%occ,eh%dm_real_den,eh%ham_real_den)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_real_den,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                    eh%dm_real_den,dm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         call elsi_init_sips(eh%ph,eh%bh)

         if(allocated(eh%ham_real_den)) then
            call elsi_deallocate(eh%bh,eh%ham_real_den,"ham_real_den")
         endif
         if(allocated(eh%ovlp_real_den)) then
            call elsi_deallocate(eh%bh,eh%ovlp_real_den,"ovlp_real_den")
         endif
         if(allocated(eh%dm_real_den)) then
            call elsi_deallocate(eh%bh,eh%dm_real_den,"dm_real_den")
         endif
         if(.not. allocated(eh%occ)) then
            call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                    eh%ph%n_kpts,"occ",caller)
         endif
         if(.not. allocated(eh%eval)) then
            call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
         endif
         if(.not. allocated(eh%evec_real)) then
            call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lcol_sp1,&
                    eh%ph%n_states,"evec_real",caller)
         endif

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,ham,&
                    ovlp,eh%eval,eh%evec_real)
            call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
            call elsi_compute_dm_sips(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%occ,dm)
            call elsi_get_energy(eh%ph,eh%bh,energy,SIPS_SOLVER)
         case(SIESTA_CSC)
            if(.not. allocated(eh%row_ind_sp1)) then
               call elsi_siesta_to_sips_hs_dim(eh%ph,eh%bh,eh%col_ptr_sp2)

               if(eh%ph%ovlp_is_unit) then
                  call elsi_allocate(eh%bh,eh%ovlp_real_csc,1,"ovlp_real_csc",&
                          caller)
               else
                  call elsi_allocate(eh%bh,eh%ovlp_real_csc,eh%bh%nnz_l_sp1,&
                          "ovlp_real_csc",caller)
               endif
               call elsi_allocate(eh%bh,eh%ham_real_csc,eh%bh%nnz_l_sp1,&
                       "ham_real_csc",caller)
               call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                       "row_ind_sp1",caller)
               call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                       "col_ptr_sp1",caller)
            endif

            call elsi_siesta_to_sips_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
                    eh%col_ptr_sp2,eh%ham_real_csc,eh%ovlp_real_csc,&
                    eh%row_ind_sp1,eh%col_ptr_sp1)

            if(.not. allocated(eh%dm_real_csc)) then
               call elsi_allocate(eh%bh,eh%dm_real_csc,eh%bh%nnz_l_sp1,&
                       "dm_real_csc",caller)
            endif
            eh%dm_real_csc = 0.0_r8

            call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                    eh%ham_real_csc,eh%ovlp_real_csc,eh%eval,eh%evec_real)
            call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
            call elsi_compute_dm_sips(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%occ,eh%dm_real_csc)
            call elsi_sips_to_siesta_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_real_csc,eh%row_ind_sp2,&
                    eh%col_ptr_sp2,dm)
            call elsi_get_energy(eh%ph,eh%bh,energy,SIPS_SOLVER)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select
      endif
   case(NTPOLY_SOLVER)
      call elsi_init_ntpoly(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_ntpoly_hs(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 ham,ovlp,eh%ph%nt_ham,eh%ph%nt_ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_ntpoly_hs(eh%ph,eh%bh,eh%row_ind_sp2,&
                 eh%col_ptr_sp2,ham,ovlp,eh%ph%nt_ham,eh%ph%nt_ovlp)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_solve_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_ovlp,eh%ph%nt_dm)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%ph%nt_dm,dm)
      case(SIESTA_CSC)
         call elsi_ntpoly_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                 eh%ph%nt_dm,dm)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(eh%ph,eh%bh,energy,NTPOLY_SOLVER)
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

   real(kind=r8)     :: t0
   character(len=29) :: dt0

   character(len=*), parameter :: caller = "elsi_dm_complex_sparse"

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

      if(.not. allocated(eh%ham_cmplx_den)) then
         call elsi_allocate(eh%bh,eh%ham_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "ham_cmplx_den",caller)
      endif
      if(.not. allocated(eh%ovlp_cmplx_den)) then
         if(.not. eh%ph%ovlp_is_unit) then
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_cmplx_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,1,1,"ovlp_cmplx_den",&
                    caller)
         endif
      endif

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 ham,ovlp,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp2,&
                 eh%col_ptr_sp2,ham,ovlp,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      endif
      if(.not. allocated(eh%evec_cmplx)) then
         call elsi_allocate(eh%bh,eh%evec_cmplx,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "evec_cmplx",caller)
      endif
      if(.not. allocated(eh%dm_cmplx_den)) then
         call elsi_allocate(eh%bh,eh%dm_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "dm_cmplx_den",caller)
      endif
      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                 eh%ph%n_kpts,"occ",caller)
      endif

      call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,eh%ham_cmplx_den,&
              eh%ovlp_cmplx_den,eh%eval,eh%evec_cmplx)
      call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
      call elsi_compute_dm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
              eh%evec_cmplx,eh%occ,eh%dm_cmplx_den,eh%ham_cmplx_den)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%dm_cmplx_den,dm)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                 eh%dm_cmplx_den,dm)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
   case(OMM_SOLVER)
      ! Set up BLACS if not done by user
      if(.not. eh%bh%blacs_ready) then
         call elsi_init_blacs(eh)
      endif

      if(.not. allocated(eh%ham_cmplx_den)) then
         call elsi_allocate(eh%bh,eh%ham_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "ham_cmplx_den",caller)
      endif
      if(.not. allocated(eh%ovlp_cmplx_den)) then
         if(.not. eh%ph%ovlp_is_unit) then
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_cmplx_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,1,1,"ovlp_cmplx_den",&
                    caller)
         endif
      endif

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 ham,ovlp,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,eh%row_ind_sp2,&
                 eh%col_ptr_sp2,ham,ovlp,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_init_elpa(eh%ph,eh%bh)

      if(eh%ph%n_calls <= eh%ph%omm_n_elpa) then
         if(eh%ph%n_calls == 1 .and. eh%ph%omm_flavor == 0) then
            ! Overlap will be destroyed by Cholesky
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_copy,eh%bh%n_lrow,&
                    eh%bh%n_lcol,"ovlp_cmplx_copy",caller)
            eh%ovlp_cmplx_copy = eh%ovlp_cmplx_den
         endif
         if(.not. allocated(eh%eval)) then
            call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
         endif
         if(.not. allocated(eh%evec_cmplx)) then
            call elsi_allocate(eh%bh,eh%evec_cmplx,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "evec_cmplx",caller)
         endif
         if(.not. allocated(eh%dm_cmplx_den)) then
            call elsi_allocate(eh%bh,eh%dm_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "dm_cmplx_den",caller)
         endif
         if(.not. allocated(eh%occ)) then
            call elsi_allocate(eh%bh,eh%occ,eh%ph%n_basis,eh%ph%n_spins,&
                    eh%ph%n_kpts,"occ",caller)
         endif

         call elsi_solve_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
                 eh%ham_cmplx_den,eh%ovlp_cmplx_den,eh%eval,eh%evec_cmplx)
         call elsi_compute_occ_elpa(eh%ph,eh%bh,eh%eval,eh%occ)
         call elsi_compute_dm_elpa(eh%ph,eh%bh,eh%row_map,eh%col_map,&
                 eh%evec_cmplx,eh%occ,eh%dm_cmplx_den,eh%ham_cmplx_den)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_cmplx_den,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                    eh%dm_cmplx_den,dm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(eh%ph,eh%bh,energy,ELPA_SOLVER)
      else ! ELPA is done
         call elsi_init_omm(eh%ph,eh%bh)

         if(allocated(eh%ovlp_cmplx_copy)) then
            ! Retrieve overlap matrix that has been destroyed by Cholesky
            eh%ovlp_cmplx_den = eh%ovlp_cmplx_copy
            call elsi_deallocate(eh%bh,eh%ovlp_cmplx_copy,"ovlp_cmplx_copy")
         endif
         if(.not. allocated(eh%dm_cmplx_den)) then
            call elsi_allocate(eh%bh,eh%dm_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                    "dm_cmplx_den",caller)
         endif
         if(.not. allocated(eh%omm_c_cmplx)) then
            call elsi_allocate(eh%bh,eh%omm_c_cmplx,eh%ph%omm_n_lrow,&
                    eh%bh%n_lcol,"omm_c_cmplx",caller)
         endif

         ! Initialize coefficient matrix with ELPA eigenvectors if possible
         if(eh%ph%omm_n_elpa > 0 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
            call pztranc(eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),&
                    eh%evec_cmplx,1,1,eh%bh%desc,(0.0_r8,0.0_r8),&
                    eh%dm_cmplx_den,1,1,eh%bh%desc)

            eh%omm_c_cmplx(1:eh%ph%omm_n_lrow,:) =&
               eh%dm_cmplx_den(1:eh%ph%omm_n_lrow,:)

            if(allocated(eh%evec_cmplx)) then
               call elsi_deallocate(eh%bh,eh%evec_cmplx,"evec_cmplx")
            endif
            if(allocated(eh%eval)) then
               call elsi_deallocate(eh%bh,eh%eval,"eval")
            endif
            if(allocated(eh%occ)) then
               call elsi_deallocate(eh%bh,eh%occ,"occ")
            endif
         endif

         call elsi_solve_omm(eh%ph,eh%bh,eh%ham_cmplx_den,eh%ovlp_cmplx_den,&
                 eh%omm_c_cmplx,eh%dm_cmplx_den)

         select case(eh%ph%matrix_format)
         case(PEXSI_CSC)
            call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                    eh%col_ptr_sp1,eh%dm_cmplx_den,dm)
         case(SIESTA_CSC)
            call elsi_blacs_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                    eh%dm_cmplx_den,dm)
         case default
            call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
         end select

         call elsi_get_energy(eh%ph,eh%bh,energy,OMM_SOLVER)
      endif
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)

      if(.not. allocated(eh%pexsi_ne_vec)) then
         call elsi_allocate(eh%bh,eh%pexsi_ne_vec,eh%ph%pexsi_options%nPoints,&
                 "pexsi_ne_vec",caller)
      endif

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%pexsi_ne_vec,ham,ovlp,dm)
      case(SIESTA_CSC)
         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_siesta_to_pexsi_hs_dim(eh%ph,eh%bh,eh%col_ptr_sp2)

            if(eh%ph%ovlp_is_unit) then
               call elsi_allocate(eh%bh,eh%ovlp_cmplx_csc,1,"ovlp_cmplx_csc",&
                       caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_cmplx_csc,eh%bh%nnz_l_sp1,&
                       "ovlp_cmplx_csc",caller)
            endif
            call elsi_allocate(eh%bh,eh%ham_cmplx_csc,eh%bh%nnz_l_sp1,&
                    "ham_cmplx_csc",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                    "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                    "col_ptr_sp1",caller)
         endif

         call elsi_siesta_to_pexsi_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
                 eh%col_ptr_sp2,eh%ham_cmplx_csc,eh%ovlp_cmplx_csc,&
                 eh%row_ind_sp1,eh%col_ptr_sp1)

         if(.not. allocated(eh%dm_cmplx_csc)) then
            call elsi_allocate(eh%bh,eh%dm_cmplx_csc,eh%bh%nnz_l_sp1,&
                    "dm_cmplx_csc",caller)
         endif
         eh%dm_cmplx_csc = (0.0_r8,0.0_r8)

         call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%pexsi_ne_vec,eh%ham_cmplx_csc,eh%ovlp_cmplx_csc,&
                 eh%dm_cmplx_csc)
         call elsi_pexsi_to_siesta_dm(eh%ph,eh%bh,eh%row_ind_sp1,&
                 eh%col_ptr_sp1,eh%dm_cmplx_csc,eh%row_ind_sp2,eh%col_ptr_sp2,&
                 dm)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(eh%ph,eh%bh,energy,PEXSI_SOLVER)
   case(NTPOLY_SOLVER)
      call elsi_init_ntpoly(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_ntpoly_hs(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 ham,ovlp,eh%ph%nt_ham,eh%ph%nt_ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_ntpoly_hs(eh%ph,eh%bh,eh%row_ind_sp2,&
                 eh%col_ptr_sp2,ham,ovlp,eh%ph%nt_ham,eh%ph%nt_ovlp)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_solve_ntpoly(eh%ph,eh%bh,eh%ph%nt_ham,eh%ph%nt_ovlp,eh%ph%nt_dm)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
                 eh%ph%nt_dm,dm)
      case(SIESTA_CSC)
         call elsi_ntpoly_to_siesta_dm(eh%bh,eh%row_ind_sp2,eh%col_ptr_sp2,&
                 eh%ph%nt_dm,dm)
      case default
         call elsi_stop(eh%bh,"Unsupported matrix format.",caller)
      end select

      call elsi_get_energy(eh%ph,eh%bh,energy,NTPOLY_SOLVER)
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

   character(len=*), parameter :: caller = "elsi_init_blacs"

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
