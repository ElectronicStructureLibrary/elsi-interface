! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Interface to ELPA.
!!
module ELSI_ELPA

   use ELSI_CONSTANT, only: LT_MAT,UT_MAT,UNSET,ELPA_SOLVER,DECISION_WIP
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI, only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_integer4,&
       mpi_comm_self
   use ELSI_OUTPUT, only: elsi_say,elsi_get_time
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTIL, only: elsi_get_nnz,elsi_get_gid,elsi_set_full_mat
   use CHECK_PD, only: elpa_check_pd_real_double,elpa_check_pd_complex_double
   use ELPA1, only: elpa_print_times,elpa_solve_evp_real_1stage_double,&
       elpa_solve_evp_complex_1stage_double,elpa_cholesky_real_double,&
       elpa_cholesky_complex_double,elpa_invert_trm_real_double,&
       elpa_invert_trm_complex_double,elpa_mult_at_b_real_double,&
       elpa_mult_ah_b_complex_double,elpa_solve_tridi_double
   use ELPA2, only: elpa_solve_evp_real_2stage_double,&
       elpa_solve_evp_complex_2stage_double

   implicit none

   private

   public :: elsi_init_elpa
   public :: elsi_cleanup_elpa
   public :: elsi_solve_elpa
   public :: elsi_update_dm_elpa
   public :: elsi_factor_ovlp_elpa
   public :: elsi_reduce_evp_elpa
   public :: elsi_back_ev_elpa
   public :: elsi_elpa_cholesky
   public :: elsi_elpa_invert
   public :: elsi_elpa_tridiag

   interface elsi_solve_elpa
      module procedure elsi_solve_elpa_real
      module procedure elsi_solve_elpa_cmplx
   end interface

   interface elsi_update_dm_elpa
      module procedure elsi_update_dm_elpa_real
      module procedure elsi_update_dm_elpa_cmplx
   end interface

   interface elsi_check_ovlp_elpa
      module procedure elsi_check_ovlp_elpa_real
      module procedure elsi_check_ovlp_elpa_cmplx
   end interface

   interface elsi_factor_ovlp_elpa
      module procedure elsi_factor_ovlp_elpa_real
      module procedure elsi_factor_ovlp_elpa_cmplx
   end interface

   interface elsi_reduce_evp_elpa
      module procedure elsi_reduce_evp_elpa_real
      module procedure elsi_reduce_evp_elpa_cmplx
   end interface

   interface elsi_back_ev_elpa
      module procedure elsi_back_ev_elpa_real
      module procedure elsi_back_ev_elpa_cmplx
   end interface

   interface elsi_elpa_evec
      module procedure elsi_elpa_evec_real
      module procedure elsi_elpa_evec_cmplx
   end interface

   interface elsi_elpa_cholesky
      module procedure elsi_elpa_cholesky_real
      module procedure elsi_elpa_cholesky_cmplx
   end interface

   interface elsi_elpa_invert
      module procedure elsi_elpa_invert_real
      module procedure elsi_elpa_invert_cmplx
   end interface

   interface elsi_elpa_multiply
      module procedure elsi_elpa_multiply_real
      module procedure elsi_elpa_multiply_cmplx
   end interface

contains

!>
!! Initialize ELPA.
!!
subroutine elsi_init_elpa(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh

   integer(kind=i4) :: ierr
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_init_elpa"

   if(.not. ph%elpa_started) then
      call MPI_Comm_split(bh%comm,bh%my_pcol,bh%my_prow,ph%elpa_comm_row,ierr)

      call elsi_check_mpi(bh,"MPI_Comm_split",ierr,caller)

      call MPI_Comm_split(bh%comm,bh%my_prow,bh%my_pcol,ph%elpa_comm_col,ierr)

      call elsi_check_mpi(bh,"MPI_Comm_split",ierr,caller)

      ph%elpa_started = .true.

      if(ph%elpa_gpu) then
         ph%elpa_gpu = .false.

         write(msg,"(A)") "No ELPA GPU acceleration available"
         call elsi_say(bh,msg)
      end if

      if(ph%elpa_n_single > 0) then
         ph%elpa_n_single = 0

         write(msg,"(A)") "No single precision ELPA available"
         call elsi_say(bh,msg)
      end if

      if(ph%elpa_autotune) then
         ph%elpa_autotune = .false.

         write(msg,"(A)") "No ELPA auto-tuning available"
         call elsi_say(bh,msg)
      end if
   end if

end subroutine

!>
!! Cholesky factorize the overlap matrix in place.
!!
subroutine elsi_factor_ovlp_elpa_real(ph,bh,ovlp)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_factor_ovlp_elpa_real"

   call elsi_get_time(t0)

   ! S = U
   call elsi_elpa_cholesky(ph,bh,ovlp)

   ! S = U^(-1)
   call elsi_elpa_invert(ph,bh,ovlp)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished Cholesky decomposition"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Transform a generalized eigenproblem to standard form using Cholesky or eigen
!! decomposition of the overlap matrix.
!!
subroutine elsi_reduce_evp_elpa_real(ph,bh,ham,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_reduce_evp_elpa_real"

   call elsi_get_time(t0)

   ! H = U^(-T) H U^(-1)
   if(ph%ill_ovlp) then
      call pdgemm("N","N",ph%n_basis,ph%n_good,ph%n_basis,1.0_r8,ham,1,1,&
           bh%desc,ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,0.0_r8,evec,1,1,&
           bh%desc)

      call pdgemm("T","N",ph%n_good,ph%n_good,ph%n_basis,1.0_r8,ovlp,1,&
           ph%n_basis-ph%n_good+1,bh%desc,evec,1,1,bh%desc,0.0_r8,ham,1,1,&
           bh%desc)
   else
      call elsi_elpa_multiply(ph,bh,"U","L",ph%n_basis,ovlp,ham,evec)

      call pdtran(ph%n_basis,ph%n_basis,1.0_r8,evec,1,1,bh%desc,0.0_r8,ham,1,1,&
           bh%desc)

      evec = ham

      call elsi_elpa_multiply(ph,bh,"U","U",ph%n_basis,ovlp,evec,ham)

      call elsi_set_full_mat(ph,bh,UT_MAT,ham)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished transformation to standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   if(ph%decision_stage == DECISION_WIP) then
      ph%decision_data(ELPA_SOLVER) = t1-t0
   end if

end subroutine

!>
!! Check the singularity of overlap matrix by computing all its eigenvalues. If
!! S is singular, it is overwritten by its eigen decomposition on exit.
!!
subroutine elsi_check_ovlp_elpa_real(ph,bh,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: i
   integer(kind=i4) :: gid
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_check_ovlp_elpa_real"

   call elsi_get_time(t0)

   ! Solve eigenvalues of S
   call elsi_elpa_evec(ph,bh,ovlp,eval,evec,.true.)

   if(ph%n_good < ph%n_basis) then ! Singular
      ph%ill_ovlp = .true.

      write(msg,"(A)") "Overlap matrix is singular"
      call elsi_say(bh,msg)
      write(msg,"(A,E12.4,A,E12.4)") "| Lowest and highest eigenvalues :",&
         eval(1),",",eval(ph%n_basis)
      call elsi_say(bh,msg)
      write(msg,"(A,I10)") "| Number of basis functions reduced to :",ph%n_good
      call elsi_say(bh,msg)

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i,gid)

         ovlp(:,i) = evec(:,i)/sqrt(eval(gid))
      end do
   else
      ph%ill_ovlp = .false.

      write(msg,"(A)") "Overlap matrix is not singular"
      call elsi_say(bh,msg)
      write(msg,"(A,E12.4,A,E12.4)") "| Lowest and highest eigenvalues :",&
         eval(1),",",eval(ph%n_basis)
      call elsi_say(bh,msg)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished singularity check of overlap matrix"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Back-transform eigenvectors in the standard form to the original generalized
!! form.
!!
subroutine elsi_back_ev_elpa_real(ph,bh,ham,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(out) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   real(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_back_ev_elpa_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)

   tmp = evec

   if(ph%ill_ovlp) then
      call pdgemm("N","N",ph%n_basis,ph%n_states_solve,ph%n_good,1.0_r8,ovlp,1,&
           ph%n_basis-ph%n_good+1,bh%desc,tmp,1,1,bh%desc,0.0_r8,evec,1,1,&
           bh%desc)
   else
      call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp,1,1,bh%desc,0.0_r8,ham,1,1,&
           bh%desc)

      call elsi_elpa_multiply(ph,bh,"L","N",ph%n_states,ham,tmp,evec)
   end if

   call elsi_deallocate(bh,tmp,"tmp")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished back-transformation of eigenvectors"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   if(ph%decision_stage == DECISION_WIP) then
      ph%decision_data(ELPA_SOLVER) = ph%decision_data(ELPA_SOLVER)+t1-t0
   end if

end subroutine

!>
!! Interface to ELPA.
!!
subroutine elsi_solve_elpa_real(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_solve_elpa_real"

   elpa_print_times = ph%elpa_output

   ! Compute sparsity
   if(bh%nnz_g == UNSET) then
      if(bh%nnz_l == UNSET) then
         call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,ham,bh%nnz_l)
      end if

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   end if

   ! Ill-conditioning check
   if(.not. ph%unit_ovlp .and. ph%elpa_first .and. ph%ill_check) then
      call elsi_check_ovlp_elpa(ph,bh,ovlp,eval,evec)
   end if

   write(msg,"(A)") "Starting ELPA eigensolver"
   call elsi_say(bh,msg)

   ! Transform to standard form
   if(.not. ph%unit_ovlp) then
      if(ph%elpa_first .and. ph%n_good == ph%n_basis) then ! Not singular
         call elsi_factor_ovlp_elpa(ph,bh,ovlp)
      end if

      call elsi_reduce_evp_elpa(ph,bh,ham,ovlp,evec)
   end if

   call elsi_get_time(t0)

   ! Solve
   call elsi_elpa_evec(ph,bh,ham,eval,evec,.false.)

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(ph%n_good < ph%n_basis) then
      eval(ph%n_good+1:ph%n_basis) = eval(ph%n_good)+10.0_r8
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished solving standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   if(ph%decision_stage == DECISION_WIP) then
      ph%decision_data(ELPA_SOLVER) = ph%decision_data(ELPA_SOLVER)+t1-t0
   end if

   ! Back-transform eigenvectors
   if(.not. ph%unit_ovlp) then
      call elsi_back_ev_elpa(ph,bh,ham,ovlp,evec)
   end if

   ph%elpa_first = .false.

end subroutine

!>
!! Extrapolate density matrix using Cholesky decomposition of the old and new
!! overlap matrices.
!!
subroutine elsi_update_dm_elpa_real(ph,bh,ovlp0,ovlp1,dm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: ovlp0(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp1(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: dm(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   real(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_update_dm_elpa_real"

   call elsi_get_time(t0)

   if(ph%elpa_first) then
      call elsi_factor_ovlp_elpa(ph,bh,ovlp1)
   end if

   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)

   ! tmp = U_1^(-T)
   call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp1,1,1,bh%desc,0.0_r8,tmp,1,1,&
        bh%desc)

   ! ovlp0 = U_0
   call elsi_elpa_invert(ph,bh,ovlp0)

   ! ovlp1 = U_1^(-1) U_0
   call elsi_elpa_multiply(ph,bh,"L","U",ph%n_basis,tmp,ovlp0,ovlp1)

   ! ovlp0 = U_0^T U_1^(-T)
   call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp1,1,1,bh%desc,0.0_r8,ovlp0,1,1,&
        bh%desc)

   ! ovlp1 = U_1^(-1) U_0 P_0
   call elsi_elpa_multiply(ph,bh,"L","U",ph%n_basis,ovlp0,dm,ovlp1)

   ! dm = P_0 U_0^T U_1^(-T)
   call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp1,1,1,bh%desc,0.0_r8,dm,1,1,&
        bh%desc)

   ! ovlp1 = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
   call elsi_elpa_multiply(ph,bh,"L","L",ph%n_basis,ovlp0,dm,ovlp1)

   call elsi_set_full_mat(ph,bh,LT_MAT,ovlp1)

   ! dm = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
   dm = ovlp1

   ! ovlp1 = U_1^(-1)
   call pdtran(ph%n_basis,ph%n_basis,1.0_r8,tmp,1,1,bh%desc,0.0_r8,ovlp1,1,1,&
        bh%desc)

   call elsi_deallocate(bh,tmp,"tmp")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished density matrix extrapolation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%elpa_first = .false.

end subroutine

!>
!! Cholesky factorize the overlap matrix in place.
!!
subroutine elsi_factor_ovlp_elpa_cmplx(ph,bh,ovlp)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_factor_ovlp_elpa_cmplx"

   call elsi_get_time(t0)

   ! S = U
   call elsi_elpa_cholesky(ph,bh,ovlp)

   ! S = U^(-1)
   call elsi_elpa_invert(ph,bh,ovlp)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished Cholesky decomposition"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Transform a generalized eigenproblem to standard form using Cholesky or eigen
!! decomposition of the overlap matrix.
!!
subroutine elsi_reduce_evp_elpa_cmplx(ph,bh,ham,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_reduce_evp_elpa_cmplx"

   call elsi_get_time(t0)

   ! H = U^(-T) H U^(-1)
   if(ph%ill_ovlp) then
      call pzgemm("N","N",ph%n_basis,ph%n_good,ph%n_basis,(1.0_r8,0.0_r8),ham,&
           1,1,bh%desc,ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,(0.0_r8,0.0_r8),&
           evec,1,1,bh%desc)

      call pzgemm("C","N",ph%n_good,ph%n_good,ph%n_basis,(1.0_r8,0.0_r8),ovlp,&
           1,ph%n_basis-ph%n_good+1,bh%desc,evec,1,1,bh%desc,(0.0_r8,0.0_r8),&
           ham,1,1,bh%desc)
   else
      call elsi_elpa_multiply(ph,bh,"U","L",ph%n_basis,ovlp,ham,evec)

      call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),evec,1,1,bh%desc,&
           (0.0_r8,0.0_r8),ham,1,1,bh%desc)

      evec = ham

      call elsi_elpa_multiply(ph,bh,"U","U",ph%n_basis,ovlp,evec,ham)

      call elsi_set_full_mat(ph,bh,UT_MAT,ham)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished transformation to standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   if(ph%decision_stage == DECISION_WIP) then
      ph%decision_data(ELPA_SOLVER) = t1-t0
   end if

end subroutine

!>
!! Check the singularity of overlap matrix by computing all its eigenvalues. If
!! S is singular, it is overwritten by its eigen decomposition on exit.
!!
subroutine elsi_check_ovlp_elpa_cmplx(ph,bh,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: i
   integer(kind=i4) :: gid
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_check_ovlp_elpa_cmplx"

   call elsi_get_time(t0)

   ! Solve eigenvalues of S
   call elsi_elpa_evec(ph,bh,ovlp,eval,evec,.true.)

   if(ph%n_good < ph%n_basis) then ! Singular
      ph%ill_ovlp = .true.

      write(msg,"(A)") "Overlap matrix is singular"
      call elsi_say(bh,msg)
      write(msg,"(A,E12.4,A,E12.4)") "| Lowest and highest eigenvalues :",&
         eval(1),",",eval(ph%n_basis)
      call elsi_say(bh,msg)
      write(msg,"(A,I10)") "| Number of basis functions reduced to :",ph%n_good
      call elsi_say(bh,msg)

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i,gid)

         ovlp(:,i) = evec(:,i)/sqrt(eval(gid))
      end do
   else
      ph%ill_ovlp = .false.

      write(msg,"(A)") "Overlap matrix is not singular"
      call elsi_say(bh,msg)
      write(msg,"(A,E12.4,A,E12.4)") "| Lowest and highest eigenvalues :",&
         eval(1),",",eval(ph%n_basis)
      call elsi_say(bh,msg)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished singularity check of overlap matrix"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Back-transform eigenvectors in the standard form to the original generalized
!! form.
!!
subroutine elsi_back_ev_elpa_cmplx(ph,bh,ham,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(out) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   complex(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_back_ev_elpa_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)

   tmp = evec

   if(ph%ill_ovlp) then
      call pzgemm("N","N",ph%n_basis,ph%n_states_solve,ph%n_good,&
           (1.0_r8,0.0_r8),ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,tmp,1,1,&
           bh%desc,(0.0_r8,0.0_r8),evec,1,1,bh%desc)
   else
      call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),ovlp,1,1,bh%desc,&
           (0.0_r8,0.0_r8),ham,1,1,bh%desc)

      call elsi_elpa_multiply(ph,bh,"L","N",ph%n_states,ham,tmp,evec)
   end if

   call elsi_deallocate(bh,tmp,"tmp")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished back-transformation of eigenvectors"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   if(ph%decision_stage == DECISION_WIP) then
      ph%decision_data(ELPA_SOLVER) = ph%decision_data(ELPA_SOLVER)+t1-t0
   end if

end subroutine

!>
!! Interface to ELPA.
!!
subroutine elsi_solve_elpa_cmplx(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_solve_elpa_cmplx"

   elpa_print_times = ph%elpa_output

   ! Compute sparsity
   if(bh%nnz_g == UNSET) then
      if(bh%nnz_l == UNSET) then
         call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,ham,bh%nnz_l)
      end if

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   end if

   ! Ill-conditioning check
   if(.not. ph%unit_ovlp .and. ph%elpa_first .and. ph%ill_check) then
      call elsi_check_ovlp_elpa(ph,bh,ovlp,eval,evec)
   end if

   write(msg,"(A)") "Starting ELPA eigensolver"
   call elsi_say(bh,msg)

   ! Transform to standard form
   if(.not. ph%unit_ovlp) then
      if(ph%elpa_first .and. ph%n_good == ph%n_basis) then ! Not singular
         call elsi_factor_ovlp_elpa(ph,bh,ovlp)
      end if

      call elsi_reduce_evp_elpa(ph,bh,ham,ovlp,evec)
   end if

   call elsi_get_time(t0)

   ! Solve
   call elsi_elpa_evec(ph,bh,ham,eval,evec,.false.)

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(ph%n_good < ph%n_basis) then
      eval(ph%n_good+1:ph%n_basis) = eval(ph%n_good)+10.0_r8
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished solving standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   if(ph%decision_stage == DECISION_WIP) then
      ph%decision_data(ELPA_SOLVER) = ph%decision_data(ELPA_SOLVER)+t1-t0
   end if

   ! Back-transform eigenvectors
   if(.not. ph%unit_ovlp) then
      call elsi_back_ev_elpa(ph,bh,ham,ovlp,evec)
   end if

   ph%elpa_first = .false.

end subroutine

!>
!! Extrapolate density matrix using Cholesky decomposition of the old and new
!! overlap matrices.
!!
subroutine elsi_update_dm_elpa_cmplx(ph,bh,ovlp0,ovlp1,dm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: ovlp0(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp1(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: dm(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   complex(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_update_dm_elpa_cmplx"

   call elsi_get_time(t0)

   if(ph%elpa_first) then
      call elsi_factor_ovlp_elpa(ph,bh,ovlp1)
   end if

   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)

   ! tmp = U_1^(-T)
   call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),ovlp1,1,1,bh%desc,&
        (0.0_r8,0.0_r8),tmp,1,1,bh%desc)

   ! ovlp0 = U_0
   call elsi_elpa_invert(ph,bh,ovlp0)

   ! ovlp1 = U_1^(-1) U_0
   call elsi_elpa_multiply(ph,bh,"L","U",ph%n_basis,tmp,ovlp0,ovlp1)

   ! ovlp0 = U_0^T U_1^(-T)
   call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),ovlp1,1,1,bh%desc,&
        (0.0_r8,0.0_r8),ovlp0,1,1,bh%desc)

   ! ovlp1 = U_1^(-1) U_0 P_0
   call elsi_elpa_multiply(ph,bh,"L","U",ph%n_basis,ovlp0,dm,ovlp1)

   ! dm = P_0 U_0^T U_1^(-T)
   call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),ovlp1,1,1,bh%desc,&
        (0.0_r8,0.0_r8),dm,1,1,bh%desc)

   ! ovlp1 = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
   call elsi_elpa_multiply(ph,bh,"L","L",ph%n_basis,ovlp0,dm,ovlp1)

   call elsi_set_full_mat(ph,bh,LT_MAT,ovlp1)

   ! dm = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
   dm = ovlp1

   ! ovlp1 = U_1^(-1)
   call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),tmp,1,1,bh%desc,&
        (0.0_r8,0.0_r8),ovlp1,1,1,bh%desc)

   call elsi_deallocate(bh,tmp,"tmp")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished density matrix extrapolation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%elpa_first = .false.

end subroutine

!>
!! Interface to ELPA eigensolver.
!!
subroutine elsi_elpa_evec_real(ph,bh,mat,eval,evec,sing_check)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: mat(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)
   logical, intent(in) :: sing_check

   logical :: ok
   character(len=200) :: msg

   real(kind=r8), allocatable :: copy(:,:)

   character(len=*), parameter :: caller = "elsi_elpa_evec_real"

   if(sing_check) then
      call elsi_allocate(bh,copy,bh%n_lrow,bh%n_lcol,"copy",caller)

      copy = mat

      ! Use modified ELPA2, which computes eigenvectors only for singular matrix
      ok = elpa_check_pd_real_double(ph%n_basis,ph%n_basis,copy,bh%n_lrow,eval,&
         evec,bh%n_lrow,bh%blk,ph%elpa_comm_row,ph%elpa_comm_col,bh%comm,&
         ph%ill_tol,ph%n_good)

      if(.not. ok) then
         write(msg,"(A)") "Singularity check failed"
         call elsi_stop(bh,msg,caller)
      end if

      call elsi_deallocate(bh,copy,"copy")

      ph%n_states_solve = min(ph%n_good,ph%n_states)
      ph%ovlp_ev_min = eval(1)
      ph%ovlp_ev_max = eval(ph%n_basis)
   else
      if(ph%elpa_solver == 2) then
         ok = elpa_solve_evp_real_2stage_double(ph%n_good,ph%n_states_solve,&
            mat,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,ph%elpa_comm_row,&
            ph%elpa_comm_col,bh%comm)
      else
         ok = elpa_solve_evp_real_1stage_double(ph%n_good,ph%n_states_solve,&
            mat,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,ph%elpa_comm_row,&
            ph%elpa_comm_col,bh%comm)
      end if

      if(.not. ok) then
         write(msg,"(A)") "ELPA eigensolver failed"
         call elsi_stop(bh,msg,caller)
      end if
   end if

end subroutine

!>
!! Interface to ELPA eigensolver.
!!
subroutine elsi_elpa_evec_cmplx(ph,bh,mat,eval,evec,sing_check)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: mat(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)
   logical, intent(in) :: sing_check

   logical :: ok
   character(len=200) :: msg

   complex(kind=r8), allocatable :: copy(:,:)

   character(len=*), parameter :: caller = "elsi_elpa_evec_cmplx"

   if(sing_check) then
      call elsi_allocate(bh,copy,bh%n_lrow,bh%n_lcol,"copy",caller)

      copy = mat

      ! Use modified ELPA2, which computes eigenvectors only for singular matrix
      ok = elpa_check_pd_complex_double(ph%n_basis,ph%n_basis,copy,bh%n_lrow,&
         eval,evec,bh%n_lrow,bh%blk,ph%elpa_comm_row,ph%elpa_comm_col,bh%comm,&
         ph%ill_tol,ph%n_good)

      if(.not. ok) then
         write(msg,"(A)") "Singularity check failed"
         call elsi_stop(bh,msg,caller)
      end if

      call elsi_deallocate(bh,copy,"copy")

      ph%n_states_solve = min(ph%n_good,ph%n_states)
      ph%ovlp_ev_min = eval(1)
      ph%ovlp_ev_max = eval(ph%n_basis)
   else
      if(ph%elpa_solver == 2) then
         ok = elpa_solve_evp_complex_2stage_double(ph%n_good,ph%n_states_solve,&
            mat,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,ph%elpa_comm_row,&
            ph%elpa_comm_col,bh%comm)
      else
         ok = elpa_solve_evp_complex_1stage_double(ph%n_good,ph%n_states_solve,&
            mat,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,ph%elpa_comm_row,&
            ph%elpa_comm_col,bh%comm)
      end if

      if(.not. ok) then
         write(msg,"(A)") "ELPA eigensolver failed"
         call elsi_stop(bh,msg,caller)
      end if
   end if

end subroutine

!>
!! Interface to ELPA Cholesky decomposition.
!!
subroutine elsi_elpa_cholesky_real(ph,bh,mat)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: mat(bh%n_lrow,bh%n_lcol)

   logical :: ok
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_elpa_cholesky_real"

   ok = elpa_cholesky_real_double(ph%n_basis,mat,bh%n_lrow,bh%blk,&
      ph%elpa_comm_row,ph%elpa_comm_col)

   if(.not. ok) then
      write(msg,"(A)") "Cholesky factorization failed"
      call elsi_stop(bh,msg,caller)
   end if

end subroutine

!>
!! Interface to ELPA Cholesky decomposition.
!!
subroutine elsi_elpa_cholesky_cmplx(ph,bh,mat)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: mat(bh%n_lrow,bh%n_lcol)

   logical :: ok
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_elpa_cholesky_cmplx"

   ok = elpa_cholesky_complex_double(ph%n_basis,mat,bh%n_lrow,bh%blk,&
      ph%elpa_comm_row,ph%elpa_comm_col)

   if(.not. ok) then
      write(msg,"(A)") "Cholesky factorization failed"
      call elsi_stop(bh,msg,caller)
   end if

end subroutine

!>
!! Interface to ELPA matrix inversion.
!!
subroutine elsi_elpa_invert_real(ph,bh,mat)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: mat(bh%n_lrow,bh%n_lcol)

   logical :: ok
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_elpa_invert_real"

   ok = elpa_invert_trm_real_double(ph%n_basis,mat,bh%n_lrow,bh%blk,&
      ph%elpa_comm_row,ph%elpa_comm_col)

   if(.not. ok) then
      write(msg,"(A)") "Matrix inversion failed"
      call elsi_stop(bh,msg,caller)
   end if

end subroutine

!>
!! Interface to ELPA matrix inversion.
!!
subroutine elsi_elpa_invert_cmplx(ph,bh,mat)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: mat(bh%n_lrow,bh%n_lcol)

   logical :: ok
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_elpa_invert_cmplx"

   ok = elpa_invert_trm_complex_double(ph%n_basis,mat,bh%n_lrow,bh%blk,&
      ph%elpa_comm_row,ph%elpa_comm_col)

   if(.not. ok) then
      write(msg,"(A)") "Matrix inversion failed"
      call elsi_stop(bh,msg,caller)
   end if

end subroutine

!>
!! Interface to ELPA matrix multiplication.
!!
subroutine elsi_elpa_multiply_real(ph,bh,uplo_a,uplo_c,n,mat_a,mat_b,mat_c)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   character, intent(in) :: uplo_a
   character, intent(in) :: uplo_c
   integer(kind=i4), intent(in) :: n
   real(kind=r8), intent(in) :: mat_a(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: mat_b(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: mat_c(bh%n_lrow,bh%n_lcol)

   logical :: ok
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_elpa_multiply_real"

   ok = elpa_mult_at_b_real_double(uplo_a,uplo_c,ph%n_basis,n,mat_a,bh%n_lrow,&
      mat_b,bh%n_lrow,bh%blk,ph%elpa_comm_row,ph%elpa_comm_col,mat_c,bh%n_lrow)

   if(.not. ok) then
      write(msg,"(A)") "Matrix multiplication failed"
      call elsi_stop(bh,msg,caller)
   end if

end subroutine

!>
!! Interface to ELPA matrix multiplication.
!!
subroutine elsi_elpa_multiply_cmplx(ph,bh,uplo_a,uplo_c,n,mat_a,mat_b,mat_c)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   character, intent(in) :: uplo_a
   character, intent(in) :: uplo_c
   integer(kind=i4), intent(in) :: n
   complex(kind=r8), intent(in) :: mat_a(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: mat_b(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: mat_c(bh%n_lrow,bh%n_lcol)

   logical :: ok
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_elpa_multiply_cmplx"

   ok = elpa_mult_ah_b_complex_double(uplo_a,uplo_c,ph%n_basis,n,mat_a,&
      bh%n_lrow,mat_b,bh%n_lrow,bh%blk,ph%elpa_comm_row,ph%elpa_comm_col,mat_c,&
      bh%n_lrow)

   if(.not. ok) then
      write(msg,"(A)") "Matrix multiplication failed"
      call elsi_stop(bh,msg,caller)
   end if

end subroutine

!>
!! Interface to ELPA tridiagonal solver.
!!
subroutine elsi_elpa_tridiag(ph,bh,d,e,q,sing_check)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: d(ph%n_basis)
   real(kind=r8), intent(inout) :: e(ph%n_good)
   real(kind=r8), intent(inout) :: q(ph%n_good,ph%n_good)
   logical, intent(in) :: sing_check

   logical :: ok
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_elpa_tridiag"

   if(sing_check) then
      ok = elpa_solve_tridi_double(ph%n_basis,ph%n_basis,d,e,q,ph%n_basis,&
         bh%blk,mpi_comm_self,mpi_comm_self)
   else
      ok = elpa_solve_tridi_double(ph%n_good,ph%n_states_solve,d,e,q,ph%n_good,&
         bh%blk,mpi_comm_self,mpi_comm_self)
   end if

   if(.not. ok) then
      write(msg,"(A)") "ELPA tridiagonal solver failed"
      call elsi_stop(bh,msg,caller)
   end if

end subroutine

!>
!! Clean up ELPA.
!!
subroutine elsi_cleanup_elpa(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_cleanup_elpa"

   if(ph%elpa_started) then
      call MPI_Comm_free(ph%elpa_comm_row,ierr)
      call MPI_Comm_free(ph%elpa_comm_col,ierr)
   end if

   ph%elpa_first = .true.
   ph%elpa_started = .false.

end subroutine

end module ELSI_ELPA
