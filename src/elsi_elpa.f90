! Copyright (c) 2015-2021, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Interface to ELPA-AEO.
!!
module ELSI_ELPA

   use ELSI_CONSTANT, only: LT_MAT,UT_MAT,UNSET
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI
   use ELSI_OUTPUT, only: elsi_say,elsi_get_time
   use ELSI_PRECISION, only: r4,r8,i4
   use ELSI_UTIL, only: elsi_check_err,elsi_get_gid,elsi_set_full_mat
   use ELPA, only: elpa_init,elpa_allocate,elpa_deallocate,&
       elpa_autotune_deallocate,ELPA_2STAGE_REAL_GPU,ELPA_2STAGE_COMPLEX_GPU,&
       ELPA_AUTOTUNE_FAST,ELPA_AUTOTUNE_MEDIUM,ELPA_AUTOTUNE_DOMAIN_REAL,&
       ELPA_AUTOTUNE_DOMAIN_COMPLEX

   implicit none

   private

   public :: elsi_init_elpa
   public :: elsi_cleanup_elpa
   public :: elsi_solve_elpa
   public :: elsi_update_dm_elpa
   public :: elsi_factor_ovlp_elpa
   public :: elsi_reduce_evp_elpa
   public :: elsi_back_ev_elpa
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

contains

!>
!! Initialize ELPA.
!!
subroutine elsi_init_elpa(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh

   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_init_elpa"

   if(.not. ph%elpa_started) then
      ierr = elpa_init(20180525)

      call elsi_check_err(bh,"ELPA initialization failed",ierr,caller)

      call MPI_Comm_split(bh%comm,bh%my_pcol,bh%my_prow,ph%elpa_comm_row,ierr)

      call elsi_check_err(bh,"MPI_Comm_split",ierr,caller)

      call MPI_Comm_split(bh%comm,bh%my_prow,bh%my_pcol,ph%elpa_comm_col,ierr)

      call elsi_check_err(bh,"MPI_Comm_split",ierr,caller)

      call elsi_elpa_setup(ph,bh,.true.)

      ph%elpa_started = .true.
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
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_factor_ovlp_elpa_real"

   call elsi_get_time(t0)

   ! S = U
   call ph%elpa_aux%cholesky(ovlp,ierr)

   call elsi_check_err(bh,"ELPA Cholesky factorization",ierr,caller)

   ! S = U^(-1)
   call ph%elpa_aux%invert_triangular(ovlp,ierr)

   call elsi_check_err(bh,"ELPA matrix inversion",ierr,caller)

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
   integer(kind=i4) :: ierr
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
      call ph%elpa_aux%hermitian_multiply("U","L",ph%n_basis,ovlp,ham,&
           bh%n_lrow,bh%n_lcol,evec,bh%n_lrow,bh%n_lcol,ierr)

      call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

      call pdtran(ph%n_basis,ph%n_basis,1.0_r8,evec,1,1,bh%desc,0.0_r8,ham,1,1,&
           bh%desc)

      evec(:,:) = ham

      call ph%elpa_aux%hermitian_multiply("U","U",ph%n_basis,ovlp,evec,&
           bh%n_lrow,bh%n_lcol,ham,bh%n_lrow,bh%n_lcol,ierr)

      call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

      call elsi_set_full_mat(ph,bh,UT_MAT,ham)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished transformation to standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

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

   if(ph%n_good < ph%n_basis) then
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
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   real(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_back_ev_elpa_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)

   tmp(:,:) = evec

   if(ph%ill_ovlp) then
      call pdgemm("N","N",ph%n_basis,ph%n_states_solve,ph%n_good,1.0_r8,ovlp,1,&
           ph%n_basis-ph%n_good+1,bh%desc,tmp,1,1,bh%desc,0.0_r8,evec,1,1,&
           bh%desc)
   else
      call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp,1,1,bh%desc,0.0_r8,ham,1,1,&
           bh%desc)

      call ph%elpa_aux%hermitian_multiply("L","N",ph%n_states,ham,tmp,&
           bh%n_lrow,bh%n_lcol,evec,bh%n_lrow,bh%n_lcol,ierr)

      call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)
   end if

   call elsi_deallocate(bh,tmp,"tmp")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished back-transformation of eigenvectors"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

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

   ! Compute sparsity
   if(bh%nnz_g == UNSET) then
      if(bh%nnz_l == UNSET) then
         bh%nnz_l = count(abs(ham) > bh%def0)
      end if

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,MPI_INTEGER4,MPI_SUM,bh%comm,ierr)

      call elsi_check_err(bh,"MPI_Allreduce",ierr,caller)
   end if

   ! Ill-conditioning check
   if(.not. ph%unit_ovlp .and. ph%elpa_first .and. ph%ill_check) then
      call elsi_check_ovlp_elpa(ph,bh,ovlp,eval,evec)
   end if

   if(ph%elpa_gpu == 1) then
      if(ph%n_calls <= ph%elpa_n_single) then
         write(msg,"(A)") "Starting ELPA eigensolver (GPU single precision)"
      else
         write(msg,"(A)") "Starting ELPA eigensolver (GPU)"
      end if
   else
      if(ph%n_calls <= ph%elpa_n_single) then
         write(msg,"(A)") "Starting ELPA eigensolver (single precision)"
      else
         write(msg,"(A)") "Starting ELPA eigensolver"
      end if
   end if
   call elsi_say(bh,msg)

   ! Transform to standard form
   if(.not. ph%unit_ovlp) then
      if(ph%elpa_first .and. ph%n_good == ph%n_basis) then
         ! Do Cholesky if not singular
         call elsi_factor_ovlp_elpa(ph,bh,ovlp)
      end if

      call elsi_reduce_evp_elpa(ph,bh,ham,ovlp,evec)
   end if

   call elsi_get_time(t0)

   ! Solve
   if(.not. associated(ph%elpa_solve)) then
      call elsi_elpa_setup(ph,bh,.false.)
   end if

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
subroutine elsi_update_dm_elpa_real(ph,bh,ovlp0,ovlp1,dm0,dm1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: ovlp0(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp1(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: dm0(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: dm1(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_update_dm_elpa_real"

   call elsi_get_time(t0)

   ! ovlp0 = U_0
   call ph%elpa_aux%cholesky(ovlp0,ierr)

   call elsi_check_err(bh,"ELPA Cholesky factorization",ierr,caller)

   ! ovlp1 = (U_1)^(-1)
   call elsi_factor_ovlp_elpa(ph,bh,ovlp1)

   ! dm1 = U_1^(-T)
   call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp1,1,1,bh%desc,0.0_r8,dm1,1,1,&
        bh%desc)

   ! ovlp1 = U_1^(-1) U_0
   call ph%elpa_aux%hermitian_multiply("L","U",ph%n_basis,dm1,ovlp0,bh%n_lrow,&
        bh%n_lcol,ovlp1,bh%n_lrow,bh%n_lcol,ierr)

   call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

   ! ovlp0 = U_0^T U_1^(-T)
   call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp1,1,1,bh%desc,0.0_r8,ovlp0,1,1,&
        bh%desc)

   ! ovlp1 = U_1^(-1) U_0 P_0
   call ph%elpa_aux%hermitian_multiply("L","U",ph%n_basis,ovlp0,dm0,bh%n_lrow,&
        bh%n_lcol,ovlp1,bh%n_lrow,bh%n_lcol,ierr)

   call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

   ! dm0 = P_0 U_0^T U_1^(-T)
   call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp1,1,1,bh%desc,0.0_r8,dm0,1,1,&
        bh%desc)

   ! ovlp1 = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
   call ph%elpa_aux%hermitian_multiply("L","L",ph%n_basis,ovlp0,dm0,bh%n_lrow,&
        bh%n_lcol,ovlp1,bh%n_lrow,bh%n_lcol,ierr)

   call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

   call elsi_set_full_mat(ph,bh,LT_MAT,ovlp1)

   ! dm1 = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
   dm1(:,:) = ovlp1

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished density matrix extrapolation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

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
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_factor_ovlp_elpa_cmplx"

   call elsi_get_time(t0)

   ! S = U
   call ph%elpa_aux%cholesky(ovlp,ierr)

   call elsi_check_err(bh,"ELPA Cholesky factorization",ierr,caller)

   ! S = U^(-1)
   call ph%elpa_aux%invert_triangular(ovlp,ierr)

   call elsi_check_err(bh,"ELPA matrix inversion",ierr,caller)

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
   integer(kind=i4) :: ierr
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
      call ph%elpa_aux%hermitian_multiply("U","L",ph%n_basis,ovlp,ham,&
           bh%n_lrow,bh%n_lcol,evec,bh%n_lrow,bh%n_lcol,ierr)

      call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

      call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),evec,1,1,bh%desc,&
           (0.0_r8,0.0_r8),ham,1,1,bh%desc)

      evec(:,:) = ham

      call ph%elpa_aux%hermitian_multiply("U","U",ph%n_basis,ovlp,evec,&
           bh%n_lrow,bh%n_lcol,ham,bh%n_lrow,bh%n_lcol,ierr)

      call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

      call elsi_set_full_mat(ph,bh,UT_MAT,ham)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished transformation to standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

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

   if(ph%n_good < ph%n_basis) then
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
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   complex(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_back_ev_elpa_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)

   tmp(:,:) = evec

   if(ph%ill_ovlp) then
      call pzgemm("N","N",ph%n_basis,ph%n_states_solve,ph%n_good,&
           (1.0_r8,0.0_r8),ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,tmp,1,1,&
           bh%desc,(0.0_r8,0.0_r8),evec,1,1,bh%desc)
   else
      call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),ovlp,1,1,bh%desc,&
           (0.0_r8,0.0_r8),ham,1,1,bh%desc)

      call ph%elpa_aux%hermitian_multiply("L","N",ph%n_states,ham,tmp,&
           bh%n_lrow,bh%n_lcol,evec,bh%n_lrow,bh%n_lcol,ierr)

      call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)
   end if

   call elsi_deallocate(bh,tmp,"tmp")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished back-transformation of eigenvectors"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

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

   ! Compute sparsity
   if(bh%nnz_g == UNSET) then
      if(bh%nnz_l == UNSET) then
         bh%nnz_l = count(abs(ham) > bh%def0)
      end if

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,MPI_INTEGER4,MPI_SUM,bh%comm,ierr)

      call elsi_check_err(bh,"MPI_Allreduce",ierr,caller)
   end if

   ! Ill-conditioning check
   if(.not. ph%unit_ovlp .and. ph%elpa_first .and. ph%ill_check) then
      call elsi_check_ovlp_elpa(ph,bh,ovlp,eval,evec)
   end if

   if(ph%elpa_gpu == 1) then
      if(ph%n_calls <= ph%elpa_n_single) then
         write(msg,"(A)") "Starting ELPA eigensolver (GPU single precision)"
      else
         write(msg,"(A)") "Starting ELPA eigensolver (GPU)"
      end if
   else
      if(ph%n_calls <= ph%elpa_n_single) then
         write(msg,"(A)") "Starting ELPA eigensolver (single precision)"
      else
         write(msg,"(A)") "Starting ELPA eigensolver"
      end if
   end if
   call elsi_say(bh,msg)

   ! Transform to standard form
   if(.not. ph%unit_ovlp) then
      if(ph%elpa_first .and. ph%n_good == ph%n_basis) then
         ! Do Cholesky if not singular
         call elsi_factor_ovlp_elpa(ph,bh,ovlp)
      end if

      call elsi_reduce_evp_elpa(ph,bh,ham,ovlp,evec)
   end if

   call elsi_get_time(t0)

   ! Solve
   if(.not. associated(ph%elpa_solve)) then
      call elsi_elpa_setup(ph,bh,.false.)
   end if

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
subroutine elsi_update_dm_elpa_cmplx(ph,bh,ovlp0,ovlp1,dm0,dm1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: ovlp0(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp1(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: dm0(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: dm1(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_update_dm_elpa_cmplx"

   call elsi_get_time(t0)

   ! ovlp0 = U_0
   call ph%elpa_aux%cholesky(ovlp0,ierr)

   call elsi_check_err(bh,"ELPA Cholesky factorization",ierr,caller)

   ! ovlp1 = (U_1)^(-1)
   call elsi_factor_ovlp_elpa(ph,bh,ovlp1)

   ! dm1 = U_1^(-T)
   call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),ovlp1,1,1,bh%desc,&
        (0.0_r8,0.0_r8),dm1,1,1,bh%desc)

   ! ovlp1 = U_1^(-1) U_0
   call ph%elpa_aux%hermitian_multiply("L","U",ph%n_basis,dm1,ovlp0,bh%n_lrow,&
        bh%n_lcol,ovlp1,bh%n_lrow,bh%n_lcol,ierr)

   call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

   ! ovlp0 = U_0^T U_1^(-T)
   call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),ovlp1,1,1,bh%desc,&
        (0.0_r8,0.0_r8),ovlp0,1,1,bh%desc)

   ! ovlp1 = U_1^(-1) U_0 P_0
   call ph%elpa_aux%hermitian_multiply("L","U",ph%n_basis,ovlp0,dm0,bh%n_lrow,&
        bh%n_lcol,ovlp1,bh%n_lrow,bh%n_lcol,ierr)

   call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

   ! dm0 = P_0 U_0^T U_1^(-T)
   call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),ovlp1,1,1,bh%desc,&
        (0.0_r8,0.0_r8),dm0,1,1,bh%desc)

   ! ovlp1 = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
   call ph%elpa_aux%hermitian_multiply("L","L",ph%n_basis,ovlp0,dm0,bh%n_lrow,&
        bh%n_lcol,ovlp1,bh%n_lrow,bh%n_lcol,ierr)

   call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

   call elsi_set_full_mat(ph,bh,LT_MAT,ovlp1)

   ! dm1 = U_1^(-1) U_0 P_0 U_0^T U_1^(-T)
   dm1(:,:) = ovlp1

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished density matrix extrapolation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Set up an instance of ELPA.
!!
subroutine elsi_elpa_setup(ph,bh,is_aux)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   logical, intent(in) :: is_aux

   integer(kind=i4) :: ierr

   integer(kind=i4), external :: numroc

   character(len=*), parameter :: caller = "elsi_elpa_setup"

   if(is_aux) then
      ph%elpa_aux => elpa_allocate(ierr)

      call ph%elpa_aux%set("na",ph%n_basis,ierr)
      call ph%elpa_aux%set("nev",ph%n_basis,ierr)
      call ph%elpa_aux%set("nblk",bh%blk,ierr)
      call ph%elpa_aux%set("local_nrows",bh%n_lrow,ierr)
      call ph%elpa_aux%set("local_ncols",bh%n_lcol,ierr)
      call ph%elpa_aux%set("mpi_comm_parent",bh%comm,ierr)
      call ph%elpa_aux%set("mpi_comm_rows",ph%elpa_comm_row,ierr)
      call ph%elpa_aux%set("mpi_comm_cols",ph%elpa_comm_col,ierr)

      ierr = ph%elpa_aux%setup()

      call elsi_check_err(bh,"ELPA setup",ierr,caller)

      if(ph%elpa_solver == UNSET .or. ph%elpa_solver == 2) then
         call ph%elpa_aux%set("solver",2,ierr)

         if(ph%elpa_gpu == UNSET .or. ph%elpa_gpu == 0) then
            call ph%elpa_aux%set("gpu",0,ierr)
         else
            call ph%elpa_aux%set("gpu",1,ierr)
            call ph%elpa_aux%set("real_kernel",ELPA_2STAGE_REAL_GPU,ierr)
            call ph%elpa_aux%set("complex_kernel",ELPA_2STAGE_COMPLEX_GPU,ierr)
         end if
      else
         call ph%elpa_aux%set("solver",1,ierr)

         if(ph%elpa_gpu == UNSET .or. ph%elpa_gpu == 0) then
            call ph%elpa_aux%set("gpu",0,ierr)
         else
            call ph%elpa_aux%set("gpu",1,ierr)
         end if
      end if
   else
      ! Dimension should be number of non-ill-conditioned basis functions
      ph%elpa_n_lrow = numroc(ph%n_good,bh%blk,bh%my_prow,0,bh%n_prow)
      ph%elpa_n_lcol = numroc(ph%n_good,bh%blk,bh%my_pcol,0,bh%n_pcol)

      ph%elpa_solve => elpa_allocate(ierr)

      call ph%elpa_solve%set("na",ph%n_good,ierr)
      call ph%elpa_solve%set("nev",ph%n_states_solve,ierr)
      call ph%elpa_solve%set("nblk",bh%blk,ierr)
      call ph%elpa_solve%set("local_nrows",ph%elpa_n_lrow,ierr)
      call ph%elpa_solve%set("local_ncols",ph%elpa_n_lcol,ierr)
      call ph%elpa_solve%set("mpi_comm_parent",bh%comm,ierr)
      call ph%elpa_solve%set("mpi_comm_rows",ph%elpa_comm_row,ierr)
      call ph%elpa_solve%set("mpi_comm_cols",ph%elpa_comm_col,ierr)

      ierr = ph%elpa_solve%setup()

      call elsi_check_err(bh,"ELPA setup",ierr,caller)

      if(ph%elpa_autotune > 0) then
         if(ph%elpa_gpu == UNSET) then
            ph%elpa_autotune = ELPA_AUTOTUNE_MEDIUM
         else
            ph%elpa_autotune = ELPA_AUTOTUNE_FAST
         end if
      else
         if(ph%elpa_gpu == UNSET) then
            ph%elpa_gpu = 0
         end if

         if(ph%elpa_solver == UNSET) then
            ph%elpa_solver = 2
         end if
      end if

      if(ph%elpa_solver /= UNSET) then
         call ph%elpa_solve%set("solver",ph%elpa_solver,ierr)
      end if

      if(ph%elpa_gpu /= UNSET) then
         call ph%elpa_solve%set("gpu",ph%elpa_gpu,ierr)
      end if

      if(ph%elpa_gpu == 1) then
         if(ph%elpa_solver == 1) then
            call ph%elpa_solve%set("gpu_tridiag",1,ierr)
            call ph%elpa_solve%set("gpu_solve_tridi",1,ierr)
            call ph%elpa_solve%set("gpu_trans_ev",1,ierr)
         else
            call ph%elpa_solve%set("gpu_bandred",1,ierr)
            call ph%elpa_solve%set("gpu_tridiag_band",1,ierr)
            call ph%elpa_solve%set("gpu_solve_tridi",1,ierr)
            call ph%elpa_solve%set("gpu_trans_ev_tridi_to_band",1,ierr)
            call ph%elpa_solve%set("gpu_trans_ev_band_to_full",1,ierr)
            call ph%elpa_solve%set("real_kernel",ELPA_2STAGE_REAL_GPU,ierr)
            call ph%elpa_solve%set("complex_kernel",ELPA_2STAGE_COMPLEX_GPU,ierr)
         end if
      end if
   end if

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

   integer(kind=i4) :: i
   integer(kind=i4) :: ierr

   real(kind=r8), allocatable :: copy(:,:)
   real(kind=r4), allocatable :: copy_r4(:,:)
   real(kind=r4), allocatable :: evec_r4(:,:)
   real(kind=r4), allocatable :: eval_r4(:)

   character(len=*), parameter :: caller = "elsi_elpa_evec_real"

   if(sing_check) then
      call elsi_allocate(bh,copy,bh%n_lrow,bh%n_lcol,"copy",caller)

      copy(:,:) = mat

      call ph%elpa_aux%set("thres_pd_double",ph%ill_tol,ierr)

      if(ierr == 0) then
         call ph%elpa_aux%set("check_pd",1,ierr)
      end if

      call ph%elpa_aux%eigenvectors(copy,eval,evec,ierr)

      call elsi_check_err(bh,"ELPA eigensolver",ierr,caller)

      call ph%elpa_aux%set("check_pd",0,ierr)

      call elsi_deallocate(bh,copy,"copy")

      do i = 1,ph%n_basis
         if(eval(i) < ph%ill_tol) then
            ph%n_good = ph%n_good-1
         end if
      end do

      ph%n_states_solve = min(ph%n_good,ph%n_states)
      ph%ovlp_ev_min = eval(1)
      ph%ovlp_ev_max = eval(ph%n_basis)
   else
      if(ph%n_calls <= ph%elpa_n_single) then
         call elsi_allocate(bh,eval_r4,ph%n_basis,"eval_r4",caller)
         call elsi_allocate(bh,evec_r4,bh%n_lrow,bh%n_lcol,"evec_r4",caller)
         call elsi_allocate(bh,copy_r4,bh%n_lrow,bh%n_lcol,"copy_r4",caller)

         copy_r4(:,:) = real(mat,kind=r4)

         call ph%elpa_solve%eigenvectors(&
              copy_r4(1:ph%elpa_n_lrow,1:ph%elpa_n_lcol),eval_r4,&
              evec_r4(1:ph%elpa_n_lrow,1:ph%elpa_n_lcol),ierr)

         call elsi_check_err(bh,"ELPA eigensolver",ierr,caller)

         eval(:) = real(eval_r4,kind=r8)
         evec(:,:) = real(evec_r4,kind=r8)

         call elsi_deallocate(bh,eval_r4,"eval_r4")
         call elsi_deallocate(bh,evec_r4,"evec_r4")
         call elsi_deallocate(bh,copy_r4,"copy_r4")
      else
         if(ph%elpa_autotune > 0) then
            if(.not. associated(ph%elpa_tune)) then
               ph%elpa_tune => ph%elpa_solve%autotune_setup(ph%elpa_autotune,&
                  ELPA_AUTOTUNE_DOMAIN_REAL,ierr)
            end if

            if(.not. ph%elpa_solve%autotune_step(ph%elpa_tune,ierr)) then
               call ph%elpa_solve%autotune_set_best(ph%elpa_tune,ierr)
               call elpa_autotune_deallocate(ph%elpa_tune,ierr)

               nullify(ph%elpa_tune)

               ph%elpa_autotune = 0
            end if
         end if

         call ph%elpa_solve%eigenvectors(&
              mat(1:ph%elpa_n_lrow,1:ph%elpa_n_lcol),eval,&
              evec(1:ph%elpa_n_lrow,1:ph%elpa_n_lcol),ierr)

         call elsi_check_err(bh,"ELPA eigensolver",ierr,caller)
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

   integer(kind=i4) :: i
   integer(kind=i4) :: ierr

   complex(kind=r8), allocatable :: copy(:,:)
   complex(kind=r4), allocatable :: copy_r4(:,:)
   complex(kind=r4), allocatable :: evec_r4(:,:)
   real(kind=r4), allocatable :: eval_r4(:)

   character(len=*), parameter :: caller = "elsi_elpa_evec_cmplx"

   if(sing_check) then
      call elsi_allocate(bh,copy,bh%n_lrow,bh%n_lcol,"copy",caller)

      copy(:,:) = mat

      call ph%elpa_aux%set("thres_pd_double",ph%ill_tol,ierr)

      if(ierr == 0) then
         call ph%elpa_aux%set("check_pd",1,ierr)
      end if

      call ph%elpa_aux%eigenvectors(copy,eval,evec,ierr)

      call elsi_check_err(bh,"ELPA eigensolver",ierr,caller)

      call ph%elpa_aux%set("check_pd",0,ierr)

      call elsi_deallocate(bh,copy,"copy")

      do i = 1,ph%n_basis
         if(eval(i) < ph%ill_tol) then
            ph%n_good = ph%n_good-1
         end if
      end do

      ph%n_states_solve = min(ph%n_good,ph%n_states)
      ph%ovlp_ev_min = eval(1)
      ph%ovlp_ev_max = eval(ph%n_basis)
   else
      if(ph%n_calls <= ph%elpa_n_single) then
         call elsi_allocate(bh,eval_r4,ph%n_basis,"eval_r4",caller)
         call elsi_allocate(bh,evec_r4,bh%n_lrow,bh%n_lcol,"evec_r4",caller)
         call elsi_allocate(bh,copy_r4,bh%n_lrow,bh%n_lcol,"copy_r4",caller)

         copy_r4(:,:) = cmplx(mat,kind=r4)

         call ph%elpa_solve%eigenvectors(&
              copy_r4(1:ph%elpa_n_lrow,1:ph%elpa_n_lcol),eval_r4,&
              evec_r4(1:ph%elpa_n_lrow,1:ph%elpa_n_lcol),ierr)

         call elsi_check_err(bh,"ELPA eigensolver",ierr,caller)

         eval(:) = real(eval_r4,kind=r8)
         evec(:,:) = cmplx(evec_r4,kind=r8)

         call elsi_deallocate(bh,eval_r4,"eval_r4")
         call elsi_deallocate(bh,evec_r4,"evec_r4")
         call elsi_deallocate(bh,copy_r4,"copy_r4")
      else
         if(ph%elpa_autotune > 0) then
            if(.not. associated(ph%elpa_tune)) then
               ph%elpa_tune => ph%elpa_solve%autotune_setup(ph%elpa_autotune,&
                  ELPA_AUTOTUNE_DOMAIN_COMPLEX,ierr)
            end if

            if(.not. ph%elpa_solve%autotune_step(ph%elpa_tune,ierr)) then
               call ph%elpa_solve%autotune_set_best(ph%elpa_tune,ierr)
               call elpa_autotune_deallocate(ph%elpa_tune,ierr)

               nullify(ph%elpa_tune)

               ph%elpa_autotune = 0
            end if
         end if

         call ph%elpa_solve%eigenvectors(&
              mat(1:ph%elpa_n_lrow,1:ph%elpa_n_lcol),eval,&
              evec(1:ph%elpa_n_lrow,1:ph%elpa_n_lcol),ierr)

         call elsi_check_err(bh,"ELPA eigensolver",ierr,caller)
      end if
   end if

end subroutine

!>
!! Interface to ELPA tridiagonal solver, to be used together with LAPACK.
!!
subroutine elsi_elpa_tridiag(ph,bh,diag,offd,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: diag(ph%n_good)
   real(kind=r8), intent(inout) :: offd(ph%n_good)
   real(kind=r8), intent(inout) :: evec(ph%n_good,ph%n_good)

   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_elpa_tridiag"

   ierr = elpa_init(20180525)

   call elsi_check_err(bh,"ELPA initialization failed",ierr,caller)

   ph%elpa_aux => elpa_allocate(ierr)

   call ph%elpa_aux%set("na",ph%n_good,ierr)
   call ph%elpa_aux%set("nev",ph%n_states_solve,ierr)
   call ph%elpa_aux%set("nblk",bh%blk,ierr)
   call ph%elpa_aux%set("local_nrows",ph%n_good,ierr)
   call ph%elpa_aux%set("local_ncols",ph%n_good,ierr)
   call ph%elpa_aux%set("mpi_comm_parent",MPI_COMM_SELF,ierr)
   call ph%elpa_aux%set("mpi_comm_rows",MPI_COMM_SELF,ierr)
   call ph%elpa_aux%set("mpi_comm_cols",MPI_COMM_SELF,ierr)

   ierr = ph%elpa_aux%setup()

   call elsi_check_err(bh,"ELPA setup",ierr,caller)

   call ph%elpa_aux%solve_tridiagonal(diag,offd,evec,ierr)

   call elsi_check_err(bh,"ELPA tridiagonal eigensolver",ierr,caller)

   call elpa_deallocate(ph%elpa_aux,ierr)

   nullify(ph%elpa_aux)

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
      if(associated(ph%elpa_solve)) then
         call elpa_deallocate(ph%elpa_solve,ierr)

         nullify(ph%elpa_solve)
      end if

      if(associated(ph%elpa_aux)) then
         call elpa_deallocate(ph%elpa_aux,ierr)

         nullify(ph%elpa_aux)
      end if

      if(associated(ph%elpa_tune)) then
         call elpa_autotune_deallocate(ph%elpa_tune,ierr)

         nullify(ph%elpa_tune)
      end if

      call MPI_Comm_free(ph%elpa_comm_row,ierr)
      call MPI_Comm_free(ph%elpa_comm_col,ierr)
   end if

   ph%elpa_first = .true.
   ph%elpa_started = .false.

end subroutine

end module ELSI_ELPA
