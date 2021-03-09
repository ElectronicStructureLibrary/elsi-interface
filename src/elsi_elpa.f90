! Copyright (c) 2015-2021, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Interface to ELPA-AEO.
!!
module ELSI_ELPA

   use ELSI_CONSTANT, only: LT_MAT,UT_MAT,UNSET,FC_BASIC,FC_PLUS_V
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI
   use ELSI_OUTPUT, only: elsi_say,elsi_get_time
   use ELSI_PRECISION, only: r4,r8,i4
   use ELSI_SORT, only: elsi_heapsort,elsi_permute
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
   public :: elsi_do_fc_elpa
   public :: elsi_undo_fc_elpa
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

   interface elsi_do_fc_elpa
      module procedure elsi_do_fc_elpa_real
      module procedure elsi_do_fc_elpa_cmplx
   end interface

   interface elsi_undo_fc_elpa
      module procedure elsi_undo_fc_elpa_real
      module procedure elsi_undo_fc_elpa_cmplx
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

   integer(kind=i4), external :: numroc

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

   ! Switch back to full dimension in case of frozen core
   if(ph%n_basis_c > 0) then
      ph%n_basis = ph%n_basis_v+ph%n_basis_c
      ph%n_states = ph%n_states+ph%n_basis_c
      ph%n_good = ph%n_good+ph%n_basis_c
      ph%n_states_solve = ph%n_states_solve+ph%n_basis_c
      bh%n_lrow = numroc(ph%n_basis,bh%blk,bh%my_prow,0,bh%n_prow)
      bh%n_lcol = numroc(ph%n_basis,bh%blk,bh%my_pcol,0,bh%n_pcol)

      call descinit(bh%desc,ph%n_basis,ph%n_basis,bh%blk,bh%blk,0,0,&
           bh%blacs_ctxt,max(1,bh%n_lrow),ierr)
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
!! Freeze core orbitals by transforming Hamiltonian and overlap.
!!
subroutine elsi_do_fc_elpa_real(ph,bh,ham,ovlp,evec,perm,ham_v,ovlp_v,evec_v)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)
   integer(kind=i4), intent(in) :: perm(ph%n_basis)
   real(kind=r8), intent(out) :: ham_v(ph%n_lrow_v,ph%n_lcol_v)
   real(kind=r8), intent(out) :: ovlp_v(ph%n_lrow_v,ph%n_lcol_v)
   real(kind=r8), intent(out) :: evec_v(ph%n_lrow_v,ph%n_lcol_v)

   integer(kind=i4) :: i
   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_do_fc_elpa_real"

   if(ph%elpa_first) then
      if(ph%fc_perm) then
         evec(:,:) = ovlp

         do i = 1,ph%n_basis
            if(i /= perm(i)) then
               call pdcopy(ph%n_basis,ovlp,1,perm(i),bh%desc,1,evec,1,i,&
                    bh%desc,1)
            end if
         end do

         ovlp(:,:) = evec

         do i = 1,ph%n_basis
            if(i /= perm(i)) then
               call pdcopy(ph%n_basis,evec,perm(i),1,bh%desc,ph%n_basis,ovlp,i,&
                    1,bh%desc,ph%n_basis)
            end if
         end do
      end if

      evec(:,:) = ovlp

      ! S_vv = S_vv - S_vc * S_cv
      call pdsyrk("U","N",ph%n_basis_v,ph%n_basis_c,-1.0_r8,evec,ph%n_basis_c+1,&
           1,bh%desc,1.0_r8,ovlp,ph%n_basis_c+1,ph%n_basis_c+1,bh%desc)

      call elsi_set_full_mat(ph,bh,UT_MAT,ovlp)

      call pdgemr2d(ph%n_basis_v,ph%n_basis_v,ovlp,ph%n_basis_c+1,&
           ph%n_basis_c+1,bh%desc,ovlp_v,1,1,ph%desc_v,bh%blacs_ctxt)
   end if

   if(ph%fc_perm) then
      evec(:,:) = ham

      do i = 1,ph%n_basis
         if(i /= perm(i)) then
            call pdcopy(ph%n_basis,ham,1,perm(i),bh%desc,1,evec,1,i,bh%desc,1)
         end if
      end do

      ham(:,:) = evec

      do i = 1,ph%n_basis
         if(i /= perm(i)) then
            call pdcopy(ph%n_basis,evec,perm(i),1,bh%desc,ph%n_basis,ham,i,1,&
                 bh%desc,ph%n_basis)
         end if
      end do
   end if

   call pdgemr2d(ph%n_basis_v,ph%n_basis_v,ham,ph%n_basis_c+1,ph%n_basis_c+1,&
        bh%desc,ham_v,1,1,ph%desc_v,bh%blacs_ctxt)

   ! Compute H_vv
   evec(:,:) = 0.0_r8

   if(ph%fc_method == FC_PLUS_V) then
      ! H_vv = H_vv + S_vc * H_cc * S_cv - H_vc * S_cv - S_vc * H_cv
      ! More accurate than H_vv = H_vv - S_vc * H_cc * S_cv
      ! H_vv = A + A^*
      ! A = 0.5*H_vv + (0.5*S_vc * H_cc - H_vc) * S_cv
      call pdgemm("N","N",ph%n_basis_v,ph%n_basis_c,ph%n_basis_c,0.5_r8,ovlp,&
           ph%n_basis_c+1,1,bh%desc,ham,1,1,bh%desc,0.0_r8,evec,ph%n_basis_c+1,&
           1,bh%desc)

      evec(:,:) = evec-ham

      call pdgemm("N","N",ph%n_basis_v,ph%n_basis_v,ph%n_basis_c,1.0_r8,evec,&
           ph%n_basis_c+1,1,bh%desc,ovlp,1,ph%n_basis_c+1,bh%desc,0.5_r8,ham_v,&
           1,1,ph%desc_v)

      call pdtran(ph%n_basis_v,ph%n_basis_v,1.0_r8,ham_v,1,1,ph%desc_v,0.0_r8,&
           evec_v,1,1,ph%desc_v)

      ham_v(:,:) = ham_v+evec_v
   else
      ! H_vv = H_vv - S_vc * H_cc * S_cv
      call pdgemm("N","N",ph%n_basis_v,ph%n_basis_c,ph%n_basis_c,1.0_r8,ovlp,&
           ph%n_basis_c+1,1,bh%desc,ham,1,1,bh%desc,0.0_r8,evec,ph%n_basis_c+1,&
           1,bh%desc)

      call pdgemm("N","N",ph%n_basis_v,ph%n_basis_v,ph%n_basis_c,-1.0_r8,evec,&
           ph%n_basis_c+1,1,bh%desc,ovlp,1,ph%n_basis_c+1,bh%desc,1.0_r8,ham_v,&
           1,1,ph%desc_v)
   end if

   ! Switch to valence dimension
   ph%n_basis = ph%n_basis-ph%n_basis_c
   ph%n_states = ph%n_states-ph%n_basis_c
   ph%n_good = ph%n_good-ph%n_basis_c
   ph%n_states_solve = ph%n_states_solve-ph%n_basis_c
   bh%n_lrow = ph%n_lrow_v
   bh%n_lcol = ph%n_lcol_v

   call descinit(bh%desc,ph%n_basis,ph%n_basis,bh%blk,bh%blk,0,0,bh%blacs_ctxt,&
        max(1,bh%n_lrow),ierr)

end subroutine

!>
!! Transforming eigenvectors back to unfrozen eigenproblem.
!!
subroutine elsi_undo_fc_elpa_real(ph,bh,ham,ovlp,evec,perm,eval_c,evec_v)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)
   integer(kind=i4), intent(in) :: perm(ph%n_basis)
   real(kind=r8), intent(out) :: eval_c(ph%n_basis_c)
   real(kind=r8), intent(in) :: evec_v(ph%n_lrow_v,ph%n_lcol_v)

   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: g_row
   integer(kind=i4) :: g_col
   integer(kind=i4) :: ierr

   real(kind=r8), allocatable :: ovlp_d(:)
   real(kind=r8), allocatable :: tmp(:)
   integer(kind=i4), allocatable :: idx(:)
   integer(kind=i4), allocatable :: tmp2(:)

   character(len=*), parameter :: caller = "elsi_undo_fc_elpa_real"

   call elsi_allocate(bh,ovlp_d,ph%n_basis_c,"ovlp_d",caller)
   call elsi_allocate(bh,tmp,ph%n_basis_c,"tmp",caller)

   eval_c(:) = 0.0_r8

   do j = 1,bh%n_lcol
      call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,j,g_col)

      if(g_col > ph%n_basis_c) then
         exit
      end if

      do i = 1,bh%n_lrow
         call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i,g_row)

         if(g_row > ph%n_basis_c) then
            exit
         end if

         if(g_row == g_col) then
            if(ph%fc_method > FC_BASIC) then
               eval_c(g_col) = ham(i,j)/ovlp(i,j)
            else
               eval_c(g_col) = ham(i,j)
            end if

            ovlp_d(g_col) = ovlp(i,j)

            exit
         end if
      end do
   end do

   call MPI_Allreduce(eval_c,tmp,ph%n_basis_c,MPI_REAL8,MPI_SUM,bh%comm,ierr)

   call elsi_check_err(bh,"MPI_Allreduce",ierr,caller)

   eval_c(:) = tmp

   call MPI_Allreduce(ovlp_d,tmp,ph%n_basis_c,MPI_REAL8,MPI_SUM,bh%comm,ierr)

   call elsi_check_err(bh,"MPI_Allreduce",ierr,caller)

   ovlp_d(:) = tmp

   call elsi_deallocate(bh,tmp,"tmp")
   call elsi_allocate(bh,idx,ph%n_basis_c,"idx",caller)
   call elsi_allocate(bh,tmp2,ph%n_basis_c,"tmp2",caller)

   do i = 1,ph%n_basis_c
      idx(i) = i
   end do

   call elsi_heapsort(ph%n_basis_c,eval_c,tmp2)
   call elsi_permute(ph%n_basis_c,tmp2,idx)
   call elsi_permute(ph%n_basis_c,tmp2,ovlp_d)

   evec(:,:) = 0.0_r8

   do j = 1,bh%n_lcol
      call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,j,g_col)

      if(g_col > ph%n_basis_c) then
         exit
      end if

      do i = 1,bh%n_lrow
         call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i,g_row)

         if(g_row > ph%n_basis_c) then
            exit
         end if

         if(g_row == idx(g_col)) then
            if(ph%fc_method > FC_BASIC) then
               evec(i,j) = 1.0_r8/sqrt(ovlp_d(g_col))
            else
               evec(i,j) = 1.0_r8
            end if

            exit
         end if
      end do
   end do

   call elsi_deallocate(bh,ovlp_d,"ovlp_d")
   call elsi_deallocate(bh,idx,"idx")
   call elsi_deallocate(bh,tmp2,"tmp2")

   call pdgemr2d(ph%n_basis_v,ph%n_basis_v,evec_v,1,1,ph%desc_v,evec,&
        ph%n_basis_c+1,ph%n_basis_c+1,bh%desc,bh%blacs_ctxt)

   ! C_cv = -S_cv * C_vv
   call pdgemm("N","N",ph%n_basis_c,ph%n_basis_v,ph%n_basis_v,-1.0_r8,ovlp,1,&
        ph%n_basis_c+1,bh%desc,evec_v,1,1,ph%desc_v,0.0_r8,evec,1,&
        ph%n_basis_c+1,bh%desc)

   if(ph%fc_perm) then
      ham(:,:) = evec

      do i = 1,ph%n_basis
         if(i /= perm(i)) then
            call pdcopy(ph%n_basis,ham,i,1,bh%desc,ph%n_basis,evec,perm(i),1,&
                 bh%desc,ph%n_basis)
         end if
      end do
   end if

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

   integer(kind=i4), external :: numroc

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

   ! Switch back to full dimension in case of frozen core
   if(ph%n_basis_c > 0) then
      ph%n_basis = ph%n_basis_v+ph%n_basis_c
      ph%n_states = ph%n_states+ph%n_basis_c
      ph%n_good = ph%n_good+ph%n_basis_c
      ph%n_states_solve = ph%n_states_solve+ph%n_basis_c
      bh%n_lrow = numroc(ph%n_basis,bh%blk,bh%my_prow,0,bh%n_prow)
      bh%n_lcol = numroc(ph%n_basis,bh%blk,bh%my_pcol,0,bh%n_pcol)

      call descinit(bh%desc,ph%n_basis,ph%n_basis,bh%blk,bh%blk,0,0,&
           bh%blacs_ctxt,max(1,bh%n_lrow),ierr)
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
!! Freeze core orbitals by transforming Hamiltonian and overlap.
!!
subroutine elsi_do_fc_elpa_cmplx(ph,bh,ham,ovlp,evec,perm,ham_v,ovlp_v,evec_v)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)
   integer(kind=i4), intent(in) :: perm(ph%n_basis)
   complex(kind=r8), intent(out) :: ham_v(ph%n_lrow_v,ph%n_lcol_v)
   complex(kind=r8), intent(out) :: ovlp_v(ph%n_lrow_v,ph%n_lcol_v)
   complex(kind=r8), intent(out) :: evec_v(ph%n_lrow_v,ph%n_lcol_v)

   integer(kind=i4) :: i
   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_do_fc_elpa_cmplx"

   if(ph%elpa_first) then
      if(ph%fc_perm) then
         evec(:,:) = ovlp

         do i = 1,ph%n_basis
            if(i /= perm(i)) then
               call pzcopy(ph%n_basis,ovlp,1,perm(i),bh%desc,1,evec,1,i,bh%desc,1)
            end if
         end do

         ovlp(:,:) = evec

         do i = 1,ph%n_basis
            if(i /= perm(i)) then
               call pzcopy(ph%n_basis,evec,perm(i),1,bh%desc,ph%n_basis,ovlp,i,1,&
                    bh%desc,ph%n_basis)
            end if
         end do
      end if

      evec(:,:) = ovlp

      ! S_vv = S_vv - S_vc * S_cv
      call pzherk("U","N",ph%n_basis_v,ph%n_basis_c,(-1.0_r8,0.0_r8),evec,&
           ph%n_basis_c+1,1,bh%desc,(1.0_r8,0.0_r8),ovlp,ph%n_basis_c+1,&
           ph%n_basis_c+1,bh%desc)

      call elsi_set_full_mat(ph,bh,UT_MAT,ovlp)

      call pzgemr2d(ph%n_basis_v,ph%n_basis_v,ovlp,ph%n_basis_c+1,&
           ph%n_basis_c+1,bh%desc,ovlp_v,1,1,ph%desc_v,bh%blacs_ctxt)
   end if

   if(ph%fc_perm) then
      evec(:,:) = ham

      do i = 1,ph%n_basis
         if(i /= perm(i)) then
            call pzcopy(ph%n_basis,ham,1,perm(i),bh%desc,1,evec,1,i,bh%desc,1)
         end if
      end do

      ham(:,:) = evec

      do i = 1,ph%n_basis
         if(i /= perm(i)) then
            call pzcopy(ph%n_basis,evec,perm(i),1,bh%desc,ph%n_basis,ham,i,1,&
                 bh%desc,ph%n_basis)
         end if
      end do
   end if

   call pzgemr2d(ph%n_basis_v,ph%n_basis_v,ham,ph%n_basis_c+1,ph%n_basis_c+1,&
        bh%desc,ham_v,1,1,ph%desc_v,bh%blacs_ctxt)

   ! Compute H_vv
   evec(:,:) = (0.0_r8,0.0_r8)

   if(ph%fc_method == FC_PLUS_V) then
      ! H_vv = H_vv + S_vc * H_cc * S_cv - H_vc * S_cv - S_vc * H_cv
      ! More accurate than H_vv = H_vv - S_vc * H_cc * S_cv
      ! H_vv = A + A^*
      ! A = 0.5*H_vv + (0.5*S_vc * H_cc - H_vc) * S_cv
      call pzgemm("N","N",ph%n_basis_v,ph%n_basis_c,ph%n_basis_c,&
           (0.5_r8,0.0_r8),ovlp,ph%n_basis_c+1,1,bh%desc,ham,1,1,bh%desc,&
           (0.0_r8,0.0_r8),evec,ph%n_basis_c+1,1,bh%desc)

      evec(:,:) = evec-ham

      call pzgemm("N","N",ph%n_basis_v,ph%n_basis_v,ph%n_basis_c,&
           (1.0_r8,0.0_r8),evec,ph%n_basis_c+1,1,bh%desc,ovlp,1,ph%n_basis_c+1,&
           bh%desc,(0.5_r8,0.0_r8),ham_v,1,1,ph%desc_v)

      call pztranc(ph%n_basis_v,ph%n_basis_v,(1.0_r8,0.0_r8),ham_v,1,1,&
           ph%desc_v,(0.0_r8,0.0_r8),evec_v,1,1,ph%desc_v)

      ham_v(:,:) = ham_v+evec_v
   else
      ! H_vv = H_vv - S_vc * H_cc * S_cv
      call pzgemm("N","N",ph%n_basis_v,ph%n_basis_c,ph%n_basis_c,&
           (1.0_r8,0.0_r8),ovlp,ph%n_basis_c+1,1,bh%desc,ham,1,1,bh%desc,&
           (0.0_r8,0.0_r8),evec,ph%n_basis_c+1,1,bh%desc)

      call pzgemm("N","N",ph%n_basis_v,ph%n_basis_v,ph%n_basis_c,&
           (-1.0_r8,0.0_r8),evec,ph%n_basis_c+1,1,bh%desc,ovlp,1,&
           ph%n_basis_c+1,bh%desc,(1.0_r8,0.0_r8),ham_v,1,1,ph%desc_v)
   end if

   ! Switch to valence dimension
   ph%n_basis = ph%n_basis-ph%n_basis_c
   ph%n_states = ph%n_states-ph%n_basis_c
   ph%n_good = ph%n_good-ph%n_basis_c
   ph%n_states_solve = ph%n_states_solve-ph%n_basis_c
   bh%n_lrow = ph%n_lrow_v
   bh%n_lcol = ph%n_lcol_v

   call descinit(bh%desc,ph%n_basis,ph%n_basis,bh%blk,bh%blk,0,0,bh%blacs_ctxt,&
        max(1,bh%n_lrow),ierr)

end subroutine

!>
!! Transforming eigenvectors back to unfrozen eigenproblem.
!!
subroutine elsi_undo_fc_elpa_cmplx(ph,bh,ham,ovlp,evec,perm,eval_c,evec_v)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)
   integer(kind=i4), intent(in) :: perm(ph%n_basis)
   real(kind=r8), intent(out) :: eval_c(ph%n_basis_c)
   complex(kind=r8), intent(in) :: evec_v(ph%n_lrow_v,ph%n_lcol_v)

   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: g_row
   integer(kind=i4) :: g_col
   integer(kind=i4) :: ierr

   real(kind=r8), allocatable :: ovlp_d(:)
   real(kind=r8), allocatable :: tmp(:)
   integer(kind=i4), allocatable :: idx(:)
   integer(kind=i4), allocatable :: tmp2(:)

   character(len=*), parameter :: caller = "elsi_undo_fc_elpa_cmplx"

   call elsi_allocate(bh,ovlp_d,ph%n_basis_c,"ovlp_d",caller)
   call elsi_allocate(bh,tmp,ph%n_basis_c,"tmp",caller)

   eval_c(:) = 0.0_r8

   do j = 1,bh%n_lcol
      call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,j,g_col)

      if(g_col > ph%n_basis_c) then
         exit
      end if

      do i = 1,bh%n_lrow
         call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i,g_row)

         if(g_row > ph%n_basis_c) then
            exit
         end if

         if(g_row == g_col) then
            if(ph%fc_method > FC_BASIC) then
               eval_c(g_col) = real(ham(i,j),kind=r8)/real(ovlp(i,j),kind=r8)
            else
               eval_c(g_col) = real(ham(i,j),kind=r8)
            end if

            ovlp_d(g_col) = real(ovlp(i,j),kind=r8)

            exit
         end if
      end do
   end do

   call MPI_Allreduce(eval_c,tmp,ph%n_basis_c,MPI_REAL8,MPI_SUM,bh%comm,ierr)

   call elsi_check_err(bh,"MPI_Allreduce",ierr,caller)

   eval_c(:) = tmp

   call MPI_Allreduce(ovlp_d,tmp,ph%n_basis_c,MPI_REAL8,MPI_SUM,bh%comm,ierr)

   call elsi_check_err(bh,"MPI_Allreduce",ierr,caller)

   ovlp_d(:) = tmp

   call elsi_deallocate(bh,tmp,"tmp")
   call elsi_allocate(bh,idx,ph%n_basis_c,"idx",caller)
   call elsi_allocate(bh,tmp2,ph%n_basis_c,"tmp2",caller)

   do i = 1,ph%n_basis_c
      idx(i) = i
   end do

   call elsi_heapsort(ph%n_basis_c,eval_c,tmp2)
   call elsi_permute(ph%n_basis_c,tmp2,idx)
   call elsi_permute(ph%n_basis_c,tmp2,ovlp_d)

   evec(:,:) = (0.0_r8,0.0_r8)

   do j = 1,bh%n_lcol
      call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,j,g_col)

      if(g_col > ph%n_basis_c) then
         exit
      end if

      do i = 1,bh%n_lrow
         call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i,g_row)

         if(g_row > ph%n_basis_c) then
            exit
         end if

         if(g_row == idx(g_col)) then
            if(ph%fc_method > FC_BASIC) then
               evec(i,j) = (1.0_r8,0.0_r8)/sqrt(ovlp_d(g_col))
            else
               evec(i,j) = (1.0_r8,0.0_r8)
            end if

            exit
         end if
      end do
   end do

   call elsi_deallocate(bh,ovlp_d,"ovlp_d")
   call elsi_deallocate(bh,idx,"idx")
   call elsi_deallocate(bh,tmp2,"tmp2")

   call pzgemr2d(ph%n_basis_v,ph%n_basis_v,evec_v,1,1,ph%desc_v,evec,&
        ph%n_basis_c+1,ph%n_basis_c+1,bh%desc,bh%blacs_ctxt)

   ! C_cv = -S_cv * C_vv
   call pzgemm("N","N",ph%n_basis_c,ph%n_basis_v,ph%n_basis_v,(-1.0_r8,0.0_r8),&
        ovlp,1,ph%n_basis_c+1,bh%desc,evec_v,1,1,ph%desc_v,(0.0_r8,0.0_r8),&
        evec,1,ph%n_basis_c+1,bh%desc)

   if(ph%fc_perm) then
      ham(:,:) = evec

      do i = 1,ph%n_basis
         if(i /= perm(i)) then
            call pzcopy(ph%n_basis,ham,i,1,bh%desc,ph%n_basis,evec,perm(i),1,&
                 bh%desc,ph%n_basis)
         end if
      end do
   end if

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
