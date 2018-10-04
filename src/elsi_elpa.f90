! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to ELPA.
!!
module ELSI_ELPA

   use ELSI_CONSTANTS, only: BLACS_DENSE,UT_MAT
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_IO, only: elsi_say,elsi_get_time
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI, only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,mpi_integer4
   use ELSI_OCC, only: elsi_mu_and_occ,elsi_entropy
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS, only: elsi_get_nnz,elsi_set_full_mat
   use CHECK_SINGULARITY, only: elpa_check_singularity_real_double,&
       elpa_check_singularity_complex_double
   use ELPA1, only: elpa_print_times,elpa_solve_evp_real_1stage_double,&
       elpa_solve_evp_complex_1stage_double,elpa_cholesky_real_double,&
       elpa_cholesky_complex_double,elpa_invert_trm_real_double,&
       elpa_invert_trm_complex_double,elpa_mult_at_b_real_double,&
       elpa_mult_ah_b_complex_double
   use ELPA2, only: elpa_solve_evp_real_2stage_double,&
       elpa_solve_evp_complex_2stage_double

   implicit none

   private

   public :: elsi_init_elpa
   public :: elsi_cleanup_elpa
   public :: elsi_compute_occ_elpa
   public :: elsi_compute_dm_elpa
   public :: elsi_compute_edm_elpa
   public :: elsi_solve_elpa

   interface elsi_compute_dm_elpa
      module procedure elsi_compute_dm_elpa_real
      module procedure elsi_compute_dm_elpa_cmplx
   end interface

   interface elsi_compute_edm_elpa
      module procedure elsi_compute_edm_elpa_real
      module procedure elsi_compute_edm_elpa_cmplx
   end interface

   interface elsi_solve_elpa
      module procedure elsi_solve_elpa_real
      module procedure elsi_solve_elpa_cmplx
   end interface

contains

!>
!! This routine initializes ELPA.
!!
subroutine elsi_init_elpa(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh

   integer(kind=i4) :: ierr
   character(len=200) :: info_str

   character(len=*), parameter :: caller = "elsi_init_elpa"

   if(.not. ph%elpa_started) then
      call MPI_Comm_split(bh%comm,bh%my_pcol,bh%my_prow,ph%elpa_comm_row,ierr)

      call elsi_check_mpi(bh,"MPI_Comm_split",ierr,caller)

      call MPI_Comm_split(bh%comm,bh%my_prow,bh%my_pcol,ph%elpa_comm_col,ierr)

      call elsi_check_mpi(bh,"MPI_Comm_split",ierr,caller)

      ph%elpa_started = .true.

      if(ph%elpa_gpu) then
         write(info_str,"(2X,A)") "No ELPA GPU acceleration available"
         call elsi_say(bh,info_str)
      end if

      if(ph%elpa_n_single > 0) then
         write(info_str,"(2X,A)") "No single precision ELPA available"
         call elsi_say(bh,info_str)
      end if

      if(ph%elpa_autotune) then
         write(info_str,"(2X,A)") "No ELPA auto-tuning available"
         call elsi_say(bh,info_str)
      end if
   end if

end subroutine

!>
!! This routine computes the chemical potential and occupation numbers.
!!
subroutine elsi_compute_occ_elpa(ph,bh,eval,occ)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: occ(ph%n_states,ph%n_spins,ph%n_kpts)

   real(kind=r8) :: mu
   real(kind=r8) :: ts
   real(kind=r8) :: n_electrons
   integer(kind=i4) :: n_states
   integer(kind=i4) :: n_spins
   integer(kind=i4) :: n_kpts
   integer(kind=i4) :: i
   integer(kind=i4) :: ierr

   real(kind=r8), allocatable :: eval_all(:,:,:)
   real(kind=r8), allocatable :: k_weight(:)
   real(kind=r8), allocatable :: tmp_real1(:)
   real(kind=r8), allocatable :: tmp_real2(:,:,:)

   character(len=*), parameter :: caller = "elsi_compute_occ_elpa"

   ! Gather eigenvalues and occupation numbers
   call elsi_allocate(bh,eval_all,ph%n_states,ph%n_spins,ph%n_kpts,"eval_all",&
           caller)
   call elsi_allocate(bh,k_weight,ph%n_kpts,"k_weight",caller)

   if(ph%n_kpts > 1) then
      call elsi_allocate(bh,tmp_real1,ph%n_kpts,"tmp_real",caller)

      if(bh%myid == 0 .and. ph%i_spin == 1) then
         tmp_real1(ph%i_kpt) = ph%i_weight
      end if

      call MPI_Allreduce(tmp_real1,k_weight,ph%n_kpts,mpi_real8,mpi_sum,&
              bh%comm_all,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,tmp_real1,"tmp_real")
   else
      k_weight = ph%i_weight
   end if

   if(ph%n_spins*ph%n_kpts > 1) then
      call elsi_allocate(bh,tmp_real2,ph%n_states,ph%n_spins,ph%n_kpts,&
              "tmp_real",caller)

      if(bh%myid == 0) then
         tmp_real2(:,ph%i_spin,ph%i_kpt) = eval(1:ph%n_states)
      end if

      call MPI_Allreduce(tmp_real2,eval_all,ph%n_states*ph%n_spins*ph%n_kpts,&
              mpi_real8,mpi_sum,bh%comm_all,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,tmp_real2,"tmp_real")
   else
      eval_all(:,ph%i_spin,ph%i_kpt) = eval(1:ph%n_states)
   end if

   ! Calculate chemical potential, occupation numbers, and electronic entropy
   n_electrons = ph%n_electrons
   n_states = ph%n_states
   n_spins = ph%n_spins
   n_kpts = ph%n_kpts

   call elsi_mu_and_occ(ph,bh,n_electrons,n_states,n_spins,n_kpts,k_weight,&
           eval_all,occ,mu)

   call elsi_entropy(ph,n_states,n_spins,n_kpts,k_weight,eval_all,occ,mu,ts)

   ph%mu = mu
   ph%ts = ts

   ! Calculate band structure energy
   ph%ebs = 0.0_r8

   do i = 1,ph%n_states_solve
      ph%ebs = ph%ebs+eval(i)*occ(i,ph%i_spin,ph%i_kpt)
   end do

   call elsi_deallocate(bh,eval_all,"eval_all")
   call elsi_deallocate(bh,k_weight,"k_weight")

end subroutine

!>
!! This routine constructs the density matrix.
!!
subroutine elsi_compute_dm_elpa_real(ph,bh,row_map,col_map,evec,occ,dm,work)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   real(kind=r8), intent(in) :: evec(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: occ(ph%n_states,ph%n_spins,ph%n_kpts)
   real(kind=r8), intent(out) :: dm(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: work(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i
   integer(kind=i4) :: max_state
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=*), parameter :: caller = "elsi_compute_dm_elpa_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,factor,ph%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,ph%n_states_solve
      if(occ(i,ph%i_spin,ph%i_kpt) > 0.0_r8) then
         factor(i) = sqrt(occ(i,ph%i_spin,ph%i_kpt))
         max_state = i
      end if
   end do

   work = evec

   do i = 1,ph%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(col_map(i) > 0) then
            work(:,col_map(i)) = work(:,col_map(i))*factor(i)
         end if
      else if(col_map(i) /= 0) then
         work(:,col_map(i)) = 0.0_r8
      end if
   end do

   dm = 0.0_r8

   ! Compute density matrix
   call pdsyrk("U","N",ph%n_basis,max_state,1.0_r8,work,1,1,bh%desc,0.0_r8,dm,&
           1,1,bh%desc)

   call elsi_deallocate(bh,factor,"factor")

   call elsi_set_full_mat(ph,bh,UT_MAT,row_map,col_map,dm)

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished density matrix calculation"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine constructs the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_elpa_real(ph,bh,row_map,col_map,eval,evec,occ,edm,&
              work)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   real(kind=r8), intent(in) :: eval(ph%n_basis)
   real(kind=r8), intent(in) :: evec(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: occ(ph%n_states,ph%n_spins,ph%n_kpts)
   real(kind=r8), intent(out) :: edm(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: work(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i
   integer(kind=i4) :: max_state
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=*), parameter :: caller = "elsi_compute_edm_elpa_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,factor,ph%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,ph%n_states_solve
      factor(i) = -occ(i,ph%i_spin,ph%i_kpt)*eval(i)
      if(factor(i) > 0.0_r8) then
         factor(i) = sqrt(factor(i))
         max_state = i
      else
         factor(i) = 0.0_r8
      end if
   end do

   work = evec

   do i = 1,ph%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(col_map(i) > 0) then
            work(:,col_map(i)) = work(:,col_map(i))*factor(i)
         end if
      else if(col_map(i) /= 0) then
         work(:,col_map(i)) = 0.0_r8
      end if
   end do

   call elsi_deallocate(bh,factor,"factor")

   edm = 0.0_r8

   ! Compute density matrix
   call pdsyrk("U","N",ph%n_basis,max_state,-1.0_r8,work,1,1,bh%desc,0.0_r8,&
           edm,1,1,bh%desc)

   call elsi_set_full_mat(ph,bh,UT_MAT,row_map,col_map,edm)

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished energy density matrix calculation"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine transforms a generalized eigenproblem to standard and returns
!! the Cholesky factor for later use.
!!
subroutine elsi_to_standard_evp_real(ph,bh,row_map,col_map,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: info_str

   character(len=*), parameter :: caller = "elsi_to_standard_evp_real"

   if(ph%n_calls == 1) then
      if(ph%check_sing) then
         call elsi_check_singularity_real(ph,bh,col_map,ovlp,eval,evec)
      end if

      if(ph%n_good == ph%n_basis) then ! Not singular
         call elsi_get_time(t0)

         ph%ovlp_is_sing = .false.

         ! S = (U^T)U, U -> S
         success = elpa_cholesky_real_double(ph%n_basis,ovlp,bh%n_lrow,bh%blk,&
                      bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,.false.)

         if(.not. success) then
            call elsi_stop(bh,"Cholesky failed.",caller)
         end if

         ! U^-1 -> S
         success = elpa_invert_trm_real_double(ph%n_basis,ovlp,bh%n_lrow,&
                      bh%blk,bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,&
                      .false.)

         if(.not. success) then
            call elsi_stop(bh,"Matrix inversion failed.",caller)
         end if

         call elsi_get_time(t1)

         write(info_str,"(2X,A)") "Finished Cholesky decomposition"
         call elsi_say(bh,info_str)
         write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
         call elsi_say(bh,info_str)
      end if
   end if

   call elsi_get_time(t0)

   if(ph%ovlp_is_sing) then ! Use scaled eigenvectors
      ! evec_real used as tmp_real
      ! tmp_real = H_real * S_real
      call pdgemm("N","N",ph%n_basis,ph%n_good,ph%n_basis,1.0_r8,ham,1,1,&
              bh%desc,ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,0.0_r8,evec,1,1,&
              bh%desc)

      ! H_real = (S_real)^T * tmp_real
      call pdgemm("T","N",ph%n_good,ph%n_good,ph%n_basis,1.0_r8,ovlp,1,&
              ph%n_basis-ph%n_good+1,bh%desc,evec,1,1,bh%desc,0.0_r8,ham,1,1,&
              bh%desc)
   else ! Use Cholesky
      success = elpa_mult_at_b_real_double("U","L",ph%n_basis,ph%n_basis,ovlp,&
                   bh%n_lrow,bh%n_lcol,ham,bh%n_lrow,bh%n_lcol,bh%blk,&
                   ph%elpa_comm_row,ph%elpa_comm_col,evec,bh%n_lrow,bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if

      call pdtran(ph%n_basis,ph%n_basis,1.0_r8,evec,1,1,bh%desc,0.0_r8,ham,1,1,&
              bh%desc)

      evec = ham

      success = elpa_mult_at_b_real_double("U","U",ph%n_basis,ph%n_basis,ovlp,&
                   bh%n_lrow,bh%n_lcol,evec,bh%n_lrow,bh%n_lcol,bh%blk,&
                   ph%elpa_comm_row,ph%elpa_comm_col,ham,bh%n_lrow,bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if

      call elsi_set_full_mat(ph,bh,UT_MAT,row_map,col_map,ham)
   end if

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished transformation to standard eigenproblem"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine checks the singularity of overlap matrix by computing all its
!! eigenvalues. On exit, S is not modified if not singular, or is overwritten by
!! scaled eigenvectors if singular, which can be used to transform the
!! generalized eigenproblem to the standard form.
!!
subroutine elsi_check_singularity_real(ph,bh,col_map,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: ev_sqrt
   integer(kind=i4) :: i
   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: copy_real(:,:)

   character(len=*), parameter :: caller = "elsi_check_singularity_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,copy_real,bh%n_lrow,bh%n_lcol,"copy_real",caller)

   ! Use copy_real to store overlap matrix, otherwise it will be destroyed by
   ! eigenvalue calculation
   copy_real = ovlp

   ! Use customized ELPA 2-stage solver to check overlap singularity
   ! Eigenvectors computed only for singular overlap matrix
   success = elpa_check_singularity_real_double(ph%n_basis,ph%n_basis,&
                copy_real,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,bh%n_lcol,&
                ph%elpa_comm_row,ph%elpa_comm_col,bh%comm,ph%sing_tol,ph%n_good)

   if(.not. success) then
      call elsi_stop(bh,"Singularity check failed.",caller)
   end if

   call elsi_deallocate(bh,copy_real,"copy_real")

   ph%n_states_solve = min(ph%n_good,ph%n_states)

   if(ph%n_good < ph%n_basis) then ! Singular
      ph%ovlp_is_sing = .true.

      write(info_str,"(2X,A)") "Overlap matrix is singular"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,E10.2,A,E10.2)")&
         "| Lowest and highest eigenvalues :",eval(ph%n_basis),",",eval(1)
      call elsi_say(bh,info_str)

      if(ph%stop_sing) then
         call elsi_stop(bh,"Overlap matrix is singular.",caller)
      end if

      write(info_str,"(2X,A,I10)") "| Number of basis functions reduced to :",&
         ph%n_good
      call elsi_say(bh,info_str)

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = ph%n_basis-ph%n_good+1,ph%n_basis
         ev_sqrt = sqrt(eval(i))

         if(col_map(i) > 0) then
            ovlp(:,col_map(i)) = evec(:,col_map(i))/ev_sqrt
         end if
      end do
   else ! Nonsingular
      ph%ovlp_is_sing = .false.

      write(info_str,"(2X,A)") "Overlap matrix is not singular"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,E10.2,A,E10.2)")&
         "| Lowest and highest eigenvalues :",eval(ph%n_basis),",",eval(1)
      call elsi_say(bh,info_str)
   end if ! Singular overlap?

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished singularity check of overlap matrix"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine back-transforms eigenvectors in the standard form to the
!! original generalized form.
!!
subroutine elsi_to_original_ev_real(ph,bh,ham,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(out) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: tmp_real(:,:)

   character(len=*), parameter :: caller = "elsi_to_original_ev_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp_real,bh%n_lrow,bh%n_lcol,"tmp_real",caller)

   tmp_real = evec

   if(ph%ovlp_is_sing) then
      call pdgemm("N","N",ph%n_basis,ph%n_states_solve,ph%n_good,1.0_r8,ovlp,1,&
              ph%n_basis-ph%n_good+1,bh%desc,tmp_real,1,1,bh%desc,0.0_r8,evec,&
              1,1,bh%desc)
   else ! Nonsingular, use Cholesky
      call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp,1,1,bh%desc,0.0_r8,ham,1,1,&
              bh%desc)

      success = elpa_mult_at_b_real_double("L","N",ph%n_basis,ph%n_states,ham,&
                   bh%n_lrow,bh%n_lcol,tmp_real,bh%n_lrow,bh%n_lcol,bh%blk,&
                   ph%elpa_comm_row,ph%elpa_comm_col,evec,bh%n_lrow,bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if
   end if

   call elsi_deallocate(bh,tmp_real,"tmp_real")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished back-transformation of eigenvectors"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_elpa_real(ph,bh,row_map,col_map,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: ierr
   logical :: success
   character(len=200) :: info_str

   character(len=*), parameter :: caller = "elsi_solve_elpa_real"

   elpa_print_times = ph%elpa_output

   ! Compute sparsity
   if(ph%n_calls == 1 .and. ph%matrix_format == BLACS_DENSE) then
      call elsi_get_nnz(bh%def0,ham,bh%n_lrow,bh%n_lcol,bh%nnz_l)

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   end if

   ! Transform to standard form
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_standard_evp_real(ph,bh,row_map,col_map,ham,ovlp,eval,evec)
   end if

   call elsi_get_time(t0)

   write(info_str,"(2X,A)") "Starting ELPA eigensolver"
   call elsi_say(bh,info_str)

   ! Solve evp, return eigenvalues and eigenvectors
   if(ph%elpa_solver == 2) then
      success = elpa_solve_evp_real_2stage_double(ph%n_good,ph%n_states_solve,&
                   ham,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,bh%n_lcol,&
                   ph%elpa_comm_row,ph%elpa_comm_col,bh%comm)
   else
      success = elpa_solve_evp_real_1stage_double(ph%n_good,ph%n_states_solve,&
                   ham,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,bh%n_lcol,&
                   ph%elpa_comm_row,ph%elpa_comm_col,bh%comm)
   end if

   if(.not. success) then
      call elsi_stop(bh,"ELPA solver failed.",caller)
   end if

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(ph%n_good < ph%n_basis) then
      eval(ph%n_good+1:ph%n_basis) = eval(ph%n_good)+10.0_r8
   end if

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished solving standard eigenproblem"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

   ! Back-transform eigenvectors
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_original_ev_real(ph,bh,ham,ovlp,evec)
   end if

end subroutine

!>
!! This routine constructs the density matrix.
!!
subroutine elsi_compute_dm_elpa_cmplx(ph,bh,row_map,col_map,evec,occ,dm,work)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   complex(kind=r8), intent(in) :: evec(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: occ(ph%n_states,ph%n_spins,ph%n_kpts)
   complex(kind=r8), intent(out) :: dm(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: work(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i
   integer(kind=i4) :: max_state
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=*), parameter :: caller = "elsi_compute_dm_elpa_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,factor,ph%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,ph%n_states_solve
      if(occ(i,ph%i_spin,ph%i_kpt) > 0.0_r8) then
         factor(i) = sqrt(occ(i,ph%i_spin,ph%i_kpt))
         max_state = i
      end if
   end do

   work = evec

   do i = 1,ph%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(col_map(i) > 0) then
            work(:,col_map(i)) = work(:,col_map(i))*factor(i)
         end if
      else if(col_map(i) /= 0) then
         work(:,col_map(i)) = (0.0_r8,0.0_r8)
      end if
   end do

   dm = (0.0_r8,0.0_r8)

   ! Compute density matrix
   call pzherk("U","N",ph%n_basis,max_state,(1.0_r8,0.0_r8),work,1,1,bh%desc,&
           (0.0_r8,0.0_r8),dm,1,1,bh%desc)

   call elsi_deallocate(bh,factor,"factor")

   call elsi_set_full_mat(ph,bh,UT_MAT,row_map,col_map,dm)

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished density matrix calculation"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine constructs the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_elpa_cmplx(ph,bh,row_map,col_map,eval,evec,occ,edm,&
              work)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   real(kind=r8), intent(in) :: eval(ph%n_basis)
   complex(kind=r8), intent(in) :: evec(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: occ(ph%n_states,ph%n_spins,ph%n_kpts)
   complex(kind=r8), intent(out) :: edm(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: work(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i
   integer(kind=i4) :: max_state
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=*), parameter :: caller = "elsi_compute_edm_elpa_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,factor,ph%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,ph%n_states_solve
      factor(i) = -occ(i,ph%i_spin,ph%i_kpt)*eval(i)
      if(factor(i) > 0.0_r8) then
         factor(i) = sqrt(factor(i))
         max_state = i
      else
         factor(i) = 0.0_r8
      end if
   end do

   work = evec

   do i = 1,ph%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(col_map(i) > 0) then
            work(:,col_map(i)) = work(:,col_map(i))*factor(i)
         end if
      else if(col_map(i) /= 0) then
         work(:,col_map(i)) = (0.0_r8,0.0_r8)
      end if
   end do

   call elsi_deallocate(bh,factor,"factor")

   edm = (0.0_r8,0.0_r8)

   ! Compute density matrix
   call pzherk("U","N",ph%n_basis,max_state,(-1.0_r8,0.0_r8),work,1,1,bh%desc,&
           (0.0_r8,0.0_r8),edm,1,1,bh%desc)

   call elsi_set_full_mat(ph,bh,UT_MAT,row_map,col_map,edm)

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished energy density matrix calculation"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine transforms a generalized eigenproblem to standard and returns
!! the Cholesky factor for later use.
!!
subroutine elsi_to_standard_evp_cmplx(ph,bh,row_map,col_map,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   complex(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: info_str

   character(len=*), parameter :: caller = "elsi_to_standard_evp_cmplx"

   if(ph%n_calls == 1) then
      if(ph%check_sing) then
         call elsi_check_singularity_cmplx(ph,bh,col_map,ovlp,eval,evec)
      end if

      if(ph%n_good == ph%n_basis) then ! Not singular
         call elsi_get_time(t0)

         ph%ovlp_is_sing = .false.

         ! S = (U^T)U, U -> S
         success = elpa_cholesky_complex_double(ph%n_basis,ovlp,bh%n_lrow,&
                      bh%blk,bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,&
                      .false.)

         if(.not. success) then
            call elsi_stop(bh,"Cholesky failed.",caller)
         end if

         ! U^-1 -> S
         success = elpa_invert_trm_complex_double(ph%n_basis,ovlp,bh%n_lrow,&
                      bh%blk,bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,&
                      .false.)

         if(.not. success) then
            call elsi_stop(bh,"Matrix inversion failed.",caller)
         end if

         call elsi_get_time(t1)

         write(info_str,"(2X,A)") "Finished Cholesky decomposition"
         call elsi_say(bh,info_str)
         write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
         call elsi_say(bh,info_str)
      end if
   end if

   call elsi_get_time(t0)

   if(ph%ovlp_is_sing) then ! Use scaled eigenvectors
      ! evec_cmplx used as tmp_cmplx
      ! tmp_cmplx = H_cmplx * S_cmplx
      call pzgemm("N","N",ph%n_basis,ph%n_good,ph%n_basis,(1.0_r8,0.0_r8),ham,&
              1,1,bh%desc,ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,&
              (0.0_r8,0.0_r8),evec,1,1,bh%desc)

      ! H_cmplx = (S_cmplx)^* * tmp_cmplx
      call pzgemm("C","N",ph%n_good,ph%n_good,ph%n_basis,(1.0_r8,0.0_r8),ovlp,&
              1,ph%n_basis-ph%n_good+1,bh%desc,evec,1,1,bh%desc,&
              (0.0_r8,0.0_r8),ham,1,1,bh%desc)
   else ! Use cholesky
      success = elpa_mult_ah_b_complex_double("U","L",ph%n_basis,ph%n_basis,&
                   ovlp,bh%n_lrow,bh%n_lcol,ham,bh%n_lrow,bh%n_lcol,bh%blk,&
                   ph%elpa_comm_row,ph%elpa_comm_col,evec,bh%n_lrow,bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if

      call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),evec,1,1,bh%desc,&
              (0.0_r8,0.0_r8),ham,1,1,bh%desc)

      evec = ham

      success = elpa_mult_ah_b_complex_double("U","U",ph%n_basis,ph%n_basis,&
                   ovlp,bh%n_lrow,bh%n_lcol,evec,bh%n_lrow,bh%n_lcol,bh%blk,&
                   ph%elpa_comm_row,ph%elpa_comm_col,ham,bh%n_lrow,bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if

      call elsi_set_full_mat(ph,bh,UT_MAT,row_map,col_map,ham)
   end if

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished transformation to standard eigenproblem"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine checks the singularity of overlap matrix by computing all its
!! eigenvalues. On exit, S is not modified if not singular, or is overwritten by
!! scaled eigenvectors if singular, which can be used to transform the
!! generalized eigenproblem to the standard form.
!!
subroutine elsi_check_singularity_cmplx(ph,bh,col_map,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: ev_sqrt
   integer(kind=i4) :: i
   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: info_str

   complex(kind=r8), allocatable :: copy_cmplx(:,:)

   character(len=*), parameter :: caller = "elsi_check_singularity_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,copy_cmplx,bh%n_lrow,bh%n_lcol,"copy_cmplx",caller)

   ! Use copy_cmplx to store overlap matrix, otherwise it will be destroyed by
   ! eigenvalue calculation
   copy_cmplx = ovlp

   ! Use customized ELPA 2-stage solver to check overlap singularity
   ! Eigenvectors computed only for singular overlap matrix
   success = elpa_check_singularity_complex_double(ph%n_basis,ph%n_basis,&
                copy_cmplx,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,bh%n_lcol,&
                ph%elpa_comm_row,ph%elpa_comm_col,bh%comm,ph%sing_tol,ph%n_good)

   if(.not. success) then
      call elsi_stop(bh,"Singularity check failed.",caller)
   end if

   call elsi_deallocate(bh,copy_cmplx,"copy_cmplx")

   ph%n_states_solve = min(ph%n_good,ph%n_states)

   if(ph%n_good < ph%n_basis) then ! Singular
      ph%ovlp_is_sing = .true.

      write(info_str,"(2X,A)") "Overlap matrix is singular"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,E10.2,A,E10.2)")&
         "| Lowest and highest eigenvalues :",eval(ph%n_basis),",",eval(1)
      call elsi_say(bh,info_str)

      if(ph%stop_sing) then
         call elsi_stop(bh,"Overlap matrix is singular.",caller)
      end if

      write(info_str,"(2X,A,I10)") "| Number of basis functions reduced to :",&
         ph%n_good
      call elsi_say(bh,info_str)

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = ph%n_basis-ph%n_good+1,ph%n_basis
         ev_sqrt = sqrt(eval(i))

         if(col_map(i) > 0) then
            ovlp(:,col_map(i)) = evec(:,col_map(i))/ev_sqrt
         end if
      end do
   else ! Nonsingular
      ph%ovlp_is_sing = .false.

      write(info_str,"(2X,A)") "Overlap matrix is not singular"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,E10.2,A,E10.2)")&
         "| Lowest and highest eigenvalues :",eval(ph%n_basis),",",eval(1)
      call elsi_say(bh,info_str)
   end if ! Singular overlap?

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished singularity check of overlap matrix"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine back-transforms eigenvectors in the standard form to the
!! original generalized form.
!!
subroutine elsi_to_original_ev_cmplx(ph,bh,ham,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(out) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: info_str

   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character(len=*), parameter :: caller = "elsi_to_original_ev_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp_cmplx,bh%n_lrow,bh%n_lcol,"tmp_cmplx",caller)

   tmp_cmplx = evec

   if(ph%ovlp_is_sing) then
      call pzgemm("N","N",ph%n_basis,ph%n_states_solve,ph%n_good,&
              (1.0_r8,0.0_r8),ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,tmp_cmplx,&
              1,1,bh%desc,(0.0_r8,0.0_r8),evec,1,1,bh%desc)
   else ! Nonsingular, use Cholesky
      call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),ovlp,1,1,bh%desc,&
              (0.0_r8,0.0_r8),ham,1,1,bh%desc)

      success = elpa_mult_ah_b_complex_double("L","N",ph%n_basis,ph%n_states,&
                   ham,bh%n_lrow,bh%n_lcol,tmp_cmplx,bh%n_lrow,bh%n_lcol,&
                   bh%blk,ph%elpa_comm_row,ph%elpa_comm_col,evec,bh%n_lrow,&
                   bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if
   end if

   call elsi_deallocate(bh,tmp_cmplx,"tmp_cmplx")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished back-transformation of eigenvectors"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_elpa_cmplx(ph,bh,row_map,col_map,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   complex(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: ierr
   logical :: success
   character(len=200) :: info_str

   character(len=*), parameter :: caller = "elsi_solve_elpa_cmplx"

   elpa_print_times = ph%elpa_output

   ! Compute sparsity
   if(ph%n_calls == 1 .and. ph%matrix_format == BLACS_DENSE) then
      call elsi_get_nnz(bh%def0,ham,bh%n_lrow,bh%n_lcol,bh%nnz_l)

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   end if

   ! Transform to standard form
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_standard_evp_cmplx(ph,bh,row_map,col_map,ham,ovlp,eval,evec)
   end if

   call elsi_get_time(t0)

   write(info_str,"(2X,A)") "Starting ELPA eigensolver"
   call elsi_say(bh,info_str)

   ! Solve evp, return eigenvalues and eigenvectors
   if(ph%elpa_solver == 2) then
      success = elpa_solve_evp_complex_2stage_double(ph%n_good,&
                   ph%n_states_solve,ham,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,&
                   bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,bh%comm)
   else
      success = elpa_solve_evp_complex_1stage_double(ph%n_good,&
                   ph%n_states_solve,ham,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,&
                   bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,bh%comm)
   end if

   if(.not. success) then
      call elsi_stop(bh,"ELPA solver failed.",caller)
   end if

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(ph%n_good < ph%n_basis) then
      eval(ph%n_good+1:ph%n_basis) = eval(ph%n_good)+10.0_r8
   end if

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished solving standard eigenproblem"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

   ! Back-transform eigenvectors
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_original_ev_cmplx(ph,bh,ham,ovlp,evec)
   end if

end subroutine

!>
!! This routine cleans up ELPA.
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

   ph%elpa_started = .false.

end subroutine

end module ELSI_ELPA
