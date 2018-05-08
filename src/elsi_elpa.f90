! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to ELPA.
!!
module ELSI_ELPA

   use ELSI_CONSTANTS,    only: BLACS_DENSE
   use ELSI_DATATYPE,     only: elsi_handle
   use ELSI_IO,           only: elsi_say,elsi_get_time
   use ELSI_MALLOC,       only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,          only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,&
                                mpi_integer4
   use ELSI_OCC,          only: elsi_compute_mu_and_occ,elsi_compute_entropy
   use ELSI_PRECISION,    only: r8,i4
   use ELSI_UTILS,        only: elsi_get_local_nnz_real,elsi_get_local_nnz_cmplx
   use CHECK_SINGULARITY, only: elpa_check_singularity_real_double,&
                                elpa_check_singularity_complex_double
   use ELPA1,             only: elpa_print_times,&
                                elpa_solve_evp_real_1stage_double,&
                                elpa_solve_evp_complex_1stage_double,&
                                elpa_cholesky_real_double,&
                                elpa_cholesky_complex_double,&
                                elpa_invert_trm_real_double,&
                                elpa_invert_trm_complex_double,&
                                elpa_mult_at_b_real_double,&
                                elpa_mult_ah_b_complex_double
   use ELPA2,             only: elpa_solve_evp_real_2stage_double,&
                                elpa_solve_evp_complex_2stage_double

   implicit none

   private

   public :: elsi_init_elpa
   public :: elsi_set_elpa_default
   public :: elsi_cleanup_elpa
   public :: elsi_compute_occ_elpa
   public :: elsi_compute_dm_elpa_real
   public :: elsi_compute_edm_elpa_real
   public :: elsi_to_standard_evp_real
   public :: elsi_solve_elpa_real
   public :: elsi_compute_dm_elpa_cmplx
   public :: elsi_compute_edm_elpa_cmplx
   public :: elsi_to_standard_evp_cmplx
   public :: elsi_solve_elpa_cmplx

contains

!>
!! This routine initializes ELPA.
!!
subroutine elsi_init_elpa(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   integer(kind=i4)   :: ierr
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_init_elpa"

   if(.not. e_h%elpa_started) then
      call MPI_Comm_split(e_h%mpi_comm,e_h%my_pcol,e_h%my_prow,&
              e_h%mpi_comm_row,ierr)
      call MPI_Comm_split(e_h%mpi_comm,e_h%my_prow,e_h%my_pcol,&
              e_h%mpi_comm_col,ierr)

      call elsi_check_mpi(e_h,"MPI_Comm_split",ierr,caller)

      e_h%elpa_started = .true.

      if(e_h%elpa_gpu) then
         write(info_str,"('  No ELPA GPU acceleration available')")
         call elsi_say(e_h%stdio,info_str)
      endif
   endif

end subroutine

!>
!! This routine computes the chemical potential and occupation numbers.
!!
subroutine elsi_compute_occ_elpa(e_h,eval)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(in)    :: eval(e_h%n_basis)

   integer(kind=i4) :: ierr

   real(kind=r8), allocatable :: tmp_real1(:)
   real(kind=r8), allocatable :: tmp_real2(:,:,:)

   character(len=40), parameter :: caller = "elsi_compute_occ_elpa"

   ! Gather eigenvalues and occupation numbers
   if(e_h%n_elsi_calls == 1) then
      call elsi_allocate(e_h,e_h%eval_all,e_h%n_states,e_h%n_spins,e_h%n_kpts,&
              "eval_all",caller)

      call elsi_allocate(e_h,e_h%occ_num,e_h%n_states,e_h%n_spins,e_h%n_kpts,&
              "occ_num",caller)

      call elsi_allocate(e_h,e_h%k_weight,e_h%n_kpts,"k_weight",caller)

      if(e_h%n_kpts > 1) then
         call elsi_allocate(e_h,tmp_real1,e_h%n_kpts,"tmp_real",caller)

         if(e_h%myid == 0 .and. e_h%i_spin == 1) then
            tmp_real1(e_h%i_kpt) = e_h%i_weight
         endif

         call MPI_Allreduce(tmp_real1,e_h%k_weight,e_h%n_kpts,mpi_real8,&
                 mpi_sum,e_h%mpi_comm_all,ierr)

         call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

         call elsi_deallocate(e_h,tmp_real1,"tmp_real")
      else
         e_h%k_weight = e_h%i_weight
      endif
   endif

   if(e_h%n_spins*e_h%n_kpts > 1) then
      call elsi_allocate(e_h,tmp_real2,e_h%n_states,e_h%n_spins,e_h%n_kpts,&
              "tmp_real",caller)

      if(e_h%myid == 0) then
         tmp_real2(:,e_h%i_spin,e_h%i_kpt) = eval(1:e_h%n_states)
      endif

      call MPI_Allreduce(tmp_real2,e_h%eval_all,&
              e_h%n_states*e_h%n_spins*e_h%n_kpts,mpi_real8,mpi_sum,&
              e_h%mpi_comm_all,ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(e_h,tmp_real2,"tmp_real")
   else
      e_h%eval_all(:,e_h%i_spin,e_h%i_kpt) = eval(1:e_h%n_states)
   endif

   ! Calculate mu and occupation numbers
   call elsi_compute_mu_and_occ(e_h,e_h%n_electrons,e_h%n_states,e_h%n_spins,&
           e_h%n_kpts,e_h%k_weight,e_h%eval_all,e_h%occ_num,e_h%mu)

   ! Calculate entropy
   call elsi_compute_entropy(e_h,e_h%n_states,e_h%n_spins,e_h%n_kpts,&
           e_h%k_weight,e_h%eval_all,e_h%occ_num,e_h%mu,e_h%ts)

end subroutine

!>
!! This routine constructs the density matrix.
!!
subroutine elsi_compute_dm_elpa_real(e_h,evec,dm,work)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(in)    :: evec(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: dm(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: work(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4)   :: i
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: max_state
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=40), parameter :: caller = "elsi_compute_dm_elpa_real"

   call elsi_get_time(t0)

   call elsi_allocate(e_h,factor,e_h%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,e_h%n_states_solve
      if(e_h%occ_num(i,e_h%i_spin,e_h%i_kpt) > 0.0_r8) then
         factor(i) = sqrt(e_h%occ_num(i,e_h%i_spin,e_h%i_kpt))
         max_state = i
      endif
   enddo

   work = evec

   do i = 1,e_h%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(e_h%loc_col(i) > 0) then
            work(:,e_h%loc_col(i)) = work(:,e_h%loc_col(i))*factor(i)
         endif
      elseif(e_h%loc_col(i) /= 0) then
         work(:,e_h%loc_col(i)) = 0.0_r8
      endif
   enddo

   dm = 0.0_r8

   ! Compute density matrix
   call pdsyrk('U','N',e_h%n_basis,max_state,1.0_r8,work,1,1,e_h%sc_desc,&
           0.0_r8,dm,1,1,e_h%sc_desc)

   call elsi_deallocate(e_h,factor,"factor")

   ! Set full matrix from upper triangle
   call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,dm,1,1,e_h%sc_desc,0.0_r8,work,1,&
           1,e_h%sc_desc)

   do i_col = 1,e_h%n_basis-1
      if(e_h%loc_col(i_col) == 0) cycle

      do i_row = i_col+1,e_h%n_basis
         if(e_h%loc_row(i_row) > 0) then
            dm(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
               work(e_h%loc_row(i_row),e_h%loc_col(i_col))
         endif
      enddo
   enddo

   call elsi_get_time(t1)

   write(info_str,"('  Finished density matrix calculation')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

end subroutine

!>
!! This routine constructs the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_elpa_real(e_h,eval,evec,edm,work)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(in)    :: eval(e_h%n_basis)
   real(kind=r8),     intent(in)    :: evec(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: edm(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: work(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4)   :: i
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: max_state
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=40), parameter :: caller = "elsi_compute_edm_elpa_real"

   call elsi_get_time(t0)

   call elsi_allocate(e_h,factor,e_h%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,e_h%n_states_solve
      factor(i) = -1.0_r8*e_h%occ_num(i,e_h%i_spin,e_h%i_kpt)*eval(i)
      if(factor(i) > 0.0_r8) then
         factor(i) = sqrt(factor(i))
         max_state = i
      else
         factor(i) = 0.0_r8
      endif
   enddo

   work = evec

   do i = 1,e_h%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(e_h%loc_col(i) > 0) then
            work(:,e_h%loc_col(i)) = work(:,e_h%loc_col(i))*factor(i)
         endif
      elseif(e_h%loc_col(i) /= 0) then
         work(:,e_h%loc_col(i)) = 0.0_r8
      endif
   enddo

   call elsi_deallocate(e_h,factor,"factor")

   edm = 0.0_r8

   ! Compute density matrix
   call pdsyrk('U','N',e_h%n_basis,max_state,1.0_r8,work,1,1,e_h%sc_desc,&
           0.0_r8,edm,1,1,e_h%sc_desc)

   edm = -1.0_r8*edm

   ! Set full matrix from upper triangle
   call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,edm,1,1,e_h%sc_desc,0.0_r8,work,&
           1,1,e_h%sc_desc)

   do i_col = 1,e_h%n_basis-1
      if(e_h%loc_col(i_col) == 0) cycle

      do i_row = i_col+1,e_h%n_basis
         if(e_h%loc_row(i_row) > 0) then
            edm(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
               work(e_h%loc_row(i_row),e_h%loc_col(i_col))
         endif
      enddo
   enddo

   call elsi_get_time(t1)

   write(info_str,"('  Finished energy density matrix calculation')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

end subroutine

!>
!! This routine transforms a generalized eigenproblem to standard and returns
!! the Cholesky factor for later use.
!!
subroutine elsi_to_standard_evp_real(e_h,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)
   real(kind=r8),     intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   logical            :: success
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_to_standard_evp_real"

   if(e_h%n_elsi_calls == 1) then
      if(e_h%check_sing) then
         call elsi_check_singularity_real(e_h,ovlp,eval,evec)
      endif

      if(e_h%n_nonsing == e_h%n_basis) then ! Not singular
         call elsi_get_time(t0)

         e_h%ovlp_is_sing = .false.

         ! S = (U^T)U, U -> S
         success = elpa_cholesky_real_double(e_h%n_basis,ovlp,e_h%n_lrow,&
                      e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                      .false.)

         if(.not. success) then
            call elsi_stop(e_h,"Cholesky failed.",caller)
         endif

         ! U^-1 -> S
         success = elpa_invert_trm_real_double(e_h%n_basis,ovlp,e_h%n_lrow,&
                      e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                      .false.)

         if(.not. success) then
            call elsi_stop(e_h,"Matrix inversion failed.",caller)
         endif

         call elsi_get_time(t1)

         write(info_str,"('  Finished Cholesky decomposition')")
         call elsi_say(e_h%stdio,info_str)
         write(info_str,"('  | Time :',F10.3,' s')") t1-t0
         call elsi_say(e_h%stdio,info_str)
      endif
   endif ! n_elsi_calls == 1

   call elsi_get_time(t0)

   if(e_h%ovlp_is_sing) then ! Use scaled eigenvectors
      ! evec_real used as tmp_real
      ! tmp_real = H_real * S_real
      call pdgemm('N','N',e_h%n_basis,e_h%n_nonsing,e_h%n_basis,1.0_r8,&
              ham,1,1,e_h%sc_desc,ovlp,1,e_h%n_basis-e_h%n_nonsing+1,&
              e_h%sc_desc,0.0_r8,evec,1,1,e_h%sc_desc)

      ! H_real = (S_real)^T * tmp_real
      call pdgemm('T','N',e_h%n_nonsing,e_h%n_nonsing,e_h%n_basis,1.0_r8,&
              ovlp,1,e_h%n_basis-e_h%n_nonsing+1,e_h%sc_desc,evec,1,1,&
              e_h%sc_desc,0.0_r8,ham,1,1,e_h%sc_desc)
   else ! Use Cholesky
      success = elpa_mult_at_b_real_double('U','L',e_h%n_basis,e_h%n_basis,&
                   ovlp,e_h%n_lrow,e_h%n_lcol,ham,e_h%n_lrow,e_h%n_lcol,&
                   e_h%blk_row,e_h%mpi_comm_row,e_h%mpi_comm_col,evec,&
                   e_h%n_lrow,e_h%n_lcol)

      if(.not. success) then
         call elsi_stop(e_h,"Matrix multiplication failed.",caller)
      endif

      call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,evec,1,1,e_h%sc_desc,0.0_r8,&
              ham,1,1,e_h%sc_desc)

      evec = ham

      success = elpa_mult_at_b_real_double('U','U',e_h%n_basis,e_h%n_basis,&
                   ovlp,e_h%n_lrow,e_h%n_lcol,evec,e_h%n_lrow,e_h%n_lcol,&
                   e_h%blk_row,e_h%mpi_comm_row,e_h%mpi_comm_col,ham,&
                   e_h%n_lrow,e_h%n_lcol)

      if(.not. success) then
         call elsi_stop(e_h,"Matrix multiplication failed.",caller)
      endif

      call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,ham,1,1,e_h%sc_desc,0.0_r8,&
              evec,1,1,e_h%sc_desc)

      ! Set the lower part from the upper
      do i_col = 1,e_h%n_basis-1
         if(e_h%loc_col(i_col) == 0) cycle

         do i_row = i_col+1,e_h%n_basis
            if(e_h%loc_row(i_row) > 0) then
               ham(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                  evec(e_h%loc_row(i_row),e_h%loc_col(i_col))
            endif
         enddo
      enddo
   endif

   call elsi_get_time(t1)

   write(info_str,"('  Finished transformation to standard eigenproblem')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

end subroutine

!>
!! This routine checks the singularity of overlap matrix by computing all its
!! eigenvalues. On exit, S is not modified if not singular, or is overwritten by
!! scaled eigenvectors if singular, which can be esed to transform the
!! generalized eigenproblem to the standard form.
!!
subroutine elsi_check_singularity_real(e_h,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)
   real(kind=r8),     intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)

   real(kind=r8)      :: ev_sqrt
   integer(kind=i4)   :: i
   logical            :: success
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: copy_real(:,:)

   character(len=40), parameter :: caller = "elsi_check_singularity_real"

   call elsi_get_time(t0)

   call elsi_allocate(e_h,copy_real,e_h%n_lrow,e_h%n_lcol,"copy_real",caller)

   ! Use copy_real to store overlap matrix, otherwise it will be destroyed by
   ! eigenvalue calculation
   copy_real = ovlp

   ! Use customized ELPA 2-stage solver to check overlap singularity
   ! Eigenvectors computed only for singular overlap matrix
   success = elpa_check_singularity_real_double(e_h%n_basis,e_h%n_basis,&
                copy_real,e_h%n_lrow,eval,evec,e_h%n_lrow,e_h%blk_row,&
                e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%mpi_comm,&
                e_h%sing_tol,e_h%n_nonsing)

   if(.not. success) then
      call elsi_stop(e_h,"Singularity check failed.",caller)
   endif

   call elsi_deallocate(e_h,copy_real,"copy_real")

   e_h%n_states_solve = min(e_h%n_nonsing,e_h%n_states)

   if(e_h%n_nonsing < e_h%n_basis) then ! Singular
      e_h%ovlp_is_sing = .true.

      write(info_str,"('  Overlap matrix is singular')")
      call elsi_say(e_h%stdio,info_str)
      write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)") eval(1)
      call elsi_say(e_h%stdio,info_str)
      write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
         eval(e_h%n_basis)
      call elsi_say(e_h%stdio,info_str)

      if(e_h%stop_sing) then
         call elsi_stop(e_h,"Overlap matrix is singular.",caller)
      endif

      write(info_str,"('  | Number of basis functions reduced to :',I10)")&
         e_h%n_nonsing
      call elsi_say(e_h%stdio,info_str)

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = e_h%n_basis-e_h%n_nonsing+1,e_h%n_basis
         ev_sqrt = sqrt(eval(i))

         if(e_h%loc_col(i) == 0) cycle

         ovlp(:,e_h%loc_col(i)) = evec(:,e_h%loc_col(i))/ev_sqrt
      enddo
   else ! Nonsingular
      e_h%ovlp_is_sing = .false.

      write(info_str,"('  Overlap matrix is not singular')")
      call elsi_say(e_h%stdio,info_str)
      write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)") eval(1)
      call elsi_say(e_h%stdio,info_str)
      write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
         eval(e_h%n_basis)
      call elsi_say(e_h%stdio,info_str)
   endif ! Singular overlap?

   call elsi_get_time(t1)

   write(info_str,"('  Finished singularity check of overlap matrix')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

end subroutine

!>
!! This routine back-transforms eigenvectors in the standard form to the
!! original generalized form.
!!
subroutine elsi_to_original_ev_real(e_h,ham,ovlp,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)

   logical            :: success
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: tmp_real(:,:)

   character(len=40), parameter :: caller = "elsi_to_original_ev_real"

   call elsi_get_time(t0)

   call elsi_allocate(e_h,tmp_real,e_h%n_lrow,e_h%n_lcol,"tmp_real",caller)
   tmp_real = evec

   if(e_h%ovlp_is_sing) then
      call pdgemm('N','N',e_h%n_basis,e_h%n_states_solve,e_h%n_nonsing,1.0_r8,&
              ovlp,1,e_h%n_basis-e_h%n_nonsing+1,e_h%sc_desc,tmp_real,1,1,&
              e_h%sc_desc,0.0_r8,evec,1,1,e_h%sc_desc)
   else ! Nonsingular, use Cholesky
      call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,ovlp,1,1,e_h%sc_desc,0.0_r8,&
              ham,1,1,e_h%sc_desc)

      success = elpa_mult_at_b_real_double('L','N',e_h%n_basis,e_h%n_states,&
                   ham,e_h%n_lrow,e_h%n_lcol,tmp_real,e_h%n_lrow,e_h%n_lcol,&
                   e_h%blk_row,e_h%mpi_comm_row,e_h%mpi_comm_col,evec,&
                   e_h%n_lrow,e_h%n_lcol)

      if(.not. success) then
         call elsi_stop(e_h,"Matrix multiplication failed.",caller)
      endif
   endif

   call elsi_deallocate(e_h,tmp_real,"tmp_real")

   call elsi_get_time(t1)

   write(info_str,"('  Finished back-transformation of eigenvectors')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_elpa_real(e_h,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)
   real(kind=r8),     intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)

   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   integer(kind=i4)   :: ierr
   logical            :: success
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_solve_elpa_real"

   elpa_print_times = e_h%elpa_output

   ! Compute sparsity
   if(e_h%n_elsi_calls == 1 .and. e_h%matrix_format == BLACS_DENSE) then
      call elsi_get_local_nnz_real(e_h,ham,e_h%n_lrow,e_h%n_lcol,e_h%nnz_l)

      call MPI_Allreduce(e_h%nnz_l,e_h%nnz_g,1,mpi_integer4,mpi_sum,&
              e_h%mpi_comm,ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)
   endif

   ! Transform to standard form
   if(.not. e_h%ovlp_is_unit) then
      call elsi_to_standard_evp_real(e_h,ham,ovlp,eval,evec)
   endif

   call elsi_get_time(t0)

   write(info_str,"('  Starting ELPA eigensolver')")
   call elsi_say(e_h%stdio,info_str)

   ! Solve evp, return eigenvalues and eigenvectors
   if(e_h%elpa_solver == 2) then
      success = elpa_solve_evp_real_2stage_double(e_h%n_nonsing,&
                   e_h%n_states_solve,ham,e_h%n_lrow,eval,evec,e_h%n_lrow,&
                   e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                   e_h%mpi_comm)
   else
      success = elpa_solve_evp_real_1stage_double(e_h%n_nonsing,&
                   e_h%n_states_solve,ham,e_h%n_lrow,eval,evec,e_h%n_lrow,&
                   e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                   e_h%mpi_comm)
   endif

   if(.not. success) then
      call elsi_stop(e_h,"ELPA solver failed.",caller)
   endif

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(e_h%n_nonsing < e_h%n_basis) then
      eval(e_h%n_nonsing+1:e_h%n_basis) = eval(e_h%n_nonsing)+10.0_r8
   endif

   call elsi_get_time(t1)

   write(info_str,"('  Finished solving standard eigenproblem')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

   ! Back-transform eigenvectors
   if(.not. e_h%ovlp_is_unit) then
      call elsi_to_original_ev_real(e_h,ham,ovlp,evec)
   endif

   call MPI_Barrier(e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Barrier",ierr,caller)

end subroutine

!>
!! This routine constructs the density matrix.
!!
subroutine elsi_compute_dm_elpa_cmplx(e_h,evec,dm,work)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(in)    :: evec(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: dm(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: work(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4)   :: i
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: max_state
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=40), parameter :: caller = "elsi_compute_dm_elpa_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(e_h,factor,e_h%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,e_h%n_states_solve
      if(e_h%occ_num(i,e_h%i_spin,e_h%i_kpt) > 0.0_r8) then
         factor(i) = sqrt(e_h%occ_num(i,e_h%i_spin,e_h%i_kpt))
         max_state = i
      endif
   enddo

   work = evec

   do i = 1,e_h%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(e_h%loc_col(i) > 0) then
            work(:,e_h%loc_col(i)) = work(:,e_h%loc_col(i))*factor(i)
         endif
      elseif(e_h%loc_col(i) /= 0) then
         work(:,e_h%loc_col(i)) = (0.0_r8,0.0_r8)
      endif
   enddo

   dm = (0.0_r8,0.0_r8)

   ! Compute density matrix
   call pzherk('U','N',e_h%n_basis,max_state,(1.0_r8,0.0_r8),work,1,1,&
           e_h%sc_desc,(0.0_r8,0.0_r8),dm,1,1,e_h%sc_desc)

   call elsi_deallocate(e_h,factor,"factor")

   ! Set full matrix from upper triangle
   call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),dm,1,1,e_h%sc_desc,&
           (0.0_r8,0.0_r8),work,1,1,e_h%sc_desc)

   do i_col = 1,e_h%n_basis-1
      if(e_h%loc_col(i_col) == 0) cycle

      do i_row = i_col+1,e_h%n_basis
         if(e_h%loc_row(i_row) > 0) then
            dm(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
               work(e_h%loc_row(i_row),e_h%loc_col(i_col))
         endif
      enddo
   enddo

   ! Make diagonal real
   do i_col = 1,e_h%n_basis
      if(e_h%loc_col(i_col) == 0 .or. e_h%loc_row(i_col) == 0) cycle

      dm(e_h%loc_row(i_col),e_h%loc_col(i_col)) = &
         real(dm(e_h%loc_row(i_col),e_h%loc_col(i_col)),kind=r8)
   enddo

   call elsi_get_time(t1)

   write(info_str,"('  Finished density matrix calculation')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

end subroutine

!>
!! This routine constructs the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_elpa_cmplx(e_h,eval,evec,edm,work)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(in)    :: eval(e_h%n_basis)
   complex(kind=r8),  intent(in)    :: evec(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: edm(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: work(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4)   :: i
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: max_state
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=40), parameter :: caller = "elsi_compute_edm_elpa_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(e_h,factor,e_h%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,e_h%n_states_solve
      factor(i) = -1.0_r8*e_h%occ_num(i,e_h%i_spin,e_h%i_kpt)*eval(i)
      if(factor(i) > 0.0_r8) then
         factor(i) = sqrt(factor(i))
         max_state = i
      else
         factor(i) = 0.0_r8
      endif
   enddo

   work = evec

   do i = 1,e_h%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(e_h%loc_col(i) > 0) then
            work(:,e_h%loc_col(i)) = work(:,e_h%loc_col(i))*factor(i)
         endif
      elseif(e_h%loc_col(i) /= 0) then
         work(:,e_h%loc_col(i)) = (0.0_r8,0.0_r8)
      endif
   enddo

   call elsi_deallocate(e_h,factor,"factor")

   edm = (0.0_r8,0.0_r8)

   ! Compute density matrix
   call pzherk('U','N',e_h%n_basis,max_state,(1.0_r8,0.0_r8),work,1,1,&
           e_h%sc_desc,(0.0_r8,0.0_r8),edm,1,1,e_h%sc_desc)

   edm = (-1.0_r8,0.0_r8)*edm

   ! Set full matrix from upper triangle
   call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),edm,1,1,e_h%sc_desc,&
           (0.0_r8,0.0_r8),work,1,1,e_h%sc_desc)

   do i_col = 1,e_h%n_basis-1
      if(e_h%loc_col(i_col) == 0) cycle

      do i_row = i_col+1,e_h%n_basis
         if(e_h%loc_row(i_row) > 0) then
            edm(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
               work(e_h%loc_row(i_row),e_h%loc_col(i_col))
         endif
      enddo
   enddo

   ! Make diagonal real
   do i_col = 1,e_h%n_basis
      if(e_h%loc_col(i_col) == 0 .or. e_h%loc_row(i_col) == 0) cycle

      edm(e_h%loc_row(i_col),e_h%loc_col(i_col)) = &
         real(edm(e_h%loc_row(i_col),e_h%loc_col(i_col)),kind=r8)
   enddo

   call elsi_get_time(t1)

   write(info_str,"('  Finished energy density matrix calculation')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

end subroutine

!>
!! This routine transforms a generalized eigenproblem to standard and returns
!! the Cholesky factor for later use.
!!
subroutine elsi_to_standard_evp_cmplx(e_h,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)
   complex(kind=r8),  intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   logical            :: success
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_to_standard_evp_cmplx"

   if(e_h%n_elsi_calls == 1) then
      if(e_h%check_sing) then
         call elsi_check_singularity_cmplx(e_h,ovlp,eval,evec)
      endif

      if(e_h%n_nonsing == e_h%n_basis) then ! Not singular
         call elsi_get_time(t0)

         e_h%ovlp_is_sing = .false.

         ! S = (U^T)U, U -> S
         success = elpa_cholesky_complex_double(e_h%n_basis,ovlp,e_h%n_lrow,&
                      e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                      .false.)

         if(.not. success) then
            call elsi_stop(e_h,"Cholesky failed.",caller)
         endif

         ! U^-1 -> S
         success = elpa_invert_trm_complex_double(e_h%n_basis,ovlp,e_h%n_lrow,&
                      e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                      .false.)

         if(.not. success) then
            call elsi_stop(e_h,"Matrix inversion failed.",caller)
         endif

         call elsi_get_time(t1)

         write(info_str,"('  Finished Cholesky decomposition')")
         call elsi_say(e_h%stdio,info_str)
         write(info_str,"('  | Time :',F10.3,' s')") t1-t0
         call elsi_say(e_h%stdio,info_str)
      endif
   endif ! n_elsi_calls == 1

   call elsi_get_time(t0)

   if(e_h%ovlp_is_sing) then ! Use scaled eigenvectors
      ! evec_cmplx used as tmp_cmplx
      ! tmp_cmplx = H_cmplx * S_cmplx
      call pzgemm('N','N',e_h%n_basis,e_h%n_nonsing,e_h%n_basis,&
              (1.0_r8,0.0_r8),ham,1,1,e_h%sc_desc,ovlp,1,&
              e_h%n_basis-e_h%n_nonsing+1,e_h%sc_desc,(0.0_r8,0.0_r8),evec,1,1,&
              e_h%sc_desc)

      ! H_cmplx = (S_cmplx)^* * tmp_cmplx
      call pzgemm('C','N',e_h%n_nonsing,e_h%n_nonsing,e_h%n_basis,&
              (1.0_r8,0.0_r8),ovlp,1,e_h%n_basis-e_h%n_nonsing+1,e_h%sc_desc,&
              evec,1,1,e_h%sc_desc,(0.0_r8,0.0_r8),ham,1,1,e_h%sc_desc)
   else ! Use cholesky
      success = elpa_mult_ah_b_complex_double('U','L',e_h%n_basis,e_h%n_basis,&
                   ovlp,e_h%n_lrow,e_h%n_lcol,ham,e_h%n_lrow,e_h%n_lcol,&
                   e_h%blk_row,e_h%mpi_comm_row,e_h%mpi_comm_col,evec,&
                   e_h%n_lrow,e_h%n_lcol)

      if(.not. success) then
         call elsi_stop(e_h,"Matrix multiplication failed.",caller)
      endif

      call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),evec,1,1,&
              e_h%sc_desc,(0.0_r8,0.0_r8),ham,1,1,e_h%sc_desc)

      evec = ham

      success = elpa_mult_ah_b_complex_double('U','U',e_h%n_basis,e_h%n_basis,&
                   ovlp,e_h%n_lrow,e_h%n_lcol,evec,e_h%n_lrow,e_h%n_lcol,&
                   e_h%blk_row,e_h%mpi_comm_row,e_h%mpi_comm_col,ham,&
                   e_h%n_lrow,e_h%n_lcol)

      if(.not. success) then
         call elsi_stop(e_h,"Matrix multiplication failed.",caller)
      endif

      call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),ham,1,1,e_h%sc_desc,&
              (0.0_r8,0.0_r8),evec,1,1,e_h%sc_desc)

      ! Set the lower part from the upper
      do i_col = 1,e_h%n_basis-1
         if(e_h%loc_col(i_col) == 0) cycle

         do i_row = i_col+1,e_h%n_basis
            if(e_h%loc_row(i_row) > 0) then
               ham(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                  evec(e_h%loc_row(i_row),e_h%loc_col(i_col))
            endif
         enddo
      enddo

      do i_col=1,e_h%n_basis
         if(e_h%loc_col(i_col) == 0 .or. e_h%loc_row(i_col) == 0) cycle

         ham(e_h%loc_row(i_col),e_h%loc_col(i_col)) = &
            real(ham(e_h%loc_row(i_col),e_h%loc_col(i_col)),kind=r8)
      enddo
   endif

   call elsi_get_time(t1)

   write(info_str,"('  Finished transformation to standard eigenproblem')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

end subroutine

!>
!! This routine checks the singularity of overlap matrix by computing all its
!! eigenvalues. On exit, S is not modified if not singular, or is overwritten by
!! scaled eigenvectors if singular, which can be esed to transform the
!! generalized eigenproblem to the standard form.
!!
subroutine elsi_check_singularity_cmplx(e_h,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)
   complex(kind=r8),  intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)

   real(kind=r8)      :: ev_sqrt
   integer(kind=i4)   :: i
   logical            :: success
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   complex(kind=r8), allocatable :: copy_cmplx(:,:)

   character(len=40), parameter :: caller = "elsi_check_singularity_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(e_h,copy_cmplx,e_h%n_lrow,e_h%n_lcol,"copy_cmplx",caller)

   ! Use copy_cmplx to store overlap matrix, otherwise it will be destroyed by
   ! eigenvalue calculation
   copy_cmplx = ovlp

   ! Use customized ELPA 2-stage solver to check overlap singularity
   ! Eigenvectors computed only for singular overlap matrix
   success = elpa_check_singularity_complex_double(e_h%n_basis,e_h%n_basis,&
                copy_cmplx,e_h%n_lrow,eval,evec,e_h%n_lrow,e_h%blk_row,&
                e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%mpi_comm,&
                e_h%sing_tol,e_h%n_nonsing)

   if(.not. success) then
      call elsi_stop(e_h,"Singularity check failed.",caller)
   endif

   call elsi_deallocate(e_h,copy_cmplx,"copy_cmplx")

   e_h%n_states_solve = min(e_h%n_nonsing,e_h%n_states)

   if(e_h%n_nonsing < e_h%n_basis) then ! Singular
      e_h%ovlp_is_sing = .true.

      write(info_str,"('  Overlap matrix is singular')")
      call elsi_say(e_h%stdio,info_str)
      write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)") eval(1)
      call elsi_say(e_h%stdio,info_str)
      write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
         eval(e_h%n_basis)
      call elsi_say(e_h%stdio,info_str)

      if(e_h%stop_sing) then
         call elsi_stop(e_h,"Overlap matrix is singular.",caller)
      endif

      write(info_str,"('  | Number of basis functions reduced to :',I10)")&
         e_h%n_nonsing
      call elsi_say(e_h%stdio,info_str)

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = e_h%n_basis-e_h%n_nonsing+1,e_h%n_basis
         ev_sqrt = sqrt(eval(i))

         if(e_h%loc_col(i) == 0) cycle

         ovlp(:,e_h%loc_col(i)) = evec(:,e_h%loc_col(i))/ev_sqrt
      enddo
   else ! Nonsingular
      e_h%ovlp_is_sing = .false.

      write(info_str,"('  Overlap matrix is not singular')")
      call elsi_say(e_h%stdio,info_str)
      write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)") eval(1)
      call elsi_say(e_h%stdio,info_str)
      write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
         eval(e_h%n_basis)
      call elsi_say(e_h%stdio,info_str)
   endif ! Singular overlap?

   call elsi_get_time(t1)

   write(info_str,"('  Finished singularity check of overlap matrix')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

end subroutine

!>
!! This routine back-transforms eigenvectors in the standard form to the
!! original generalized form.
!!
subroutine elsi_to_original_ev_cmplx(e_h,ham,ovlp,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)

   logical            :: success
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character(len=40), parameter :: caller = "elsi_to_original_ev_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(e_h,tmp_cmplx,e_h%n_lrow,e_h%n_lcol,"tmp_cmplx",caller)
   tmp_cmplx = evec

   if(e_h%ovlp_is_sing) then
      call pzgemm('N','N',e_h%n_basis,e_h%n_states_solve,e_h%n_nonsing,&
              (1.0_r8,0.0_r8),ovlp,1,e_h%n_basis-e_h%n_nonsing+1,e_h%sc_desc,&
              tmp_cmplx,1,1,e_h%sc_desc,(0.0_r8,0.0_r8),evec,1,1,e_h%sc_desc)
   else ! Nonsingular, use Cholesky
      call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),ovlp,1,1,&
              e_h%sc_desc,(0.0_r8,0.0_r8),ham,1,1,e_h%sc_desc)

      success = elpa_mult_ah_b_complex_double('L','N',e_h%n_basis,e_h%n_states,&
                   ham,e_h%n_lrow,e_h%n_lcol,tmp_cmplx,e_h%n_lrow,e_h%n_lcol,&
                   e_h%blk_row,e_h%mpi_comm_row,e_h%mpi_comm_col,evec,&
                   e_h%n_lrow,e_h%n_lcol)

      if(.not. success) then
         call elsi_stop(e_h,"Matrix multiplication failed.",caller)
      endif
   endif

   call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")

   call elsi_get_time(t1)

   write(info_str,"('  Finished back-transformation of eigenvectors')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_elpa_cmplx(e_h,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)
   complex(kind=r8),  intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)

   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   integer(kind=i4)   :: ierr
   logical            :: success
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_solve_elpa_cmplx"

   elpa_print_times = e_h%elpa_output

   ! Compute sparsity
   if(e_h%n_elsi_calls == 1 .and. e_h%matrix_format == BLACS_DENSE) then
      call elsi_get_local_nnz_cmplx(e_h,ham,e_h%n_lrow,e_h%n_lcol,e_h%nnz_l)

      call MPI_Allreduce(e_h%nnz_l,e_h%nnz_g,1,mpi_integer4,mpi_sum,&
              e_h%mpi_comm,ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)
   endif

   ! Transform to standard form
   if(.not. e_h%ovlp_is_unit) then
      call elsi_to_standard_evp_cmplx(e_h,ham,ovlp,eval,evec)
   endif

   call elsi_get_time(t0)

   write(info_str,"('  Starting ELPA eigensolver')")
   call elsi_say(e_h%stdio,info_str)

   ! Solve evp, return eigenvalues and eigenvectors
   if(e_h%elpa_solver == 2) then
      success = elpa_solve_evp_complex_2stage_double(e_h%n_nonsing,&
                   e_h%n_states_solve,ham,e_h%n_lrow,eval,evec,e_h%n_lrow,&
                   e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                   e_h%mpi_comm)
   else
      success = elpa_solve_evp_complex_1stage_double(e_h%n_nonsing,&
                   e_h%n_states_solve,ham,e_h%n_lrow,eval,evec,e_h%n_lrow,&
                   e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                   e_h%mpi_comm)
   endif

   if(.not. success) then
      call elsi_stop(e_h,"ELPA solver failed.",caller)
   endif

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(e_h%n_nonsing < e_h%n_basis) then
      eval(e_h%n_nonsing+1:e_h%n_basis) = eval(e_h%n_nonsing)+10.0_r8
   endif

   call elsi_get_time(t1)

   write(info_str,"('  Finished solving standard eigenproblem')")
   call elsi_say(e_h%stdio,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h%stdio,info_str)

   ! Back-transform eigenvectors
   if(.not. e_h%ovlp_is_unit) then
      call elsi_to_original_ev_cmplx(e_h,ham,ovlp,evec)
   endif

   call MPI_Barrier(e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Barrier",ierr,caller)

end subroutine

!>
!! This routine sets default ELPA parameters.
!!
subroutine elsi_set_elpa_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   character(len=40), parameter :: caller = "elsi_set_elpa_default"

   ! ELPA solver
   e_h%elpa_solver = 2

   ! Number of SCF steps with single precision
   e_h%elpa_n_single = 0

   ! ELPA output?
   e_h%elpa_output = .false.

   ! Use GPU acceleration?
   e_h%elpa_gpu = .false.

   ! Use GPU kernels?
   e_h%elpa_gpu_kernels = .false.

end subroutine

!>
!! This routine cleans up ELPA.
!!
subroutine elsi_cleanup_elpa(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   integer(kind=i4) :: ierr

   character(len=40), parameter :: caller = "elsi_cleanup_elpa"

   ! ELPA
   if(allocated(e_h%ham_real_elpa)) then
      call elsi_deallocate(e_h,e_h%ham_real_elpa,"ham_real_elpa")
   endif
   if(allocated(e_h%ham_cmplx_elpa)) then
      call elsi_deallocate(e_h,e_h%ham_cmplx_elpa,"ham_cmplx_elpa")
   endif
   if(allocated(e_h%ovlp_real_elpa)) then
      call elsi_deallocate(e_h,e_h%ovlp_real_elpa,"ovlp_real_elpa")
   endif
   if(allocated(e_h%ovlp_cmplx_elpa)) then
      call elsi_deallocate(e_h,e_h%ovlp_cmplx_elpa,"ovlp_cmplx_elpa")
   endif
   if(allocated(e_h%eval_elpa)) then
      call elsi_deallocate(e_h,e_h%eval_elpa,"eval_elpa")
   endif
   if(allocated(e_h%evec_real_elpa)) then
      call elsi_deallocate(e_h,e_h%evec_real_elpa,"evec_real_elpa")
   endif
   if(allocated(e_h%evec_cmplx_elpa)) then
      call elsi_deallocate(e_h,e_h%evec_cmplx_elpa,"evec_cmplx_elpa")
   endif
   if(allocated(e_h%dm_real_elpa)) then
      call elsi_deallocate(e_h,e_h%dm_real_elpa,"dm_real_elpa")
   endif
   if(allocated(e_h%dm_cmplx_elpa)) then
      call elsi_deallocate(e_h,e_h%dm_cmplx_elpa,"dm_cmplx_elpa")
   endif
   if(allocated(e_h%occ_num)) then
      call elsi_deallocate(e_h,e_h%occ_num,"occ_num")
   endif
   if(allocated(e_h%k_weight)) then
      call elsi_deallocate(e_h,e_h%k_weight,"k_weight")
   endif
   if(allocated(e_h%eval_all)) then
      call elsi_deallocate(e_h,e_h%eval_all,"eval_all")
   endif

   if(e_h%elpa_started) then
      call MPI_Comm_free(e_h%mpi_comm_row,ierr)
      call MPI_Comm_free(e_h%mpi_comm_col,ierr)
   endif

   e_h%elpa_started = .false.

end subroutine

end module ELSI_ELPA