! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to ELPA-AEO.
!!
module ELSI_ELPA

   use ELSI_CONSTANTS, only: BLACS_DENSE
   use ELSI_DATATYPE,  only: elsi_handle
   use ELSI_IO,        only: elsi_say,elsi_get_time
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,&
                             mpi_integer4
   use ELSI_OCC,       only: elsi_compute_mu_and_occ,elsi_compute_entropy
   use ELSI_PRECISION, only: r4,r8,i4
   use ELSI_UTILS,     only: elsi_get_local_nnz_real,elsi_get_local_nnz_cmplx
   use ELPA,           only: elpa_t,elpa_init,elpa_uninit,elpa_allocate,&
                             elpa_deallocate,elpa_autotune_deallocate,&
                             ELPA_SOLVER_1STAGE,ELPA_SOLVER_2STAGE,&
                             ELPA_AUTOTUNE_FAST,ELPA_AUTOTUNE_DOMAIN_REAL,&
                             ELPA_AUTOTUNE_DOMAIN_COMPLEX,ELPA_2STAGE_REAL_GPU,&
                             ELPA_2STAGE_COMPLEX_GPU,ELPA_2STAGE_REAL_DEFAULT,&
                             ELPA_2STAGE_COMPLEX_DEFAULT

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

   interface elsi_elpa_evec
      module procedure elsi_elpa_evec_cmplx,&
                       elsi_elpa_evec_real
   end interface

   interface elsi_elpa_mult
      module procedure elsi_elpa_mult_cmplx,&
                       elsi_elpa_mult_real
   end interface

   interface elsi_elpa_chol
      module procedure elsi_elpa_chol_cmplx,&
                       elsi_elpa_chol_real
   end interface

   interface elsi_elpa_invt
      module procedure elsi_elpa_invt_cmplx,&
                       elsi_elpa_invt_real
   end interface

contains

!>
!! This routine initializes ELPA.
!!
subroutine elsi_init_elpa(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   integer(kind=i4) :: ierr

   character(len=40), parameter :: caller = "elsi_init_elpa"

   if(.not. e_h%elpa_started) then
      ierr = elpa_init(20171201)

      if(ierr /= 0) then
         call elsi_stop(e_h,"Initialization failed.",caller)
      endif

      call MPI_Comm_split(e_h%mpi_comm,e_h%my_pcol,e_h%my_prow,&
              e_h%mpi_comm_row,ierr)
      call MPI_Comm_split(e_h%mpi_comm,e_h%my_prow,e_h%my_pcol,&
              e_h%mpi_comm_col,ierr)

      call elsi_check_mpi(e_h,"MPI_Comm_split",ierr,caller)

      e_h%elpa_started = .true.
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

         if(e_h%myid == 0) then
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
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_to_standard_evp_real"

   call elsi_get_time(t0)

   if(e_h%n_elsi_calls == 1) then
      if(e_h%check_sing) then
         call elsi_check_singularity_real(e_h,ovlp,eval,evec)
      endif

      if(e_h%n_nonsing == e_h%n_basis) then ! Not singular

         e_h%ovlp_is_sing = .false.

         ! S = (U^T)U, U -> S
         call elsi_elpa_chol(e_h,ovlp)

         ! U^-1 -> S
         call elsi_elpa_invt(e_h,ovlp)

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
      call elsi_elpa_mult(e_h,'U','L',ovlp,ham,evec)

      call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,evec,1,1,e_h%sc_desc,0.0_r8,&
              ham,1,1,e_h%sc_desc)

      evec = ham

      call elsi_elpa_mult(e_h,'U','U',ovlp,evec,ham)

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

   integer(kind=i4)   :: i
   real(kind=r8)      :: ev_sqrt
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_check_singularity_real"

   call elsi_get_time(t0)

   ! Use ELPA to check overlap singularity
   call elsi_elpa_evec(e_h,ovlp,eval,evec,.true.)

   do i = 1,e_h%n_basis
      if(eval(i) < e_h%sing_tol) then
         e_h%n_nonsing = e_h%n_nonsing-1
      endif
   enddo

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

      call elsi_elpa_mult(e_h,'L','N',ham,tmp_real,evec)
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

   integer(kind=i4)   :: ierr
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_solve_elpa_real"

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

   ! Solve evp, return eigenvalues and eigenvectors
   call elsi_elpa_evec(e_h,ham,eval,evec,.false.)

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
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_to_standard_evp_cmplx"

   call elsi_get_time(t0)

   if(e_h%n_elsi_calls == 1) then
      if(e_h%check_sing) then
         call elsi_check_singularity_cmplx(e_h,ovlp,eval,evec)
      endif

      if(e_h%n_nonsing == e_h%n_basis) then ! Not singular
         e_h%ovlp_is_sing = .false.

         ! S = (U^T)U, U -> S
         call elsi_elpa_chol(e_h,ovlp)

         ! U^-1 -> S
         call elsi_elpa_invt(e_h,ovlp)

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
      call elsi_elpa_mult(e_h,'U','L',ovlp,ham,evec)

      call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),evec,1,1,&
              e_h%sc_desc,(0.0_r8,0.0_r8),ham,1,1,e_h%sc_desc)

      evec = ham

      call elsi_elpa_mult(e_h,'U','U',ovlp,evec,ham)

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

   integer(kind=i4)   :: i
   real(kind=r8)      :: ev_sqrt
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_check_singularity_cmplx"

   call elsi_get_time(t0)

   ! Use ELPA to check overlap singularity
   call elsi_elpa_evec(e_h,ovlp,eval,evec,.true.)

   do i = 1,e_h%n_basis
      if(eval(i) < e_h%sing_tol) then
         e_h%n_nonsing = e_h%n_nonsing-1
      endif
   enddo

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

      call elsi_elpa_mult(e_h,'L','N',ham,tmp_cmplx,evec)
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

   integer(kind=i4)   :: ierr
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_solve_elpa_cmplx"

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

   ! Solve evp, return eigenvalues and eigenvectors
   call elsi_elpa_evec(e_h,ham,eval,evec,.false.)

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

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character(len=40), parameter :: caller = "elsi_set_elpa_default"

   ! ELPA solver
   e_h%elpa_solver = 2

   ! How many single precision steps?
   e_h%elpa_n_single = 0

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
      nullify(e_h%elpa_main)

      call elpa_uninit()

      call MPI_Comm_free(e_h%mpi_comm_row,ierr)
      call MPI_Comm_free(e_h%mpi_comm_col,ierr)
   endif

   e_h%elpa_started = .false.

end subroutine

!>
!! This routine calls ELPA real eigensolver.
!!
subroutine elsi_elpa_evec_real(e_h,ham,eval,evec,sing_check)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)
   real(kind=r8),     intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)
   logical,           intent(in)    :: sing_check

   integer(kind=i4)   :: ierr
   character(len=200) :: info_str

   class(elpa_t), pointer :: elpa_main

   real(kind=r8), allocatable :: copy_ham(:,:)
   real(kind=r4), allocatable :: copy_ham_single(:,:)
   real(kind=r4), allocatable :: eval_single(:)
   real(kind=r4), allocatable :: evec_single(:,:)

   character(len=40), parameter :: caller = "elsi_elpa_evec_real"

   if(sing_check) then
      call elsi_elpa_setup(e_h,elpa_main,e_h%n_basis,e_h%n_basis)
      ! TODO: Define ill-conditioning tolerance (not yet available in ELPA)
!      call elpa_main%set("check_pd",1,ierr)
      call elpa_main%set("solver",ELPA_SOLVER_2STAGE,ierr)

      call elsi_allocate(e_h,copy_ham,e_h%n_lrow,e_h%n_lcol,"copy_ham",caller)

      copy_ham = ham

      call elpa_main%eigenvectors(copy_ham,eval,evec,ierr)
      call elpa_deallocate(elpa_main)

      nullify(elpa_main)

      call elsi_deallocate(e_h,copy_ham,"copy_ham")
   else
      call elsi_elpa_setup(e_h,e_h%elpa_main,e_h%n_nonsing,e_h%n_states_solve)
!      call elsi_elpa_autotuning(e_h,"real")

      elpa_main => e_h%elpa_main

      if(e_h%n_elsi_calls <= e_h%elpa_n_single) then
         write(info_str,"('  Starting ELPA eigensolver (single precision)')")
         call elsi_say(e_h%stdio,info_str)

         call elsi_allocate(e_h,eval_single,e_h%n_basis,"eval_single",caller)
         call elsi_allocate(e_h,evec_single,e_h%n_lrow,e_h%n_lcol,&
                 "evec_single",caller)
         call elsi_allocate(e_h,copy_ham,e_h%n_lrow,e_h%n_lcol,&
                 "copy_ham_single",caller)

         copy_ham_single=real(ham,kind=r4)

         call elpa_main%eigenvectors(copy_ham_single,eval_single,evec_single,&
                 ierr)

         eval = real(eval_single,kind=r8)
         evec = real(evec_single,kind=r8)

         call elsi_deallocate(e_h,eval_single,"eval_single")
         call elsi_deallocate(e_h,evec_single,"evec_single")
         call elsi_deallocate(e_h,copy_ham_single,"copy_ham_single")
      else
         write(info_str,"('  Starting ELPA eigensolver')")
         call elsi_say(e_h%stdio,info_str)

         call elpa_main%eigenvectors(ham,eval,evec,ierr)
      endif
   endif

   if(ierr /= 0) then
      call elsi_stop(e_h,"ELPA eigensolver failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA complex eigensolver.
!!
subroutine elsi_elpa_evec_cmplx(e_h,ham,eval,evec,sing_check)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)
   complex(kind=r8),  intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)
   logical,           intent(in)    :: sing_check

   integer(kind=i4)   :: ierr
   character(len=200) :: info_str

   class(elpa_t), pointer :: elpa_main

   complex(kind=r8), allocatable :: copy_ham(:,:)
   complex(kind=r4), allocatable :: copy_ham_single(:,:)
   real(kind=r4),    allocatable :: eval_single(:)
   complex(kind=r4), allocatable :: evec_single(:,:)

   character(len=40), parameter :: caller = "elsi_elpa_evec_cmplx"

   if(sing_check) then
      call elsi_elpa_setup(e_h,elpa_main,e_h%n_basis,e_h%n_basis)
      ! TODO: Define ill-conditioning tolerance (not yet available in ELPA)
!      call elpa_main%set("check_pd",1,ierr)
      call elpa_main%set("solver",ELPA_SOLVER_2STAGE,ierr)

      call elsi_allocate(e_h,copy_ham,e_h%n_lrow,e_h%n_lcol,"copy_ham",caller)

      copy_ham = ham

      call elpa_main%eigenvectors(copy_ham,eval,evec,ierr)
      call elpa_deallocate(elpa_main)

      nullify(elpa_main)

      call elsi_deallocate(e_h,copy_ham,"copy_ham")
   else
      call elsi_elpa_setup(e_h,e_h%elpa_main,e_h%n_nonsing,e_h%n_states_solve)
!      call elsi_elpa_autotuning(e_h,"complex")

      elpa_main => e_h%elpa_main

      if(e_h%n_elsi_calls <= e_h%elpa_n_single) then
         write(info_str,"('  Starting ELPA eigensolver (single precision)')")
         call elsi_say(e_h%stdio,info_str)

         call elsi_allocate(e_h,eval_single,e_h%n_basis,"eval_single",caller)
         call elsi_allocate(e_h,evec_single,e_h%n_lrow,e_h%n_lcol,&
                 "evec_single",caller)
         call elsi_allocate(e_h,copy_ham,e_h%n_lrow,e_h%n_lcol,&
                 "copy_ham_single",caller)

         copy_ham_single=cmplx(ham,kind=r4)

         call elpa_main%eigenvectors(copy_ham_single,eval_single,evec_single,&
                 ierr)

         eval = real(eval_single,kind=r8)
         evec = cmplx(evec_single,kind=r8)

         call elsi_deallocate(e_h,eval_single,"eval_single")
         call elsi_deallocate(e_h,evec_single,"evec_single")
         call elsi_deallocate(e_h,copy_ham_single,"copy_ham_single")
      else
         write(info_str,"('  Starting ELPA eigensolver')")
         call elsi_say(e_h%stdio,info_str)

         call elpa_main%eigenvectors(ham,eval,evec,ierr)
      endif
   endif

   if(ierr /= 0) then
      call elsi_stop(e_h,"ELPA eigensolver failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA real hermitian_multiply.
!!
subroutine elsi_elpa_mult_real(e_h,uplo,uplo2,a,b,c)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   character,         intent(in)    :: uplo
   character,         intent(in)    :: uplo2
   real(kind=r8),     intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: b(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: c(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_mult_real"

   call elsi_elpa_setup(e_h,elpa_main,e_h%n_basis,e_h%n_basis)
   call elpa_main%hermitian_multiply(uplo,uplo2,e_h%n_basis,a,b,e_h%n_lrow,&
           e_h%n_lcol,c,e_h%n_lrow,e_h%n_lcol,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(e_h,"Matrix multiplication failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA complex hermitian_multiply.
!!
subroutine elsi_elpa_mult_cmplx(e_h,uplo,uplo2,a,b,c)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle
   character,         intent(in)    :: uplo
   character,         intent(in)    :: uplo2
   complex(kind=r8),  intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: b(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: c(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_mult_cmplx"

   call elsi_elpa_setup(e_h,elpa_main,e_h%n_basis,e_h%n_basis)
   call elpa_main%hermitian_multiply(uplo,uplo2,e_h%n_basis,a,b,e_h%n_lrow,&
           e_h%n_lcol,c,e_h%n_lrow,e_h%n_lcol,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(e_h,"Matrix multiplication failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA real cholesky.
!!
subroutine elsi_elpa_chol_real(e_h,a)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_chol_real"

   call elsi_elpa_setup(e_h,elpa_main,e_h%n_basis,e_h%n_basis)
   call elpa_main%cholesky(a,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(e_h,"Cholesky decomposition failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA complex cholesky.
!!
subroutine elsi_elpa_chol_cmplx(e_h,a)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_chol_cmplx"

   call elsi_elpa_setup(e_h,elpa_main,e_h%n_basis,e_h%n_basis)
   call elpa_main%cholesky(a,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(e_h,"Cholesky decomposition failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA real invert_triangular.
!!
subroutine elsi_elpa_invt_real(e_h,a)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_invt_real"

   call elsi_elpa_setup(e_h,elpa_main,e_h%n_nonsing,e_h%n_nonsing)
   call elpa_main%invert_triangular(a,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(e_h,"Matrix inversion failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA complex invert_triangular.
!!
subroutine elsi_elpa_invt_cmplx(e_h,a)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_invt_cmplx"

   call elsi_elpa_setup(e_h,elpa_main,e_h%n_nonsing,e_h%n_nonsing)
   call elpa_main%invert_triangular(a,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(e_h,"Matrix inversion failed.",caller)
   endif

end subroutine

!>
!! This routine sets ELPA-AEO parameters.
!!
subroutine elsi_elpa_setup(e_h,elpa_i,na,nev)

   implicit none

   type(elsi_handle), intent(inout)          :: e_h
   class(elpa_t),     intent(inout), pointer :: elpa_i
   integer,           intent(in)             :: na
   integer,           intent(in)             :: nev

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: ierr2
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_elpa_setup"

   if(e_h%n_elsi_calls == e_h%elpa_n_single .and. associated(elpa_i)) then
      call elpa_deallocate(elpa_i)

      nullify(elpa_i)
   endif

   if(.not. associated(elpa_i)) then
      elpa_i => elpa_allocate()

      call elpa_i%set("na",na,ierr)
      call elpa_i%set("nev",nev,ierr)
      call elpa_i%set("local_nrows",e_h%n_lrow,ierr)
      call elpa_i%set("local_ncols",e_h%n_lcol,ierr)
      call elpa_i%set("nblk",e_h%blk_row,ierr)
      call elpa_i%set("mpi_comm_parent",e_h%mpi_comm,ierr)
      call elpa_i%set("process_row",e_h%my_prow,ierr)
      call elpa_i%set("process_col",e_h%my_pcol,ierr)

      ierr = elpa_i%setup()

      if(ierr /= 0) then
         call elsi_stop(e_h,"ELPA setup failed.",caller)
      endif

      if(e_h%elpa_solver == 1) then
         call elpa_i%set("solver",ELPA_SOLVER_1STAGE,ierr)
      else
         call elpa_i%set("solver",ELPA_SOLVER_2STAGE,ierr)
      endif

      ! Try to enable ELPA GPU acceleration
      if(e_h%elpa_gpu) then
         call elpa_i%set("gpu",1,ierr)

         if(ierr /= 0) then ! Failed
            call elpa_i%set("gpu",0,ierr)

            e_h%elpa_gpu         = .false.
            e_h%elpa_gpu_kernels = .false.

            write(info_str,"('  No ELPA GPU acceleration available')")
            call elsi_say(e_h%stdio,info_str)
         else
            write(info_str,"('  ELPA GPU acceleration activated')")
            call elsi_say(e_h%stdio,info_str)
         endif
      endif

      ! Try to enable ELPA2 GPU kernels
      if(e_h%elpa_gpu_kernels) then
         if(e_h%elpa_solver == 2) then
            call elpa_i%set("real_kernel",ELPA_2STAGE_REAL_GPU,ierr)
            call elpa_i%set("complex_kernel",ELPA_2STAGE_COMPLEX_GPU,ierr2)

            if(ierr /= 0 .or. ierr2 /= 0) then
               call elpa_i%set("real_kernel",ELPA_2STAGE_REAL_DEFAULT,ierr)
               call elpa_i%set("complex_kernel",ELPA_2STAGE_COMPLEX_DEFAULT,&
                       ierr)

               e_h%elpa_gpu_kernels = .false.

               write(info_str,"('  ELPA GPU kernels not available')")
               call elsi_say(e_h%stdio,info_str)
            else
               write(info_str,"('  ELPA GPU kernels will be used')")
               call elsi_say(e_h%stdio,info_str)
            endif
         else ! GPU kernels currently only make sense with ELPA2
            write(info_str,"('  No GPU kernels available with 1-stage ELPA')")
            call elsi_say(e_h%stdio,info_str)
         endif
      endif
   endif

end subroutine

!>
!! This routine sets up ELPA AEO auto-tuning.
!!
subroutine elsi_elpa_autotuning(e_h,real_cmplx)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   character(len=*),  intent(in)    :: real_cmplx

   integer(kind=i4) :: ierr

   character(len=40), parameter :: caller = "elsi_elpa_autotuning"

   if(e_h%n_elsi_calls == e_h%elpa_n_single) then
      call elpa_autotune_deallocate(e_h%elpa_tune)
   endif

   if(e_h%n_elsi_calls == 1 .or. e_h%n_elsi_calls == e_h%elpa_n_single) then
      if(real_cmplx == "complex") then
         e_h%elpa_tune => e_h%elpa_main%autotune_setup(ELPA_AUTOTUNE_FAST,&
                             ELPA_AUTOTUNE_DOMAIN_COMPLEX,ierr)
      elseif(real_cmplx == "real") then
         e_h%elpa_tune => e_h%elpa_main%autotune_setup(ELPA_AUTOTUNE_FAST,&
                             ELPA_AUTOTUNE_DOMAIN_REAL,ierr)
      endif

      if(ierr /= 0) then
         call elsi_stop(e_h,"ELPA auto-tuning failed.",caller)
      endif
   endif

   if(associated(e_h%elpa_tune)) then
      if(.not. e_h%elpa_main%autotune_step(e_h%elpa_tune)) then
         call e_h%elpa_main%autotune_set_best(e_h%elpa_tune)
         call elpa_autotune_deallocate(e_h%elpa_tune)

         nullify(e_h%elpa_tune)
      endif
   endif

end subroutine

end module ELSI_ELPA