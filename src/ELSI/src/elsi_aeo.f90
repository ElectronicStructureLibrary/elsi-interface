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
   use ELSI_IO,        only: elsi_say
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,&
                             mpi_integer4
   use ELSI_MU,        only: elsi_compute_mu_and_occ,elsi_compute_entropy
   use ELSI_PRECISION, only: r4,r8,i4
   use ELSI_TIMINGS,   only: elsi_get_time
   use ELSI_UTILS,     only: elsi_get_local_nnz_real,elsi_get_local_nnz_cmplx,&
                             elsi_trace_mat_mat_real,elsi_trace_mat_mat_cmplx
!   use ELPA
   use ELPA1,             only: elpa_get_communicators
   use elsi_elpa_api
   implicit none

   private

   public :: elsi_get_elpa_comms
   public :: elsi_set_elpa_default
   public :: elsi_compute_occ_elpa
   public :: elsi_compute_dm_elpa_real
   public :: elsi_compute_edm_elpa_real
   public :: elsi_normalize_dm_elpa_real
   public :: elsi_to_standard_evp_real
   public :: elsi_solve_evp_elpa_real
   public :: elsi_compute_dm_elpa_cmplx
   public :: elsi_compute_edm_elpa_cmplx
   public :: elsi_normalize_dm_elpa_cmplx
   public :: elsi_to_standard_evp_cmplx
   public :: elsi_solve_evp_elpa_cmplx

contains

!>
!! This routine gets the row and column communicators for ELPA.
!!
subroutine elsi_get_elpa_comms(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   integer(kind=i4) :: success

   character(len=40), parameter :: caller = "elsi_get_elpa_comms"

   success = elpa_get_communicators(e_h%mpi_comm,e_h%my_prow,e_h%my_pcol,&
                e_h%mpi_comm_row,e_h%mpi_comm_col)

   if(success /= 0) then
      call elsi_stop(" Failed to get MPI communicators.",e_h,caller)
   endif

   e_h%elpa_started = .true.

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
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

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
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine normalizes the density matrix to the exact number of electrons.
!!
subroutine elsi_normalize_dm_elpa_real(e_h,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: dm(e_h%n_lrow,e_h%n_lcol)

   real(kind=r8)      :: l_ne ! Local number of electrons
   real(kind=r8)      :: g_ne ! Global number of electrons
   real(kind=r8)      :: factor ! Normalization factor
   complex(kind=r8)   :: tmp_cmplx
   integer(kind=i4)   :: ierr
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_normalize_dm_elpa_real"

   call elsi_trace_mat_mat_real(e_h,dm,ovlp,g_ne)

   if(e_h%n_spins*e_h%n_kpts > 1) then
      if(e_h%myid == 0) then
         l_ne = g_ne
      else
         l_ne = 0.0_r8
      endif

      call MPI_Allreduce(l_ne,g_ne,1,mpi_real8,mpi_sum,e_h%mpi_comm_all,ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)
   endif

   ! Scale density matrix
   if(abs(g_ne-e_h%n_electrons) > e_h%occ_tolerance) then
      factor = e_h%n_electrons/g_ne

      write(info_str,"('  | Scaled density matrix by factor',F12.8)") factor
      call elsi_say(e_h,info_str)

      dm = dm*factor
   endif

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
   integer(kind=i4)   :: ierr
   logical            :: success
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
         call elpa_cholesky(e_h,ovlp)
  
         ! U^-1 -> S
         call elpa_invert_triangular(e_h,ovlp)

         call elsi_get_time(t1)

         write(info_str,"('  Finished Cholesky decomposition')")
         call elsi_say(e_h,info_str)
         write(info_str,"('  | Time :',F10.3,' s')") t1-t0
         call elsi_say(e_h,info_str)
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
      call elpa_hermitian_multiply(e_h,'U','L',ovlp,ham,evec)

      call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,evec,1,1,e_h%sc_desc,0.0_r8,&
              ham,1,1,e_h%sc_desc)

      evec = ham

      call elpa_hermitian_multiply(e_h,'U','U',ovlp,evec,ham)

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
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

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
   integer(kind=i4)   :: ierr
   real(kind=r8)      :: ev_sqrt
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str


   character(len=40), parameter :: caller = "elsi_check_singularity_real"

   call elsi_get_time(t0)

   ! TODO: ill-conditioning tolerance should be set here

   ! Use ELPA to check overlap singularity
   call elpa_eigenvectors(e_h,ovlp,eval,evec,2)
 
  do i = 1,e_h%n_basis
      if(eval(i) < e_h%sing_tol) then
         e_h%n_nonsing = e_h%n_nonsing-1
      endif
   enddo

   if(ierr /= 0) then
      call elsi_stop(" Singularity check failed.",e_h,caller)
   endif


   e_h%n_states_solve = min(e_h%n_nonsing,e_h%n_states)

   if(e_h%n_nonsing < e_h%n_basis) then ! Singular
      e_h%ovlp_is_sing = .true.

      call elsi_say(e_h,"  Overlap matrix is singular")
      write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)") eval(1)
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
         eval(e_h%n_basis)
      call elsi_say(e_h,info_str)

      if(e_h%stop_sing) then
         call elsi_stop(" Overlap matrix is singular.",e_h,caller)
      endif

      write(info_str,"('  | Number of basis functions reduced to :',I10)")&
         e_h%n_nonsing
      call elsi_say(e_h,info_str)

      call elsi_say(e_h,"  Using scaled eigenvectors of overlap matrix"//&
              " for transformation")

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = e_h%n_basis-e_h%n_nonsing+1,e_h%n_basis
         ev_sqrt = sqrt(eval(i))

         if(e_h%loc_col(i) == 0) cycle

         ovlp(:,e_h%loc_col(i)) = evec(:,e_h%loc_col(i))/ev_sqrt
      enddo
   else ! Nonsingular
      e_h%ovlp_is_sing = .false.
      call elsi_say(e_h,"  Overlap matrix is nonsingular")
      write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)") eval(1)
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
         eval(e_h%n_basis)
      call elsi_say(e_h,info_str)
   endif ! Singular overlap?

   call elpa_uninit()

   call elsi_get_time(t1)

   write(info_str,"('  Finished singularity check of overlap matrix')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)
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
   integer(kind=i4)   :: ierr
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
      call elpa_hermitian_multiply(e_h,'L','N',ham,tmp_real,evec)

   endif

   call elsi_deallocate(e_h,tmp_real,"tmp_real")

   call elsi_get_time(t1)

   write(info_str,"('  Finished back-transformation of eigenvectors')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_evp_elpa_real(e_h,ham,ovlp,eval,evec)

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


   character(len=40), parameter :: caller = "elsi_solve_evp_elpa_real"


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
   if(e_h%n_elsi_calls <= e_h%elpa_n_single) then
      call elsi_say(e_h,"  Starting ELPA eigensolver (single precision)")


      ! Solve with single precision
      call elpa_eigenvectors(e_h,real(ham,kind=r4),eval,evec)

   else ! No single precision
      call elsi_say(e_h,"  Starting ELPA eigensolver")

      ! Use copy_cmplx to store overlap matrix, otherwise it will
      ! be destroyed by eigenvalue calculation     
      call elpa_eigenvectors(e_h,ham,eval,evec)

   endif

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(e_h%n_nonsing < e_h%n_basis) then
      eval(e_h%n_nonsing+1:e_h%n_basis) = eval(e_h%n_nonsing)+10.0_r8
   endif

   call elsi_get_time(t1)

   write(info_str,"('  Finished solving standard eigenproblem')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

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
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

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
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine normalizes the density matrix to the exact number of electrons.
!!
subroutine elsi_normalize_dm_elpa_cmplx(e_h,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: dm(e_h%n_lrow,e_h%n_lcol)

   real(kind=r8)      :: l_ne ! Local number of electrons
   real(kind=r8)      :: g_ne ! Global number of electrons
   real(kind=r8)      :: factor ! Normalization factor
   complex(kind=r8)   :: tmp_cmplx
   integer(kind=i4)   :: ierr
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_normalize_dm_elpa_cmplx"

   call elsi_trace_mat_mat_cmplx(e_h,dm,ovlp,tmp_cmplx)

   g_ne = real(tmp_cmplx,kind=r8)

   if(e_h%n_spins*e_h%n_kpts > 1) then
      if(e_h%myid == 0) then
         l_ne = g_ne
      else
         l_ne = 0.0_r8
      endif

      call MPI_Allreduce(l_ne,g_ne,1,mpi_real8,mpi_sum,e_h%mpi_comm_all,ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)
   endif

   ! Scale density matrix
   if(abs(g_ne-e_h%n_electrons) > e_h%occ_tolerance) then
      factor = e_h%n_electrons/g_ne

      write(info_str,"('  | Scaled density matrix by factor',F12.8)") factor
      call elsi_say(e_h,info_str)

      dm = dm*factor
   endif

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
   integer(kind=i4)   :: ierr
   logical            :: success
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
         call elpa_cholesky(e_h,ovlp)

         ! U^-1 -> S
         call elpa_invert_triangular(e_h,ovlp)

         call elsi_get_time(t1)

         write(info_str,"('  Finished Cholesky decomposition')")
         call elsi_say(e_h,info_str)
         write(info_str,"('  | Time :',F10.3,' s')") t1-t0
         call elsi_say(e_h,info_str)
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
      call elpa_hermitian_multiply(e_h,'U','L',ovlp,ham,evec)

      call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),evec,1,1,&
              e_h%sc_desc,(0.0_r8,0.0_r8),ham,1,1,e_h%sc_desc)
      evec = ham

      call elpa_hermitian_multiply(e_h,'U','U',ovlp,evec,ham)

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
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

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
   integer(kind=i4)   :: ierr
   real(kind=r8)      :: ev_sqrt
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   class(elpa_t), pointer :: elpa_t

   complex(kind=r8), allocatable :: copy_cmplx(:,:)

   character(len=40), parameter :: caller = "elsi_check_singularity_cmplx"

   call elsi_get_time(t0)

   ! TODO: ill-conditioning tolerance should be set here


   ! Use ELPA to check overlap singularity
   call elpa_eigenvectors(e_h,ovlp,eval,evec)
   do i = 1,e_h%n_basis
      if(eval(i) < e_h%sing_tol) then
         e_h%n_nonsing = e_h%n_nonsing-1
      endif
   enddo

   if(ierr /= 0) then
      call elsi_stop(" Singularity check failed.",e_h,caller)
   endif


   e_h%n_states_solve = min(e_h%n_nonsing,e_h%n_states)

   if(e_h%n_nonsing < e_h%n_basis) then ! Singular
      e_h%ovlp_is_sing = .true.

      call elsi_say(e_h,"  Overlap matrix is singular")
      write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)") eval(1)
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
         eval(e_h%n_basis)
      call elsi_say(e_h,info_str)

      if(e_h%stop_sing) then
         call elsi_stop(" Overlap matrix is singular.",e_h,caller)
      endif

      write(info_str,"('  | Number of basis functions reduced to :',I10)")&
         e_h%n_nonsing
      call elsi_say(e_h,info_str)

      call elsi_say(e_h,"  Using scaled eigenvectors of overlap matrix"//&
              " for transformation")

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = e_h%n_basis-e_h%n_nonsing+1,e_h%n_basis
         ev_sqrt = sqrt(eval(i))

         if(e_h%loc_col(i) == 0) cycle

         ovlp(:,e_h%loc_col(i)) = evec(:,e_h%loc_col(i))/ev_sqrt
      enddo
   else ! Nonsingular
      e_h%ovlp_is_sing = .false.
      call elsi_say(e_h,"  Overlap matrix is nonsingular")
      write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)") eval(1)
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
         eval(e_h%n_basis)
      call elsi_say(e_h,info_str)
   endif ! Singular overlap?

   call elpa_uninit()

   call elsi_get_time(t1)

   write(info_str,"('  Finished singularity check of overlap matrix')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)
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
   integer(kind=i4)   :: ierr

   class(elpa_t), pointer :: elpa_t

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

      call elpa_hermitian_multiply(e_h,'L','N',ham,tmp_cmplx,evec)

      if(ierr /= 0) then
         call elsi_stop(" Matrix multiplication failed.",e_h,caller)
      endif
   endif

   call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")

   call elsi_get_time(t1)

   write(info_str,"('  Finished back-transformation of eigenvectors')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_evp_elpa_cmplx(e_h,ham,ovlp,eval,evec)

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

   character(len=40), parameter :: caller = "elsi_solve_evp_elpa_cmplx"


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
   if(e_h%n_elsi_calls <= e_h%elpa_n_single) then
      call elsi_say(e_h,"  Starting ELPA eigensolver (single precision)")

      ! Solve with single precision
      call elpa_eigenvectors(e_h,cmplx(ham,kind=r4),eval,evec)

   else ! No single precision
      call elsi_say(e_h,"  Starting ELPA eigensolver")

      call elpa_eigenvectors(e_h,ham,eval,evec)

   endif

!   if(ierr /= 0) then
!      call elsi_stop(" ELPA solver failed.",e_h,caller)
!   endif

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(e_h%n_nonsing < e_h%n_basis) then
      eval(e_h%n_nonsing+1:e_h%n_basis) = eval(e_h%n_nonsing)+10.0_r8
   endif

   call elsi_get_time(t1)

   write(info_str,"('  Finished solving standard eigenproblem')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

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

   if(e_h%handle_ready) then
      e_h%handle_changed = .true.
   endif

   ! ELPA solver
   e_h%elpa_solver = 2

   ! How many single precision steps?
   e_h%elpa_n_single = 0

end subroutine
end module ELSI_ELPA
