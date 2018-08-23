! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to ELPA-AEO.
!!
module ELSI_ELPA

   use ELSI_CONSTANTS, only: BLACS_DENSE,UT_MAT
   use ELSI_DATATYPE,  only: elsi_param_t,elsi_basic_t
   use ELSI_IO,        only: elsi_say,elsi_get_time
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,&
                             mpi_integer4
   use ELSI_OCC,       only: elsi_mu_and_occ,elsi_entropy
   use ELSI_PRECISION, only: r4,r8,i4
   use ELSI_UTILS,     only: elsi_get_nnz,elsi_set_full_mat
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

   interface elsi_elpa_evec
      module procedure elsi_elpa_evec_cmplx
      module procedure elsi_elpa_evec_real
   end interface

   interface elsi_elpa_mult
      module procedure elsi_elpa_mult_cmplx
      module procedure elsi_elpa_mult_real
   end interface

   interface elsi_elpa_chol
      module procedure elsi_elpa_chol_cmplx
      module procedure elsi_elpa_chol_real
   end interface

   interface elsi_elpa_invt
      module procedure elsi_elpa_invt_cmplx
      module procedure elsi_elpa_invt_real
   end interface

contains

!>
!! This routine initializes ELPA.
!!
subroutine elsi_init_elpa(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh

   integer(kind=i4) :: ierr

   character(len=40), parameter :: caller = "elsi_init_elpa"

   if(.not. ph%elpa_started) then
      ierr = elpa_init(20171201)

      if(ierr /= 0) then
         call elsi_stop(bh,"Initialization failed.",caller)
      endif

      call MPI_Comm_split(bh%comm,bh%my_pcol,bh%my_prow,ph%elpa_comm_row,ierr)

      call elsi_check_mpi(bh,"MPI_Comm_split",ierr,caller)

      call MPI_Comm_split(bh%comm,bh%my_prow,bh%my_pcol,ph%elpa_comm_col,ierr)

      call elsi_check_mpi(bh,"MPI_Comm_split",ierr,caller)

      ph%elpa_started = .true.
   endif

end subroutine

!>
!! This routine computes the chemical potential and occupation numbers.
!!
subroutine elsi_compute_occ_elpa(ph,bh,eval,occ)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   real(kind=r8),      intent(in)    :: eval(ph%n_basis)
   real(kind=r8),      intent(out)   :: occ(ph%n_states,ph%n_spins,ph%n_kpts)

   real(kind=r8)    :: mu
   real(kind=r8)    :: ts
   real(kind=r8)    :: n_electrons
   integer(kind=i4) :: n_states
   integer(kind=i4) :: n_spins
   integer(kind=i4) :: n_kpts
   integer(kind=i4) :: i
   integer(kind=i4) :: ierr

   real(kind=r8), allocatable :: eval_all(:,:,:)
   real(kind=r8), allocatable :: k_weight(:)
   real(kind=r8), allocatable :: tmp_real1(:)
   real(kind=r8), allocatable :: tmp_real2(:,:,:)

   character(len=40), parameter :: caller = "elsi_compute_occ_elpa"

   ! Gather eigenvalues and occupation numbers
   call elsi_allocate(bh,eval_all,ph%n_states,ph%n_spins,ph%n_kpts,"eval_all",&
           caller)
   call elsi_allocate(bh,k_weight,ph%n_kpts,"k_weight",caller)

   if(ph%n_kpts > 1) then
      call elsi_allocate(bh,tmp_real1,ph%n_kpts,"tmp_real",caller)

      if(bh%myid == 0 .and. ph%i_spin == 1) then
         tmp_real1(ph%i_kpt) = ph%i_weight
      endif

      call MPI_Allreduce(tmp_real1,k_weight,ph%n_kpts,mpi_real8,mpi_sum,&
              bh%comm_all,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,tmp_real1,"tmp_real")
   else
      k_weight = ph%i_weight
   endif

   if(ph%n_spins*ph%n_kpts > 1) then
      call elsi_allocate(bh,tmp_real2,ph%n_states,ph%n_spins,ph%n_kpts,&
              "tmp_real",caller)

      if(bh%myid == 0) then
         tmp_real2(:,ph%i_spin,ph%i_kpt) = eval(1:ph%n_states)
      endif

      call MPI_Allreduce(tmp_real2,eval_all,ph%n_states*ph%n_spins*ph%n_kpts,&
              mpi_real8,mpi_sum,bh%comm_all,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,tmp_real2,"tmp_real")
   else
      eval_all(:,ph%i_spin,ph%i_kpt) = eval(1:ph%n_states)
   endif

   ! Calculate chemical potential, occupation numbers, and electronic entropy
   n_electrons = ph%n_electrons
   n_states    = ph%n_states
   n_spins     = ph%n_spins
   n_kpts      = ph%n_kpts

   call elsi_mu_and_occ(ph,bh,n_electrons,n_states,n_spins,n_kpts,k_weight,&
           eval_all,occ,mu)

   call elsi_entropy(ph,n_states,n_spins,n_kpts,k_weight,eval_all,occ,mu,ts)

   ph%mu = mu
   ph%ts = ts

   ! Calculate band structure energy
   ph%ebs = 0.0_r8

   do i = 1,ph%n_states_solve
      ph%ebs = ph%ebs+eval(i)*occ(i,ph%i_spin,ph%i_kpt)
   enddo

   call elsi_deallocate(bh,eval_all,"eval_all")
   call elsi_deallocate(bh,k_weight,"k_weight")

end subroutine

!>
!! This routine constructs the density matrix.
!!
subroutine elsi_compute_dm_elpa_real(ph,bh,row_map,col_map,evec,occ,dm,work)

   implicit none

   type(elsi_param_t), intent(in)  :: ph
   type(elsi_basic_t), intent(in)  :: bh
   integer(kind=i4),   intent(in)  :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)  :: col_map(ph%n_basis)
   real(kind=r8),      intent(in)  :: evec(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)  :: occ(ph%n_states,ph%n_spins,ph%n_kpts)
   real(kind=r8),      intent(out) :: dm(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out) :: work(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: i
   integer(kind=i4)   :: max_state
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=40), parameter :: caller = "elsi_compute_dm_elpa_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,factor,ph%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,ph%n_states_solve
      if(occ(i,ph%i_spin,ph%i_kpt) > 0.0_r8) then
         factor(i) = sqrt(occ(i,ph%i_spin,ph%i_kpt))
         max_state = i
      endif
   enddo

   work = evec

   do i = 1,ph%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(col_map(i) > 0) then
            work(:,col_map(i)) = work(:,col_map(i))*factor(i)
         endif
      elseif(col_map(i) /= 0) then
         work(:,col_map(i)) = 0.0_r8
      endif
   enddo

   dm = 0.0_r8

   ! Compute density matrix
   call pdsyrk('U','N',ph%n_basis,max_state,1.0_r8,work,1,1,bh%desc,0.0_r8,dm,&
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

   type(elsi_param_t), intent(in)  :: ph
   type(elsi_basic_t), intent(in)  :: bh
   integer(kind=i4),   intent(in)  :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)  :: col_map(ph%n_basis)
   real(kind=r8),      intent(in)  :: eval(ph%n_basis)
   real(kind=r8),      intent(in)  :: evec(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)  :: occ(ph%n_states,ph%n_spins,ph%n_kpts)
   real(kind=r8),      intent(out) :: edm(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out) :: work(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: i
   integer(kind=i4)   :: max_state
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=40), parameter :: caller = "elsi_compute_edm_elpa_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,factor,ph%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,ph%n_states_solve
      factor(i) = -1.0_r8*occ(i,ph%i_spin,ph%i_kpt)*eval(i)
      if(factor(i) > 0.0_r8) then
         factor(i) = sqrt(factor(i))
         max_state = i
      else
         factor(i) = 0.0_r8
      endif
   enddo

   work = evec

   do i = 1,ph%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(col_map(i) > 0) then
            work(:,col_map(i)) = work(:,col_map(i))*factor(i)
         endif
      elseif(col_map(i) /= 0) then
         work(:,col_map(i)) = 0.0_r8
      endif
   enddo

   call elsi_deallocate(bh,factor,"factor")

   edm = 0.0_r8

   ! Compute density matrix
   call pdsyrk('U','N',ph%n_basis,max_state,1.0_r8,work,1,1,bh%desc,0.0_r8,edm,&
           1,1,bh%desc)

   edm = -1.0_r8*edm

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
   type(elsi_basic_t), intent(in)    :: bh
   integer(kind=i4),   intent(in)    :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)    :: col_map(ph%n_basis)
   real(kind=r8),      intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: eval(ph%n_basis)
   real(kind=r8),      intent(out)   :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_to_standard_evp_real"

   call elsi_get_time(t0)

   if(ph%n_calls == 1) then
      if(ph%check_sing) then
         call elsi_check_singularity_real(ph,bh,col_map,ovlp,eval,evec)
      endif

      if(ph%n_good == ph%n_basis) then ! Not singular

         ph%ovlp_is_sing = .false.

         ! S = (U^T)U, U -> S
         call elsi_elpa_chol(ph,bh,ovlp)

         ! U^-1 -> S
         call elsi_elpa_invt(ph,bh,ovlp)

         call elsi_get_time(t1)

         write(info_str,"(2X,A)") "Finished Cholesky decomposition"
         call elsi_say(bh,info_str)
         write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
         call elsi_say(bh,info_str)
      endif
   endif

   call elsi_get_time(t0)

   if(ph%ovlp_is_sing) then ! Use scaled eigenvectors
      ! evec_real used as tmp_real
      ! tmp_real = H_real * S_real
      call pdgemm('N','N',ph%n_basis,ph%n_good,ph%n_basis,1.0_r8,ham,1,1,&
              bh%desc,ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,0.0_r8,evec,1,1,&
              bh%desc)

      ! H_real = (S_real)^T * tmp_real
      call pdgemm('T','N',ph%n_good,ph%n_good,ph%n_basis,1.0_r8,ovlp,1,&
              ph%n_basis-ph%n_good+1,bh%desc,evec,1,1,bh%desc,0.0_r8,ham,1,1,&
              bh%desc)
   else ! Use Cholesky
      call elsi_elpa_mult(ph,bh,'U','L',ovlp,ham,evec)

      call pdtran(ph%n_basis,ph%n_basis,1.0_r8,evec,1,1,bh%desc,0.0_r8,ham,1,1,&
              bh%desc)

      evec = ham

      call elsi_elpa_mult(ph,bh,'U','U',ovlp,evec,ham)

      call elsi_set_full_mat(ph,bh,UT_MAT,row_map,col_map,ham)
   endif

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
   type(elsi_basic_t), intent(in)    :: bh
   integer(kind=i4),   intent(in)    :: col_map(ph%n_basis)
   real(kind=r8),      intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: eval(ph%n_basis)
   real(kind=r8),      intent(out)   :: evec(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: i
   real(kind=r8)      :: ev_sqrt
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_check_singularity_real"

   call elsi_get_time(t0)

   ! Use ELPA to check overlap singularity
   call elsi_elpa_evec(ph,bh,ovlp,eval,evec,.true.)

   do i = 1,ph%n_basis
      if(eval(i) < ph%sing_tol) then
         ph%n_good = ph%n_good-1
      endif
   enddo

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
      endif

      write(info_str,"(2X,A,I10)") "| Number of basis functions reduced to :",&
         ph%n_good
      call elsi_say(bh,info_str)

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = ph%n_basis-ph%n_good+1,ph%n_basis
         ev_sqrt = sqrt(eval(i))

         if(col_map(i) > 0) then
            ovlp(:,col_map(i)) = evec(:,col_map(i))/ev_sqrt
         endif
      enddo
   else ! Nonsingular
      ph%ovlp_is_sing = .false.

      write(info_str,"(2X,A)") "Overlap matrix is not singular"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,E10.2,A,E10.2)")&
         "| Lowest and highest eigenvalues :",eval(ph%n_basis),",",eval(1)
      call elsi_say(bh,info_str)
   endif ! Singular overlap?

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

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   real(kind=r8),      intent(out)   :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)    :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: tmp_real(:,:)

   character(len=40), parameter :: caller = "elsi_to_original_ev_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp_real,bh%n_lrow,bh%n_lcol,"tmp_real",caller)
   tmp_real = evec

   if(ph%ovlp_is_sing) then
      call pdgemm('N','N',ph%n_basis,ph%n_states_solve,ph%n_good,1.0_r8,ovlp,1,&
              ph%n_basis-ph%n_good+1,bh%desc,tmp_real,1,1,bh%desc,0.0_r8,evec,&
              1,1,bh%desc)
   else ! Nonsingular, use Cholesky
      call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp,1,1,bh%desc,0.0_r8,ham,1,1,&
              bh%desc)

      call elsi_elpa_mult(ph,bh,'L','N',ham,tmp_real,evec)
   endif

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
   integer(kind=i4),   intent(in)    :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)    :: col_map(ph%n_basis)
   real(kind=r8),      intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: eval(ph%n_basis)
   real(kind=r8),      intent(out)   :: evec(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: ierr
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_solve_elpa_real"

   ! Compute sparsity
   if(ph%n_calls == 1 .and. ph%matrix_format == BLACS_DENSE) then
      call elsi_get_nnz(bh%def0,ham,bh%n_lrow,bh%n_lcol,bh%nnz_l)

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   endif

   ! Transform to standard form
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_standard_evp_real(ph,bh,row_map,col_map,ham,ovlp,eval,evec)
   endif

   call elsi_get_time(t0)

   ! Solve evp, return eigenvalues and eigenvectors
   call elsi_elpa_evec(ph,bh,ham,eval,evec,.false.)

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(ph%n_good < ph%n_basis) then
      eval(ph%n_good+1:ph%n_basis) = eval(ph%n_good)+10.0_r8
   endif

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished solving standard eigenproblem"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

   ! Back-transform eigenvectors
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_original_ev_real(ph,bh,ham,ovlp,evec)
   endif

end subroutine

!>
!! This routine constructs the density matrix.
!!
subroutine elsi_compute_dm_elpa_cmplx(ph,bh,row_map,col_map,evec,occ,dm,work)

   implicit none

   type(elsi_param_t), intent(in)  :: ph
   type(elsi_basic_t), intent(in)  :: bh
   integer(kind=i4),   intent(in)  :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)  :: col_map(ph%n_basis)
   complex(kind=r8),   intent(in)  :: evec(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)  :: occ(ph%n_states,ph%n_spins,ph%n_kpts)
   complex(kind=r8),   intent(out) :: dm(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(out) :: work(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: i
   integer(kind=i4)   :: max_state
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=40), parameter :: caller = "elsi_compute_dm_elpa_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,factor,ph%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,ph%n_states_solve
      if(occ(i,ph%i_spin,ph%i_kpt) > 0.0_r8) then
         factor(i) = sqrt(occ(i,ph%i_spin,ph%i_kpt))
         max_state = i
      endif
   enddo

   work = evec

   do i = 1,ph%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(col_map(i) > 0) then
            work(:,col_map(i)) = work(:,col_map(i))*factor(i)
         endif
      elseif(col_map(i) /= 0) then
         work(:,col_map(i)) = (0.0_r8,0.0_r8)
      endif
   enddo

   dm = (0.0_r8,0.0_r8)

   ! Compute density matrix
   call pzherk('U','N',ph%n_basis,max_state,(1.0_r8,0.0_r8),work,1,1,bh%desc,&
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

   type(elsi_param_t), intent(in)  :: ph
   type(elsi_basic_t), intent(in)  :: bh
   integer(kind=i4),   intent(in)  :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)  :: col_map(ph%n_basis)
   real(kind=r8),      intent(in)  :: eval(ph%n_basis)
   complex(kind=r8),   intent(in)  :: evec(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)  :: occ(ph%n_states,ph%n_spins,ph%n_kpts)
   complex(kind=r8),   intent(out) :: edm(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(out) :: work(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: i
   integer(kind=i4)   :: max_state
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: factor(:)

   character(len=40), parameter :: caller = "elsi_compute_edm_elpa_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,factor,ph%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,ph%n_states_solve
      factor(i) = -1.0_r8*occ(i,ph%i_spin,ph%i_kpt)*eval(i)
      if(factor(i) > 0.0_r8) then
         factor(i) = sqrt(factor(i))
         max_state = i
      else
         factor(i) = 0.0_r8
      endif
   enddo

   work = evec

   do i = 1,ph%n_states_solve
      if(factor(i) > 0.0_r8) then
         if(col_map(i) > 0) then
            work(:,col_map(i)) = work(:,col_map(i))*factor(i)
         endif
      elseif(col_map(i) /= 0) then
         work(:,col_map(i)) = (0.0_r8,0.0_r8)
      endif
   enddo

   call elsi_deallocate(bh,factor,"factor")

   edm = (0.0_r8,0.0_r8)

   ! Compute density matrix
   call pzherk('U','N',ph%n_basis,max_state,(1.0_r8,0.0_r8),work,1,1,bh%desc,&
           (0.0_r8,0.0_r8),edm,1,1,bh%desc)

   edm = (-1.0_r8,0.0_r8)*edm

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
   type(elsi_basic_t), intent(in)    :: bh
   integer(kind=i4),   intent(in)    :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)    :: col_map(ph%n_basis)
   complex(kind=r8),   intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: eval(ph%n_basis)
   complex(kind=r8),   intent(out)   :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_to_standard_evp_cmplx"

   call elsi_get_time(t0)

   if(ph%n_calls == 1) then
      if(ph%check_sing) then
         call elsi_check_singularity_cmplx(ph,bh,col_map,ovlp,eval,evec)
      endif

      if(ph%n_good == ph%n_basis) then ! Not singular
         ph%ovlp_is_sing = .false.

         ! S = (U^T)U, U -> S
         call elsi_elpa_chol(ph,bh,ovlp)

         ! U^-1 -> S
         call elsi_elpa_invt(ph,bh,ovlp)

         call elsi_get_time(t1)

         write(info_str,"(2X,A)") "Finished Cholesky decomposition"
         call elsi_say(bh,info_str)
         write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
         call elsi_say(bh,info_str)
      endif
   endif

   call elsi_get_time(t0)

   if(ph%ovlp_is_sing) then ! Use scaled eigenvectors
      ! evec_cmplx used as tmp_cmplx
      ! tmp_cmplx = H_cmplx * S_cmplx
      call pzgemm('N','N',ph%n_basis,ph%n_good,ph%n_basis,(1.0_r8,0.0_r8),ham,&
              1,1,bh%desc,ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,&
              (0.0_r8,0.0_r8),evec,1,1,bh%desc)

      ! H_cmplx = (S_cmplx)^* * tmp_cmplx
      call pzgemm('C','N',ph%n_good,ph%n_good,ph%n_basis,(1.0_r8,0.0_r8),ovlp,&
              1,ph%n_basis-ph%n_good+1,bh%desc,evec,1,1,bh%desc,&
              (0.0_r8,0.0_r8),ham,1,1,bh%desc)
   else ! Use cholesky
      call elsi_elpa_mult(ph,bh,'U','L',ovlp,ham,evec)

      call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),evec,1,1,bh%desc,&
              (0.0_r8,0.0_r8),ham,1,1,bh%desc)

      evec = ham

      call elsi_elpa_mult(ph,bh,'U','U',ovlp,evec,ham)

      call elsi_set_full_mat(ph,bh,UT_MAT,row_map,col_map,ham)
   endif

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
   type(elsi_basic_t), intent(in)    :: bh
   integer(kind=i4),   intent(in)    :: col_map(ph%n_basis)
   complex(kind=r8),   intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: eval(ph%n_basis)
   complex(kind=r8),   intent(out)   :: evec(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: i
   real(kind=r8)      :: ev_sqrt
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_check_singularity_cmplx"

   call elsi_get_time(t0)

   ! Use ELPA to check overlap singularity
   call elsi_elpa_evec(ph,bh,ovlp,eval,evec,.true.)

   do i = 1,ph%n_basis
      if(eval(i) < ph%sing_tol) then
         ph%n_good = ph%n_good-1
      endif
   enddo

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
      endif

      write(info_str,"(2X,A,I10)") "| Number of basis functions reduced to :",&
         ph%n_good
      call elsi_say(bh,info_str)

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = ph%n_basis-ph%n_good+1,ph%n_basis
         ev_sqrt = sqrt(eval(i))

         if(col_map(i) > 0) then
            ovlp(:,col_map(i)) = evec(:,col_map(i))/ev_sqrt
         endif
      enddo
   else ! Nonsingular
      ph%ovlp_is_sing = .false.

      write(info_str,"(2X,A)") "Overlap matrix is not singular"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,E10.2,A,E10.2)")&
         "| Lowest and highest eigenvalues :",eval(ph%n_basis),",",eval(1)
      call elsi_say(bh,info_str)
   endif ! Singular overlap?

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

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   complex(kind=r8),   intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character(len=40), parameter :: caller = "elsi_to_original_ev_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp_cmplx,bh%n_lrow,bh%n_lcol,"tmp_cmplx",caller)
   tmp_cmplx = evec

   if(ph%ovlp_is_sing) then
      call pzgemm('N','N',ph%n_basis,ph%n_states_solve,ph%n_good,&
              (1.0_r8,0.0_r8),ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,tmp_cmplx,&
              1,1,bh%desc,(0.0_r8,0.0_r8),evec,1,1,bh%desc)
   else ! Nonsingular, use Cholesky
      call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),ovlp,1,1,bh%desc,&
              (0.0_r8,0.0_r8),ham,1,1,bh%desc)

      call elsi_elpa_mult(ph,bh,'L','N',ham,tmp_cmplx,evec)
   endif

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
   integer(kind=i4),   intent(in)    :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)    :: col_map(ph%n_basis)
   complex(kind=r8),   intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: eval(ph%n_basis)
   complex(kind=r8),   intent(out)   :: evec(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: ierr
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_solve_elpa_cmplx"

   ! Compute sparsity
   if(ph%n_calls == 1 .and. ph%matrix_format == BLACS_DENSE) then
      call elsi_get_nnz(bh%def0,ham,bh%n_lrow,bh%n_lcol,bh%nnz_l)

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   endif

   ! Transform to standard form
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_standard_evp_cmplx(ph,bh,row_map,col_map,ham,ovlp,eval,evec)
   endif

   call elsi_get_time(t0)

   ! Solve evp, return eigenvalues and eigenvectors
   call elsi_elpa_evec(ph,bh,ham,eval,evec,.false.)

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(ph%n_good < ph%n_basis) then
      eval(ph%n_good+1:ph%n_basis) = eval(ph%n_good)+10.0_r8
   endif

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished solving standard eigenproblem"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

   ! Back-transform eigenvectors
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_original_ev_cmplx(ph,bh,ham,ovlp,evec)
   endif

end subroutine

!>
!! This routine cleans up ELPA.
!!
subroutine elsi_cleanup_elpa(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   integer(kind=i4) :: ierr

   character(len=40), parameter :: caller = "elsi_cleanup_elpa"

   if(ph%elpa_started) then
      nullify(ph%elpa_main)

      call elpa_uninit()

      call MPI_Comm_free(ph%elpa_comm_row,ierr)
      call MPI_Comm_free(ph%elpa_comm_col,ierr)
   endif

   ph%elpa_started = .false.

end subroutine

!>
!! This routine calls ELPA real eigensolver.
!!
subroutine elsi_elpa_evec_real(ph,bh,ham,eval,evec,sing_check)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   real(kind=r8),      intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(inout) :: eval(ph%n_basis)
   real(kind=r8),      intent(inout) :: evec(bh%n_lrow,bh%n_lcol)
   logical,            intent(in)    :: sing_check

   integer(kind=i4)   :: ierr
   character(len=200) :: info_str

   class(elpa_t), pointer :: elpa_main

   real(kind=r8), allocatable :: copy_ham(:,:)
   real(kind=r4), allocatable :: copy_ham_single(:,:)
   real(kind=r4), allocatable :: eval_single(:)
   real(kind=r4), allocatable :: evec_single(:,:)

   character(len=40), parameter :: caller = "elsi_elpa_evec_real"

   if(sing_check) then
      call elsi_elpa_setup(ph,bh,elpa_main,ph%n_basis,ph%n_basis)
      ! TODO: Define ill-conditioning tolerance (not yet available in ELPA)
!      call elpa_main%set("check_pd",1,ierr)
      call elpa_main%set("solver",ELPA_SOLVER_2STAGE,ierr)

      call elsi_allocate(bh,copy_ham,bh%n_lrow,bh%n_lcol,"copy_ham",caller)

      copy_ham = ham

      call elpa_main%eigenvectors(copy_ham,eval,evec,ierr)
      call elpa_deallocate(elpa_main)

      nullify(elpa_main)

      call elsi_deallocate(bh,copy_ham,"copy_ham")
   else
      call elsi_elpa_setup(ph,bh,ph%elpa_main,ph%n_good,ph%n_states_solve)
      call elsi_elpa_autotuning(ph,bh,"real")

      elpa_main => ph%elpa_main

      if(ph%n_calls <= ph%elpa_n_single) then
         write(info_str,"(2X,A)") "Starting ELPA eigensolver (single precision)"
         call elsi_say(bh,info_str)

         call elsi_allocate(bh,eval_single,ph%n_basis,"eval_single",caller)
         call elsi_allocate(bh,evec_single,bh%n_lrow,bh%n_lcol,"evec_single",&
                 caller)
         call elsi_allocate(bh,copy_ham,bh%n_lrow,bh%n_lcol,"copy_ham_single",&
                 caller)

         copy_ham_single=real(ham,kind=r4)

         call elpa_main%eigenvectors(copy_ham_single,eval_single,evec_single,&
                 ierr)

         eval = real(eval_single,kind=r8)
         evec = real(evec_single,kind=r8)

         call elsi_deallocate(bh,eval_single,"eval_single")
         call elsi_deallocate(bh,evec_single,"evec_single")
         call elsi_deallocate(bh,copy_ham_single,"copy_ham_single")
      else
         write(info_str,"(2X,A)") "Starting ELPA eigensolver"
         call elsi_say(bh,info_str)

         call elpa_main%eigenvectors(ham,eval,evec,ierr)
      endif
   endif

   if(ierr /= 0) then
      call elsi_stop(bh,"ELPA eigensolver failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA complex eigensolver.
!!
subroutine elsi_elpa_evec_cmplx(ph,bh,ham,eval,evec,sing_check)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   complex(kind=r8),   intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(inout) :: eval(ph%n_basis)
   complex(kind=r8),   intent(inout) :: evec(bh%n_lrow,bh%n_lcol)
   logical,            intent(in)    :: sing_check

   integer(kind=i4)   :: ierr
   character(len=200) :: info_str

   class(elpa_t), pointer :: elpa_main

   complex(kind=r8), allocatable :: copy_ham(:,:)
   complex(kind=r4), allocatable :: copy_ham_single(:,:)
   real(kind=r4),    allocatable :: eval_single(:)
   complex(kind=r4), allocatable :: evec_single(:,:)

   character(len=40), parameter :: caller = "elsi_elpa_evec_cmplx"

   if(sing_check) then
      call elsi_elpa_setup(ph,bh,elpa_main,ph%n_basis,ph%n_basis)
      ! TODO: Define ill-conditioning tolerance (not yet available in ELPA)
!      call elpa_main%set("check_pd",1,ierr)
      call elpa_main%set("solver",ELPA_SOLVER_2STAGE,ierr)

      call elsi_allocate(bh,copy_ham,bh%n_lrow,bh%n_lcol,"copy_ham",caller)

      copy_ham = ham

      call elpa_main%eigenvectors(copy_ham,eval,evec,ierr)
      call elpa_deallocate(elpa_main)

      nullify(elpa_main)

      call elsi_deallocate(bh,copy_ham,"copy_ham")
   else
      call elsi_elpa_setup(ph,bh,ph%elpa_main,ph%n_good,ph%n_states_solve)
      call elsi_elpa_autotuning(ph,bh,"complex")

      elpa_main => ph%elpa_main

      if(ph%n_calls <= ph%elpa_n_single) then
         write(info_str,"(2X,A)") "Starting ELPA eigensolver (single precision)"
         call elsi_say(bh,info_str)

         call elsi_allocate(bh,eval_single,ph%n_basis,"eval_single",caller)
         call elsi_allocate(bh,evec_single,bh%n_lrow,bh%n_lcol,"evec_single",&
                 caller)
         call elsi_allocate(bh,copy_ham,bh%n_lrow,bh%n_lcol,"copy_ham_single",&
                 caller)

         copy_ham_single=cmplx(ham,kind=r4)

         call elpa_main%eigenvectors(copy_ham_single,eval_single,evec_single,&
                 ierr)

         eval = real(eval_single,kind=r8)
         evec = cmplx(evec_single,kind=r8)

         call elsi_deallocate(bh,eval_single,"eval_single")
         call elsi_deallocate(bh,evec_single,"evec_single")
         call elsi_deallocate(bh,copy_ham_single,"copy_ham_single")
      else
         write(info_str,"(2X,A)") "Starting ELPA eigensolver"
         call elsi_say(bh,info_str)

         call elpa_main%eigenvectors(ham,eval,evec,ierr)
      endif
   endif

   if(ierr /= 0) then
      call elsi_stop(bh,"ELPA eigensolver failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA real hermitian_multiply.
!!
subroutine elsi_elpa_mult_real(ph,bh,uplo,uplo2,a,b,c)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   character,         intent(in)     :: uplo
   character,         intent(in)     :: uplo2
   real(kind=r8),     intent(inout)  :: a(bh%n_lrow,bh%n_lcol)
   real(kind=r8),     intent(inout)  :: b(bh%n_lrow,bh%n_lcol)
   real(kind=r8),     intent(inout)  :: c(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_mult_real"

   call elsi_elpa_setup(ph,bh,elpa_main,ph%n_basis,ph%n_basis)
   call elpa_main%hermitian_multiply(uplo,uplo2,ph%n_basis,a,b,bh%n_lrow,&
           bh%n_lcol,c,bh%n_lrow,bh%n_lcol,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(bh,"Matrix multiplication failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA complex hermitian_multiply.
!!
subroutine elsi_elpa_mult_cmplx(ph,bh,uplo,uplo2,a,b,c)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   character,          intent(in)    :: uplo
   character,          intent(in)    :: uplo2
   complex(kind=r8),   intent(inout) :: a(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(inout) :: b(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(inout) :: c(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_mult_cmplx"

   call elsi_elpa_setup(ph,bh,elpa_main,ph%n_basis,ph%n_basis)
   call elpa_main%hermitian_multiply(uplo,uplo2,ph%n_basis,a,b,bh%n_lrow,&
           bh%n_lcol,c,bh%n_lrow,bh%n_lcol,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(bh,"Matrix multiplication failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA real cholesky.
!!
subroutine elsi_elpa_chol_real(ph,bh,a)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   real(kind=r8),      intent(inout) :: a(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_chol_real"

   call elsi_elpa_setup(ph,bh,elpa_main,ph%n_basis,ph%n_basis)
   call elpa_main%cholesky(a,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(bh,"Cholesky decomposition failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA complex cholesky.
!!
subroutine elsi_elpa_chol_cmplx(ph,bh,a)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   complex(kind=r8),   intent(inout) :: a(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_chol_cmplx"

   call elsi_elpa_setup(ph,bh,elpa_main,ph%n_basis,ph%n_basis)
   call elpa_main%cholesky(a,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(bh,"Cholesky decomposition failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA real invert_triangular.
!!
subroutine elsi_elpa_invt_real(ph,bh,a)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   real(kind=r8),      intent(inout) :: a(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_invt_real"

   call elsi_elpa_setup(ph,bh,elpa_main,ph%n_good,ph%n_good)
   call elpa_main%invert_triangular(a,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(bh,"Matrix inversion failed.",caller)
   endif

end subroutine

!>
!! This routine calls ELPA complex invert_triangular.
!!
subroutine elsi_elpa_invt_cmplx(ph,bh,a)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   complex(kind=r8),   intent(inout) :: a(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr

   class(elpa_t), pointer :: elpa_main

   character(len=40), parameter :: caller = "elsi_elpa_invt_cmplx"

   call elsi_elpa_setup(ph,bh,elpa_main,ph%n_good,ph%n_good)
   call elpa_main%invert_triangular(a,ierr)
   call elpa_deallocate(elpa_main)

   nullify(elpa_main)

   if(ierr /= 0) then
      call elsi_stop(bh,"Matrix inversion failed.",caller)
   endif

end subroutine

!>
!! This routine sets ELPA-AEO parameters.
!!
subroutine elsi_elpa_setup(ph,bh,elpa_i,na,nev)

   implicit none

   type(elsi_param_t), intent(inout)          :: ph
   type(elsi_basic_t), intent(in)             :: bh
   class(elpa_t),      intent(inout), pointer :: elpa_i
   integer,            intent(in)             :: na
   integer,            intent(in)             :: nev

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: ierr2
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_elpa_setup"

   if(ph%n_calls == ph%elpa_n_single .and. associated(elpa_i)) then
      call elpa_deallocate(elpa_i)

      nullify(elpa_i)
   endif

   if(.not. associated(elpa_i)) then
      elpa_i => elpa_allocate()

      call elpa_i%set("na",na,ierr)
      call elpa_i%set("nev",nev,ierr)
      call elpa_i%set("local_nrows",bh%n_lrow,ierr)
      call elpa_i%set("local_ncols",bh%n_lcol,ierr)
      call elpa_i%set("nblk",bh%blk,ierr)
      call elpa_i%set("mpi_comm_parent",bh%comm,ierr)
      call elpa_i%set("process_row",bh%my_prow,ierr)
      call elpa_i%set("process_col",bh%my_pcol,ierr)

      ierr = elpa_i%setup()

      if(ierr /= 0) then
         call elsi_stop(bh,"ELPA setup failed.",caller)
      endif

      if(ph%elpa_solver == 1) then
         call elpa_i%set("solver",ELPA_SOLVER_1STAGE,ierr)
      else
         call elpa_i%set("solver",ELPA_SOLVER_2STAGE,ierr)
      endif

      ! Try to enable ELPA GPU acceleration
      if(ph%elpa_gpu) then
         call elpa_i%set("gpu",1,ierr)

         if(ierr /= 0) then ! Failed
            call elpa_i%set("gpu",0,ierr)

            ph%elpa_gpu         = .false.
            ph%elpa_gpu_kernels = .false.

            write(info_str,"(2X,A)") "No ELPA GPU acceleration available"
            call elsi_say(bh,info_str)
         else
            write(info_str,"(2X,A)") "ELPA GPU acceleration activated"
            call elsi_say(bh,info_str)
         endif
      endif

      ! Try to enable ELPA2 GPU kernels
      if(ph%elpa_gpu_kernels) then
         if(ph%elpa_solver == 2) then
            call elpa_i%set("real_kernel",ELPA_2STAGE_REAL_GPU,ierr)
            call elpa_i%set("complex_kernel",ELPA_2STAGE_COMPLEX_GPU,ierr2)

            if(ierr /= 0 .or. ierr2 /= 0) then
               call elpa_i%set("real_kernel",ELPA_2STAGE_REAL_DEFAULT,ierr)
               call elpa_i%set("complex_kernel",ELPA_2STAGE_COMPLEX_DEFAULT,&
                       ierr)

               ph%elpa_gpu_kernels = .false.

               write(info_str,"(2X,A)") "ELPA GPU kernels not available"
               call elsi_say(bh,info_str)
            else
               write(info_str,"(2X,A)") "ELPA GPU kernels will be used"
               call elsi_say(bh,info_str)
            endif
         else ! GPU kernels currently only make sense with ELPA2
            write(info_str,"(2X,A)")&
               "No GPU kernels available with 1-stage ELPA"
            call elsi_say(bh,info_str)
         endif
      endif
   endif

end subroutine

!>
!! This routine sets up ELPA AEO auto-tuning.
!!
subroutine elsi_elpa_autotuning(ph,bh,real_cmplx)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   character(len=*),   intent(in)    :: real_cmplx

   integer(kind=i4) :: ierr

   character(len=40), parameter :: caller = "elsi_elpa_autotuning"

   if(ph%elpa_autotune) then
      if(ph%n_calls == ph%elpa_n_single) then
         call elpa_autotune_deallocate(ph%elpa_tune)
      endif

      if(ph%n_calls == 1 .or. ph%n_calls == ph%elpa_n_single) then
         if(real_cmplx == "complex") then
            ph%elpa_tune => ph%elpa_main%autotune_setup(ELPA_AUTOTUNE_FAST,&
                               ELPA_AUTOTUNE_DOMAIN_COMPLEX,ierr)
         elseif(real_cmplx == "real") then
            ph%elpa_tune => ph%elpa_main%autotune_setup(ELPA_AUTOTUNE_FAST,&
                               ELPA_AUTOTUNE_DOMAIN_REAL,ierr)
         endif

         if(ierr /= 0) then
            call elsi_stop(bh,"ELPA auto-tuning failed.",caller)
         endif
      endif

      if(associated(ph%elpa_tune)) then
         if(.not. ph%elpa_main%autotune_step(ph%elpa_tune)) then
            call ph%elpa_main%autotune_set_best(ph%elpa_tune)
            call elpa_autotune_deallocate(ph%elpa_tune)

            nullify(ph%elpa_tune)
         endif
      endif
   endif

end subroutine

end module ELSI_ELPA
