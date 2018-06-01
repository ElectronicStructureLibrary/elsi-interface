! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module implements the density matrix purification algorithm.
!!
module ELSI_DMP

   use ELSI_CONSTANTS, only: BLACS_DENSE,UT_MAT,LT_MAT,CANONICAL,&
                             TRACE_CORRECTING
   use ELSI_DATATYPE,  only: elsi_param_t,elsi_basic_t
   use ELSI_ELPA,      only: elsi_to_standard_evp_real
   use ELSI_IO,        only: elsi_say,elsi_get_time
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,&
                             mpi_integer4
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS,     only: elsi_get_nnz_real,elsi_trace_mat_real,&
                             elsi_trace_mat_mat_real

   implicit none

   private

   public :: elsi_init_dmp
   public :: elsi_cleanup_dmp
   public :: elsi_solve_dmp_real

contains

!>
!! This routine initialzes density matrix purification.
!!
subroutine elsi_init_dmp(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   character(len=40), parameter :: caller = "elsi_init_dmp"

   if(.not. ph%dmp_started) then
      ! No initialization needed
      ph%dmp_started = .true.
   endif

end subroutine

!>
!! This routine computes the density matrix using the density matrix
!! purification algorithm.
!!
subroutine elsi_solve_dmp_real(ph,bh,row_map,col_map,ham,ovlp,ham_copy,&
              ovlp_copy,ovlp_inv,vec1,vec2,dm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)    :: col_map(ph%n_basis)
   real(kind=r8),      intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)    :: ham_copy(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)    :: ovlp_copy(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(inout) :: ovlp_inv(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(inout) :: vec1(bh%n_lrow)
   real(kind=r8),      intent(inout) :: vec2(bh%n_lrow)
   real(kind=r8),      intent(out)   :: dm(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: i
   integer(kind=i4)   :: j
   integer(kind=i4)   :: i_iter
   integer(kind=i4)   :: ierr
   logical            :: ev_min_found
   logical            :: dmp_conv
   real(kind=r8)      :: dmp_ne
   real(kind=r8)      :: this_ev
   real(kind=r8)      :: prev_ev
   real(kind=r8)      :: nrm2
   real(kind=r8)      :: diff
   real(kind=r8)      :: mu
   real(kind=r8)      :: lambda
   real(kind=r8)      :: c1
   real(kind=r8)      :: c2
   real(kind=r8)      :: c
   real(kind=r8)      :: tmp
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   real(kind=r8), allocatable :: tmp_real1(:)
   real(kind=r8), allocatable :: tmp_real2(:)
   real(kind=r8), allocatable :: tmp_real3(:,:)
   real(kind=r8), allocatable :: dsd(:,:)
   real(kind=r8), allocatable :: dsdsd(:,:)

   character(len=40), parameter :: caller = "elsi_solve_dmp_real"

   ! Compute sparsity
   if(ph%n_calls == 1 .and. ph%matrix_format == BLACS_DENSE) then
      call elsi_get_nnz_real(bh%zero_def,ham,bh%n_lrow,bh%n_lcol,bh%nnz_l)

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   endif

   call elsi_get_time(t0)

   write(info_str,"(2X,A)") "Starting density matrix purification"
   call elsi_say(bh,info_str)

   ! Transform the generalized evp to the standard form
   ph%check_sing = .false.

   call elsi_allocate(bh,tmp_real1,ph%n_basis,"tmp_real1",caller)
   call elsi_to_standard_evp_real(ph,bh,row_map,col_map,ham,ovlp,tmp_real1,dm)
   call elsi_deallocate(bh,tmp_real1,"tmp_real1")

   ! Compute inverse of overlap
   if(ph%n_calls == 1) then
      ovlp_inv = ovlp_copy

      call pdpotrf('U',ph%n_basis,ovlp_inv,1,1,bh%desc,ierr)
      call pdpotri('U',ph%n_basis,ovlp_inv,1,1,bh%desc,ierr)

      call elsi_allocate(bh,tmp_real3,bh%n_lrow,bh%n_lcol,"tmp_real3",caller)

      call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp_inv,1,1,bh%desc,0.0_r8,&
              tmp_real3,1,1,bh%desc)

      do j = 1,ph%n_basis-1
         if(col_map(j) > 0) then
            do i = j+1,ph%n_basis
               if(row_map(i) > 0) then
                  ovlp_inv(row_map(i),col_map(j)) = &
                     tmp_real3(row_map(i),col_map(j))
               endif
            enddo
         endif
      enddo

      call elsi_deallocate(bh,tmp_real3,"tmp_real3")
   endif

   ! Use power iteration to find the largest in magnitude eigenvalue
   ! Usually this is the smallest, which is a better estimate of dmp_ev_ham_min
   ev_min_found = .false.

   call elsi_allocate(bh,tmp_real1,bh%n_lrow,"tmp_real1",caller)
   call elsi_allocate(bh,tmp_real2,bh%n_lrow,"tmp_real2",caller)

   if(ph%n_calls == 1) then
      tmp_real1 = 1.0_r8/sqrt(real(ph%n_basis,kind=r8))
   else
      tmp_real1 = vec1
   endif

   prev_ev = 0.0_r8

   do i_iter = 1,ph%dmp_max_power
      call pdgemv('N',ph%n_basis,ph%n_basis,1.0_r8,ham,1,1,bh%desc,tmp_real1,1,&
              1,bh%desc,1,0.0_r8,tmp_real2,1,1,bh%desc,1)

      this_ev = 0.0_r8
      nrm2    = 0.0_r8

      do i = 1,ph%n_basis
         if(row_map(i) > 0) then
            this_ev = this_ev+tmp_real1(row_map(i))*tmp_real2(row_map(i))
            nrm2    = nrm2+tmp_real2(row_map(i))*tmp_real2(row_map(i))
         endif
      enddo

      tmp = this_ev

      call MPI_Allreduce(tmp,this_ev,1,mpi_real8,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      tmp = nrm2

      call MPI_Allreduce(tmp,nrm2,1,mpi_real8,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      tmp_real1 = tmp_real2/sqrt(nrm2)

      if(abs(this_ev-prev_ev) < 1.0e-4_r8) then
         exit
      endif

      prev_ev = this_ev
   enddo

   if(this_ev > 0) then
      ph%dmp_ev_ham_max = this_ev
   else
      ph%dmp_ev_ham_min = this_ev
      ev_min_found      = .true.
   endif

   vec1 = tmp_real1

   if(ph%n_calls == 1) then
      tmp_real1 = 1.0_r8/sqrt(real(ph%n_basis,kind=r8))
   else
      tmp_real1 = vec2
   endif

   ! Shift H and use power iteration to find a better estimate of ev_ham_max
   dm = ham

   do i = 1,ph%n_basis
      if(row_map(i) > 0 .and. col_map(i) > 0) then
         ham(row_map(i),col_map(i)) = ham(row_map(i),col_map(i))+&
                                         abs(ph%dmp_ev_ham_min)
      endif
   enddo

   prev_ev = 0.0_r8

   do i_iter = 1,ph%dmp_max_power
      call pdgemv('N',ph%n_basis,ph%n_basis,1.0_r8,ham,1,1,bh%desc,tmp_real1,1,&
              1,bh%desc,1,0.0_r8,tmp_real2,1,1,bh%desc,1)

      this_ev = 0.0_r8
      nrm2    = 0.0_r8

      do i = 1,ph%n_basis
         if(row_map(i) > 0) then
            this_ev = this_ev+tmp_real1(row_map(i))*tmp_real2(row_map(i))
            nrm2    = nrm2+tmp_real2(row_map(i))*tmp_real2(row_map(i))
         endif
      enddo

      tmp = this_ev

      call MPI_Allreduce(tmp,this_ev,1,mpi_real8,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      tmp = nrm2

      call MPI_Allreduce(tmp,nrm2,1,mpi_real8,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      tmp_real1 = tmp_real2/sqrt(nrm2)

      if(abs(this_ev-prev_ev) < 1.0e-4_r8) then
         exit
      endif

      prev_ev = this_ev
   enddo

   if(ev_min_found .and. this_ev > 0) then
      ph%dmp_ev_ham_max = this_ev-abs(ph%dmp_ev_ham_min)
   elseif(this_ev < 0) then
      ph%dmp_ev_ham_min = this_ev+abs(ph%dmp_ev_ham_max)
   endif

   vec2 = tmp_real1
   ham  = dm

   call elsi_deallocate(bh,tmp_real1,"tmp_real1")
   call elsi_deallocate(bh,tmp_real2,"tmp_real2")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished power iteration"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

   ! Initialization
   call elsi_get_time(t0)

   call elsi_trace_mat_real(ph,bh,row_map,col_map,ham,mu)

   mu     = mu/ph%n_basis
   lambda = min(ph%dmp_n_states/(ph%dmp_ev_ham_max-mu),&
               (ph%n_basis-ph%dmp_n_states)/(mu-ph%dmp_ev_ham_min))

   call elsi_allocate(bh,dsd,bh%n_lrow,bh%n_lcol,"dsd",caller)
   if(ph%dmp_method == CANONICAL) then
      call elsi_allocate(bh,dsdsd,bh%n_lrow,bh%n_lcol,"dsdsd",caller)
   endif

   ! ham_real used as tmp after this point
   call pdgemm('N','N',ph%n_basis,ph%n_basis,ph%n_basis,1.0_r8,ham_copy,1,1,&
           bh%desc,ovlp_inv,1,1,bh%desc,0.0_r8,ham,1,1,bh%desc)
   call pdgemm('N','N',ph%n_basis,ph%n_basis,ph%n_basis,1.0_r8,ovlp_inv,1,1,&
           bh%desc,ham,1,1,bh%desc,0.0_r8,dm,1,1,bh%desc)

   dm = (mu*ovlp_inv-dm)*lambda/ph%n_basis
   dm = dm+ovlp_inv*ph%dmp_n_states/ph%n_basis

   ! Start main density matrix purification loop
   dmp_conv = .false.

   do i_iter = 1,ph%dmp_max_iter
      ham = dm

      select case(ph%dmp_method)
      case(TRACE_CORRECTING)
         call pdgemm('N','N',ph%n_basis,ph%n_basis,ph%n_basis,1.0_r8,ham,1,1,&
                 bh%desc,ovlp_copy,1,1,bh%desc,0.0_r8,dm,1,1,bh%desc)

         call elsi_trace_mat_real(ph,bh,row_map,col_map,dm,c1)

         call pdgemm('N','N',ph%n_basis,ph%n_basis,ph%n_basis,1.0_r8,dm,1,1,&
                 bh%desc,ham,1,1,bh%desc,0.0_r8,dsd,1,1,bh%desc)

         if(ph%dmp_n_states-c1 > 0.0_r8) then
            dm = 2*ham-dsd
         else
            dm = dsd
         endif
      case(CANONICAL)
         call pdgemm('N','N',ph%n_basis,ph%n_basis,ph%n_basis,1.0_r8,ovlp_copy,&
                 1,1,bh%desc,ham,1,1,bh%desc,0.0_r8,dm,1,1,bh%desc)
         call pdgemm('N','N',ph%n_basis,ph%n_basis,ph%n_basis,1.0_r8,ham,1,1,&
                 bh%desc,dm,1,1,bh%desc,0.0_r8,dsd,1,1,bh%desc)
         call pdgemm('N','N',ph%n_basis,ph%n_basis,ph%n_basis,1.0_r8,ovlp_copy,&
                 1,1,bh%desc,dsd,1,1,bh%desc,0.0_r8,dm,1,1,bh%desc)
         call pdgemm('N','N',ph%n_basis,ph%n_basis,ph%n_basis,1.0_r8,ham,1,1,&
                 bh%desc,dm,1,1,bh%desc,0.0_r8,dsdsd,1,1,bh%desc)

         dm = dsd-dsdsd

         call elsi_trace_mat_mat_real(bh,dm,ovlp_copy,c1)

         dm = ham-dsd

         call elsi_trace_mat_mat_real(bh,dm,ovlp_copy,c2)

         c = c1/c2

         if(c <= 0.5) then
            dm = ((1-2*c)*ham+(1+c)*dsd-dsdsd)/(1-c)
         else
            dm = ((1+c)*dsd-dsdsd)/c
         endif
      end select

      diff = 0.0_r8

      do j = 1,bh%n_lcol
         do i = 1,bh%n_lrow
            diff = diff+(dm(i,j)-ham(i,j))**2
         enddo
      enddo

      tmp = sqrt(diff)

      call MPI_Allreduce(tmp,diff,1,mpi_real8,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      if(diff < ph%dmp_tol) then
         dmp_conv = .true.
         exit
      endif
   enddo

   call elsi_deallocate(bh,dsd,"dsd")
   if(ph%dmp_method == CANONICAL) then
      call elsi_deallocate(bh,dsdsd,"dsdsd")
   endif

   if(dmp_conv) then
      ! E = Trace(H * DM)
      call elsi_trace_mat_mat_real(bh,ham_copy,dm,ph%ebs)

      ! n_electrons = Trace(S * DM)
      call elsi_trace_mat_mat_real(bh,ovlp_copy,dm,dmp_ne)
      dmp_ne = ph%spin_degen*dmp_ne

      write(info_str,"(2X,A,I10,A)") &
         "Density matrix purification converged in",i_iter," iterations"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,F10.3)") "| Number of electrons :",dmp_ne
      call elsi_say(bh,info_str)
   else
      call elsi_stop(bh,"Density matrix purification failed.",caller)
   endif

   dm = ph%spin_degen*dm

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished density matrix purification"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine cleans up DMP.
!!
subroutine elsi_cleanup_dmp(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   character(len=40), parameter :: caller = "elsi_cleanup_dmp"

   ph%dmp_started = .false.

end subroutine

end module ELSI_DMP
