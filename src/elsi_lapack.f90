! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module interfaces to LAPACK routines that solve a generalized
!! eigenproblem in serial.
!!
module ELSI_LAPACK

   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_IO, only: elsi_say,elsi_get_time
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI, only: elsi_stop,mpi_comm_self
   use ELSI_PRECISION, only: r8,i4
   use ELPA1, only: elpa_solve_tridi_double

   implicit none

   private

   public :: elsi_solve_lapack

   interface elsi_solve_lapack
      module procedure elsi_solve_lapack_real
      module procedure elsi_solve_lapack_cmplx
   end interface

contains

!>
!! This routine transforms a generalized eigenproblem to standard and returns
!! the Cholesky factor for later use.
!!
subroutine elsi_to_standard_evp_sp_real(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: nwork
   integer(kind=i4) :: n
   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: ierr
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   integer(kind=i4), parameter :: nblk = 128
   character(len=*), parameter :: caller = "elsi_to_standard_evp_sp_real"

   if(ph%check_sing) then
      call elsi_check_singularity_sp_real(ph,bh,ovlp,eval,evec)
   end if

   if(ph%n_good == ph%n_basis) then ! Not singular
      call elsi_get_time(t0)

      ph%ovlp_is_sing = .false.

      ! Erase the lower triangle
      do i = 1,ph%n_basis
         do j= 1,i-1
            ovlp(i,j) = 0.0_r8
         end do
      end do

      ! S = U
      call dpotrf("U",ph%n_basis,ovlp,ph%n_basis,ierr)

      ! S = U^(-1)
      call dtrtri("U","N",ph%n_basis,ovlp,ph%n_basis,ierr)

      call elsi_get_time(t1)

      write(msg,"(2X,A)") "Finished Cholesky decomposition"
      call elsi_say(bh,msg)
      write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,msg)
   end if

   call elsi_get_time(t0)

   ! H = U^(-T) H U^(-1)
   if(ph%ovlp_is_sing) then
      call dgemm("N","N",ph%n_basis,ph%n_good,ph%n_basis,1.0_r8,ham,ph%n_basis,&
              ovlp,ph%n_basis,0.0_r8,evec,ph%n_basis)

      call dgemm("T","N",ph%n_good,ph%n_good,ph%n_basis,1.0_r8,ovlp,ph%n_basis,&
              evec,ph%n_basis,0.0_r8,ham,ph%n_basis)
   else ! Use Cholesky
      do n = 1,ph%n_basis,nblk
         nwork = nblk

         if(n+nwork-1 > ph%n_basis) then
            nwork = ph%n_basis-n+1
         end if

         call dgemm("N","N",n+nwork-1,nwork,n+nwork-1,1.0_r8,ham,ph%n_basis,&
                 ovlp(1,n),ph%n_basis,0.0_r8,evec(1,n),ph%n_basis)
      end do

      do n = 1,ph%n_basis,nblk
         nwork = nblk

         if(n+nwork-1 > ph%n_basis) then
            nwork = ph%n_basis-n+1
         end if

         call dgemm("T","N",nwork,ph%n_basis-n+1,n+nwork-1,1.0_r8,ovlp(1,n),&
                 ph%n_basis,evec(1,n),ph%n_basis,0.0_r8,ham(n,n),ph%n_basis)
      end do
   end if

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished transformation to standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! This routine back-transforms eigenvectors in the standard form to the
!! original generalized form.
!!
subroutine elsi_to_original_ev_sp_real(ph,bh,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   real(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_to_original_ev_sp_real"

   call elsi_get_time(t0)

   if(ph%ovlp_is_sing) then
      call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)
      tmp = evec

      call dgemm("N","N",ph%n_basis,ph%n_states_solve,ph%n_good,1.0_r8,ovlp,&
              ph%n_basis,tmp,ph%n_basis,0.0_r8,evec,ph%n_basis)

      call elsi_deallocate(bh,tmp,"tmp")
   else ! Nonsingular, use Cholesky
      call dtrmm("L","U","N","N",ph%n_basis,ph%n_states,1.0_r8,ovlp,ph%n_basis,&
              evec,ph%n_basis)
   end if

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished back-transformation of eigenvectors"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! This routine interfaces to LAPACK.
!!
subroutine elsi_solve_lapack_real(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   logical :: success
   integer(kind=i4) :: ierr
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   real(kind=r8), allocatable :: off_diag(:)
   real(kind=r8), allocatable :: tau(:)
   real(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_solve_lapack_real"

   ! Transform to standard form
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_standard_evp_sp_real(ph,bh,ham,ovlp,eval,evec)
   end if

   ! Solve evp, return eigenvalues and eigenvectors
   call elsi_get_time(t0)

   write(msg,"(2X,A)") "Starting LAPACK eigensolver"
   call elsi_say(bh,msg)

   call elsi_allocate(bh,off_diag,ph%n_good,"off_diag",caller)
   call elsi_allocate(bh,tau,ph%n_good,"tau",caller)
   call elsi_allocate(bh,tmp,ph%n_good,ph%n_good,"tmp",caller)

   call dsytrd("U",ph%n_good,ham,ph%n_basis,eval,off_diag,tau,tmp,&
           ph%n_good*ph%n_good,ierr)

   success = elpa_solve_tridi_double(ph%n_good,ph%n_states_solve,eval,off_diag,&
                tmp,ph%n_good,64,ph%n_good,mpi_comm_self,mpi_comm_self,.false.)

   if(.not. success) then
      call elsi_stop(bh,"ELPA tridiagonal solver failed.",caller)
   end if

   evec(1:ph%n_good,1:ph%n_states_solve) = tmp(1:ph%n_good,1:ph%n_states_solve)

   call dormtr("L","U","N",ph%n_good,ph%n_states_solve,ham,ph%n_basis,tau,evec,&
           ph%n_basis,tmp,ph%n_good*ph%n_good,ierr)

   call elsi_deallocate(bh,off_diag,"off_diag")
   call elsi_deallocate(bh,tau,"tau")
   call elsi_deallocate(bh,tmp,"tmp")

   ! Overwrite zero eigenvalues
   if(ph%n_good < ph%n_basis) then
      eval(ph%n_good+1:ph%n_basis) = eval(ph%n_good)+10.0_r8
   end if

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished solving standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ! Back-transform eigenvectors
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_original_ev_sp_real(ph,bh,ovlp,evec)
   end if

end subroutine

!>
!! This routine checks the singularity of overlap matrix by computing all its
!! eigenvalues. On exit, S is not modified if not singular, or is overwritten by
!! scaled eigenvectors if singular, which can be used to transform the
!! generalized eigenproblem to the standard form.
!!
subroutine elsi_check_singularity_sp_real(ph,bh,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i
   integer(kind=i4) :: ierr
   real(kind=r8) :: ev_sqrt
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   logical :: success
   character(len=200) :: msg

   real(kind=r8), allocatable :: off_diag(:)
   real(kind=r8), allocatable :: tau(:)
   real(kind=r8), allocatable :: tmp(:,:)
   real(kind=r8), allocatable :: copy(:,:)

   character(len=*), parameter :: caller = "elsi_check_singularity_sp_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,copy,bh%n_lrow,bh%n_lcol,"copy",caller)

   ! Overlap will be destroyed by eigenvalue calculation
   copy = -ovlp

   call elsi_allocate(bh,off_diag,ph%n_basis,"off_diag",caller)
   call elsi_allocate(bh,tau,ph%n_basis,"tau",caller)
   call elsi_allocate(bh,tmp,ph%n_basis,ph%n_basis,"tmp",caller)

   call dsytrd("U",ph%n_basis,copy,ph%n_basis,eval,off_diag,tau,tmp,&
           ph%n_basis*ph%n_basis,ierr)

   success = elpa_solve_tridi_double(ph%n_basis,ph%n_basis,eval,off_diag,&
                tmp,ph%n_basis,64,ph%n_basis,mpi_comm_self,mpi_comm_self,&
                .false.)

   if(.not. success) then
      call elsi_stop(bh,"ELPA tridiagonal solver failed.",caller)
   end if

   ! Get the number of nonsingular eigenvalues
   eval = -eval

   do i = 1,ph%n_basis
      if(eval(i) < ph%sing_tol) then
         exit
      end if
   end do

   ph%n_good = i-1

   ! Eigenvectors computed only for singular overlap matrix
   if(ph%n_good < ph%n_basis) then
      evec = tmp

      call dormtr("L","U","N",ph%n_basis,ph%n_basis,copy,ph%n_basis,tau,evec,&
              ph%n_basis,tmp,ph%n_basis*ph%n_basis,ierr)
   end if

   call elsi_deallocate(bh,off_diag,"off_diag")
   call elsi_deallocate(bh,tau,"tau")
   call elsi_deallocate(bh,tmp,"tmp")
   call elsi_deallocate(bh,copy,"copy")

   ph%n_states_solve = min(ph%n_good,ph%n_states)

   if(ph%n_good < ph%n_basis) then ! Singular
      ph%ovlp_is_sing = .true.

      write(msg,"(2X,A)") "Overlap matrix is singular"
      call elsi_say(bh,msg)
      write(msg,"(2X,A,E10.2,A,E10.2)") "| Lowest and highest eigenvalues :",&
         eval(ph%n_basis),",",eval(1)
      call elsi_say(bh,msg)

      if(ph%stop_sing) then
         call elsi_stop(bh,"Overlap matrix is singular.",caller)
      end if

      write(msg,"(2X,A,I10)") "| Number of basis functions reduced to :",&
         ph%n_good
      call elsi_say(bh,msg)

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = 1,ph%n_good
         ev_sqrt = sqrt(eval(i))
         ovlp(:,i) = evec(:,i)/ev_sqrt
      end do
   else ! Nonsingular
      ph%ovlp_is_sing = .false.

      write(msg,"(2X,A)") "Overlap matrix is not singular"
      call elsi_say(bh,msg)
      write(msg,"(2X,A,E10.2,A,E10.2)") "| Lowest and highest eigenvalues :",&
         eval(ph%n_basis),",",eval(1)
      call elsi_say(bh,msg)
   end if

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished singularity check of overlap matrix"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! This routine transforms a generalized eigenproblem to standard and returns
!! the Cholesky factor for later use.
!!
subroutine elsi_to_standard_evp_sp_cmplx(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: nwork
   integer(kind=i4) :: n
   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: ierr
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   integer(kind=i4), parameter :: nblk = 128
   character(len=*), parameter :: caller = "elsi_to_standard_evp_sp_cmplx"

   if(ph%check_sing) then
      call elsi_check_singularity_sp_cmplx(ph,bh,ovlp,eval,evec)
   end if

   if(ph%n_good == ph%n_basis) then ! Not singular
      call elsi_get_time(t0)

      ph%ovlp_is_sing = .false.

      ! Erase the lower triangle
      do i = 1,ph%n_basis
         do j= 1,i-1
            ovlp(i,j) = (0.0_r8,0.0_r8)
         end do
      end do

      ! S = U
      call zpotrf("U",ph%n_basis,ovlp,ph%n_basis,ierr)

      ! S = U^(-1)
      call ztrtri("U","N",ph%n_basis,ovlp,ph%n_basis,ierr)

      call elsi_get_time(t1)

      write(msg,"(2X,A)") "Finished Cholesky decomposition"
      call elsi_say(bh,msg)
      write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,msg)
   end if

   call elsi_get_time(t0)

   ! H = U^(-T) H U^(-1)
   if(ph%ovlp_is_sing) then
      call zgemm("N","N",ph%n_basis,ph%n_good,ph%n_basis,(1.0_r8,0.0_r8),ham,&
              ph%n_basis,ovlp,ph%n_basis,(0.0_r8,0.0_r8),evec,ph%n_basis)

      call zgemm("C","N",ph%n_good,ph%n_good,ph%n_basis,(1.0_r8,0.0_r8),ovlp,&
              ph%n_basis,evec,ph%n_basis,(0.0_r8,0.0_r8),ham,ph%n_basis)
   else ! Use cholesky
      do n = 1,ph%n_basis,nblk
         nwork = nblk

         if(n+nwork-1 > ph%n_basis) then
            nwork = ph%n_basis-n+1
         end if

         call zgemm("N","N",n+nwork-1,nwork,n+nwork-1,(1.0_r8,0.0_r8),ham,&
                 ph%n_basis,ovlp(1,n),ph%n_basis,(0.0_r8,0.0_r8),evec(1,n),&
                 ph%n_basis)
      end do

      do n = 1,ph%n_basis,nblk
         nwork = nblk

         if(n+nwork-1 > ph%n_basis) then
            nwork = ph%n_basis-n+1
         end if

         call zgemm("C","N",nwork,ph%n_basis-n+1,n+nwork-1,(1.0_r8,0.0_r8),&
                 ovlp(1,n),ph%n_basis,evec(1,n),ph%n_basis,(0.0_r8,0.0_r8),&
                 ham(n,n),ph%n_basis)
      end do
   end if

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished transformation to standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! This routine back-transforms eigenvectors in the standard form to the
!! original generalized form.
!!
subroutine elsi_to_original_ev_sp_cmplx(ph,bh,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   complex(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_to_original_ev_sp_cmplx"

   call elsi_get_time(t0)

   if(ph%ovlp_is_sing) then
      call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)
      tmp = evec

      call zgemm("N","N",ph%n_basis,ph%n_states_solve,ph%n_good,&
              (1.0_r8,0.0_r8),ovlp,ph%n_basis,tmp,ph%n_basis,(0.0_r8,0.0_r8),&
              evec,ph%n_basis)

      call elsi_deallocate(bh,tmp,"tmp")
   else ! Nonsingular, use Cholesky
      call ztrmm("L","U","N","N",ph%n_basis,ph%n_states,(1.0_r8,0.0_r8),ovlp,&
              ph%n_basis,evec,ph%n_basis)
   end if

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished back-transformation of eigenvectors"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! This routine interfaces to LAPACK.
!!
subroutine elsi_solve_lapack_cmplx(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   logical :: success
   integer(kind=i4) :: ierr
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   real(kind=r8), allocatable :: off_diag(:)
   complex(kind=r8), allocatable :: tau(:)
   real(kind=r8), allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character(len=*), parameter :: caller = "elsi_solve_lapack_cmplx"

   ! Transform to standard form
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_standard_evp_sp_cmplx(ph,bh,ham,ovlp,eval,evec)
   end if

   ! Solve evp, return eigenvalues and eigenvectors
   call elsi_get_time(t0)

   write(msg,"(2X,A)") "Starting LAPACK eigensolver"
   call elsi_say(bh,msg)

   call elsi_allocate(bh,off_diag,ph%n_good,"off_diag",caller)
   call elsi_allocate(bh,tau,ph%n_good,"tau",caller)
   call elsi_allocate(bh,tmp_real,ph%n_good,ph%n_good,"tmp_real",caller)
   call elsi_allocate(bh,tmp_cmplx,ph%n_good,ph%n_good,"tmp_cmplx",caller)

   call zhetrd("U",ph%n_good,ham,ph%n_basis,eval,off_diag,tau,tmp_cmplx,&
           ph%n_good*ph%n_good,ierr)

   success = elpa_solve_tridi_double(ph%n_good,ph%n_states_solve,eval,off_diag,&
                tmp_real,ph%n_good,64,ph%n_good,mpi_comm_self,mpi_comm_self,&
                .false.)

   if(.not. success) then
      call elsi_stop(bh,"ELPA tridiagonal solver failed.",caller)
   end if

   evec(1:ph%n_good,1:ph%n_states_solve) =&
      tmp_real(1:ph%n_good,1:ph%n_states_solve)

   call zunmtr("L","U","N",ph%n_good,ph%n_states_solve,ham,ph%n_basis,tau,evec,&
           ph%n_basis,tmp_cmplx,ph%n_good*ph%n_good,ierr)

   call elsi_deallocate(bh,off_diag,"off_diag")
   call elsi_deallocate(bh,tau,"tau")
   call elsi_deallocate(bh,tmp_real,"tmp_real")
   call elsi_deallocate(bh,tmp_cmplx,"tmp_cmplx")

   ! Overwrite zero eigenvalues
   if(ph%n_good < ph%n_basis) then
      eval(ph%n_good+1:ph%n_basis) = eval(ph%n_good)+10.0_r8
   end if

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished solving standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ! Back-transform eigenvectors
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_original_ev_sp_cmplx(ph,bh,ovlp,evec)
   end if

end subroutine

!>
!! This routine checks the singularity of overlap matrix by computing all its
!! eigenvalues. On exit, S is not modified if not singular, or is overwritten by
!! scaled eigenvectors if singular, which can be used to transform the
!! generalized eigenproblem to the standard form.
!!
subroutine elsi_check_singularity_sp_cmplx(ph,bh,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i
   integer(kind=i4) :: ierr
   real(kind=r8) :: ev_sqrt
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   logical :: success
   character(len=200) :: msg

   real(kind=r8), allocatable :: off_diag(:)
   complex(kind=r8), allocatable :: tau(:)
   real(kind=r8), allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_cmplx(:,:)
   complex(kind=r8), allocatable :: copy(:,:)

   character(len=*), parameter :: caller = "elsi_check_singularity_sp_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,copy,bh%n_lrow,bh%n_lcol,"copy",caller)

   ! Overlap will be destroyed by eigenvalue calculation
   copy = -ovlp

   call elsi_allocate(bh,off_diag,ph%n_basis,"off_diag",caller)
   call elsi_allocate(bh,tau,ph%n_basis,"tau",caller)
   call elsi_allocate(bh,tmp_real,ph%n_basis,ph%n_basis,"tmp_real",caller)
   call elsi_allocate(bh,tmp_cmplx,ph%n_basis,ph%n_basis,"tmp_cmplx",caller)

   call zhetrd("U",ph%n_basis,copy,ph%n_basis,eval,off_diag,tau,tmp_cmplx,&
           ph%n_basis*ph%n_basis,ierr)

   success = elpa_solve_tridi_double(ph%n_basis,ph%n_basis,eval,off_diag,&
                tmp_real,ph%n_basis,64,ph%n_basis,mpi_comm_self,mpi_comm_self,&
                .false.)

   if(.not. success) then
      call elsi_stop(bh,"ELPA tridiagonal solver failed.",caller)
   end if

   ! Get the number of nonsingular eigenvalues
   eval = -eval

   do i = 1,ph%n_basis
      if(eval(i) < ph%sing_tol) then
         exit
      end if
   end do

   ph%n_good = i-1

   ! Eigenvectors computed only for singular overlap matrix
   if(ph%n_good < ph%n_basis) then
      evec = tmp_real

      call zunmtr("L","U","N",ph%n_basis,ph%n_basis,copy,ph%n_basis,tau,evec,&
              ph%n_basis,tmp_cmplx,ph%n_basis*ph%n_basis,ierr)
   end if

   call elsi_deallocate(bh,off_diag,"off_diag")
   call elsi_deallocate(bh,tau,"tau")
   call elsi_deallocate(bh,tmp_real,"tmp_real")
   call elsi_deallocate(bh,tmp_cmplx,"tmp_cmplx")
   call elsi_deallocate(bh,copy,"copy")

   ph%n_states_solve = min(ph%n_good,ph%n_states)

   if(ph%n_good < ph%n_basis) then ! Singular
      ph%ovlp_is_sing = .true.

      write(msg,"(2X,A)") "Overlap matrix is singular"
      call elsi_say(bh,msg)
      write(msg,"(2X,A,E10.2,A,E10.2)") "| Lowest and highest eigenvalues :",&
         eval(ph%n_basis),",",eval(1)
      call elsi_say(bh,msg)

      if(ph%stop_sing) then
         call elsi_stop(bh,"Overlap matrix is singular.",caller)
      end if

      write(msg,"(2X,A,I10)") "| Number of basis functions reduced to :",&
         ph%n_good
      call elsi_say(bh,msg)

      ! Overlap matrix is overwritten with scaled eigenvectors
      do i = 1,ph%n_good
         ev_sqrt = sqrt(eval(i))
         ovlp(:,i) = evec(:,i)/ev_sqrt
      end do
   else ! Nonsingular
      ph%ovlp_is_sing = .false.

      write(msg,"(2X,A)") "Overlap matrix is not singular"
      call elsi_say(bh,msg)
      write(msg,"(2X,A,E10.2,A,E10.2)") "| Lowest and highest eigenvalues :",&
         eval(ph%n_basis),",",eval(1)
      call elsi_say(bh,msg)
   end if

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished singularity check of overlap matrix"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

end module ELSI_LAPACK
