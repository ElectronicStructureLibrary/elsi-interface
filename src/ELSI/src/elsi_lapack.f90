! Copyright (c) 2015-2017, the ELSI team. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!  * Neither the name of the "ELectronic Structure Infrastructure" project nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!>
!! This module interfaces to LAPACK routines that solve a generalized
!! eigenproblem in serial.
!!
module ELSI_LAPACK

   use ELSI_CONSTANTS, only: REAL_VALUES,COMPLEX_VALUES
   use ELSI_DATATYPE
   use ELSI_MALLOC
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS
   use ELPA1,          only: elpa_solve_tridi_double

   implicit none

   private

   public :: elsi_solve_evp_lapack

contains

!>
!! This routine transforms a generalized eigenproblem to standard and returns
!! the Cholesky factor for later use.
!!
subroutine elsi_to_standard_evp_sp(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: nwork
   integer(kind=i4) :: n
   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: ierr
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), parameter :: nblk = 128
   character*40,     parameter :: caller = "elsi_to_standard_evp_sp"

   select case(e_h%data_type)
   case(COMPLEX_VALUES)
      if(e_h%check_sing) then
         call elsi_check_singularity_sp(e_h)
      endif

      if(e_h%n_nonsing == e_h%n_basis) then ! Not singular
         call elsi_get_time(e_h,t0)

         e_h%ovlp_is_sing = .false.

         ! Erase the lower triangle
         do i = 1,e_h%n_basis
            do j= 1,i-1
               e_h%ovlp_cmplx(i,j) = (0.0_r8,0.0_r8)
            enddo
         enddo

         ! Compute S = (U^H)U, U -> S
         call zpotrf('U',e_h%n_basis,e_h%ovlp_cmplx,e_h%n_basis,ierr)

         ! compute U^-1 -> S
         call ztrtri('U','N',e_h%n_basis,e_h%ovlp_cmplx,e_h%n_basis,ierr)

         call elsi_get_time(e_h,t1)

         write(info_str,"('  Finished Cholesky decomposition')")
         call elsi_say(e_h,info_str)
         write(info_str,"('  | Time :',F10.3,' s')") t1-t0
         call elsi_say(e_h,info_str)
      endif

      call elsi_get_time(e_h,t0)

      if(e_h%ovlp_is_sing) then ! Use scaled eigenvectors
         ! evec_cmplx used as tmp_cmplx
         ! tmp_cmplx = H_cmplx * S_cmplx
         call zgemm('N','N',e_h%n_basis,e_h%n_nonsing,e_h%n_basis,&
                 (1.0_r8,0.0_r8),e_h%ham_cmplx(1,1),e_h%n_basis,&
                 e_h%ovlp_cmplx(1,1),e_h%n_basis,(0.0_r8,0.0_r8),&
                 e_h%evec_cmplx(1,1),e_h%n_basis)

         ! H_cmplx = (S_cmplx)^* * tmp_cmplx
         call zgemm('C','N',e_h%n_nonsing,e_h%n_nonsing,e_h%n_basis,&
                 (1.0_r8,0.0_r8),e_h%ovlp_cmplx(1,1),e_h%n_basis,&
                 e_h%evec_cmplx(1,1),e_h%n_basis,(0.0_r8,0.0_r8),&
                 e_h%ham_cmplx(1,1),e_h%n_basis)
      else ! Use cholesky
         ! tmp_cmplx = H_cmplx * S_cmplx
         do n = 1,e_h%n_basis,nblk
            nwork = nblk

            if(n+nwork-1 > e_h%n_basis) then
               nwork = e_h%n_basis-n+1
            endif

            call zgemm('N','N',n+nwork-1,nwork,n+nwork-1,(1.0_r8,0.0_r8),&
                    e_h%ham_cmplx(1,1),e_h%n_basis,e_h%ovlp_cmplx(1,n),&
                    e_h%n_basis,(0.0_r8,0.0_r8),e_h%evec_cmplx(1,n),e_h%n_basis)
         enddo

         ! H_cmplx = (tmp_cmplx)^* * S_cmplx
         do n = 1,e_h%n_basis,nblk
            nwork = nblk

            if(n+nwork-1 > e_h%n_basis) then
               nwork = e_h%n_basis-n+1
            endif

            call zgemm('C','N',nwork,e_h%n_basis-n+1,n+nwork-1,(1.0_r8,0.0_r8),&
                    e_h%ovlp_cmplx(1,n),e_h%n_basis,e_h%evec_cmplx(1,n),&
                    e_h%n_basis,(0.0_r8,0.0_r8),e_h%ham_cmplx(n,n),e_h%n_basis)
         enddo
      endif

      call elsi_get_time(e_h,t1)

      write(info_str,"('  Finished transformation to standard eigenproblem')")
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_say(e_h,info_str)
   case(REAL_VALUES)
      if(e_h%check_sing) then
         call elsi_check_singularity_sp(e_h)
      endif

      if(e_h%n_nonsing == e_h%n_basis) then ! Not singular
         call elsi_get_time(e_h,t0)

         e_h%ovlp_is_sing = .false.

         ! Erase the lower triangle
         do i = 1,e_h%n_basis
            do j= 1,i-1
               e_h%ovlp_real(i,j) = 0.0_r8
            enddo
         enddo

         ! Compute S = (U^T)U, U -> S
         call dpotrf('U',e_h%n_basis,e_h%ovlp_real,e_h%n_basis,ierr)

         ! compute U^-1 -> S
         call dtrtri('U','N',e_h%n_basis,e_h%ovlp_real,e_h%n_basis,ierr)

         call elsi_get_time(e_h,t1)

         write(info_str,"('  Finished Cholesky decomposition')")
         call elsi_say(e_h,info_str)
         write(info_str,"('  | Time :',F10.3,' s')") t1-t0
         call elsi_say(e_h,info_str)
      endif

      call elsi_get_time(e_h,t0)

      if(e_h%ovlp_is_sing) then ! Use scaled eigenvectors
         ! evec_real used as tmp_real
         ! tmp_real = H_real * S_real
         call dgemm('N','N',e_h%n_basis,e_h%n_nonsing,e_h%n_basis,1.0_r8,&
                 e_h%ham_real(1,1),e_h%n_basis,e_h%ovlp_real(1,1),e_h%n_basis,&
                 0.0_r8,e_h%evec_real(1,1),e_h%n_basis)

         ! H_real = (S_real)^T * tmp_real
         call dgemm('T','N',e_h%n_nonsing,e_h%n_nonsing,e_h%n_basis,1.0_r8,&
                 e_h%ovlp_real(1,1),e_h%n_basis,e_h%evec_real(1,1),e_h%n_basis,&
                 0.0_r8,e_h%ham_real(1,1),e_h%n_basis)
      else ! Use Cholesky
         ! tmp_real = H_real * S_real
         do n = 1,e_h%n_basis,nblk
            nwork = nblk

            if(n+nwork-1 > e_h%n_basis) then
               nwork = e_h%n_basis-n+1
            endif

            call dgemm('N','N',n+nwork-1,nwork,n+nwork-1,1.0_r8,&
                    e_h%ham_real(1,1),e_h%n_basis,e_h%ovlp_real(1,n),&
                    e_h%n_basis,0.0_r8,e_h%evec_real(1,n),e_h%n_basis)
         enddo

         ! H_real = (tmp_real)*T * S_real
         do n = 1,e_h%n_basis,nblk
            nwork = nblk

            if(n+nwork-1 > e_h%n_basis) then
               nwork = e_h%n_basis-n+1
            endif

            call dgemm('T','N',nwork,e_h%n_basis-n+1,n+nwork-1,1.0_r8,&
                    e_h%ovlp_real(1,n),e_h%n_basis,e_h%evec_real(1,n),&
                    e_h%n_basis,0.0_r8,e_h%ham_real(n,n),e_h%n_basis)
         enddo
      endif

      call elsi_get_time(e_h,t1)

      write(info_str,"('  Finished transformation to standard eigenproblem')")
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_say(e_h,info_str)
   end select

end subroutine

!>
!! This routine back-transforms eigenvectors in the standard form to the
!! original generalized form.
!!
subroutine elsi_to_original_ev_sp(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character*200 :: info_str

   real(kind=r8),    allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character*40, parameter :: caller = "elsi_to_original_ev_sp"

   call elsi_get_time(e_h,t0)

   select case(e_h%data_type)
   case(COMPLEX_VALUES)
      if(e_h%ovlp_is_sing) then
         call elsi_allocate(e_h,tmp_cmplx,e_h%n_lrow,e_h%n_lcol,"tmp_cmplx",&
                 caller)
         tmp_cmplx = e_h%evec_cmplx

         ! Transform matrix is stored in S_cmplx after elsi_to_standard_evp
         call zgemm('N','N',e_h%n_basis,e_h%n_states_solve,e_h%n_nonsing,&
                 (1.0_r8,0.0_r8),e_h%ovlp_cmplx(1,1),e_h%n_basis,&
                 tmp_cmplx(1,1),e_h%n_basis,(0.0_r8,0.0_r8),&
                 e_h%evec_cmplx(1,1),e_h%n_basis)

         call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")
      else ! Nonsingular, use Cholesky
         ! (U^-1) is stored in S_cmplx after elsi_to_standard_evp
         call ztrmm('L','U','N','N',e_h%n_basis,e_h%n_states,(1.0_r8,0.0_r8),&
                 e_h%ovlp_cmplx(1,1),e_h%n_basis,e_h%evec_cmplx(1,1),&
                 e_h%n_basis)
      endif
   case(REAL_VALUES)
      if(e_h%ovlp_is_sing) then
         call elsi_allocate(e_h,tmp_real,e_h%n_lrow,e_h%n_lcol,"tmp_real",&
                 caller)
         tmp_real = e_h%evec_real

         ! Transform matrix is stored in S_real after elsi_to_standard_evp
         call dgemm('N','N',e_h%n_basis,e_h%n_states_solve,e_h%n_nonsing,&
                 1.0_r8,e_h%ovlp_real(1,1),e_h%n_basis,tmp_real(1,1),&
                 e_h%n_basis,0.0_r8,e_h%evec_real(1,1),e_h%n_basis)

         call elsi_deallocate(e_h,tmp_real,"tmp_real")
      else ! Nonsingular, use Cholesky
         ! (U^-1) is stored in S_real after elsi_to_standard_evp
         call dtrmm('L','U','N','N',e_h%n_basis,e_h%n_states,1.0_r8,&
                 e_h%ovlp_real(1,1),e_h%n_basis,e_h%evec_real(1,1),e_h%n_basis)
      endif
   end select

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished back-transformation of eigenvectors')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine interfaces to LAPACK.
!!
subroutine elsi_solve_evp_lapack(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   real(kind=r8),    allocatable :: off_diag(:)
   real(kind=r8),    allocatable :: tau_real(:)
   complex(kind=r8), allocatable :: tau_cmplx(:)
   real(kind=r8),    allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   logical          :: success
   integer(kind=i4) :: ierr
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_solve_evp_lapack"

   ! Transform to standard form
   if(.not. e_h%ovlp_is_unit) then
      call elsi_to_standard_evp_sp(e_h)
   endif

   ! Solve evp, return eigenvalues and eigenvectors
   call elsi_say(e_h,"  Starting LAPACK eigensolver")
   call elsi_get_time(e_h,t0)

   select case(e_h%data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(e_h,off_diag,e_h%n_nonsing,"off_diag",caller)
      call elsi_allocate(e_h,tau_cmplx,e_h%n_nonsing,"tau_cmplx",caller)
      call elsi_allocate(e_h,tmp_real,e_h%n_nonsing,e_h%n_nonsing,"tmp_real",&
              caller)
      call elsi_allocate(e_h,tmp_cmplx,e_h%n_nonsing,e_h%n_nonsing,"tmp_cmplx",&
              caller)

      call zhetrd('U',e_h%n_nonsing,e_h%ham_cmplx,e_h%n_basis,e_h%eval,&
              off_diag,tau_cmplx,tmp_cmplx,e_h%n_nonsing*e_h%n_nonsing,ierr)

      success = elpa_solve_tridi_double(e_h%n_nonsing,e_h%n_states_solve,&
                   e_h%eval,off_diag,tmp_real,e_h%n_nonsing,64,e_h%n_nonsing,&
                   mpi_comm_self,mpi_comm_self,.false.)

      if(.not. success) then
         call elsi_stop(" ELPA tridiagonal solver failed.",e_h,caller)
      endif

      e_h%evec_cmplx(1:e_h%n_nonsing,1:e_h%n_states_solve) = &
         tmp_real(1:e_h%n_nonsing,1:e_h%n_states_solve)

      call zunmtr('L','U','N',e_h%n_nonsing,e_h%n_states_solve,e_h%ham_cmplx,&
              e_h%n_basis,tau_cmplx,e_h%evec_cmplx,e_h%n_basis,tmp_cmplx,&
              e_h%n_nonsing*e_h%n_nonsing,ierr)

      call elsi_deallocate(e_h,off_diag,"off_diag")
      call elsi_deallocate(e_h,tau_cmplx,"tau_cmplx")
      call elsi_deallocate(e_h,tmp_real,"tmp_real")
      call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")
   case(REAL_VALUES)
      call elsi_allocate(e_h,off_diag,e_h%n_nonsing,"off_diag",caller)
      call elsi_allocate(e_h,tau_real,e_h%n_nonsing,"tau_real",caller)
      call elsi_allocate(e_h,tmp_real,e_h%n_nonsing,e_h%n_nonsing,"tmp_real",&
              caller)

      call dsytrd('U',e_h%n_nonsing,e_h%ham_real,e_h%n_basis,e_h%eval,off_diag,&
              tau_real,tmp_real,e_h%n_nonsing*e_h%n_nonsing,ierr)

      success = elpa_solve_tridi_double(e_h%n_nonsing,e_h%n_states_solve,&
                   e_h%eval,off_diag,tmp_real,e_h%n_nonsing,64,e_h%n_nonsing,&
                   mpi_comm_self,mpi_comm_self,.false.)

      if(.not. success) then
         call elsi_stop(" ELPA tridiagonal solver failed.",e_h,caller)
      endif

      e_h%evec_real(1:e_h%n_nonsing,1:e_h%n_states_solve) = &
         tmp_real(1:e_h%n_nonsing,1:e_h%n_states_solve)

      call dormtr('L','U','N',e_h%n_nonsing,e_h%n_states_solve,e_h%ham_real,&
              e_h%n_basis,tau_real,e_h%evec_real,e_h%n_basis,tmp_real,&
              e_h%n_nonsing*e_h%n_nonsing,ierr)

      call elsi_deallocate(e_h,off_diag,"off_diag")
      call elsi_deallocate(e_h,tau_real,"tau_real")
      call elsi_deallocate(e_h,tmp_real,"tmp_real")
   end select

   ! Overwrite zero eigenvalues
   if(e_h%n_nonsing < e_h%n_basis) then
      e_h%eval(e_h%n_nonsing+1:e_h%n_basis) = e_h%eval(e_h%n_nonsing)+10.0_r8
   endif

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished solving standard eigenproblem')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

   ! Back-transform eigenvectors
   if(.not. e_h%ovlp_is_unit) then
      call elsi_to_original_ev_sp(e_h)
   endif

end subroutine

!>
!! This routine checks the singularity of overlap matrix by computing all its
!! eigenvalues. On exit, S is not modified if not singular, or is overwritten by
!! scaled eigenvectors if singular, which can be esed to transform the
!! generalized eigenproblem to the standard form.
!!
subroutine elsi_check_singularity_sp(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: i
   integer(kind=i4) :: ierr
   real(kind=r8)    :: ev_sqrt
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   logical          :: success
   character*200    :: info_str

   real(kind=r8),    allocatable :: off_diag(:)
   real(kind=r8),    allocatable :: tau_real(:)
   complex(kind=r8), allocatable :: tau_cmplx(:)
   real(kind=r8),    allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_cmplx(:,:)
   real(kind=r8),    allocatable :: copy_real(:,:)
   complex(kind=r8), allocatable :: copy_cmplx(:,:)

   character*40, parameter :: caller = "elsi_check_singularity_sp"

   call elsi_get_time(e_h,t0)

   select case(e_h%data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(e_h,copy_cmplx,e_h%n_lrow,e_h%n_lcol,"copy_cmplx",&
              caller)

      ! Use copy_cmplx to store overlap matrix, otherwise it will
      ! be destroyed by eigenvalue calculation
      copy_cmplx = -e_h%ovlp_cmplx

      call elsi_allocate(e_h,off_diag,e_h%n_basis,"off_diag",caller)
      call elsi_allocate(e_h,tau_cmplx,e_h%n_basis,"tau_cmplx",caller)
      call elsi_allocate(e_h,tmp_real,e_h%n_basis,e_h%n_basis,"tmp_real",caller)
      call elsi_allocate(e_h,tmp_cmplx,e_h%n_basis,e_h%n_basis,"tmp_cmplx",&
              caller)

      call zhetrd('U',e_h%n_basis,copy_cmplx,e_h%n_basis,e_h%eval,off_diag,&
              tau_cmplx,tmp_cmplx,e_h%n_basis*e_h%n_basis,ierr)

      success = elpa_solve_tridi_double(e_h%n_basis,e_h%n_basis,e_h%eval,&
                   off_diag,tmp_real,e_h%n_basis,64,e_h%n_basis,mpi_comm_self,&
                   mpi_comm_self,.false.)

      if(.not. success) then
         call elsi_stop(" ELPA tridiagonal solver failed.",e_h,caller)
      endif

      ! Get the number of nonsingular eigenvalues
      e_h%eval = -e_h%eval

      do i = 1,e_h%n_basis
         if(e_h%eval(i) < e_h%sing_tol) exit
      enddo

      e_h%n_nonsing = i-1

      ! Eigenvectors computed only for singular overlap matrix
      if(e_h%n_nonsing < e_h%n_basis) then
         e_h%evec_cmplx = tmp_real

         call zunmtr('L','U','N',e_h%n_basis,e_h%n_basis,copy_cmplx,&
                 e_h%n_basis,tau_cmplx,e_h%evec_cmplx,e_h%n_basis,tmp_cmplx,&
                 e_h%n_basis*e_h%n_basis,ierr)
      endif

      call elsi_deallocate(e_h,off_diag,"off_diag")
      call elsi_deallocate(e_h,tau_cmplx,"tau_cmplx")
      call elsi_deallocate(e_h,tmp_real,"tmp_real")
      call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")
      call elsi_deallocate(e_h,copy_cmplx,"copy_cmplx")

      e_h%n_states_solve = min(e_h%n_nonsing,e_h%n_states)

      if(e_h%n_nonsing < e_h%n_basis) then ! Singular
         e_h%ovlp_is_sing = .true.

         call elsi_say(e_h,"  Overlap matrix is singular")
         write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)")&
            e_h%eval(e_h%n_basis)
         call elsi_say(e_h,info_str)
         write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
            e_h%eval(1)
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
         do i = 1,e_h%n_nonsing
            ev_sqrt = sqrt(e_h%eval(i))
            e_h%ovlp_cmplx(:,i) = e_h%evec_cmplx(:,i)/ev_sqrt
         enddo
      else ! Nonsingular
         e_h%ovlp_is_sing = .false.
         call elsi_say(e_h,"  Overlap matrix is nonsingular")
         write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)")&
            e_h%eval(e_h%n_basis)
         call elsi_say(e_h,info_str)
         write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
            e_h%eval(1)
         call elsi_say(e_h,info_str)
      endif ! Singular overlap?
   case(REAL_VALUES)
      call elsi_allocate(e_h,copy_real,e_h%n_lrow,e_h%n_lcol,"copy_real",caller)

      ! Use copy_real to store overlap matrix, otherwise it will be
      ! destroyed by eigenvalue calculation
      copy_real = -e_h%ovlp_real

      call elsi_allocate(e_h,off_diag,e_h%n_basis,"off_diag",caller)
      call elsi_allocate(e_h,tau_real,e_h%n_basis,"tau_real",caller)
      call elsi_allocate(e_h,tmp_real,e_h%n_basis,e_h%n_basis,"tmp_real",caller)

      call dsytrd('U',e_h%n_basis,copy_real,e_h%n_basis,e_h%eval,off_diag,&
              tau_real,tmp_real,e_h%n_basis*e_h%n_basis,ierr)

      success = elpa_solve_tridi_double(e_h%n_basis,e_h%n_basis,e_h%eval,&
                   off_diag,tmp_real,e_h%n_basis,64,e_h%n_basis,mpi_comm_self,&
                   mpi_comm_self,.false.)

      if(.not. success) then
         call elsi_stop(" ELPA tridiagonal solver failed.",e_h,caller)
      endif

      ! Get the number of nonsingular eigenvalues
      e_h%eval = -e_h%eval

      do i = 1,e_h%n_basis
         if(e_h%eval(i) < e_h%sing_tol) exit
      enddo

      e_h%n_nonsing = i-1

      ! Eigenvectors computed only for singular overlap matrix
      if(e_h%n_nonsing < e_h%n_basis) then
         e_h%evec_real = tmp_real

         call dormtr('L','U','N',e_h%n_basis,e_h%n_basis,copy_real,e_h%n_basis,&
                 tau_real,e_h%evec_real,e_h%n_basis,tmp_real,&
                 e_h%n_basis*e_h%n_basis,ierr)
      endif

      call elsi_deallocate(e_h,off_diag,"off_diag")
      call elsi_deallocate(e_h,tau_real,"tau_real")
      call elsi_deallocate(e_h,tmp_real,"tmp_real")
      call elsi_deallocate(e_h,copy_real,"copy_real")

      e_h%n_states_solve = min(e_h%n_nonsing,e_h%n_states)

      if(e_h%n_nonsing < e_h%n_basis) then ! Singular
         e_h%ovlp_is_sing = .true.

         call elsi_say(e_h,"  Overlap matrix is singular")
         write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)")&
            e_h%eval(e_h%n_basis)
         call elsi_say(e_h,info_str)
         write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
            e_h%eval(1)
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
         do i = 1,e_h%n_nonsing
            ev_sqrt = sqrt(e_h%eval(i))
            e_h%ovlp_real(:,i) = e_h%evec_real(:,i)/ev_sqrt
         enddo
      else ! Nonsingular
         e_h%ovlp_is_sing = .false.
         call elsi_say(e_h,"  Overlap matrix is nonsingular")
         write(info_str,"('  | Lowest eigenvalue of overlap  :',E10.2)")&
            e_h%eval(e_h%n_basis)
         call elsi_say(e_h,info_str)
         write(info_str,"('  | Highest eigenvalue of overlap :',E10.2)")&
            e_h%eval(1)
         call elsi_say(e_h,info_str)
      endif ! Singular overlap?
   end select

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished singularity check of overlap matrix')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

end module ELSI_LAPACK
