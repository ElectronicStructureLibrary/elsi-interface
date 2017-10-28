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
!! This module implements the density matrix purification algorithm.
!!
module ELSI_DMP

   use ELSI_CONSTANTS
   use ELSI_DATATYPE
   use ELSI_ELPA,      only: elsi_to_standard_evp
   use ELSI_MALLOC
   use ELSI_MATRICES,  only: elsi_set_full_mat
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS

   implicit none

   private

   public :: elsi_solve_evp_dmp
   public :: elsi_set_dmp_default

contains

!>
!! This routine computes the density matrix using the density matrix
!! purification algorithm.
!!
subroutine elsi_solve_evp_dmp(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: i_iter
   integer(kind=i4) :: uplo_save
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: ierr
   logical          :: ev_min_found
   logical          :: dmp_conv
   real(kind=r8)    :: this_ev
   real(kind=r8)    :: prev_ev
   real(kind=r8)    :: nrm2
   real(kind=r8)    :: diff
   real(kind=r8)    :: mu
   real(kind=r8)    :: lambda
   real(kind=r8)    :: c1
   real(kind=r8)    :: c2
   real(kind=r8)    :: c
   real(kind=r8)    :: tmp
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   real(kind=r8), allocatable :: row_sum(:)
   real(kind=r8), allocatable :: diag(:)
   real(kind=r8), allocatable :: tmp_real1(:)
   real(kind=r8), allocatable :: tmp_real2(:)
   real(kind=r8), allocatable :: dsd(:,:)
   real(kind=r8), allocatable :: dsdsd(:,:)

   character*40, parameter :: caller = "elsi_solve_evp_dmp"

   ! Compute sparsity
   if(e_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz_real(e_h,e_h%ham_real,e_h%n_lrow,e_h%n_lcol,&
              e_h%nnz_l)

      call MPI_Allreduce(e_h%nnz_l,e_h%nnz_g,1,mpi_integer4,mpi_sum,&
              e_h%mpi_comm,mpierr)
   endif

   call elsi_get_time(e_h,t0)

   call elsi_say(e_h,"  Starting density matrix purification")

   ! Transform the generalized evp to the standard form
   e_h%evec_real => e_h%dm_real

   call elsi_to_standard_evp(e_h)

   nullify(e_h%evec_real)

   ! Compute inverse of overlap
   if(e_h%n_elsi_calls == 1) then
      call elsi_allocate(e_h,e_h%ovlp_real_inv,e_h%n_lrow,e_h%n_lcol,&
              "ovlp_real_inv",caller)

      e_h%ovlp_real_inv = e_h%ovlp_real_copy

      call pdpotrf('U',e_h%n_basis,e_h%ovlp_real_inv,1,1,e_h%sc_desc,ierr)
      call pdpotri('U',e_h%n_basis,e_h%ovlp_real_inv,1,1,e_h%sc_desc,ierr)

      uplo_save = e_h%uplo
      e_h%uplo  = UT_MAT

      call elsi_set_full_mat(e_h,e_h%ovlp_real_inv)

      e_h%uplo = uplo_save
   endif

   ! Use Gershgorin's theorem to find the bounds of the spectrum
   if(e_h%n_elsi_calls == 1) then
      call elsi_allocate(e_h,row_sum,e_h%n_basis,"row_sum",caller)
      call elsi_allocate(e_h,diag,e_h%n_basis,"diag",caller)
      call elsi_allocate(e_h,tmp_real1,e_h%n_basis,"tmp_real1",caller)

      do i = 1,e_h%n_basis
         do j = 1,e_h%n_basis
            if(e_h%loc_row(i) > 0 .and. e_h%loc_col(j) > 0) then
               if(i /= j) then
                  row_sum(i) = row_sum(i)+abs(e_h%ham_real(e_h%loc_row(i),&
                                  e_h%loc_col(j)))
               else
                  diag(i) = e_h%ham_real(e_h%loc_row(i),e_h%loc_col(j))
               endif
            endif
         enddo
      enddo

      tmp_real1 = row_sum

      call MPI_Allreduce(tmp_real1,row_sum,e_h%n_basis,mpi_real8,mpi_sum,&
              e_h%mpi_comm,mpierr)

      tmp_real1 = diag

      call MPI_Allreduce(tmp_real1,diag,e_h%n_basis,mpi_real8,mpi_sum,&
              e_h%mpi_comm,mpierr)

      e_h%ev_ham_min = minval(diag-row_sum)
      e_h%ev_ham_max = maxval(diag+row_sum)

      call elsi_deallocate(e_h,row_sum,"row_sum")
      call elsi_deallocate(e_h,diag,"diag")
      call elsi_deallocate(e_h,tmp_real1,"tmp_real1")
   endif

   ! Use power iteration to find the largest in magnitude eigenvalue
   ! Usually this is the smallest, which is a better estimate of ev_ham_min
   ev_min_found = .false.

   call elsi_allocate(e_h,tmp_real1,e_h%n_lrow,"tmp_real1",caller)
   call elsi_allocate(e_h,tmp_real2,e_h%n_lrow,"tmp_real2",caller)

   if(e_h%n_elsi_calls == 1) then
      call elsi_allocate(e_h,e_h%evec1,e_h%n_lrow,"e_h%evec1",caller)
      call elsi_allocate(e_h,e_h%evec2,e_h%n_lrow,"e_h%evec2",caller)

      tmp_real1 = 1.0_r8/sqrt(real(e_h%n_basis,kind=r8))
   else
      tmp_real1 = e_h%evec1
   endif

   prev_ev = 0.0_r8

   do i_iter = 1,e_h%max_power_iter
      call pdgemv('N',e_h%n_basis,e_h%n_basis,1.0_r8,e_h%ham_real,1,1,&
              e_h%sc_desc,tmp_real1,1,1,e_h%sc_desc,1,0.0_r8,tmp_real2,1,1,&
              e_h%sc_desc,1)

      this_ev = 0.0_r8
      nrm2    = 0.0_r8

      do i = 1,e_h%n_basis
         if(e_h%loc_row(i) > 0) then
            this_ev = this_ev+tmp_real1(e_h%loc_row(i))*&
                         tmp_real2(e_h%loc_row(i))
            nrm2    = nrm2+tmp_real2(e_h%loc_row(i))*tmp_real2(e_h%loc_row(i))
         endif
      enddo

      tmp = this_ev

      call MPI_Allreduce(tmp,this_ev,1,mpi_real8,mpi_sum,e_h%mpi_comm,mpierr)

      tmp = nrm2

      call MPI_Allreduce(tmp,nrm2,1,mpi_real8,mpi_sum,e_h%mpi_comm,mpierr)

      tmp_real1 = tmp_real2/sqrt(nrm2)

      if(abs(this_ev-prev_ev) < 1.0e-4_r8) then
         exit
      endif

      prev_ev = this_ev
   enddo

   if(this_ev > 0) then
      e_h%ev_ham_max = this_ev
   else
      e_h%ev_ham_min = this_ev
      ev_min_found   = .true.
   endif

   e_h%evec1 = tmp_real1

   if(e_h%n_elsi_calls == 1) then
      tmp_real1 = 1.0_r8/sqrt(real(e_h%n_basis,kind=r8))
   else
      tmp_real1 = e_h%evec2
   endif

   ! Shift H and use power iteration to find a better estimate of ev_ham_max
   e_h%dm_real = e_h%ham_real

   do i = 1,e_h%n_basis
      if(e_h%loc_row(i) > 0 .and. e_h%loc_col(i) > 0) then
         e_h%ham_real(e_h%loc_row(i),e_h%loc_col(i)) = &
            e_h%ham_real(e_h%loc_row(i),e_h%loc_col(i))+abs(e_h%ev_ham_min)
      endif
   enddo

   prev_ev = 0.0_r8

   do i_iter = 1,e_h%max_power_iter
      call pdgemv('N',e_h%n_basis,e_h%n_basis,1.0_r8,e_h%ham_real,1,1,&
              e_h%sc_desc,tmp_real1,1,1,e_h%sc_desc,1,0.0_r8,tmp_real2,1,1,&
              e_h%sc_desc,1)

      this_ev = 0.0_r8
      nrm2    = 0.0_r8

      do i = 1,e_h%n_basis
         if(e_h%loc_row(i) > 0) then
            this_ev = this_ev+tmp_real1(e_h%loc_row(i))*tmp_real2(e_h%loc_row(i))
            nrm2    = nrm2+tmp_real2(e_h%loc_row(i))*tmp_real2(e_h%loc_row(i))
         endif
      enddo

      tmp = this_ev

      call MPI_Allreduce(tmp,this_ev,1,mpi_real8,mpi_sum,e_h%mpi_comm,mpierr)

      tmp = nrm2

      call MPI_Allreduce(tmp,nrm2,1,mpi_real8,mpi_sum,e_h%mpi_comm,mpierr)

      tmp_real1 = tmp_real2/sqrt(nrm2)

      if(abs(this_ev-prev_ev) < 1.0e-4_r8) then
         exit
      endif

      prev_ev = this_ev
   enddo

   if(ev_min_found .and. this_ev > 0) then
      e_h%ev_ham_max = this_ev-abs(e_h%ev_ham_min)
   elseif(this_ev < 0) then
      e_h%ev_ham_min = this_ev+abs(e_h%ev_ham_max)
   endif

   e_h%evec2    = tmp_real1
   e_h%ham_real = e_h%dm_real

   call elsi_deallocate(e_h,tmp_real1,"tmp_real1")
   call elsi_deallocate(e_h,tmp_real2,"tmp_real2")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished power iteration')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

   ! Initialization
   call elsi_get_time(e_h,t0)

   call elsi_trace_mat(e_h,e_h%ham_real,mu)

   mu     = mu/e_h%n_basis
   lambda = min(e_h%n_states_dmp/(e_h%ev_ham_max-mu),&
               (e_h%n_basis-e_h%n_states_dmp)/(mu-e_h%ev_ham_min))

   call elsi_allocate(e_h,dsd,e_h%n_lrow,e_h%n_lcol,"dsd",caller)
   if(e_h%dmp_method == CANONICAL) then
      call elsi_allocate(e_h,dsdsd,e_h%n_lrow,e_h%n_lcol,"dsdsd",caller)
   endif

   ! ham_real used as tmp after this point
   call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
           e_h%ham_real_copy,1,1,e_h%sc_desc,e_h%ovlp_real_inv,1,1,e_h%sc_desc,&
           0.0_r8,e_h%ham_real,1,1,e_h%sc_desc)
   call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
           e_h%ovlp_real_inv,1,1,e_h%sc_desc,e_h%ham_real,1,1,e_h%sc_desc,&
           0.0_r8,e_h%dm_real,1,1,e_h%sc_desc)

   e_h%dm_real = (mu*e_h%ovlp_real_inv-e_h%dm_real)*lambda/e_h%n_basis
   e_h%dm_real = e_h%dm_real+e_h%ovlp_real_inv*e_h%n_states_dmp/e_h%n_basis

   ! Start main density matrix purification loop
   dmp_conv = .false.

   do i_iter = 1,e_h%max_dmp_iter
      e_h%ham_real = e_h%dm_real

      select case(e_h%dmp_method)
      case(TRACE_CORRECTING)
         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 e_h%ham_real,1,1,e_h%sc_desc,e_h%ovlp_real_copy,1,1,&
                 e_h%sc_desc,0.0_r8,e_h%dm_real,1,1,e_h%sc_desc)

         call elsi_trace_mat(e_h,e_h%dm_real,c1)

         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 e_h%dm_real,1,1,e_h%sc_desc,e_h%ham_real,1,1,e_h%sc_desc,&
                 0.0_r8,dsd,1,1,e_h%sc_desc)

         if(e_h%n_states_dmp-c1 > 0.0_r8) then
            e_h%dm_real = 2*e_h%ham_real-dsd
         else
            e_h%dm_real = dsd
         endif
      case(CANONICAL)
         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 e_h%ovlp_real_copy,1,1,e_h%sc_desc,e_h%ham_real,1,1,&
                 e_h%sc_desc,0.0_r8,e_h%dm_real,1,1,e_h%sc_desc)
         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 e_h%ham_real,1,1,e_h%sc_desc,e_h%dm_real,1,1,e_h%sc_desc,&
                 0.0_r8,dsd,1,1,e_h%sc_desc)
         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 e_h%ovlp_real_copy,1,1,e_h%sc_desc,dsd,1,1,e_h%sc_desc,0.0_r8,&
                 e_h%dm_real,1,1,e_h%sc_desc)
         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 e_h%ham_real,1,1,e_h%sc_desc,e_h%dm_real,1,1,e_h%sc_desc,&
                 0.0_r8,dsdsd,1,1,e_h%sc_desc)

         e_h%dm_real = dsd-dsdsd

         call elsi_trace_mat_mat(e_h,e_h%dm_real,e_h%ovlp_real_copy,c1)

         e_h%dm_real = e_h%ham_real-dsd

         call elsi_trace_mat_mat(e_h,e_h%dm_real,e_h%ovlp_real_copy,c2)

         c = c1/c2

         if(c <= 0.5) then
            e_h%dm_real = ((1-2*c)*e_h%ham_real+(1+c)*dsd-dsdsd)/(1-c)
         else
            e_h%dm_real = ((1+c)*dsd-dsdsd)/c
         endif
      end select

      diff = 0.0_r8

      do j = 1,e_h%n_lcol
         do i = 1,e_h%n_lrow
            diff = diff+(e_h%dm_real(i,j)-e_h%ham_real(i,j))**2
         enddo
      enddo

      tmp = sqrt(diff)

      call MPI_Allreduce(tmp,diff,1,mpi_real8,mpi_sum,e_h%mpi_comm,mpierr)

      if(diff < e_h%dmp_tol) then
         dmp_conv = .true.
         exit
      endif
   enddo

   call elsi_deallocate(e_h,dsd,"dsd")
   if(e_h%dmp_method == CANONICAL) then
      call elsi_deallocate(e_h,dsdsd,"dsdsd")
   endif

   if(dmp_conv) then
      ! E = Trace(H * DM)
      call elsi_trace_mat_mat(e_h,e_h%ham_real_copy,e_h%dm_real,e_h%energy_hdm)

      ! n_electrons = Trace(S * DM)
      call elsi_trace_mat_mat(e_h,e_h%ovlp_real_copy,e_h%dm_real,e_h%ne_dmp)
      e_h%ne_dmp = e_h%spin_degen*e_h%ne_dmp

      call elsi_say(e_h,"  Density matrix purification converged")
      write(info_str,"('  | Number of iterations :',I10)") i_iter
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Number of electrons  :',F10.3)") e_h%ne_dmp
      call elsi_say(e_h,info_str)
   else
      call elsi_stop(" Density matrix purification failed to converge.",e_h,&
              caller)
   endif

   call MPI_Barrier(e_h%mpi_comm,mpierr)

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished density matrix purification')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine sets default DMP parameters.
!!
subroutine elsi_set_dmp_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*40, parameter :: caller = "elsi_set_dmp_default"

   ! Purification method
   e_h%dmp_method = 0

   ! Maximum number of power iterations
   e_h%max_power_iter = 50

   ! Maximum number of purification steps
   e_h%max_dmp_iter = 500

   ! Tolerance for purification
   e_h%dmp_tol = 1.0e-10_r8

end subroutine

end module ELSI_DMP
