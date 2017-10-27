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
   use ELSI_MALLOC
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS
   use ELSI_ELPA,      only: elsi_to_standard_evp

   implicit none

   private

   public :: elsi_solve_evp_dmp
   public :: elsi_set_dmp_default

! TODO
!  - Declarations
!  - H, S, S^-1 need to be stored
!  - dm = spin * dm

contains

!>
!! This routine computes the density matrix using the density matrix
!! purification algorithm.
!!
subroutine elsi_solve_evp_dmp(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: success

   character*40, parameter :: caller = "elsi_solve_evp_dmp"

   call elsi_allocate(e_h,row_sum,e_h%n_basis,"row_sum",caller)
   call elsi_allocate(e_h,diag,e_h%n_basis,"diag",caller)

   ! Transform the generalized evp to the standard form
   call elsi_to_standard_evp(e_h)

   ! Use Gershgorin's theorem to find the bounds of the spectrum
   if(e_h%n_elsi_calls == 1) then
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

      tmp_real = row_sum

      call MPI_Allreduce(tmp_real,row_sum,e_h%n_basis,mpi_real8,mpi_sum,&
              e_h%mpi_comm,mpierr)

      tmp_real = diag

      call MPI_Allreduce(tmp_real,diag,e_h%n_basis,mpi_real8,mpi_sum,&
              e_h%mpi_comm,mpierr)

      e_h%ev_ham_min = minval(diag-row_sum)
      e_h%ev_ham_max = maxval(diag+row_sum)
   endif

   call elsi_deallocate(e_h,row_sum,"row_sum")
   call elsi_deallocate(e_h,diag,"diag")

   ! Use power iteration to find the largest in magnitude eigenvalue
   ! Usually this is the smallest, which is a better estimate of ev_ham_min
   ev_min_found = .false.

   if(.not. allocated(evec1)) then
      call elsi_allocate(e_h,evec1,e_h%n_lrow,"evec1",caller)
      vec1 = 1.0_r8/sqrt(real(e_h%n_basis,kind=r8))
   else
      vec1 = evec1
   endif

   prev_ev = 0.0_r8

   do i_iter = 1,e_h%max_power_iter
      call pdgemv('N',e_h%n_basis,e_h%n_basis,1.0_r8,e_h%ham_real,1,1,&
              e_h%sc_desc,vec1,1,1,e_h%sc_desc,1,0.0_r8,vec2,1,1,e_h%sc_desc,1)

      this_ev = 0.0_r8
      nrm2    = 0.0_r8

      do i = 1,e_h%n_basis
         if(e_h%loc_row(i) > 0) then
            this_ev = this_ev+vec1(e_h%loc_row(i))*vec2(e_h%loc_row(i))
            nrm2    = nrm2+vec2(e_h%loc_row(i))*vec2(e_h%loc_row(i))
         endif
      enddo

      tmp_real = this_ev

      call MPI_Allreduce(tmp_real,this_ev,e_h%n_basis,mpi_real8,mpi_sum,&
              e_h%mpi_comm,mpierr)

      tmp_real = nrm2

      call MPI_Allreduce(tmp_real,nrm2,e_h%n_basis,mpi_real8,mpi_sum,&
              e_h%mpi_comm,mpierr)

      vec1 = vec2/sqrt(nrm2)

      if(abs(this_ev-prev_ev) < 1.0e-4_r8) then
         exit
      endif

      prev_ev = this_ev
   enddo

   if(this_ev > 0) then
      e_h%ev_ham_max = this_ev
   else
      e_h%ev_ham_min = this_ev
      ev_min_found = .true.
   endif

   evec1 = vec1

   if(.not. allocated(dm_temp)) then
      call elsi_allocate(e_h,dm_temp,e_h%n_lrow,e_h%n_lcol,"dm_temp",caller)
   endif

   if(.not. allocated(evec2)) then
      call elsi_allocate(e_h,evec2,e_h%n_lrow,"evec2",caller)
      vec1 = 1.0_r8/sqrt(dble(e_h%n_basis))
   else
      vec1 = evec2
   endif

   ! Shift H and use power iteration to find a better estimate of ev_ham_max
   dm_temp = e_h%ham_real

   do i = 1,e_h%n_basis
      if(e_h%loc_row(i) > 0 .and. e_h%loc_col(i) > 0) then
         dm_temp(e_h%loc_row(i),e_h%loc_col(i)) = &
            dm_temp(e_h%loc_row(i),e_h%loc_col(i))+abs(e_h%ev_ham_min)
      endif
   enddo

   prev_ev = 0.0_r8

   do i_iter = 1,e_h%max_power_iter
      call pdgemv('N',e_h%n_basis,e_h%n_basis,1.0_r8,dm_temp,1,1,e_h%sc_desc,&
              vec1,1,1,e_h%sc_desc,1,0.0_r8,vec2,1,1,e_h%sc_desc,1)

      this_ev = 0.0_r8
      nrm2    = 0.0_r8

      do i = 1,e_h%n_basis
         if(e_h%loc_row(i) > 0) then
            this_ev = this_ev+vec1(e_h%loc_row(i))*vec2(e_h%loc_row(i))
            nrm2    = nrm2+vec2(e_h%loc_row(i))*vec2(e_h%loc_row(i))
         endif
      enddo

      tmp_real = this_ev

      call MPI_Allreduce(tmp_real,this_ev,e_h%n_basis,mpi_real8,mpi_sum,&
              e_h%mpi_comm,mpierr)

      tmp_real = nrm2

      call MPI_Allreduce(tmp_real,nrm2,e_h%n_basis,mpi_real8,mpi_sum,&
              e_h%mpi_comm,mpierr)

      vec1 = vec2/sqrt(nrm2)

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

   evec2 = vec1

   ! Purification parameters
   mu     = mat_trace(e_h%ham_real)/e_h%n_basis
   lambda = min(e_h%n_states_dmp/(e_h%ev_ham_max-mu),&
               (e_h%n_basis-e_h%n_states_dmp)/(mu-e_h%ev_ham_min))

   if(.not. allocated(dm2)) then
      call elsi_allocate(e_h,dm2,e_h%n_lrow,e_h%n_lcol,"dm2",caller)
   endif
   if(.not. allocated(dm3)) then
      call elsi_allocate(e_h,dm3,e_h%n_lrow,e_h%n_lcol,"dm3",caller)
   endif
   if(.not. allocated(dm_prev)) then
      call elsi_allocate(e_h,dm_prev,e_h%n_lrow,e_h%n_lcol,"dm_prev",caller)
   endif
   if(.not. allocated(dm_temp)) then
      call elsi_allocate(e_h,dm_temp,e_h%n_lrow,e_h%n_lcol,"dm_temp",caller)
   endif
   if(.not. allocated(dm_temp2)) then
      call elsi_allocate(e_h,dm_temp2,e_h%n_lrow,e_h%n_lcol,"dm_temp2",caller)
   endif

   ! Initialize density matrix
   call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
           e_h%ham_real_copy,1,1,e_h%sc_desc,e_h%ovlp_real_inv,1,1,e_h%sc_desc,&
           0.0_r8,dm_temp,1,1,e_h%sc_desc)
   call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
           e_h%ovlp_real_inv,1,1,e_h%sc_desc,dm_temp,1,1,e_h%sc_desc,0.0_r8,&
           e_h%dm_real,1,1,e_h%sc_desc)

   e_h%dm_real = (mu*e_h%ovlp_real_inv-e_h%dm_real)*lambda/e_h%n_basis
   e_h%dm_real = e_h%dm_real+e_h%ovlp_real_inv*e_h%n_states_dmp/e_h%n_basis

   ! Start main density matrix purification loop
   dmp_conv = .false.

   do i_iter = 1,e_h%max_dmp_iter
      dm_prev = e_h%dm_real

      select case(e_h%dmp_method)
      case(TRACE_CORRECTING) ! J. Chem. Phys. 123, 044107 (2005)
         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 e_h%dm_real,1,1,e_h%sc_desc,e_h%ovlp_real_copy,1,1,&
                 e_h%sc_desc,0.0_r8,dm_temp,1,1,e_h%sc_desc)

         c1 = trace(dm_temp)

         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 dm_temp,1,1,e_h%sc_desc,e_h%dm_real,1,1,e_h%sc_desc,0.0_r8,&
                 dm_temp2,1,1,e_h%sc_desc)

         if(e_h%n_states_dmp-c1 > 0) then
            e_h%dm_real = 2.0*e_h%dm_real-dm_temp2
         else
            e_h%dm_real = dm_temp2
         endif

      case(CANONICAL) ! Phys. Rev. B 58, 12704 (1998)
         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 e_h%ovlp_real_copy,1,1,e_h%sc_desc,e_h%dm_real,1,1,&
                 e_h%sc_desc,0.0_r8,dm_temp,1,1,e_h%sc_desc)
         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 e_h%dm_real,1,1,e_h%sc_desc,dm_temp,1,1,e_h%sc_desc,0.0_r8,&
                 dm2,1,1,e_h%sc_desc)
         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 e_h%ovlp_real_copy,1,1,e_h%sc_desc,dm2,1,1,e_h%sc_desc,&
                 0.0_r8,dm_temp,1,1,e_h%sc_desc)
         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 e_h%dm_real,1,1,e_h%sc_desc,dm_temp,1,1,e_h%sc_desc,0.0_r8,&
                 dm3,1,1,e_h%sc_desc)

         dm_temp = dm2-dm3

         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 dm_temp,1,1,e_h%sc_desc,e_h%ovlp_real_copy,1,1,e_h%sc_desc,&
                 0.0_r8,dm_temp2,1,1,e_h%sc_desc)

         c1 = trace_scalapack(dm_temp2)

         dm_temp = e_h%dm_real-dm2

         call pdgemm('N','N',e_h%n_basis,e_h%n_basis,e_h%n_basis,1.0_r8,&
                 dm_temp,1,1,e_h%sc_desc,e_h%ovlp_real_copy,1,1,e_h%sc_desc,&
                 0.0_r8,dm_temp2,1,1,e_h%sc_desc)

         c2 = trace_scalapack(dm_temp2)

         c = c1/c2
         if(c <= 0.5) then
            e_h%dm_real = ((1-2*c)*e_h%dm_real+(1+c)*dm2-dm3)/(1-c)
         else
            e_h%dm_real = ((1+c)*dm2-dm3)/c
         endif

      end select

      conv = 0.0_r8

      do i = 1,e_h%n_basis
         do j = 1,e_h%n_basis
            if(e_h%loc_row(i) > 0 .and. e_h%loc_col(j) > 0) then
               conv = conv+(e_h%dm_real(e_h%loc_row(i),e_h%loc_col(j))-&
                         dm_prev(e_h%loc_row(i),e_h%loc_col(j)))**2
            endif
         enddo
      enddo

      call sync_real_number(conv)
      conv = sqrt(conv)

      if(conv < e_h%dmp_tol) then
         dmp_conv = .true.
         exit
      endif
   enddo

   if(dmp_conv) then
      ! Energy = Trace(H * DM)
      tmp_real = ddot(e_h%n_lrow*e_h%n_lcol,e_h%ham_real_copy,1,e_h%dm_real,1)

      call MPI_Allreduce(tmp_real,e_h%energy_hdm,1,mpi_real8,mpi_sum,0,&
              e_h%mpi_comm,mpierr)

      ! n_electrons = Trace(S * DM)
      tmp_real = ddot(e_h%n_lrow*e_h%n_lcol,e_h%ovlp_real_copy,1,e_h%dm_real,1)

      call MPI_Allreduce(tmp_real,e_h%ne_dmp,1,mpi_real8,mpi_sum,0,&
              e_h%mpi_comm,mpierr)

      call elsi_say(e_h,"  Density matrix purification converged")
      write(info_str,"('  | Number of iterations :',I10)") i_iter
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Number of electrons  :',F10.3)") e_h%ne_dmp
      call elsi_say(e_h,info_str)
   else
      call elsi_stop(" Density matrix purification failed to converge.",e_h,&
              caller)
   endif

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
   max_power_iter = 50

   ! Maximum number of purification steps
   max_dmp_iter = 200

   ! Tolerance for purification
   dmp_tol = 1.0e-10_r8

end subroutine

end module ELSI_DMP
