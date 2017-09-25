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
!! This module provides interfaces to ELPA.
!!
module ELSI_ELPA

   use ELSI_CONSTANTS,    only: REAL_VALUES,COMPLEX_VALUES
   use ELSI_DATATYPE
   use ELSI_MALLOC
   use ELSI_MU,           only: elsi_compute_mu_and_occ
   use ELSI_PRECISION,    only: r8,i4
   use ELSI_UTILS
   use CHECK_SINGULARITY, only: elpa_check_singularity_real_double,&
                                elpa_check_singularity_complex_double
   use ELPA1,             only: elpa_print_times,elpa_get_communicators,&
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

   public :: elsi_get_elpa_comms
   public :: elsi_set_elpa_default
   public :: elsi_compute_occ_elpa
   public :: elsi_compute_dm_elpa
   public :: elsi_compute_edm_elpa
   public :: elsi_solve_evp_elpa

contains

!>
!! This routine gets the row and column communicators for ELPA.
!!
subroutine elsi_get_elpa_comms(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: success

   character*40, parameter :: caller = "elsi_get_elpa_comms"

   success = elpa_get_communicators(e_h%mpi_comm,e_h%my_prow,e_h%my_pcol,&
                e_h%mpi_comm_row,e_h%mpi_comm_col)

   if(success /= 0) then
      call elsi_stop(" Failed to get MPI communicators.",e_h,caller)
   endif

end subroutine

!>
!! This routine computes the chemical potential and occupation numbers.
!!
subroutine elsi_compute_occ_elpa(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   real(kind=r8), allocatable :: tmp_real1(:)
   real(kind=r8), allocatable :: tmp_real2(:,:,:)

   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_compute_occ_elpa"

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
                 mpi_sum,e_h%mpi_comm_all,mpierr)

         call elsi_deallocate(e_h,tmp_real1,"tmp_real")
      else
         e_h%k_weight = e_h%i_weight
      endif
   endif

   if(e_h%n_spins*e_h%n_kpts > 1) then
      call elsi_allocate(e_h,tmp_real2,e_h%n_states,e_h%n_spins,e_h%n_kpts,&
              "tmp_real",caller)

      if(e_h%myid == 0) then
         tmp_real2(:,e_h%i_spin,e_h%i_kpt) = e_h%eval_elpa(1:e_h%n_states)
      endif

      call MPI_Allreduce(tmp_real2,e_h%eval_all,&
              e_h%n_states*e_h%n_spins*e_h%n_kpts,mpi_real8,mpi_sum,&
              e_h%mpi_comm_all,mpierr)

      call elsi_deallocate(e_h,tmp_real2,"tmp_real")
   else
      e_h%eval_all(:,e_h%i_spin,e_h%i_kpt) = e_h%eval_elpa(1:e_h%n_states)
   endif

   ! Calculate mu and occupation numbers
   call elsi_compute_mu_and_occ(e_h,e_h%n_electrons,e_h%n_states,e_h%n_spins,&
           e_h%n_kpts,e_h%k_weight,e_h%eval_all,e_h%occ_num,e_h%mu)

end subroutine

!>
!! This routine constructs the density matrix.
!!
subroutine elsi_compute_dm_elpa(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: i
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_row
   integer(kind=i4) :: max_state
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   real(kind=r8), allocatable :: factor(:)

   character*40, parameter :: caller = "elsi_compute_dm_elpa"

   call elsi_get_time(e_h,t0)

   call elsi_allocate(e_h,factor,e_h%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,e_h%n_states_solve
      if(e_h%occ_num(i,e_h%i_spin,e_h%i_kpt) > 0.0_r8) then
         factor(i) = sqrt(e_h%occ_num(i,e_h%i_spin,e_h%i_kpt))
         max_state = i
      endif
   enddo

   select case(e_h%data_type)
   case(REAL_VALUES)
      e_h%ham_real = e_h%evec_real

      do i = 1,e_h%n_states_solve
         if(factor(i) > 0.0_r8) then
            if(e_h%loc_col(i) > 0) then
               e_h%ham_real(:,e_h%loc_col(i)) = e_h%ham_real(:,e_h%loc_col(i))*&
                                                   factor(i)
            endif
         elseif(e_h%loc_col(i) .ne. 0) then
            e_h%ham_real(:,e_h%loc_col(i)) = 0.0_r8
         endif
      enddo

      e_h%dm_real = 0.0_r8

      ! Compute density matrix
      call pdsyrk('U','N',e_h%n_basis,max_state,1.0_r8,e_h%ham_real,1,1,&
              e_h%sc_desc,0.0_r8,e_h%dm_real,1,1,e_h%sc_desc)

      call elsi_deallocate(e_h,factor,"factor")

      ! Set full matrix from upper triangle
      call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,e_h%dm_real,1,1,e_h%sc_desc,&
              0.0_r8,e_h%ham_real,1,1,e_h%sc_desc)

      do i_col = 1,e_h%n_basis-1
         if(e_h%loc_col(i_col) == 0) cycle

         do i_row = i_col+1,e_h%n_basis
            if(e_h%loc_row(i_row) > 0) then
               e_h%dm_real(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                  e_h%ham_real(e_h%loc_row(i_row),e_h%loc_col(i_col))
            endif
         enddo
      enddo
   case(COMPLEX_VALUES)
      e_h%ham_cmplx = e_h%evec_cmplx

      do i = 1,e_h%n_states_solve
         if(factor(i) > 0.0_r8) then
            if(e_h%loc_col(i) > 0) then
               e_h%ham_cmplx(:,e_h%loc_col(i)) = &
                  e_h%ham_cmplx(:,e_h%loc_col(i))*factor(i)
            endif
         elseif(e_h%loc_col(i) .ne. 0) then
            e_h%ham_cmplx(:,e_h%loc_col(i)) = (0.0_r8,0.0_r8)
         endif
      enddo

      e_h%dm_cmplx = (0.0_r8,0.0_r8)

      ! Compute density matrix
      call pzherk('U','N',e_h%n_basis,max_state,(1.0_r8,0.0_r8),e_h%ham_cmplx,&
              1,1,e_h%sc_desc,(0.0_r8,0.0_r8),e_h%dm_cmplx,1,1,e_h%sc_desc)

      call elsi_deallocate(e_h,factor,"factor")

      ! Set full matrix from upper triangle
      call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),e_h%dm_cmplx,1,1,&
              e_h%sc_desc,(0.0_r8,0.0_r8),e_h%ham_cmplx,1,1,e_h%sc_desc)

      do i_col = 1,e_h%n_basis-1
         if(e_h%loc_col(i_col) == 0) cycle

         do i_row = i_col+1,e_h%n_basis
            if(e_h%loc_row(i_row) > 0) then
               e_h%dm_cmplx(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                  e_h%ham_cmplx(e_h%loc_row(i_row),e_h%loc_col(i_col))
            endif
         enddo
      enddo

      ! Make diagonal real
      do i_col = 1,e_h%n_basis
         if(e_h%loc_col(i_col) == 0 .or. e_h%loc_row(i_col) == 0) cycle

         e_h%dm_cmplx(e_h%loc_row(i_col),e_h%loc_col(i_col)) = &
            real(e_h%dm_cmplx(e_h%loc_row(i_col),e_h%loc_col(i_col)),kind=r8)
      enddo
   end select

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished density matrix calculation')")
   call elsi_say(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,e_h)

end subroutine

!>
!! This routine constructs the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_elpa(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: i
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_row
   integer(kind=i4) :: max_state
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   real(kind=r8),    allocatable :: factor(:)
   real(kind=r8),    allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character*40, parameter :: caller = "elsi_compute_edm_elpa"

   call elsi_get_time(e_h,t0)

   call elsi_allocate(e_h,factor,e_h%n_states_solve,"factor",caller)

   max_state = 0

   do i = 1,e_h%n_states_solve
      factor(i) = -1.0_r8*e_h%occ_num(i,e_h%i_spin,e_h%i_kpt)*e_h%eval(i)
      if(factor(i) > 0.0_r8) then
         factor(i) = sqrt(factor(i))
         max_state = i
      else
         factor(i) = 0.0_r8
      endif
   enddo

   select case(e_h%data_type)
   case(REAL_VALUES)
      call elsi_allocate(e_h,tmp_real,e_h%n_lrow,e_h%n_lcol,"tmp_real",&
              caller)
      tmp_real = e_h%evec_real

      do i = 1,e_h%n_states_solve
         if(factor(i) > 0.0_r8) then
            if(e_h%loc_col(i) > 0) then
               tmp_real(:,e_h%loc_col(i)) = tmp_real(:,e_h%loc_col(i))*factor(i)
            endif
         elseif(e_h%loc_col(i) .ne. 0) then
            tmp_real(:,e_h%loc_col(i)) = 0.0_r8
         endif
      enddo

      call elsi_deallocate(e_h,factor,"factor")

      e_h%dm_real = 0.0_r8

      ! Compute density matrix
      call pdsyrk('U','N',e_h%n_basis,max_state,1.0_r8,tmp_real,1,1,&
              e_h%sc_desc,0.0_r8,e_h%dm_real,1,1,e_h%sc_desc)

      e_h%dm_real = -1.0_r8*e_h%dm_real

      ! Set full matrix from upper triangle
      call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,e_h%dm_real,1,1,e_h%sc_desc,&
              0.0_r8,tmp_real,1,1,e_h%sc_desc)

      do i_col = 1,e_h%n_basis-1
         if(e_h%loc_col(i_col) == 0) cycle

         do i_row = i_col+1,e_h%n_basis
            if(e_h%loc_row(i_row) > 0) then
               e_h%dm_real(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                  tmp_real(e_h%loc_row(i_row),e_h%loc_col(i_col))
            endif
         enddo
      enddo

      call elsi_deallocate(e_h,tmp_real,"tmp_real")
   case(COMPLEX_VALUES)
      call elsi_allocate(e_h,tmp_cmplx,e_h%n_lrow,e_h%n_lcol,"tmp_cmplx",&
              caller)
      tmp_cmplx = e_h%evec_cmplx

      do i = 1,e_h%n_states_solve
         if(factor(i) > 0.0_r8) then
            if(e_h%loc_col(i) > 0) then
               tmp_cmplx(:,e_h%loc_col(i)) = tmp_cmplx(:,e_h%loc_col(i))*&
                                                factor(i)
            endif
         elseif(e_h%loc_col(i) .ne. 0) then
            tmp_cmplx(:,e_h%loc_col(i)) = (0.0_r8,0.0_r8)
         endif
      enddo

      call elsi_deallocate(e_h,factor,"factor")

      e_h%dm_cmplx = (0.0_r8,0.0_r8)

      ! Compute density matrix
      call pzherk('U','N',e_h%n_basis,max_state,(1.0_r8,0.0_r8),tmp_cmplx,1,1,&
              e_h%sc_desc,(0.0_r8,0.0_r8),e_h%dm_cmplx,1,1,e_h%sc_desc)

      e_h%dm_cmplx = (-1.0_r8,0.0_r8)*e_h%dm_cmplx

      ! Set full matrix from upper triangle
      call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),e_h%dm_cmplx,1,1,&
              e_h%sc_desc,(0.0_r8,0.0_r8),tmp_cmplx,1,1,e_h%sc_desc)

      do i_col = 1,e_h%n_basis-1
         if(e_h%loc_col(i_col) == 0) cycle

         do i_row = i_col+1,e_h%n_basis
            if(e_h%loc_row(i_row) > 0) then
               e_h%dm_cmplx(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                  tmp_cmplx(e_h%loc_row(i_row),e_h%loc_col(i_col))
            endif
         enddo
      enddo

      call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")

      ! Make diagonal real
      do i_col = 1,e_h%n_basis
         if(e_h%loc_col(i_col) == 0 .or. e_h%loc_row(i_col) == 0) cycle

         e_h%dm_cmplx(e_h%loc_row(i_col),e_h%loc_col(i_col)) = &
            real(e_h%dm_cmplx(e_h%loc_row(i_col),e_h%loc_col(i_col)),kind=r8)
      enddo
   end select

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished energy density matrix calculation')")
   call elsi_say(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,e_h)

end subroutine

!>
!! This routine transforms a generalized eigenproblem to standard and returns
!! the Cholesky factor for later use.
!!
subroutine elsi_to_standard_evp(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   logical          :: success
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_to_standard_evp"

   select case(e_h%data_type)
   case(COMPLEX_VALUES)
      if(e_h%n_elsi_calls == 1) then
         if(.not. e_h%no_sing_check) then
            call elsi_check_singularity(e_h)
         endif

         if(e_h%n_nonsing == e_h%n_basis) then ! Not singular
            call elsi_get_time(e_h,t0)

            e_h%ovlp_is_sing = .false.

            ! Compute S = (U^T)U, U -> S
            success = elpa_cholesky_complex_double(e_h%n_basis,e_h%ovlp_cmplx,&
                         e_h%n_lrow,e_h%blk_row,e_h%n_lcol,&
                         e_h%mpi_comm_row,e_h%mpi_comm_col,.false.)

            if(.not. success) then
               call elsi_stop(" Cholesky failed.",e_h,caller)
            endif

            ! compute U^-1 -> S
            success = elpa_invert_trm_complex_double(e_h%n_basis,&
                         e_h%ovlp_cmplx,e_h%n_lrow,e_h%blk_row,e_h%n_lcol,&
                         e_h%mpi_comm_row,e_h%mpi_comm_col,.false.)

            if(.not. success) then
               call elsi_stop(" Matrix inversion failed.",e_h,caller)
            endif

            call elsi_get_time(e_h,t1)

            write(info_str,"('  Finished Cholesky decomposition')")
            call elsi_say(info_str,e_h)
            write(info_str,"('  | Time :',F10.3,' s')") t1-t0
            call elsi_say(info_str,e_h)
         endif
      endif ! n_elsi_calls == 1

      call elsi_get_time(e_h,t0)

      if(e_h%ovlp_is_sing) then ! Use scaled eigenvectors
         ! evec_cmplx used as tmp_cmplx
         ! tmp_cmplx = H_cmplx * S_cmplx
         call pzgemm('N','N',e_h%n_basis,e_h%n_nonsing,e_h%n_basis,&
                 (1.0_r8,0.0_r8),e_h%ham_cmplx,1,1,e_h%sc_desc,e_h%ovlp_cmplx,&
                 1,e_h%n_basis-e_h%n_nonsing+1,e_h%sc_desc,(0.0_r8,0.0_r8),&
                 e_h%evec_cmplx,1,1,e_h%sc_desc)

         ! H_cmplx = (S_cmplx)^* * tmp_cmplx
         call pzgemm('C','N',e_h%n_nonsing,e_h%n_nonsing,e_h%n_basis,&
                 (1.0_r8,0.0_r8),e_h%ovlp_cmplx,1,e_h%n_basis-e_h%n_nonsing+1,&
                 e_h%sc_desc,e_h%evec_cmplx,1,1,e_h%sc_desc,(0.0_r8,0.0_r8),&
                 e_h%ham_cmplx,1,1,e_h%sc_desc)
      else ! Use cholesky
         success = elpa_mult_ah_b_complex_double('U','L',e_h%n_basis,&
                      e_h%n_basis,e_h%ovlp_cmplx,e_h%n_lrow,e_h%n_lcol,&
                      e_h%ham_cmplx,e_h%n_lrow,e_h%n_lcol,e_h%blk_row,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%evec_cmplx,&
                      e_h%n_lrow,e_h%n_lcol)

         if(.not. success) then
            call elsi_stop(" Matrix multiplication failed.",e_h,caller)
         endif

         call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),e_h%evec_cmplx,1,&
                 1,e_h%sc_desc,(0.0_r8,0.0_r8),e_h%ham_cmplx,1,1,e_h%sc_desc)

         e_h%evec_cmplx = e_h%ham_cmplx

         success = elpa_mult_ah_b_complex_double('U','U',e_h%n_basis,&
                      e_h%n_basis,e_h%ovlp_cmplx,e_h%n_lrow,e_h%n_lcol,&
                      e_h%evec_cmplx,e_h%n_lrow,e_h%n_lcol,e_h%blk_row,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%ham_cmplx,&
                      e_h%n_lrow,e_h%n_lcol)

         if(.not. success) then
            call elsi_stop(" Matrix multiplication failed.",e_h,caller)
         endif

         call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),e_h%ham_cmplx,1,&
                 1,e_h%sc_desc,(0.0_r8,0.0_r8),e_h%evec_cmplx,1,1,e_h%sc_desc)

         ! Set the lower part from the upper
         do i_col = 1,e_h%n_basis-1
            if(e_h%loc_col(i_col) == 0) cycle

            do i_row = i_col+1,e_h%n_basis
               if(e_h%loc_row(i_row) > 0) then
                  e_h%ham_cmplx(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                     e_h%evec_cmplx(e_h%loc_row(i_row),e_h%loc_col(i_col))
               endif
            enddo
         enddo

         do i_col=1,e_h%n_basis
            if(e_h%loc_col(i_col) == 0 .or. e_h%loc_row(i_col) == 0) cycle

            e_h%ham_cmplx(e_h%loc_row(i_col),e_h%loc_col(i_col)) = &
               real(e_h%ham_cmplx(e_h%loc_row(i_col),e_h%loc_col(i_col)),&
               kind=r8)
         enddo
      endif

      call elsi_get_time(e_h,t1)

      write(info_str,"('  Finished transformation to standard eigenproblem')")
      call elsi_say(info_str,e_h)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_say(info_str,e_h)
   case(REAL_VALUES)
      if(e_h%n_elsi_calls == 1) then
         if(.not. e_h%no_sing_check) then
            call elsi_check_singularity(e_h)
         endif

         if(e_h%n_nonsing == e_h%n_basis) then ! Not singular
            call elsi_get_time(e_h,t0)

            e_h%ovlp_is_sing = .false.

            ! Compute S = (U^T)U, U -> S
            success = elpa_cholesky_real_double(e_h%n_basis,e_h%ovlp_real,&
                         e_h%n_lrow,e_h%blk_row,e_h%n_lcol,&
                         e_h%mpi_comm_row,e_h%mpi_comm_col,.false.)

            if(.not. success) then
               call elsi_stop(" Cholesky failed.",e_h,caller)
            endif

            ! compute U^-1 -> S
            success = elpa_invert_trm_real_double(e_h%n_basis,e_h%ovlp_real,&
                         e_h%n_lrow,e_h%blk_row,e_h%n_lcol,&
                         e_h%mpi_comm_row,e_h%mpi_comm_col,.false.)

            if(.not. success) then
               call elsi_stop(" Matrix inversion failed.",e_h,caller)
            endif

            call elsi_get_time(e_h,t1)

            write(info_str,"('  Finished Cholesky decomposition')")
            call elsi_say(info_str,e_h)
            write(info_str,"('  | Time :',F10.3,' s')") t1-t0
            call elsi_say(info_str,e_h)
         endif
      endif ! n_elsi_calls == 1

      call elsi_get_time(e_h,t0)

      if(e_h%ovlp_is_sing) then ! Use scaled eigenvectors
         ! evec_real used as tmp_real
         ! tmp_real = H_real * S_real
         call pdgemm('N','N',e_h%n_basis,e_h%n_nonsing,e_h%n_basis,1.0_r8,&
                 e_h%ham_real,1,1,e_h%sc_desc,e_h%ovlp_real,1,&
                 e_h%n_basis-e_h%n_nonsing+1,e_h%sc_desc,0.0_r8,e_h%evec_real,&
                 1,1,e_h%sc_desc)

         ! H_real = (S_real)^T * tmp_real
         call pdgemm('T','N',e_h%n_nonsing,e_h%n_nonsing,e_h%n_basis,1.0_r8,&
                 e_h%ovlp_real,1,e_h%n_basis-e_h%n_nonsing+1,e_h%sc_desc,&
                 e_h%evec_real,1,1,e_h%sc_desc,0.0_r8,e_h%ham_real,1,1,&
                 e_h%sc_desc)
      else ! Use Cholesky
         success = elpa_mult_at_b_real_double('U','L',e_h%n_basis,e_h%n_basis,&
                      e_h%ovlp_real,e_h%n_lrow,e_h%n_lcol,e_h%ham_real,&
                      e_h%n_lrow,e_h%n_lcol,e_h%blk_row,e_h%mpi_comm_row,&
                      e_h%mpi_comm_col,e_h%evec_real,e_h%n_lrow,e_h%n_lcol)

         if(.not. success) then
            call elsi_stop(" Matrix multiplication failed.",e_h,caller)
         endif

         call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,e_h%evec_real,1,1,&
                 e_h%sc_desc,0.0_r8,e_h%ham_real,1,1,e_h%sc_desc)

         e_h%evec_real = e_h%ham_real

         success = elpa_mult_at_b_real_double('U','U',e_h%n_basis,e_h%n_basis,&
                      e_h%ovlp_real,e_h%n_lrow,e_h%n_lcol,e_h%evec_real,&
                      e_h%n_lrow,e_h%n_lcol,e_h%blk_row,e_h%mpi_comm_row,&
                      e_h%mpi_comm_col,e_h%ham_real,e_h%n_lrow,e_h%n_lcol)

         if(.not. success) then
            call elsi_stop(" Matrix multiplication failed.",e_h,caller)
         endif

         call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,e_h%ham_real,1,1,&
                 e_h%sc_desc,0.0_r8,e_h%evec_real,1,1,e_h%sc_desc)

         ! Set the lower part from the upper
         do i_col = 1,e_h%n_basis-1
            if(e_h%loc_col(i_col) == 0) cycle

            do i_row = i_col+1,e_h%n_basis
               if(e_h%loc_row(i_row) > 0) then
                  e_h%ham_real(e_h%loc_row(i_row),e_h%loc_col(i_col)) = &
                     e_h%evec_real(e_h%loc_row(i_row),e_h%loc_col(i_col))
               endif
            enddo
         enddo
      endif

      call elsi_get_time(e_h,t1)

      write(info_str,"('  Finished transformation to standard eigenproblem')")
      call elsi_say(info_str,e_h)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_say(info_str,e_h)
   end select

end subroutine

!>
!! This routine checks the singularity of overlap matrix by computing all its
!! eigenvalues. On exit, S is not modified if not singular, or is overwritten by
!! scaled eigenvectors if singular, which can be esed to transform the
!! generalized eigenproblem to the standard form.
!!
subroutine elsi_check_singularity(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   real(kind=r8)    :: ev_sqrt
   integer(kind=i4) :: i
   logical          :: success
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   real(kind=r8),    allocatable :: copy_real(:,:)
   complex(kind=r8), allocatable :: copy_cmplx(:,:)

   character*40, parameter :: caller = "elsi_check_singularity"

   call elsi_get_time(e_h,t0)

   select case(e_h%data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(e_h,copy_cmplx,e_h%n_lrow,e_h%n_lcol,"copy_cmplx",&
              caller)

      ! Use copy_cmplx to store overlap matrix, otherwise it will
      ! be destroyed by eigenvalue calculation
      copy_cmplx = e_h%ovlp_cmplx

      ! Use customized ELPA 2-stage solver to check overlap singularity
      ! Eigenvectors computed only for singular overlap matrix
      success = elpa_check_singularity_complex_double(e_h%n_basis,e_h%n_basis,&
                   copy_cmplx,e_h%n_lrow,e_h%eval,e_h%evec_cmplx,&
                   e_h%n_lrow,e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,&
                   e_h%mpi_comm_col,e_h%mpi_comm,e_h%sing_tol,e_h%n_nonsing)

      if(.not. success) then
         call elsi_stop(" Singularity check failed.",e_h,caller)
      endif

      call elsi_deallocate(e_h,copy_cmplx,"copy_cmplx")

      e_h%n_states_solve = min(e_h%n_nonsing,e_h%n_states)

      if(e_h%n_nonsing < e_h%n_basis) then ! Singular
         e_h%ovlp_is_sing = .true.

         if(e_h%stop_sing) then
            call elsi_stop(" Overlap matrix is singular.",e_h,caller)
         endif

         call elsi_say("  Overlap matrix is singular.",e_h)

         write(info_str,"('  | Number of basis functions reduced to: ',I13)")&
            e_h%n_nonsing
         call elsi_say(info_str,e_h)

         call elsi_say("  Using scaled eigenvectors of overlap matrix for"//&
                 " transformation",e_h)

         ! Overlap matrix is overwritten with scaled eigenvectors
         do i = e_h%n_basis-e_h%n_nonsing+1,e_h%n_basis
            ev_sqrt = sqrt(e_h%eval(i))

            if(e_h%loc_col(i) == 0) cycle

            e_h%ovlp_cmplx(:,e_h%loc_col(i)) = &
               e_h%evec_cmplx(:,e_h%loc_col(i))/ev_sqrt
         enddo
      else ! Nonsingular
         e_h%ovlp_is_sing = .false.
         call elsi_say("  Overlap matrix is nonsingular",e_h)
      endif ! Singular overlap?
   case(REAL_VALUES)
      call elsi_allocate(e_h,copy_real,e_h%n_lrow,e_h%n_lcol,"copy_real",&
              caller)

      ! Use copy_real to store overlap matrix, otherwise it will be
      ! destroyed by eigenvalue calculation
      copy_real = e_h%ovlp_real

      ! Use customized ELPA 2-stage solver to check overlap singularity
      ! Eigenvectors computed only for singular overlap matrix
      success = elpa_check_singularity_real_double(e_h%n_basis,e_h%n_basis,&
                   copy_real,e_h%n_lrow,e_h%eval,e_h%evec_real,e_h%n_lrow,&
                   e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                   e_h%mpi_comm,e_h%sing_tol,e_h%n_nonsing)

      if(.not. success) then
         call elsi_stop(" Singularity check failed.",e_h,caller)
      endif

      call elsi_deallocate(e_h,copy_real,"copy_real")

      e_h%n_states_solve = min(e_h%n_nonsing,e_h%n_states)

      if(e_h%n_nonsing < e_h%n_basis) then ! Singular
         e_h%ovlp_is_sing = .true.

         if(e_h%stop_sing) then
            call elsi_stop(" Overlap matrix is singular.",e_h,caller)
         endif

         call elsi_say("  Overlap matrix is singular.",e_h)

         write(info_str,"('  | Number of basis functions reduced to: ',I13)")&
            e_h%n_nonsing
         call elsi_say(info_str,e_h)

         call elsi_say("  Using scaled eigenvectors of overlap matrix for"//&
                 " transformation",e_h)

         ! Overlap matrix is overwritten with scaled eigenvectors
         do i = e_h%n_basis-e_h%n_nonsing+1,e_h%n_basis
            ev_sqrt = sqrt(e_h%eval(i))

            if(e_h%loc_col(i) == 0) cycle

            e_h%ovlp_real(:,e_h%loc_col(i)) = e_h%evec_real(:,e_h%loc_col(i))/&
                                                 ev_sqrt
         enddo
      else ! Nonsingular
         e_h%ovlp_is_sing = .false.
         call elsi_say("  Overlap matrix is nonsingular",e_h)
      endif ! Singular overlap?
   end select

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished singularity check of overlap matrix')")
   call elsi_say(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,e_h)

end subroutine

!>
!! This routine back-transforms eigenvectors in the standard form to the
!! original generalized form.
!!
subroutine elsi_to_original_ev(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   logical       :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character*200 :: info_str

   real(kind=r8),    allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character*40, parameter :: caller = "elsi_to_original_ev"

   call elsi_get_time(e_h,t0)

   select case(e_h%data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(e_h,tmp_cmplx,e_h%n_lrow,e_h%n_lcol,"tmp_cmplx",&
              caller)
      tmp_cmplx = e_h%evec_cmplx

      if(e_h%ovlp_is_sing) then
         ! Transform matrix is stored in S_cmplx after elsi_to_standard_evp
         call pzgemm('N','N',e_h%n_basis,e_h%n_states_solve,e_h%n_nonsing,&
                 (1.0_r8,0.0_r8),e_h%ovlp_cmplx,1,e_h%n_basis-e_h%n_nonsing+1,&
                 e_h%sc_desc,tmp_cmplx,1,1,e_h%sc_desc,(0.0_r8,0.0_r8),&
                 e_h%evec_cmplx,1,1,e_h%sc_desc)
      else ! Nonsingular, use Cholesky
         ! (U^-1) is stored in S_cmplx after elsi_to_standard_evp
         ! C_cmplx = S_cmplx * C_cmplx = S_cmplx * tmp_cmplx
         call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),e_h%ovlp_cmplx,1,&
                 1,e_h%sc_desc,(0.0_r8,0.0_r8),e_h%ham_cmplx,1,1,e_h%sc_desc)

         success = elpa_mult_ah_b_complex_double('L','N',e_h%n_basis,&
                      e_h%n_states,e_h%ham_cmplx,e_h%n_lrow,e_h%n_lcol,&
                      tmp_cmplx,e_h%n_lrow,e_h%n_lcol,e_h%blk_row,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%evec_cmplx,&
                      e_h%n_lrow,e_h%n_lcol)

         if(.not. success) then
            call elsi_stop(" Matrix multiplication failed.",e_h,caller)
         endif
      endif

      call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")
   case(REAL_VALUES)
      call elsi_allocate(e_h,tmp_real,e_h%n_lrow,e_h%n_lcol,"tmp_real",&
              caller)
      tmp_real = e_h%evec_real

      if(e_h%ovlp_is_sing) then
         ! Transform matrix is stored in S_real after elsi_to_standard_evp
         call pdgemm('N','N',e_h%n_basis,e_h%n_states_solve,e_h%n_nonsing,&
                 1.0_r8,e_h%ovlp_real,1,e_h%n_basis-e_h%n_nonsing+1,&
                 e_h%sc_desc,tmp_real,1,1,e_h%sc_desc,0.0_r8,e_h%evec_real,1,1,&
                 e_h%sc_desc)
      else ! Nonsingular, use Cholesky
         ! (U^-1) is stored in S_real after elsi_to_standard_evp
         ! C_real = S_real * C_real = S_real * tmp_real
         call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,e_h%ovlp_real,1,1,&
                 e_h%sc_desc,0.0_r8,e_h%ham_real,1,1,e_h%sc_desc)

         success = elpa_mult_at_b_real_double('L','N',e_h%n_basis,e_h%n_states,&
                      e_h%ham_real,e_h%n_lrow,e_h%n_lcol,tmp_real,&
                      e_h%n_lrow,e_h%n_lcol,e_h%blk_row,e_h%mpi_comm_row,&
                      e_h%mpi_comm_col,e_h%evec_real,e_h%n_lrow,e_h%n_lcol)

         if(.not. success) then
            call elsi_stop(" Matrix multiplication failed.",e_h,caller)
         endif
      endif

      call elsi_deallocate(e_h,tmp_real,"tmp_real")
   end select

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished back-transformation of eigenvectors')")
   call elsi_say(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,e_h)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_evp_elpa(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   integer(kind=i4) :: mpierr
   logical          :: success
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_solve_evp_elpa"

   elpa_print_times = e_h%elpa_output

   ! Compute sparsity
   if(e_h%n_elsi_calls == 1) then
      select case(e_h%data_type)
      case(COMPLEX_VALUES)
         call elsi_get_local_nnz_complex(e_h,e_h%ham_cmplx,e_h%n_lrow,&
                 e_h%n_lcol,e_h%nnz_l)

         call MPI_Allreduce(e_h%nnz_l,e_h%nnz_g,1,mpi_integer4,mpi_sum,&
                 e_h%mpi_comm,mpierr)
      case(REAL_VALUES)
         call elsi_get_local_nnz_real(e_h,e_h%ham_real,e_h%n_lrow,&
                 e_h%n_lcol,e_h%nnz_l)

         call MPI_Allreduce(e_h%nnz_l,e_h%nnz_g,1,mpi_integer4,mpi_sum,&
                 e_h%mpi_comm,mpierr)
      end select
   endif

   ! Transform to standard form
   if(.not. e_h%ovlp_is_unit) then
      call elsi_to_standard_evp(e_h)
   endif

   call elsi_get_time(e_h,t0)

   call elsi_say("  Starting ELPA eigensolver",e_h)

   ! Solve evp, return eigenvalues and eigenvectors
   select case(e_h%data_type)
   case(COMPLEX_VALUES)
      if(e_h%elpa_solver == 2) then
         success = elpa_solve_evp_complex_2stage_double(e_h%n_nonsing,&
                      e_h%n_states_solve,e_h%ham_cmplx,e_h%n_lrow,e_h%eval,&
                      e_h%evec_cmplx,e_h%n_lrow,e_h%blk_row,e_h%n_lcol,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%mpi_comm)
      else
         success = elpa_solve_evp_complex_1stage_double(e_h%n_nonsing,&
                      e_h%n_states_solve,e_h%ham_cmplx,e_h%n_lrow,e_h%eval,&
                      e_h%evec_cmplx,e_h%n_lrow,e_h%blk_row,e_h%n_lcol,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%mpi_comm)
      endif
   case(REAL_VALUES)
      if(e_h%elpa_solver == 2) then
         success = elpa_solve_evp_real_2stage_double(e_h%n_nonsing,&
                      e_h%n_states_solve,e_h%ham_real,e_h%n_lrow,e_h%eval,&
                      e_h%evec_real,e_h%n_lrow,e_h%blk_row,e_h%n_lcol,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%mpi_comm)
      else
         success = elpa_solve_evp_real_1stage_double(e_h%n_nonsing,&
                      e_h%n_states_solve,e_h%ham_real,e_h%n_lrow,e_h%eval,&
                      e_h%evec_real,e_h%n_lrow,e_h%blk_row,e_h%n_lcol,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%mpi_comm)
      endif
   end select

   if(.not. success) then
      call elsi_stop(" ELPA solver failed.",e_h,caller)
   endif

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(e_h%n_nonsing < e_h%n_basis) then
      e_h%eval(e_h%n_nonsing+1:e_h%n_basis) = e_h%eval(e_h%n_nonsing)+10.0_r8
   endif

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished solving standard eigenproblem')")
   call elsi_say(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,e_h)

   ! Back-transform eigenvectors
   if(.not. e_h%ovlp_is_unit) then
      call elsi_to_original_ev(e_h)
   endif

   call MPI_Barrier(e_h%mpi_comm,mpierr)

end subroutine

!>
!! This routine sets default ELPA parameters.
!!
subroutine elsi_set_elpa_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*40, parameter :: caller = "elsi_set_elpa_default"

   ! ELPA solver
   if(e_h%n_basis < 250) then
      e_h%elpa_solver = 1
   else
      e_h%elpa_solver = 2
   endif

   ! How many single precision steps?
   e_h%n_single_steps = 0

   ! ELPA output?
   e_h%elpa_output = .false.

end subroutine

end module ELSI_ELPA
