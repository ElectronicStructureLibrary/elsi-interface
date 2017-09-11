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

   use ELSI_CONSTANTS, only: REAL_VALUES,COMPLEX_VALUES,GAUSSIAN,FERMI,&
                             METHFESSEL_PAXTON_0,METHFESSEL_PAXTON_1,&
                             INVERT_SQRT_PI

   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_MU, only: elsi_compute_mu_and_occ
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS
   use CHECK_SINGULARITY, only: elpa_check_singularity_real_double,&
                                elpa_check_singularity_complex_double
   use ELPA1
   use ELPA2

   implicit none

   private

   public :: elsi_get_elpa_comms
   public :: elsi_set_elpa_default
   public :: elsi_compute_occ_elpa
   public :: elsi_compute_dm_elpa
   public :: elsi_compute_edm_elpa
   public :: elsi_solve_evp_elpa
   public :: elsi_solve_evp_elpa_sp

contains

!>
!! This routine gets the row and column communicators for ELPA.
!!
subroutine elsi_get_elpa_comms(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: success

   character*40, parameter :: caller = "elsi_get_elpa_comms"

   success = elpa_get_communicators(e_h%mpi_comm,e_h%my_p_row,e_h%my_p_col,&
                e_h%mpi_comm_row,e_h%mpi_comm_col)

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

   select case(e_h%matrix_data_type)
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
            dble(e_h%dm_cmplx(e_h%loc_row(i_col),e_h%loc_col(i_col)))
      enddo
   end select

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished density matrix calculation')")
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

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

   select case(e_h%matrix_data_type)
   case(REAL_VALUES)
      call elsi_allocate(e_h,tmp_real,e_h%n_l_rows,e_h%n_l_cols,"tmp_real",caller)
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
      call elsi_allocate(e_h,tmp_cmplx,e_h%n_l_rows,e_h%n_l_cols,"tmp_cmplx",&
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
            dble(e_h%dm_cmplx(e_h%loc_row(i_col),e_h%loc_col(i_col)))
      enddo
   end select

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished energy density matrix calculation')")
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

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

   select case(e_h%matrix_data_type)
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
                         e_h%n_l_rows,e_h%n_b_rows,e_h%n_l_cols,&
                         e_h%mpi_comm_row,e_h%mpi_comm_col,.false.)

            ! compute U^-1 -> S
            success = elpa_invert_trm_complex_double(e_h%n_basis,&
                         e_h%ovlp_cmplx,e_h%n_l_rows,e_h%n_b_rows,e_h%n_l_cols,&
                         e_h%mpi_comm_row,e_h%mpi_comm_col,.false.)

            call elsi_get_time(e_h,t1)

            write(info_str,"('  Finished Cholesky decomposition')")
            call elsi_statement_print(info_str,e_h)
            write(info_str,"('  | Time :',F10.3,' s')") t1-t0
            call elsi_statement_print(info_str,e_h)
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
                      e_h%n_basis,e_h%ovlp_cmplx,e_h%n_l_rows,e_h%n_l_cols,&
                      e_h%ham_cmplx,e_h%n_l_rows,e_h%n_l_cols,e_h%n_b_rows,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%evec_cmplx,&
                      e_h%n_l_rows,e_h%n_l_cols)

         call pztranc(e_h%n_basis,e_h%n_basis,(1.0_r8,0.0_r8),e_h%evec_cmplx,1,&
                 1,e_h%sc_desc,(0.0_r8,0.0_r8),e_h%ham_cmplx,1,1,e_h%sc_desc)

         e_h%evec_cmplx = e_h%ham_cmplx

         success = elpa_mult_ah_b_complex_double('U','U',e_h%n_basis,&
                      e_h%n_basis,e_h%ovlp_cmplx,e_h%n_l_rows,e_h%n_l_cols,&
                      e_h%evec_cmplx,e_h%n_l_rows,e_h%n_l_cols,e_h%n_b_rows,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%ham_cmplx,&
                      e_h%n_l_rows,e_h%n_l_cols)

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
               dble(e_h%ham_cmplx(e_h%loc_row(i_col),e_h%loc_col(i_col)))
         enddo
      endif

      call elsi_get_time(e_h,t1)

      write(info_str,"('  Finished transformation to standard eigenproblem')")
      call elsi_statement_print(info_str,e_h)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_statement_print(info_str,e_h)
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
                         e_h%n_l_rows,e_h%n_b_rows,e_h%n_l_cols,&
                         e_h%mpi_comm_row,e_h%mpi_comm_col,.false.)

            ! compute U^-1 -> S
            success = elpa_invert_trm_real_double(e_h%n_basis,e_h%ovlp_real,&
                         e_h%n_l_rows,e_h%n_b_rows,e_h%n_l_cols,&
                         e_h%mpi_comm_row,e_h%mpi_comm_col,.false.)

            call elsi_get_time(e_h,t1)

            write(info_str,"('  Finished Cholesky decomposition')")
            call elsi_statement_print(info_str,e_h)
            write(info_str,"('  | Time :',F10.3,' s')") t1-t0
            call elsi_statement_print(info_str,e_h)
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
                      e_h%ovlp_real,e_h%n_l_rows,e_h%n_l_cols,e_h%ham_real,&
                      e_h%n_l_rows,e_h%n_l_cols,e_h%n_b_rows,e_h%mpi_comm_row,&
                      e_h%mpi_comm_col,e_h%evec_real,e_h%n_l_rows,e_h%n_l_cols)

         call pdtran(e_h%n_basis,e_h%n_basis,1.0_r8,e_h%evec_real,1,1,&
                 e_h%sc_desc,0.0_r8,e_h%ham_real,1,1,e_h%sc_desc)

         e_h%evec_real = e_h%ham_real

         success = elpa_mult_at_b_real_double('U','U',e_h%n_basis,e_h%n_basis,&
                      e_h%ovlp_real,e_h%n_l_rows,e_h%n_l_cols,e_h%evec_real,&
                      e_h%n_l_rows,e_h%n_l_cols,e_h%n_b_rows,e_h%mpi_comm_row,&
                      e_h%mpi_comm_col,e_h%ham_real,e_h%n_l_rows,e_h%n_l_cols)

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
      call elsi_statement_print(info_str,e_h)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_statement_print(info_str,e_h)
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

   select case(e_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(e_h,copy_cmplx,e_h%n_l_rows,e_h%n_l_cols,"copy_cmplx",&
              caller)

      ! Use copy_cmplx to store overlap matrix, otherwise it will
      ! be destroyed by eigenvalue calculation
      copy_cmplx = e_h%ovlp_cmplx

      ! Use customized ELPA 2-stage solver to check overlap singularity
      ! Eigenvectors computed only for singular overlap matrix
      success = elpa_check_singularity_complex_double(e_h%n_basis,e_h%n_basis,&
                   copy_cmplx,e_h%n_l_rows,e_h%eval,e_h%evec_cmplx,&
                   e_h%n_l_rows,e_h%n_b_rows,e_h%n_l_cols,e_h%mpi_comm_row,&
                   e_h%mpi_comm_col,e_h%mpi_comm,e_h%sing_tol,e_h%n_nonsing)

      call elsi_deallocate(e_h,copy_cmplx,"copy_cmplx")

      e_h%n_states_solve = min(e_h%n_nonsing,e_h%n_states)

      if(e_h%n_nonsing < e_h%n_basis) then ! Singular
         e_h%ovlp_is_sing = .true.

         if(e_h%stop_sing) then
            call elsi_stop(" Overlap matrix is singular. Running with a near"//&
                    "-singular basis set may lead to completely wrong"//&
                    " numerical results. Exiting...",e_h,caller)
         endif

         call elsi_statement_print("  Overlap matrix is singular. A large"//&
                 " basis set may be in use. Note that running with a naer"//&
                 " -singular basis set may lead to completely wrong"//&
                 " numerical results.",e_h)

         write(info_str,"(A,I13)") "  | Basis functions reduced to: ",&
            e_h%n_nonsing
         call elsi_statement_print(info_str,e_h)

         call elsi_statement_print("  Using scaled eigenvectors of overlap"//&
                 " matrix for transformation",e_h)

         ! Overlap matrix is overwritten with scaled eigenvectors
         do i = e_h%n_basis-e_h%n_nonsing+1,e_h%n_basis
            ev_sqrt = sqrt(e_h%eval(i))

            if(e_h%loc_col(i) == 0) cycle

            e_h%ovlp_cmplx(:,e_h%loc_col(i)) = &
               e_h%evec_cmplx(:,e_h%loc_col(i))/ev_sqrt
         enddo
      else ! Nonsingular
         e_h%ovlp_is_sing = .false.
         call elsi_statement_print("  Overlap matrix is nonsingular",e_h)
      endif ! Singular overlap?
   case(REAL_VALUES)
      call elsi_allocate(e_h,copy_real,e_h%n_l_rows,e_h%n_l_cols,"copy_real",&
              caller)

      ! Use copy_real to store overlap matrix, otherwise it will be
      ! destroyed by eigenvalue calculation
      copy_real = e_h%ovlp_real

      ! Use customized ELPA 2-stage solver to check overlap singularity
      ! Eigenvectors computed only for singular overlap matrix
      success = elpa_check_singularity_real_double(e_h%n_basis,e_h%n_basis,&
                   copy_real,e_h%n_l_rows,e_h%eval,e_h%evec_real,e_h%n_l_rows,&
                   e_h%n_b_rows,e_h%n_l_cols,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                   e_h%mpi_comm,e_h%sing_tol,e_h%n_nonsing)

      call elsi_deallocate(e_h,copy_real,"copy_real")

      e_h%n_states_solve = min(e_h%n_nonsing,e_h%n_states)

      if(e_h%n_nonsing < e_h%n_basis) then ! Singular
         e_h%ovlp_is_sing = .true.

         if(e_h%stop_sing) then
            call elsi_stop(" Overlap matrix is singular. Running with a near"//&
                    "-singular basis set may lead to completely wrong"//&
                    " numerical results. Exiting...",e_h,caller)
         endif

         call elsi_statement_print("  Overlap matrix is singular. Note that"//&
                 " running with a near-singular basis set may lead to"//&
                 " completely wrong numerical results.",e_h)

         write(info_str,"(A,I13)") "  | Basis functions reduced to: ",&
            e_h%n_nonsing
         call elsi_statement_print(info_str,e_h)

         call elsi_statement_print("  Using scaled eigenvectors of overlap"//&
                 " matrix for transformation",e_h)

         ! Overlap matrix is overwritten with scaled eigenvectors
         do i = e_h%n_basis-e_h%n_nonsing+1,e_h%n_basis
            ev_sqrt = sqrt(e_h%eval(i))

            if(e_h%loc_col(i) == 0) cycle

            e_h%ovlp_real(:,e_h%loc_col(i)) = e_h%evec_real(:,e_h%loc_col(i))/&
                                                 ev_sqrt
         enddo
      else ! Nonsingular
         e_h%ovlp_is_sing = .false.
         call elsi_statement_print("  Overlap matrix is nonsingular",e_h)
      endif ! Singular overlap?
   end select ! select matrix_data_type

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished singularity check of overlap matrix')")
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

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

   select case(e_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(e_h,tmp_cmplx,e_h%n_l_rows,e_h%n_l_cols,"tmp_cmplx",&
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
                      e_h%n_states,e_h%ham_cmplx,e_h%n_l_rows,e_h%n_l_cols,&
                      tmp_cmplx,e_h%n_l_rows,e_h%n_l_cols,e_h%n_b_rows,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%evec_cmplx,&
                      e_h%n_l_rows,e_h%n_l_cols)
      endif

      call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")
   case(REAL_VALUES)
      call elsi_allocate(e_h,tmp_real,e_h%n_l_rows,e_h%n_l_cols,"tmp_real",&
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
                      e_h%ham_real,e_h%n_l_rows,e_h%n_l_cols,tmp_real,&
                      e_h%n_l_rows,e_h%n_l_cols,e_h%n_b_rows,e_h%mpi_comm_row,&
                      e_h%mpi_comm_col,e_h%evec_real,e_h%n_l_rows,e_h%n_l_cols)
      endif

      call elsi_deallocate(e_h,tmp_real,"tmp_real")
   end select

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished back-transformation of eigenvectors')")
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

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

   ! Transform to standard form
   if(.not. e_h%ovlp_is_unit) then
      call elsi_to_standard_evp(e_h)
   endif

   call elsi_get_time(e_h,t0)

   call elsi_statement_print("  Starting ELPA eigensolver",e_h)

   ! Solve evp, return eigenvalues and eigenvectors
   select case(e_h%matrix_data_type)
   case(COMPLEX_VALUES)
      if(e_h%elpa_solver == 2) then
         success = elpa_solve_evp_complex_2stage_double(e_h%n_nonsing,&
                      e_h%n_states_solve,e_h%ham_cmplx,e_h%n_l_rows,e_h%eval,&
                      e_h%evec_cmplx,e_h%n_l_rows,e_h%n_b_rows,e_h%n_l_cols,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%mpi_comm)
      else
         success = elpa_solve_evp_complex_1stage_double(e_h%n_nonsing,&
                      e_h%n_states_solve,e_h%ham_cmplx,e_h%n_l_rows,e_h%eval,&
                      e_h%evec_cmplx,e_h%n_l_rows,e_h%n_b_rows,e_h%n_l_cols,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%mpi_comm)
      endif
   case(REAL_VALUES)
      if(e_h%elpa_solver == 2) then
         success = elpa_solve_evp_real_2stage_double(e_h%n_nonsing,&
                      e_h%n_states_solve,e_h%ham_real,e_h%n_l_rows,e_h%eval,&
                      e_h%evec_real,e_h%n_l_rows,e_h%n_b_rows,e_h%n_l_cols,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%mpi_comm)
      else
         success = elpa_solve_evp_real_1stage_double(e_h%n_nonsing,&
                      e_h%n_states_solve,e_h%ham_real,e_h%n_l_rows,e_h%eval,&
                      e_h%evec_real,e_h%n_l_rows,e_h%n_b_rows,e_h%n_l_cols,&
                      e_h%mpi_comm_row,e_h%mpi_comm_col,e_h%mpi_comm)
      endif
   end select

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
   if(e_h%n_nonsing < e_h%n_basis) then
      e_h%eval(e_h%n_nonsing+1:e_h%n_basis) = e_h%eval(e_h%n_nonsing)+10.0_r8
   endif

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished solving standard eigenproblem')")
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

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

   ! ELPA output?
   e_h%elpa_output = .false.

end subroutine

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

   select case(e_h%matrix_data_type)
   case(COMPLEX_VALUES)
      if(e_h%n_elsi_calls == 1) then
         if(.not. e_h%no_sing_check) then
            call elsi_check_singularity_sp(e_h)
         endif
      endif ! n_elsi_calls == 1

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
         call elsi_statement_print(info_str,e_h)
         write(info_str,"('  | Time :',F10.3,' s')") t1-t0
         call elsi_statement_print(info_str,e_h)
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
      call elsi_statement_print(info_str,e_h)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_statement_print(info_str,e_h)
   case(REAL_VALUES)
      if(e_h%n_elsi_calls == 1) then
         if(.not. e_h%no_sing_check) then
            call elsi_check_singularity_sp(e_h)
         endif
      endif ! n_elsi_calls == 1

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
         call elsi_statement_print(info_str,e_h)
         write(info_str,"('  | Time :',F10.3,' s')") t1-t0
         call elsi_statement_print(info_str,e_h)
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
      call elsi_statement_print(info_str,e_h)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_statement_print(info_str,e_h)
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

   select case(e_h%matrix_data_type)
   case(COMPLEX_VALUES)
      if(e_h%ovlp_is_sing) then
         call elsi_allocate(e_h,tmp_cmplx,e_h%n_l_rows,e_h%n_l_cols,&
                 "tmp_cmplx",caller)
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
         call elsi_allocate(e_h,tmp_real,e_h%n_l_rows,e_h%n_l_cols,"tmp_real",&
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
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_evp_elpa_sp(e_h)

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

   character*40, parameter :: caller = "elsi_solve_evp_elpa_sp"

   ! Transform to standard form
   if(.not. e_h%ovlp_is_unit) then
      call elsi_to_standard_evp_sp(e_h)
   endif

   ! Solve evp, return eigenvalues and eigenvectors
   call elsi_statement_print("  Starting ELPA eigensolver",e_h)
   call elsi_get_time(e_h,t0)

   select case(e_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(e_h,off_diag,e_h%n_nonsing,"off_diag",caller)
      call elsi_allocate(e_h,tau_cmplx,e_h%n_nonsing,"tau_cmplx",caller)
      call elsi_allocate(e_h,tmp_real,e_h%n_nonsing,e_h%n_nonsing,"tmp_real",&
              caller)
      call elsi_allocate(e_h,tmp_cmplx,e_h%n_nonsing,e_h%n_nonsing,"tmp_cmplx",&
              caller)

      call zhetrd('U',e_h%n_nonsing,e_h%ham_cmplx,e_h%n_nonsing,e_h%eval,&
              off_diag,tau_cmplx,tmp_cmplx,e_h%n_nonsing*e_h%n_nonsing,ierr)

      success = elpa_solve_tridi_double(e_h%n_nonsing,e_h%n_states_solve,&
                   e_h%eval,off_diag,tmp_real,e_h%n_nonsing,64,e_h%n_nonsing,&
                   mpi_comm_self,mpi_comm_self,.false.)

      e_h%evec_cmplx(1:e_h%n_nonsing,1:e_h%n_states_solve) = &
         tmp_real(1:e_h%n_nonsing,1:e_h%n_states_solve)

      call zunmtr('L','U','N',e_h%n_nonsing,e_h%n_states_solve,e_h%ham_cmplx,&
              e_h%n_nonsing,tau_cmplx,e_h%evec_cmplx,e_h%n_nonsing,tmp_cmplx,&
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

      call dsytrd('U',e_h%n_nonsing,e_h%ham_real,e_h%n_nonsing,e_h%eval,&
              off_diag,tau_real,tmp_real,e_h%n_nonsing*e_h%n_nonsing,ierr)

      success = elpa_solve_tridi_double(e_h%n_nonsing,e_h%n_states_solve,&
                   e_h%eval,off_diag,tmp_real,e_h%n_nonsing,64,e_h%n_nonsing,&
                   mpi_comm_self,mpi_comm_self,.false.)

      e_h%evec_real(1:e_h%n_nonsing,1:e_h%n_states_solve) = &
         tmp_real(1:e_h%n_nonsing,1:e_h%n_states_solve)

      call dormtr('L','U','N',e_h%n_nonsing,e_h%n_states_solve,e_h%ham_real,&
              e_h%n_nonsing,tau_real,e_h%evec_real,e_h%n_nonsing,tmp_real,&
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
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

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

   select case(e_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(e_h,copy_cmplx,e_h%n_l_rows,e_h%n_l_cols,"copy_cmplx",&
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

         if(e_h%stop_sing) then
            call elsi_stop(" Overlap matrix is singular. Running with a near"//&
                    "-singular basis set may lead to completely wrong"//&
                    " numerical results. Exiting...",e_h,caller)
         endif

         call elsi_statement_print("  Overlap matrix is singular. Note that"//&
                 " running with a near-singular basis set may lead to"//&
                 " completely wrong numerical results.",e_h)

         write(info_str,"(A,I13)") "  | Basis functions reduced to: ",&
            e_h%n_nonsing
         call elsi_statement_print(info_str,e_h)

         call elsi_statement_print("  Using scaled eigenvectors of overlap"//&
                 " matrix for transformation",e_h)

         ! Overlap matrix is overwritten with scaled eigenvectors
         do i = 1,e_h%n_nonsing
            ev_sqrt = sqrt(e_h%eval(i))
            e_h%ovlp_cmplx(:,i) = e_h%evec_cmplx(:,i)/ev_sqrt
         enddo
      else ! Nonsingular
         e_h%ovlp_is_sing = .false.
         call elsi_statement_print("  Overlap matrix is nonsingular",e_h)
      endif ! Singular overlap?
   case(REAL_VALUES)
      call elsi_allocate(e_h,copy_real,e_h%n_l_rows,e_h%n_l_cols,"copy_real",&
              caller)

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

         if(e_h%stop_sing) then
            call elsi_stop(" Overlap matrix is singular. Running with a near"//&
                    "-singular basis set may lead to completely wrong"//&
                    " numerical results. Exiting...",e_h,caller)
         endif

         call elsi_statement_print("  Overlap matrix is singular. Note that"//&
                 " running with a near-singular basis set may lead to"//&
                 " completely wrong numerical results.",e_h)

         write(info_str,"(A,I13)") "  | Basis functions reduced to: ",&
            e_h%n_nonsing
         call elsi_statement_print(info_str,e_h)

         call elsi_statement_print("  Using scaled eigenvectors of overlap"//&
                 " overlap matrix for transformation",e_h)

         ! Overlap matrix is overwritten with scaled eigenvectors
         do i = 1,e_h%n_nonsing
            ev_sqrt = sqrt(e_h%eval(i))
            e_h%ovlp_real(:,i) = e_h%evec_real(:,i)/ev_sqrt
         enddo
      else ! Nonsingular
         e_h%ovlp_is_sing = .false.
         call elsi_statement_print("  Overlap matrix is nonsingular",e_h)
      endif ! Singular overlap?
   end select ! select matrix_data_type

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished singularity check of overlap matrix')")
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

end subroutine

end module ELSI_ELPA
