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

   use ELSI_CONSTANTS, only: REAL_VALUES,COMPLEX_VALUES,UNSET,&
                             GAUSSIAN,FERMI,METHFESSEL_PAXTON_0,&
                             METHFESSEL_PAXTON_1,INVERT_SQRT_PI

   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4
   use ELSI_TIMERS
   use ELSI_UTILS
   use CHECK_SINGULARITY, only: elpa_check_singularity_real_double,&
                                elpa_check_singularity_complex_double
   use ELPA1
   use ELPA2

   implicit none
   private

   public :: elsi_get_elpa_comms
   public :: elsi_compute_occ_elpa
   public :: elsi_compute_dm_elpa
   public :: elsi_compute_edm_elpa
   public :: elsi_solve_evp_elpa
   public :: elsi_solve_evp_elpa_sp
   public :: elsi_set_elpa_default

contains

!>
!! This routine gets the row and column communicators for ELPA.
!!
subroutine elsi_get_elpa_comms(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: success

   character*40, parameter :: caller = "elsi_get_elpa_comms"

   success = elpa_get_communicators(elsi_h%mpi_comm,elsi_h%my_p_row,&
      elsi_h%my_p_col,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col)

end subroutine

!>
!! This routine constructs the density matrix.
!!
subroutine elsi_compute_dm_elpa(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: i
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_row
   integer(kind=i4) :: max_state

   real(kind=r8), allocatable :: factor(:)

   character*40, parameter :: caller = "elsi_compute_dm_elpa"

   call elsi_start_density_matrix_time(elsi_h)

   call elsi_allocate(elsi_h,factor,elsi_h%n_states,"factor",caller)

   max_state = 0

   do i = 1,elsi_h%n_states
      if(elsi_h%occ_num(i,elsi_h%i_spin,elsi_h%i_kpt) > 0.0_r8) then
         factor(i) = sqrt(elsi_h%occ_num(i,elsi_h%i_spin,elsi_h%i_kpt))
         max_state = i
      endif
   enddo

   select case(elsi_h%matrix_data_type)
   case(REAL_VALUES)
      elsi_h%ham_real = elsi_h%evec_real

      do i = 1,elsi_h%n_states
         if(factor(i) > 0.0_r8) then
            if(elsi_h%local_col(i) > 0) then
               elsi_h%ham_real(:,elsi_h%local_col(i)) = &
                  elsi_h%ham_real(:,elsi_h%local_col(i))*factor(i)
            endif
         elseif(elsi_h%local_col(i) .ne. 0) then
            elsi_h%ham_real(:,elsi_h%local_col(i)) = 0.0_r8
         endif
      enddo

      elsi_h%dm_real = 0.0_r8

      ! Compute density matrix
      call pdsyrk('U','N',elsi_h%n_basis,max_state,1.0_r8,elsi_h%ham_real,1,1,&
              elsi_h%sc_desc,0.0_r8,elsi_h%dm_real,1,1,elsi_h%sc_desc)

      call elsi_deallocate(elsi_h,factor,"factor")

      ! Set full matrix from upper triangle
      call pdtran(elsi_h%n_basis,elsi_h%n_basis,1.0_r8,elsi_h%dm_real,&
              1,1,elsi_h%sc_desc,0.0_r8,elsi_h%ham_real,1,1,elsi_h%sc_desc)

      do i_col = 1,elsi_h%n_basis-1
         if(elsi_h%local_col(i_col) == 0) cycle
         do i_row = i_col+1,elsi_h%n_basis
            if(elsi_h%local_row(i_row) > 0) then
               elsi_h%dm_real(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                  elsi_h%ham_real(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
            endif
         enddo
      enddo

   case(COMPLEX_VALUES)
      elsi_h%ham_complex = elsi_h%evec_complex

      do i = 1,elsi_h%n_states
         if(factor(i) > 0.0_r8) then
            if(elsi_h%local_col(i) > 0) then
               elsi_h%ham_complex(:,elsi_h%local_col(i)) = &
                  elsi_h%ham_complex(:,elsi_h%local_col(i))*factor(i)
            endif
         elseif(elsi_h%local_col(i) .ne. 0) then
            elsi_h%ham_complex(:,elsi_h%local_col(i)) = (0.0_r8,0.0_r8)
         endif
      enddo

      elsi_h%dm_complex = (0.0_r8,0.0_r8)

      ! Compute density matrix
      call pzherk('U','N',elsi_h%n_basis,max_state,(1.0_r8,0.0_r8),elsi_h%ham_complex,&
              1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),elsi_h%dm_complex,1,1,elsi_h%sc_desc)

      call elsi_deallocate(elsi_h,factor,"factor")

      ! Set full matrix from upper triangle
      call pztranc(elsi_h%n_basis,elsi_h%n_basis,(1.0_r8,0.0_r8),elsi_h%dm_complex,&
              1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),elsi_h%ham_complex,1,1,elsi_h%sc_desc)

      do i_col = 1,elsi_h%n_basis-1
         if(elsi_h%local_col(i_col) == 0) cycle
         do i_row = i_col+1,elsi_h%n_basis
            if(elsi_h%local_row(i_row) > 0) then
               elsi_h%dm_complex(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                  elsi_h%ham_complex(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
            endif
         enddo
      enddo

      ! Make diagonal real
      do i_col = 1,elsi_h%n_basis
         if((elsi_h%local_col(i_col) == 0) .or. (elsi_h%local_row(i_col) == 0)) cycle

         elsi_h%dm_complex(elsi_h%local_row(i_col),elsi_h%local_col(i_col)) = &
            dble(elsi_h%dm_complex(elsi_h%local_row(i_col),elsi_h%local_col(i_col)))
      enddo
   end select

   call elsi_stop_density_matrix_time(elsi_h)

end subroutine

!>
!! This routine constructs the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_elpa(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: i
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_row
   integer(kind=i4) :: max_state

   real(kind=r8),    allocatable :: factor(:)
   real(kind=r8),    allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_complex(:,:)

   character*40, parameter :: caller = "elsi_compute_edm_elpa"

   call elsi_allocate(elsi_h,factor,elsi_h%n_states,"factor",caller)

   max_state = 0

   do i = 1,elsi_h%n_states
      factor(i) = -1.0_r8*elsi_h%occ_num(i,elsi_h%i_spin,elsi_h%i_kpt)*&
                     elsi_h%eval(i)
      if(factor(i) > 0.0_r8) then
         factor(i) = sqrt(factor(i))
         max_state = i
      else
         factor(i) = 0.0_r8
      endif
   enddo

   select case(elsi_h%matrix_data_type)
   case(REAL_VALUES)
      call elsi_allocate(elsi_h,tmp_real,elsi_h%n_l_rows,elsi_h%n_l_cols,&
              "tmp_real",caller)
      tmp_real = elsi_h%evec_real

      do i = 1,elsi_h%n_states
         if(factor(i) > 0.0_r8) then
            if(elsi_h%local_col(i) > 0) then
               tmp_real(:,elsi_h%local_col(i)) = &
                  tmp_real(:,elsi_h%local_col(i))*factor(i)
            endif
         elseif(elsi_h%local_col(i) .ne. 0) then
            tmp_real(:,elsi_h%local_col(i)) = 0.0_r8
         endif
      enddo

      call elsi_deallocate(elsi_h,factor,"factor")

      elsi_h%dm_real = 0.0_r8

      ! Compute density matrix
      call pdsyrk('U','N',elsi_h%n_basis,max_state,1.0_r8,tmp_real,1,1,&
              elsi_h%sc_desc,0.0_r8,elsi_h%dm_real,1,1,elsi_h%sc_desc)

      elsi_h%dm_real = -1.0_r8*elsi_h%dm_real

      ! Set full matrix from upper triangle
      call pdtran(elsi_h%n_basis,elsi_h%n_basis,1.0_r8,elsi_h%dm_real,&
              1,1,elsi_h%sc_desc,0.0_r8,tmp_real,1,1,elsi_h%sc_desc)

      do i_col = 1,elsi_h%n_basis-1
         if(elsi_h%local_col(i_col) == 0) cycle
         do i_row = i_col+1,elsi_h%n_basis
            if(elsi_h%local_row(i_row) > 0) then
               elsi_h%dm_real(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                  tmp_real(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
            endif
         enddo
      enddo

      call elsi_deallocate(elsi_h,tmp_real,"tmp_real")

   case(COMPLEX_VALUES)
      call elsi_allocate(elsi_h,tmp_real,elsi_h%n_l_rows,elsi_h%n_l_cols,&
              "tmp_complex",caller)
      tmp_complex = elsi_h%evec_complex

      do i = 1,elsi_h%n_states
         if(factor(i) > 0.0_r8) then
            if(elsi_h%local_col(i) > 0) then
               tmp_complex(:,elsi_h%local_col(i)) = &
                  tmp_complex(:,elsi_h%local_col(i))*factor(i)
            endif
         elseif(elsi_h%local_col(i) .ne. 0) then
            tmp_complex(:,elsi_h%local_col(i)) = (0.0_r8,0.0_r8)
         endif
      enddo

      call elsi_deallocate(elsi_h,factor,"factor")

      elsi_h%dm_complex = (0.0_r8,0.0_r8)

      ! Compute density matrix
      call pzherk('U','N',elsi_h%n_basis,max_state,(1.0_r8,0.0_r8),tmp_complex,&
              1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),elsi_h%dm_complex,1,1,elsi_h%sc_desc)

      elsi_h%dm_complex = (-1.0_r8,0.0_r8)*elsi_h%dm_complex

      ! Set full matrix from upper triangle
      call pztranc(elsi_h%n_basis,elsi_h%n_basis,(1.0_r8,0.0_r8),elsi_h%dm_complex,&
              1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),tmp_complex,1,1,elsi_h%sc_desc)

      do i_col = 1,elsi_h%n_basis-1
         if(elsi_h%local_col(i_col) == 0) cycle
         do i_row = i_col+1,elsi_h%n_basis
            if(elsi_h%local_row(i_row) > 0) then
               elsi_h%dm_complex(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                  tmp_complex(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
            endif
         enddo
      enddo

      call elsi_deallocate(elsi_h,tmp_complex,"tmp_complex")

      ! Make diagonal real
      do i_col = 1,elsi_h%n_basis
         if((elsi_h%local_col(i_col) == 0) .or. (elsi_h%local_row(i_col) == 0)) cycle

         elsi_h%dm_complex(elsi_h%local_row(i_col),elsi_h%local_col(i_col)) = &
            dble(elsi_h%dm_complex(elsi_h%local_row(i_col),elsi_h%local_col(i_col)))
      enddo
   end select

end subroutine

!> 
!! This routine transforms a generalized eigenproblem (Hv = eSv) to
!! the standard form (H'v' = ev').
!!
!! First perform a Cholesky decomposition: S = (U^T)U, so that Hv = e(U^T)Uv.
!!
!! Then the new standard eigenproblem is
!! H(U^-1)(Uv) = e(U^T)(Uv) => ((U^-1)^T)H(U^-1)(Uv) = e(Uv).
!!
!! On exit, (U^-1) is stored in S to be reused later.
!!
subroutine elsi_to_standard_evp(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   logical          :: success

   character*40, parameter :: caller = "elsi_to_standard_evp"

   select case(elsi_h%matrix_data_type)
   case(COMPLEX_VALUES)
      if(elsi_h%n_elsi_calls == 1) then
         if(.not. elsi_h%no_sing_check) then
            call elsi_check_singularity(elsi_h)
         endif

         if(elsi_h%n_nonsing == elsi_h%n_basis) then ! Not singular
            call elsi_start_cholesky_time(elsi_h)

            elsi_h%ovlp_is_sing = .false.

            ! Compute S = (U^T)U, U -> S
            success = elpa_cholesky_complex_double(elsi_h%n_basis,&
                         elsi_h%ovlp_complex,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                         elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                         elsi_h%mpi_comm_col,.false.)

            ! compute U^-1 -> S
            success = elpa_invert_trm_complex_double(elsi_h%n_basis,&
                         elsi_h%ovlp_complex,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                         elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                         elsi_h%mpi_comm_col,.false.)

            call elsi_stop_cholesky_time(elsi_h)
         endif
      endif ! n_elsi_calls == 1

      call elsi_start_transform_evp_time(elsi_h)

      if(elsi_h%ovlp_is_sing) then ! Use scaled eigenvectors
         ! evec_complex used as tmp_complex
         ! tmp_complex = H_complex * S_complex
         call pzgemm('N','N',elsi_h%n_basis,elsi_h%n_nonsing,elsi_h%n_basis,&
                 (1.0_r8,0.0_r8),elsi_h%ham_complex,1,1,elsi_h%sc_desc,&
                 elsi_h%ovlp_complex,1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),&
                 elsi_h%evec_complex,1,1,elsi_h%sc_desc)

         ! H_complex = (S_complex)^* * tmp_complex
         call pzgemm('C','N',elsi_h%n_nonsing,elsi_h%n_nonsing,elsi_h%n_basis,&
                 (1.0_r8,0.0_r8),elsi_h%ovlp_complex,1,1,elsi_h%sc_desc,&
                 elsi_h%evec_complex,1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),&
                 elsi_h%ham_complex,1,1,elsi_h%sc_desc)

      else ! Use cholesky
         success = elpa_mult_ah_b_complex_double('U','L',elsi_h%n_basis,&
                      elsi_h%n_basis,elsi_h%ovlp_complex,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols,elsi_h%ham_complex,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols,elsi_h%n_b_rows,elsi_h%mpi_comm_row,&
                      elsi_h%mpi_comm_col,elsi_h%evec_complex,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols)

         call pztranc(elsi_h%n_basis,elsi_h%n_basis,(1.0_r8,0.0_r8),&
                 elsi_h%evec_complex,1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),&
                 elsi_h%ham_complex,1,1,elsi_h%sc_desc)

         elsi_h%evec_complex = elsi_h%ham_complex

         success = elpa_mult_ah_b_complex_double('U','U',elsi_h%n_basis,&
                      elsi_h%n_basis,elsi_h%ovlp_complex,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols,elsi_h%evec_complex,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols,elsi_h%n_b_rows,elsi_h%mpi_comm_row,&
                      elsi_h%mpi_comm_col,elsi_h%ham_complex,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols)

         call pztranc(elsi_h%n_basis,elsi_h%n_basis,(1.0_r8,0.0_r8),&
                 elsi_h%ham_complex,1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),&
                 elsi_h%evec_complex,1,1,elsi_h%sc_desc)

         ! Set the lower part from the upper
         do i_col = 1,elsi_h%n_basis-1
            if(elsi_h%local_col(i_col) == 0) cycle
            do i_row = i_col+1,elsi_h%n_basis
               if(elsi_h%local_row(i_row) > 0) then
                  elsi_h%ham_complex(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                     elsi_h%evec_complex(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
               endif
            enddo
         enddo

         do i_col=1,elsi_h%n_basis
            if(elsi_h%local_col(i_col) == 0 .or. elsi_h%local_row(i_col) == 0) cycle
            elsi_h%ham_complex(elsi_h%local_row(i_col),elsi_h%local_col(i_col)) = &
               dble(elsi_h%ham_complex(elsi_h%local_row(i_col),elsi_h%local_col(i_col)))
         enddo
      endif

      call elsi_stop_transform_evp_time(elsi_h)

   case(REAL_VALUES)
      if(elsi_h%n_elsi_calls == 1) then
         if(.not. elsi_h%no_sing_check) then
            call elsi_check_singularity(elsi_h)
         endif

         if(elsi_h%n_nonsing == elsi_h%n_basis) then ! Not singular
            call elsi_start_cholesky_time(elsi_h)

            elsi_h%ovlp_is_sing = .false.

            ! Compute S = (U^T)U, U -> S
            success = elpa_cholesky_real_double(elsi_h%n_basis,elsi_h%ovlp_real,&
                         elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,&
                         elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,.false.)

            ! compute U^-1 -> S
            success = elpa_invert_trm_real_double(elsi_h%n_basis,elsi_h%ovlp_real,&
                         elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,&
                         elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,.false.)

            call elsi_stop_cholesky_time(elsi_h)
         endif
      endif ! n_elsi_calls == 1

      call elsi_start_transform_evp_time(elsi_h)

      if(elsi_h%ovlp_is_sing) then ! Use scaled eigenvectors
         ! evec_real used as tmp_real
         ! tmp_real = H_real * S_real
         call pdgemm('N','N',elsi_h%n_basis,elsi_h%n_nonsing,elsi_h%n_basis,&
                 1.0_r8,elsi_h%ham_real,1,1,elsi_h%sc_desc,elsi_h%ovlp_real,&
                 1,1,elsi_h%sc_desc,0.0_r8,elsi_h%evec_real,1,1,elsi_h%sc_desc)

         ! H_real = (S_real)^T * tmp_real
         call pdgemm('T','N',elsi_h%n_nonsing,elsi_h%n_nonsing,elsi_h%n_basis,&
                 1.0_r8,elsi_h%ovlp_real,1,1,elsi_h%sc_desc,elsi_h%evec_real,&
                 1,1,elsi_h%sc_desc,0.0_r8,elsi_h%ham_real,1,1,elsi_h%sc_desc)

      else ! Use Cholesky
         success = elpa_mult_at_b_real_double('U','L',elsi_h%n_basis,&
                      elsi_h%n_basis,elsi_h%ovlp_real,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols,elsi_h%ham_real,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols,elsi_h%n_b_rows,elsi_h%mpi_comm_row,&
                      elsi_h%mpi_comm_col,elsi_h%evec_real,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols)

         call pdtran(elsi_h%n_basis,elsi_h%n_basis,1.0_r8,elsi_h%evec_real,1,1,&
                 elsi_h%sc_desc,0.0_r8,elsi_h%ham_real,1,1,elsi_h%sc_desc)

         elsi_h%evec_real = elsi_h%ham_real

         success = elpa_mult_at_b_real_double('U','U',elsi_h%n_basis,&
                      elsi_h%n_basis,elsi_h%ovlp_real,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols,elsi_h%evec_real,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols,elsi_h%n_b_rows,elsi_h%mpi_comm_row,&
                      elsi_h%mpi_comm_col,elsi_h%ham_real,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols)

         call pdtran(elsi_h%n_basis,elsi_h%n_basis,1.0_r8,elsi_h%ham_real,&
                 1,1,elsi_h%sc_desc,0.0_r8,elsi_h%evec_real,1,1,elsi_h%sc_desc)

         ! Set the lower part from the upper
         do i_col = 1,elsi_h%n_basis-1
            if(elsi_h%local_col(i_col) == 0) cycle
            do i_row = i_col+1,elsi_h%n_basis
               if(elsi_h%local_row(i_row) > 0) then
                  elsi_h%ham_real(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                     elsi_h%evec_real(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
               endif
            enddo
         enddo
      endif

      call elsi_stop_transform_evp_time(elsi_h)

   end select

end subroutine

!> 
!! This routine checks the singularity of overlap matrix by computing all
!! its eigenvalues.
!!
!! On exit, S is not modified if not singular, while overwritten by scaled
!! eigenvectors if singular, which is used to transform the generalized
!! eigenproblem to the standard form without using Cholesky.
!!
subroutine elsi_check_singularity(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   real(kind=r8)    :: ev_sqrt
   integer(kind=i4) :: i
   logical          :: success
   character*200    :: info_str

   real(kind=r8),    allocatable :: copy_real(:,:)
   complex(kind=r8), allocatable :: copy_complex(:,:)

   character*40, parameter :: caller = "elsi_check_singularity"

   call elsi_start_singularity_check_time(elsi_h)

   select case(elsi_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(elsi_h,copy_complex,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,"copy_complex",caller)

      ! Use copy_complex to store overlap matrix, otherwise it will
      ! be destroyed by eigenvalue calculation
      copy_complex = -elsi_h%ovlp_complex

      ! Use customized ELPA 2-stage solver to check overlap singularity
      ! Eigenvectors computed only for singular overlap matrix
      success = elpa_check_singularity_complex_double(elsi_h%n_basis,&
                   elsi_h%n_basis,copy_complex,elsi_h%n_l_rows,elsi_h%eval,&
                   elsi_h%evec_complex,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                   elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                   elsi_h%mpi_comm,elsi_h%sing_tol,elsi_h%n_nonsing)

      call elsi_deallocate(elsi_h,copy_complex,"copy_complex")

      ! Stop if n_states is larger than n_nonsing
      if(elsi_h%n_nonsing < elsi_h%n_states) then ! Too singular to continue
         call elsi_stop(" Overlap matrix is singular. The number of"//&
                 " basis functions after removing singularity is smaller"//&
                 " than the number of states. Exiting...",elsi_h,caller)
      elseif(elsi_h%n_nonsing < elsi_h%n_basis) then ! Singular
         elsi_h%ovlp_is_sing = .true.

         if(elsi_h%stop_sing) then
            call elsi_stop(" Overlap matrix is singular. Running with a"//&
                    " a near-singular basis set may lead to completely"//&
                    " wrong numerical results. Exiting...",elsi_h,caller)
         endif

         call elsi_statement_print("  Overlap matrix is singular. A large"//&
                 " basis sest may be in use. Note that running with a"//&
                 " near-singular basis set may lead to completely wrong"//&
                 " numerical results.",elsi_h)

         write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
            elsi_h%n_nonsing
         call elsi_statement_print(info_str,elsi_h)

         call elsi_statement_print("  Using scaled eigenvectors of overlap"//&
                 " matrix for transformation",elsi_h)

         ! Overlap matrix is overwritten with scaled eigenvectors
         do i = 1,elsi_h%n_nonsing
            ev_sqrt = sqrt(elsi_h%eval(i))
            if(elsi_h%local_col(i) == 0) cycle
            elsi_h%ovlp_complex(:,elsi_h%local_col(i)) = &
               elsi_h%evec_complex(:,elsi_h%local_col(i))/ev_sqrt
         enddo

      else ! Nonsingular
         elsi_h%ovlp_is_sing = .false.
         call elsi_statement_print("  Overlap matrix is nonsingular",elsi_h)
      endif ! Singular overlap?

   case(REAL_VALUES)
      call elsi_allocate(elsi_h,copy_real,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,"copy_real",caller)

      ! Use copy_real to store overlap matrix, otherwise it will be
      ! destroyed by eigenvalue calculation
      copy_real = -elsi_h%ovlp_real

      ! Use customized ELPA 2-stage solver to check overlap singularity
      ! Eigenvectors computed only for singular overlap matrix
      success = elpa_check_singularity_real_double(elsi_h%n_basis,&
                   elsi_h%n_basis,copy_real,elsi_h%n_l_rows,elsi_h%eval,&
                   elsi_h%evec_real,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                   elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                   elsi_h%mpi_comm,elsi_h%sing_tol,elsi_h%n_nonsing)

      call elsi_deallocate(elsi_h,copy_real,"copy_real")

      ! Stop if n_states is larger than n_nonsing
      if(elsi_h%n_nonsing < elsi_h%n_states) then ! Too singular to continue
         call elsi_stop(" Overlap matrix is singular. The number of basis"//&
                 " functions after removing singularity is smaller than the"//&
                 " number of states. Exiting...",elsi_h,caller)
      elseif(elsi_h%n_nonsing < elsi_h%n_basis) then ! Singular
         elsi_h%ovlp_is_sing = .true.

         if(elsi_h%stop_sing) then
            call elsi_stop(" Overlap matrix is singular. Running with a"//&
                    " near-singular basis set may lead to completely"//&
                    " wrong numerical results. Exiting...",elsi_h,caller)
         endif

         call elsi_statement_print("  Overlap matrix is singular. Note that"//&
                 " running with a near-singular basis set may lead to"//&
                 " completely wrong numerical results.",elsi_h)

         write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
            elsi_h%n_nonsing
         call elsi_statement_print(info_str,elsi_h)

         call elsi_statement_print("  Using scaled eigenvectors of"//&
                 " overlap matrix for transformation",elsi_h)

         ! Overlap matrix is overwritten with scaled eigenvectors
         do i = 1,elsi_h%n_nonsing
            ev_sqrt = sqrt(elsi_h%eval(i))
            if(elsi_h%local_col(i) == 0) cycle
            elsi_h%ovlp_real(:,elsi_h%local_col(i)) = &
               elsi_h%evec_real(:,elsi_h%local_col(i))/ev_sqrt
         enddo

      else ! Nonsingular
         elsi_h%ovlp_is_sing = .false.
         call elsi_statement_print("  Overlap matrix is nonsingular",elsi_h)
      endif ! Singular overlap?

   end select ! select matrix_data_type

   call elsi_stop_singularity_check_time(elsi_h)

end subroutine

!> 
!! This routine back-transforms eigenvectors in the standard form (H'v' = ev')
!! to the original generalized form (Hv = eSv), by computing v = (U^-1)v'.
!!
subroutine elsi_to_original_ev(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   logical :: success

   real(kind=r8),    allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_complex(:,:)

   character*40, parameter :: caller = "elsi_to_original_ev"

   call elsi_start_back_transform_ev_time(elsi_h)

   select case(elsi_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(elsi_h,tmp_complex,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,"tmp_complex",caller)
      tmp_complex = elsi_h%evec_complex

      if(elsi_h%ovlp_is_sing) then
         ! Transform matrix is stored in S_complex after elsi_to_standard_evp
         call pzgemm('N','N',elsi_h%n_basis,elsi_h%n_states,elsi_h%n_nonsing,&
                 (1.0_r8,0.0_r8),elsi_h%ovlp_complex,1,1,elsi_h%sc_desc,&
                 tmp_complex,1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),&
                 elsi_h%evec_complex,1,1,elsi_h%sc_desc)
      else ! Nonsingular, use Cholesky
         ! (U^-1) is stored in S_complex after elsi_to_standard_evp
         ! C_complex = S_complex * C_complex = S_complex * tmp_complex
         call pztranc(elsi_h%n_basis,elsi_h%n_basis,(1.0_r8,0.0_r8),&
                 elsi_h%ovlp_complex,1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),&
                 elsi_h%ham_complex,1,1,elsi_h%sc_desc)

         success = elpa_mult_ah_b_complex_double('L','N',elsi_h%n_basis,&
                      elsi_h%n_states,elsi_h%ham_complex,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols,tmp_complex,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols,elsi_h%n_b_rows,elsi_h%mpi_comm_row,&
                      elsi_h%mpi_comm_col,elsi_h%evec_complex,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols)
      endif

   case(REAL_VALUES)
      call elsi_allocate(elsi_h,tmp_real,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,"tmp_real",caller)
      tmp_real = elsi_h%evec_real

      if(elsi_h%ovlp_is_sing) then
         ! Transform matrix is stored in S_real after elsi_to_standard_evp
         call pdgemm('N','N',elsi_h%n_basis,elsi_h%n_states,elsi_h%n_nonsing,&
                 1.0_r8,elsi_h%ovlp_real,1,1,elsi_h%sc_desc,tmp_real,1,1,&
                 elsi_h%sc_desc,0.0_r8,elsi_h%evec_real,1,1,elsi_h%sc_desc)
      else ! Nonsingular, use Cholesky
         ! (U^-1) is stored in S_real after elsi_to_standard_evp
         ! C_real = S_real * C_real = S_real * tmp_real
         call pdtran(elsi_h%n_basis,elsi_h%n_basis,1.0_r8,elsi_h%ovlp_real,&
                 1,1,elsi_h%sc_desc,0.0_r8,elsi_h%ham_real,1,1,elsi_h%sc_desc)

         success = elpa_mult_at_b_real_double('L','N',elsi_h%n_basis,&
                      elsi_h%n_states,elsi_h%ham_real,elsi_h%n_l_rows,&
                      elsi_h%n_l_cols,tmp_real,elsi_h%n_l_rows,elsi_h%n_l_cols,&
                      elsi_h%n_b_rows,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                      elsi_h%evec_real,elsi_h%n_l_rows,elsi_h%n_l_cols)
      endif

   end select

   if(allocated(tmp_real))    call elsi_deallocate(elsi_h,tmp_real,"tmp_real")
   if(allocated(tmp_complex)) call elsi_deallocate(elsi_h,tmp_complex,"tmp_complex")

   call elsi_stop_back_transform_ev_time(elsi_h)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_evp_elpa(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: mpierr
   logical          :: success

   character*40, parameter :: caller = "elsi_solve_evp_elpa"

   elpa_print_times = elsi_h%elpa_output

   ! Transform to standard form
   if(.not. elsi_h%ovlp_is_unit) then
      call elsi_to_standard_evp(elsi_h)
   endif

   call elsi_start_standard_evp_time(elsi_h)

   ! Solve evp, return eigenvalues and eigenvectors
   if(elsi_h%elpa_solver == 2) then ! 2-stage solver
      call elsi_statement_print("  Starting ELPA 2-stage eigensolver",elsi_h)
      select case(elsi_h%matrix_data_type)
      case(COMPLEX_VALUES)
         success = elpa_solve_evp_complex_2stage_double(elsi_h%n_nonsing,&
                      elsi_h%n_states,elsi_h%ham_complex,elsi_h%n_l_rows,&
                      elsi_h%eval,elsi_h%evec_complex,elsi_h%n_l_rows,&
                      elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                      elsi_h%mpi_comm_col,elsi_h%mpi_comm)
      case(REAL_VALUES)
         success = elpa_solve_evp_real_2stage_double(elsi_h%n_nonsing,&
                      elsi_h%n_states,elsi_h%ham_real,elsi_h%n_l_rows,&
                      elsi_h%eval,elsi_h%evec_real,elsi_h%n_l_rows,&
                      elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                      elsi_h%mpi_comm_col,elsi_h%mpi_comm)
      end select
   else ! 1-stage solver
      call elsi_statement_print("  Starting ELPA 1-stage eigensolver",elsi_h)
      select case(elsi_h%matrix_data_type)
      case(COMPLEX_VALUES)
         success = elpa_solve_evp_complex_1stage_double(elsi_h%n_nonsing,&
                      elsi_h%n_states,elsi_h%ham_complex,elsi_h%n_l_rows,&
                      elsi_h%eval,elsi_h%evec_complex,elsi_h%n_l_rows,&
                      elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                      elsi_h%mpi_comm_col,elsi_h%mpi_comm)
      case(REAL_VALUES)
         success = elpa_solve_evp_real_1stage_double(elsi_h%n_nonsing,&
                      elsi_h%n_states,elsi_h%ham_real,elsi_h%n_l_rows,&
                      elsi_h%eval,elsi_h%evec_real,elsi_h%n_l_rows,&
                      elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                      elsi_h%mpi_comm_col,elsi_h%mpi_comm)
      end select
   endif

   call elsi_stop_standard_evp_time(elsi_h)

   ! Back-transform eigenvectors
   if(.not. elsi_h%ovlp_is_unit) then
      call elsi_to_original_ev(elsi_h)
   endif

   call MPI_Barrier(elsi_h%mpi_comm,mpierr)

end subroutine

!> 
!! This routine sets default ELPA parameters.
!! 
subroutine elsi_set_elpa_default(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_set_elpa_default"

   ! ELPA solver
   if(elsi_h%n_basis < 250) then
      elsi_h%elpa_solver = 1
   else
      elsi_h%elpa_solver = 2
   endif

   ! ELPA output?
   elsi_h%elpa_output = .false.

end subroutine

!> 
!! This routine transforms a generalized eigenproblem (Hv = eSv) to
!! the standard form (H'v' = ev').
!!
!! First perform a Cholesky decomposition: S = (U^T)U, so that Hv = e(U^T)Uv.
!!
!! Then the new standard eigenproblem is
!! H(U^-1)(Uv) = e(U^T)(Uv) => ((U^-1)^T)H(U^-1)(Uv) = e(Uv).
!!
!! On exit, (U^-1) is stored in S to be reused later.
!!
subroutine elsi_to_standard_evp_sp(elsi_h)

   implicit none
   include 'mpif.h'

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: nwork
   integer(kind=i4) :: n
   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: ierr
   logical          :: success

   integer(kind=i4), parameter :: nblk = 128
   character*40,     parameter :: caller = "elsi_to_standard_evp_sp"

   select case(elsi_h%matrix_data_type)
   case(COMPLEX_VALUES)
      if(elsi_h%n_elsi_calls == 1) then
         if(.not. elsi_h%no_sing_check) then
            call elsi_check_singularity_sp(elsi_h)
         endif
      endif ! n_elsi_calls == 1

      if(elsi_h%n_nonsing == elsi_h%n_basis) then ! Not singular
         call elsi_start_cholesky_time(elsi_h)

         elsi_h%ovlp_is_sing = .false.

         ! Erase the lower triangle
         do i = 1,elsi_h%n_basis
            do j= 1,i-1
               elsi_h%ovlp_complex(i,j) = (0.0_r8,0.0_r8)
            enddo
         enddo

         ! Compute S = (U^H)U, U -> S
         call zpotrf('U',elsi_h%n_basis,elsi_h%ovlp_complex,elsi_h%n_basis,ierr)

         ! compute U^-1 -> S
         call ztrtri('U','N',elsi_h%n_basis,elsi_h%ovlp_complex,elsi_h%n_basis,ierr)

         call elsi_stop_cholesky_time(elsi_h)
      endif

      call elsi_start_transform_evp_time(elsi_h)

      if(elsi_h%ovlp_is_sing) then ! Use scaled eigenvectors
         ! evec_complex used as tmp_complex
         ! tmp_complex = H_complex * S_complex
         call zgemm('N','N',elsi_h%n_basis,elsi_h%n_nonsing,elsi_h%n_basis,&
                 (1.0_r8,0.0_r8),elsi_h%ham_complex(1,1),elsi_h%n_basis,&
                 elsi_h%ovlp_complex(1,1),elsi_h%n_basis,(0.0_r8,0.0_r8),&
                 elsi_h%evec_complex(1,1),elsi_h%n_basis)

         ! H_complex = (S_complex)^* * tmp_complex
         call zgemm('C','N',elsi_h%n_nonsing,elsi_h%n_nonsing,elsi_h%n_basis,&
                 (1.0_r8,0.0_r8),elsi_h%ovlp_complex(1,1),elsi_h%n_basis,&
                 elsi_h%evec_complex(1,1),elsi_h%n_basis,(0.0_r8,0.0_r8),&
                 elsi_h%ham_complex(1,1),elsi_h%n_basis)

      else ! Use cholesky
         ! tmp_complex = H_complex * S_complex
         do n = 1,elsi_h%n_basis,nblk
            nwork = nblk

            if(n+nwork-1 > elsi_h%n_basis) then
               nwork = elsi_h%n_basis-n+1
            endif

            call zgemm('N','N',n+nwork-1,nwork,n+nwork-1,(1.0_r8,0.0_r8),&
                    elsi_h%ham_complex(1,1),elsi_h%n_basis,elsi_h%ovlp_complex(1,n),&
                    elsi_h%n_basis,(0.0_r8,0.0_r8),elsi_h%evec_complex(1,n),&
                    elsi_h%n_basis)
         enddo

         ! H_complex = (tmp_complex)^* * S_complex
         do n = 1,elsi_h%n_basis,nblk
            nwork = nblk

            if(n+nwork-1 > elsi_h%n_basis) then
               nwork = elsi_h%n_basis-n+1
            endif

            call zgemm('C','N',nwork,elsi_h%n_basis-n+1,n+nwork-1,(1.0_r8,0.0_r8),&
                    elsi_h%ovlp_complex(1,n),elsi_h%n_basis,elsi_h%evec_complex(1,n),&
                    elsi_h%n_basis,(0.0_r8,0.0_r8),elsi_h%ham_complex(n,n),&
                    elsi_h%n_basis)
         enddo
      endif

      call elsi_stop_transform_evp_time(elsi_h)

   case(REAL_VALUES)
      if(elsi_h%n_elsi_calls == 1) then
         if(.not. elsi_h%no_sing_check) then
            call elsi_check_singularity_sp(elsi_h)
         endif
      endif ! n_elsi_calls == 1

      if(elsi_h%n_nonsing == elsi_h%n_basis) then ! Not singular
         call elsi_start_cholesky_time(elsi_h)

         elsi_h%ovlp_is_sing = .false.

         ! Erase the lower triangle
         do i = 1,elsi_h%n_basis
            do j= 1,i-1
               elsi_h%ovlp_real(i,j) = 0.0_r8
            enddo
         enddo

         ! Compute S = (U^T)U, U -> S
         call dpotrf('U',elsi_h%n_basis,elsi_h%ovlp_real,elsi_h%n_basis,ierr)

         ! compute U^-1 -> S
         call dtrtri('U','N',elsi_h%n_basis,elsi_h%ovlp_real,elsi_h%n_basis,ierr)

         call elsi_stop_cholesky_time(elsi_h)
      endif

      call elsi_start_transform_evp_time(elsi_h)

      if(elsi_h%ovlp_is_sing) then ! Use scaled eigenvectors
         ! evec_real used as tmp_real
         ! tmp_real = H_real * S_real
         call dgemm('N','N',elsi_h%n_basis,elsi_h%n_nonsing,elsi_h%n_basis,&
                 1.0_r8,elsi_h%ham_real(1,1),elsi_h%n_basis,elsi_h%ovlp_real(1,1),&
                 elsi_h%n_basis,0.0_r8,elsi_h%evec_real(1,1),elsi_h%n_basis)

         ! H_real = (S_real)^T * tmp_real
         call dgemm('T','N',elsi_h%n_nonsing,elsi_h%n_nonsing,elsi_h%n_basis,&
                 1.0_r8,elsi_h%ovlp_real(1,1),elsi_h%n_basis,elsi_h%evec_real(1,1),&
                 elsi_h%n_basis,0.0_r8,elsi_h%ham_real(1,1),elsi_h%n_basis)

      else ! Use Cholesky
        ! tmp_real = H_real * S_real
        do n = 1,elsi_h%n_basis,nblk
           nwork = nblk

           if(n+nwork-1 > elsi_h%n_basis) then
              nwork = elsi_h%n_basis-n+1
           endif

           call dgemm('N','N',n+nwork-1,nwork,n+nwork-1,1.0_r8,elsi_h%ham_real(1,1),&
                   elsi_h%n_basis,elsi_h%ovlp_real(1,n),elsi_h%n_basis,0.0_r8,&
                   elsi_h%evec_real(1,n),elsi_h%n_basis) 
        enddo

        ! H_real = (tmp_real)*T * S_real
        do n = 1,elsi_h%n_basis,nblk
           nwork = nblk

           if(n+nwork-1 > elsi_h%n_basis) then
              nwork = elsi_h%n_basis-n+1
           endif

           call dgemm('T','N',nwork,elsi_h%n_basis-n+1,n+nwork-1,1.0_r8,&
                   elsi_h%ovlp_real(1,n),elsi_h%n_basis,elsi_h%evec_real(1,n),&
                   elsi_h%n_basis,0.0_r8,elsi_h%ham_real(n,n),elsi_h%n_basis)
        enddo
     endif

     call elsi_stop_transform_evp_time(elsi_h)

   end select

end subroutine

!> 
!! This routine back-transforms eigenvectors in the standard form (H'v' = ev')
!! to the original generalized form (Hv = eSv), by computing v = (U^-1)v'.
!!
subroutine elsi_to_original_ev_sp(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   real(kind=r8),    allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_complex(:,:)

   character*40, parameter :: caller = "elsi_to_original_ev_sp"

   call elsi_start_back_transform_ev_time(elsi_h)

   select case(elsi_h%matrix_data_type)
   case(COMPLEX_VALUES)
      if(elsi_h%ovlp_is_sing) then
         call elsi_allocate(elsi_h,tmp_complex,elsi_h%n_l_rows,&
                 elsi_h%n_l_cols,"tmp_complex",caller)
         tmp_complex = elsi_h%evec_complex

         ! Transform matrix is stored in S_complex after elsi_to_standard_evp
         call zgemm('N','N',elsi_h%n_basis,elsi_h%n_states,elsi_h%n_nonsing,&
                 (1.0_r8,0.0_r8),elsi_h%ovlp_complex(1,1),elsi_h%n_basis,&
                 tmp_complex(1,1),elsi_h%n_basis,(0.0_r8,0.0_r8),&
                 elsi_h%evec_complex(1,1),elsi_h%n_basis)

         call elsi_deallocate(elsi_h,tmp_complex,"tmp_complex")
      else ! Nonsingular, use Cholesky
         ! (U^-1) is stored in S_complex after elsi_to_standard_evp
         call ztrmm('L','U','N','N',elsi_h%n_basis,elsi_h%n_states,(1.0_r8,0.0_r8),&
                 elsi_h%ovlp_complex(1,1),elsi_h%n_basis,elsi_h%evec_complex(1,1),&
                 elsi_h%n_basis)
      endif

   case(REAL_VALUES)
      if(elsi_h%ovlp_is_sing) then
         call elsi_allocate(elsi_h,tmp_real,elsi_h%n_l_rows,&
                 elsi_h%n_l_cols,"tmp_real",caller)
         tmp_real = elsi_h%evec_real

         ! Transform matrix is stored in S_real after elsi_to_standard_evp
         call dgemm('N','N',elsi_h%n_basis,elsi_h%n_states,elsi_h%n_nonsing,1.0_r8,&
                 elsi_h%ovlp_real(1,1),elsi_h%n_basis,tmp_real(1,1),elsi_h%n_basis,&
                 0.0_r8,elsi_h%evec_real(1,1),elsi_h%n_basis)

         call elsi_deallocate(elsi_h,tmp_real,"tmp_real")
      else ! Nonsingular, use Cholesky
         ! (U^-1) is stored in S_real after elsi_to_standard_evp
         call dtrmm('L','U','N','N',elsi_h%n_basis,elsi_h%n_states,1.0_r8,&
                 elsi_h%ovlp_real(1,1),elsi_h%n_basis,elsi_h%evec_real(1,1),&
                 elsi_h%n_basis)
      endif

   end select

   call elsi_stop_back_transform_ev_time(elsi_h)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_evp_elpa_sp(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   real(kind=r8),    allocatable :: off_diag(:)
   real(kind=r8),    allocatable :: tau_real(:)
   complex(kind=r8), allocatable :: tau_complex(:)
   real(kind=r8),    allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_complex(:,:)

   logical          :: success
   integer(kind=i4) :: ierr

   character*40, parameter :: caller = "elsi_solve_evp_elpa_sp"

   ! Transform to standard form
   if(.not. elsi_h%ovlp_is_unit) then
      call elsi_to_standard_evp_sp(elsi_h)
   endif

   ! Solve evp, return eigenvalues and eigenvectors
   call elsi_statement_print("  Starting ELPA eigensolver",elsi_h)
   call elsi_start_standard_evp_time(elsi_h)

   select case(elsi_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(elsi_h,off_diag,elsi_h%n_basis,"off_diag",caller)
      call elsi_allocate(elsi_h,tau_complex,elsi_h%n_basis,"tau_complex",caller)
      call elsi_allocate(elsi_h,tmp_real,elsi_h%n_basis,elsi_h%n_basis,"tmp_real",caller)
      call elsi_allocate(elsi_h,tmp_complex,elsi_h%n_basis,elsi_h%n_basis,"tmp_complex",caller)

      call zhetrd('U',elsi_h%n_basis,elsi_h%ham_complex,elsi_h%n_basis,elsi_h%eval,&
              off_diag,tau_complex,tmp_complex,elsi_h%n_basis*elsi_h%n_basis,ierr)

      success = elpa_solve_tridi_double(elsi_h%n_basis,elsi_h%n_states,elsi_h%eval,&
                   off_diag,tmp_real,elsi_h%n_basis,64,elsi_h%n_basis,mpi_comm_self,&
                   mpi_comm_self,.false.)

      elsi_h%evec_complex(1:elsi_h%n_basis,1:elsi_h%n_states) = &
         tmp_real(1:elsi_h%n_basis,1:elsi_h%n_states)

      call zunmtr('L','U','N',elsi_h%n_basis,elsi_h%n_states,elsi_h%ham_complex,&
              elsi_h%n_basis,tau_complex,elsi_h%evec_complex,elsi_h%n_basis,&
              tmp_complex,elsi_h%n_basis*elsi_h%n_basis,ierr)

      call elsi_deallocate(elsi_h,off_diag,"off_diag")
      call elsi_deallocate(elsi_h,tau_complex,"tau_complex")
      call elsi_deallocate(elsi_h,tmp_real,"tmp_real")
      call elsi_deallocate(elsi_h,tmp_complex,"tmp_complex")

   case(REAL_VALUES)
      call elsi_allocate(elsi_h,off_diag,elsi_h%n_basis,"off_diag",caller)
      call elsi_allocate(elsi_h,tau_real,elsi_h%n_basis,"tau_real",caller)
      call elsi_allocate(elsi_h,tmp_real,elsi_h%n_basis,elsi_h%n_basis,"tmp_real",caller)

      call dsytrd('U',elsi_h%n_basis,elsi_h%ham_real,elsi_h%n_basis,elsi_h%eval,&
              off_diag,tau_real,tmp_real,elsi_h%n_basis*elsi_h%n_basis,ierr)

      success = elpa_solve_tridi_double(elsi_h%n_basis,elsi_h%n_states,elsi_h%eval,&
                   off_diag,tmp_real,elsi_h%n_basis,64,elsi_h%n_basis,mpi_comm_self,&
                   mpi_comm_self,.false.)

      elsi_h%evec_real(1:elsi_h%n_basis,1:elsi_h%n_states) = &
         tmp_real(1:elsi_h%n_basis,1:elsi_h%n_states)

      call dormtr('L','U','N',elsi_h%n_basis,elsi_h%n_states,elsi_h%ham_real,&
              elsi_h%n_basis,tau_real,elsi_h%evec_real,elsi_h%n_basis,tmp_real,&
              elsi_h%n_basis*elsi_h%n_basis,ierr)

      call elsi_deallocate(elsi_h,off_diag,"off_diag")
      call elsi_deallocate(elsi_h,tau_real,"tau_real")
      call elsi_deallocate(elsi_h,tmp_real,"tmp_real")
   end select

   call elsi_stop_standard_evp_time(elsi_h)

   ! Back-transform eigenvectors
   if(.not. elsi_h%ovlp_is_unit) then
      call elsi_to_original_ev_sp(elsi_h)
   endif

end subroutine

!> 
!! This routine checks the singularity of overlap matrix by computing all
!! its eigenvalues.
!!
!! On exit, S is not modified if not singular, while overwritten by scaled
!! eigenvectors if singular, which is used to transform the generalized
!! eigenvalue problem to standard form without using Cholesky.
!!
subroutine elsi_check_singularity_sp(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   real(kind=r8)    :: ev_sqrt
   integer(kind=i4) :: i
   integer(kind=i4) :: ierr
   logical          :: success
   character*200    :: info_str

   real(kind=r8),    allocatable :: off_diag(:)
   real(kind=r8),    allocatable :: tau_real(:)
   complex(kind=r8), allocatable :: tau_complex(:)
   real(kind=r8),    allocatable :: tmp_real(:,:)
   complex(kind=r8), allocatable :: tmp_complex(:,:)
   real(kind=r8),    allocatable :: copy_real(:,:)
   complex(kind=r8), allocatable :: copy_complex(:,:)

   character*40, parameter :: caller = "elsi_check_singularity_sp"

   call elsi_start_singularity_check_time(elsi_h)

   select case(elsi_h%matrix_data_type)
   case(COMPLEX_VALUES)
      call elsi_allocate(elsi_h,copy_complex,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,"copy_complex",caller)

      ! Use copy_complex to store overlap matrix, otherwise it will
      ! be destroyed by eigenvalue calculation
      copy_complex = -elsi_h%ovlp_complex

      call elsi_allocate(elsi_h,off_diag,elsi_h%n_basis,"off_diag",caller)
      call elsi_allocate(elsi_h,tau_complex,elsi_h%n_basis,"tau_complex",caller)
      call elsi_allocate(elsi_h,tmp_real,elsi_h%n_basis,elsi_h%n_basis,"tmp_real",caller)
      call elsi_allocate(elsi_h,tmp_complex,elsi_h%n_basis,elsi_h%n_basis,"tmp_complex",caller)

      call zhetrd('U',elsi_h%n_basis,copy_complex,elsi_h%n_basis,elsi_h%eval,off_diag,&
              tau_complex,tmp_complex,elsi_h%n_basis*elsi_h%n_basis,ierr)

      success = elpa_solve_tridi_double(elsi_h%n_basis,elsi_h%n_basis,elsi_h%eval,&
                   off_diag,tmp_real,elsi_h%n_basis,64,elsi_h%n_basis,mpi_comm_self,&
                   mpi_comm_self,.false.)

      ! Get the number of nonsingular eigenvalues
      elsi_h%eval = -elsi_h%eval

      do i = 1,elsi_h%n_basis
         if(elsi_h%eval(i) < elsi_h%sing_tol) exit
      enddo

      elsi_h%n_nonsing = i-1

      ! Eigenvectors computed only for singular overlap matrix
      if(elsi_h%n_nonsing < elsi_h%n_basis) then
         elsi_h%evec_complex = tmp_real

         call zunmtr('L','U','N',elsi_h%n_basis,elsi_h%n_basis,copy_complex,&
                 elsi_h%n_basis,tau_complex,elsi_h%evec_complex,elsi_h%n_basis,&
                 tmp_complex,elsi_h%n_basis*elsi_h%n_basis,ierr)
      endif

      call elsi_deallocate(elsi_h,off_diag,"off_diag")
      call elsi_deallocate(elsi_h,tau_complex,"tau_complex")
      call elsi_deallocate(elsi_h,tmp_real,"tmp_real")
      call elsi_deallocate(elsi_h,tmp_complex,"tmp_complex")
      call elsi_deallocate(elsi_h,copy_complex,"copy_complex")

      ! Stop if n_states is larger than n_nonsing
      if(elsi_h%n_nonsing < elsi_h%n_states) then ! Too singular to continue
         call elsi_stop(" Overlap matrix is singular. The number of"//&
                 " basis functions after removing singularity is smaller"//&
                 " than the number of states. Exiting...",elsi_h,caller)
      elseif(elsi_h%n_nonsing < elsi_h%n_basis) then ! Singular
         elsi_h%ovlp_is_sing = .true.

         if(elsi_h%stop_sing) then
            call elsi_stop(" Overlap matrix is singular. Running with a"//&
                    " near-singular basis set may lead to completely"//&
                    " wrong numerical results. Exiting...",elsi_h,caller)
         endif

         call elsi_statement_print("  Overlap matrix is singular. Note"//&
                 " that running with a near-singular basis set may lead"//&
                 " to completely wrong numerical results.",elsi_h)

         write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
            elsi_h%n_nonsing
         call elsi_statement_print(info_str,elsi_h)

         call elsi_statement_print("  Using scaled eigenvectors of"//&
                 " overlap matrix for transformation",elsi_h)

         ! Overlap matrix is overwritten with scaled eigenvectors
         do i = 1,elsi_h%n_nonsing
            ev_sqrt = sqrt(elsi_h%eval(i))
            elsi_h%ovlp_complex(:,i) = elsi_h%evec_complex(:,i)/ev_sqrt
         enddo

      else ! Nonsingular
         elsi_h%ovlp_is_sing = .false.
         call elsi_statement_print("  Overlap matrix is nonsingular",elsi_h)
      endif ! Singular overlap?

   case(REAL_VALUES)
      call elsi_allocate(elsi_h,copy_real,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,"copy_real",caller)

      ! Use copy_real to store overlap matrix, otherwise it will be
      ! destroyed by eigenvalue calculation
      copy_real = -elsi_h%ovlp_real

      call elsi_allocate(elsi_h,off_diag,elsi_h%n_basis,"off_diag",caller)
      call elsi_allocate(elsi_h,tau_real,elsi_h%n_basis,"tau_real",caller)
      call elsi_allocate(elsi_h,tmp_real,elsi_h%n_basis,elsi_h%n_basis,"tmp_real",caller)

      call dsytrd('U',elsi_h%n_basis,copy_real,elsi_h%n_basis,elsi_h%eval,off_diag,&
              tau_real,tmp_real,elsi_h%n_basis*elsi_h%n_basis,ierr)

      success = elpa_solve_tridi_double(elsi_h%n_basis,elsi_h%n_basis,elsi_h%eval,&
                   off_diag,tmp_real,elsi_h%n_basis,64,elsi_h%n_basis,mpi_comm_self,&
                   mpi_comm_self,.false.)

      ! Get the number of nonsingular eigenvalues
      elsi_h%eval = -elsi_h%eval

      do i = 1,elsi_h%n_basis
         if(elsi_h%eval(i) < elsi_h%sing_tol) exit
      enddo

      elsi_h%n_nonsing = i-1

      ! Eigenvectors computed only for singular overlap matrix
      if(elsi_h%n_nonsing < elsi_h%n_basis) then
         elsi_h%evec_real = tmp_real

         call dormtr('L','U','N',elsi_h%n_basis,elsi_h%n_basis,copy_real,&
                 elsi_h%n_basis,tau_real,elsi_h%evec_real,elsi_h%n_basis,tmp_real,&
                 elsi_h%n_basis*elsi_h%n_basis,ierr)
      endif

      call elsi_deallocate(elsi_h,off_diag,"off_diag")
      call elsi_deallocate(elsi_h,tau_real,"tau_real")
      call elsi_deallocate(elsi_h,tmp_real,"tmp_real")
      call elsi_deallocate(elsi_h,copy_real,"copy_real")

      ! Stop if n_states is larger than n_nonsing
      if(elsi_h%n_nonsing < elsi_h%n_states) then ! Too singular to continue
         call elsi_stop(" Overlap matrix is singular. The number of"//&
                 " basis functions after removing singularity is smaller"//&
                 " than the number of states. Exiting...",elsi_h,caller)
      elseif(elsi_h%n_nonsing < elsi_h%n_basis) then ! Singular
         elsi_h%ovlp_is_sing = .true.

         if(elsi_h%stop_sing) then
            call elsi_stop(" Overlap matrix is singular. Running with a"//&
                    " near-singular basis set may lead to completely wrong"//&
                    " numerical results. Exiting...",elsi_h,caller)
         endif

         call elsi_statement_print("  Overlap matrix is singular. Note"//&
                 " that running with a near-singular basis set may lead to"//&
                 " completely wrong numerical results.",elsi_h)

         write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
            elsi_h%n_nonsing
         call elsi_statement_print(info_str,elsi_h)

         call elsi_statement_print("  Using scaled eigenvectors of"//&
                 " overlap matrix for transformation",elsi_h)

         ! Overlap matrix is overwritten with scaled eigenvectors
         do i = 1,elsi_h%n_nonsing
            ev_sqrt = sqrt(elsi_h%eval(i))
            elsi_h%ovlp_real(:,i) = elsi_h%evec_real(:,i)/ev_sqrt
         enddo

      else ! Nonsingular
         elsi_h%ovlp_is_sing = .false.
         call elsi_statement_print("  Overlap matrix is nonsingular",elsi_h)
      endif ! Singular overlap?

   end select ! select matrix_data_type

   call elsi_stop_singularity_check_time(elsi_h)

end subroutine

!>
!! This routine computes the chemical potential and occupation numbers.
!!
subroutine elsi_compute_occ_elpa(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   real(kind=r8)    :: e_low         ! Lowest eigenvalue
   real(kind=r8)    :: e_high        ! Highest eigenvalue
   real(kind=r8)    :: mu_lower      ! Lower bound of chemical potential
   real(kind=r8)    :: mu_upper      ! Upper bound of chemical potential
   real(kind=r8)    :: diff_ne_lower ! Difference in number of electrons on lower bound
   real(kind=r8)    :: diff_ne_upper ! Difference in number of electrons on upper bound
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: n_steps
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_compute_occ_elpa"

   call elsi_get_eval_all(elsi_h)

   ! Determine the smallest and largest eivenvalues
   e_low = elsi_h%eval_all(1,1,1)
   e_high = elsi_h%eval_all(elsi_h%n_states,1,1)

   do i_kpt = 1,elsi_h%n_kpts
      do i_spin = 1,elsi_h%n_spins
         do i_state = 1,elsi_h%n_states
            if(elsi_h%eval_all(i_state,i_spin,i_kpt) < e_low) then
               e_low = elsi_h%eval_all(i_state,i_spin,i_kpt)
            endif
            if(elsi_h%eval_all(i_state,i_spin,i_kpt) > e_high) then
               e_high = elsi_h%eval_all(i_state,i_spin,i_kpt)
            endif
         enddo
      enddo
   enddo

   ! Determine the upper and lower bounds for the chemical potential search
   mu_lower = e_low

   if(e_low == e_high) then
      mu_upper = 0.0_r8
   else
      mu_upper = e_high
   endif

   elsi_h%occ_num = 0.0_r8

   ! Compute the error in electron count
   call elsi_check_electrons(elsi_h,mu_lower,diff_ne_lower)
   call elsi_check_electrons(elsi_h,mu_upper,diff_ne_upper)

   ! If diff_ne_lower*diff_ne_upper > 0, the solution is not in this interval.
   ! Enlarge the interval towards both sides, then recheck.
   n_steps = 0
   do while(diff_ne_lower*diff_ne_upper > 0)
      n_steps = n_steps+1
      if(n_steps > elsi_h%max_mu_steps) then
         write(info_str,"(A,I13,A)") " Chemical potential not found in ",&
            elsi_h%max_mu_steps," iterations! Exiting..."
         call elsi_stop(info_str,elsi_h,caller)
      endif

      mu_lower = mu_lower-0.5_r8*abs(e_high-e_low)
      mu_upper = mu_upper+0.5_r8*abs(e_high-e_low)

      call elsi_check_electrons(elsi_h,mu_lower,diff_ne_lower)
      call elsi_check_electrons(elsi_h,mu_upper,diff_ne_upper)
   enddo

   ! Now the solution should lie in the interval. Use bisection to find it.
   call elsi_find_mu(elsi_h,mu_lower,mu_upper)

end subroutine

!>
!! This routine computes the number of electrons for a given chemical potential,
!! and returns the error in the number of electrons. The occupation numbers will
!! be updated as well.
!!
subroutine elsi_check_electrons(elsi_h,mu_in,diff_ne_out)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h      !< Handle
   real(kind=r8),     intent(in)    :: mu_in       !< Input chemical potential
   real(kind=r8),     intent(out)   :: diff_ne_out !< Difference in number of electrons

   real(kind=r8) :: invert_width ! 1/broaden_width
   real(kind=r8) :: max_exp ! Maximum possible exponent
   real(kind=r8) :: this_exp
   real(kind=r8) :: this_hermite

   integer(kind=i4) :: n_steps
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin

   character*40, parameter :: caller = "elsi_check_electrons"

   if(elsi_h%broadening_width .le. 0.0_r8) then
      call elsi_stop(" Broadening width in chemical potential determination must"//&
              " be a positive number. Exiting...",elsi_h,caller)
   endif

   invert_width = 1.0_r8/elsi_h%broadening_width
   diff_ne_out = 0.0_r8

   if(elsi_h%n_spins == 2) then
      elsi_h%spin_degen = 1.0_r8
   else
      elsi_h%spin_degen = 2.0_r8
   endif

   select case(elsi_h%broadening_scheme)
   case(GAUSSIAN)
      do i_kpt = 1,elsi_h%n_kpts
         do i_spin = 1,elsi_h%n_spins
            do i_state = 1,elsi_h%n_states
               elsi_h%occ_num(i_state,i_spin,i_kpt) = elsi_h%spin_degen*0.5_r8*&
                  (1.0_r8-erf((elsi_h%eval_all(i_state,i_spin,i_kpt)-mu_in)*invert_width))

               diff_ne_out = elsi_h%occ_num(i_state,i_spin,i_kpt)*elsi_h%k_weight(i_kpt)+&
                                diff_ne_out
            enddo
         enddo
      enddo

   case(FERMI)
      max_exp = maxexponent(mu_in)*log(2.0_r8)

      do i_kpt = 1,elsi_h%n_kpts
         do i_spin = 1,elsi_h%n_spins
            do i_state = 1,elsi_h%n_states
               this_exp = (elsi_h%eval_all(i_state,i_spin,i_kpt)-mu_in)*invert_width

               if(this_exp < max_exp) then
                  elsi_h%occ_num(i_state,i_spin,i_kpt) = elsi_h%spin_degen/(1.0_r8+exp(this_exp))

                  diff_ne_out = elsi_h%occ_num(i_state,i_spin,i_kpt)*elsi_h%k_weight(i_kpt)+&
                                   diff_ne_out
               else ! Exponent in this step is larger than the largest possible exponent
                  elsi_h%occ_num(i_state,i_spin,i_kpt) = 0.0_r8
               endif
            enddo
         enddo
      enddo

   case(METHFESSEL_PAXTON_0)
      do i_kpt = 1,elsi_h%n_kpts
         do i_spin = 1,elsi_h%n_spins
            do i_state = 1,elsi_h%n_states
               elsi_h%occ_num(i_state,i_spin,i_kpt) = elsi_h%spin_degen*0.5_r8*&
                  (1.0_r8-erf((elsi_h%eval_all(i_state,i_spin,i_kpt)-mu_in)*invert_width))

               diff_ne_out = elsi_h%occ_num(i_state,i_spin,i_kpt)*elsi_h%k_weight(i_kpt)+&
                                diff_ne_out
            enddo
         enddo
      enddo

   case(METHFESSEL_PAXTON_1)
      do i_kpt = 1,elsi_h%n_kpts
         do i_spin = 1,elsi_h%n_spins
            do i_state = 1,elsi_h%n_states
               this_hermite = (elsi_h%eval_all(i_state,i_spin,i_kpt)-mu_in)*invert_width

               elsi_h%occ_num(i_state,i_spin,i_kpt) = elsi_h%spin_degen*0.5_r8*&
                  (1.0_r8-erf(this_hermite))-0.5_r8*invert_sqrt_pi*this_hermite*&
                  exp(-this_hermite*this_hermite)

               diff_ne_out = elsi_h%occ_num(i_state,i_spin,i_kpt)*elsi_h%k_weight(i_kpt)+&
                                diff_ne_out
            enddo
         enddo
      enddo
   end select

   diff_ne_out = diff_ne_out-elsi_h%n_electrons

end subroutine

!>
!! This routine computes the chemical potential using a bisection algorithm.
!!
subroutine elsi_find_mu(elsi_h,mu_lower_in,mu_upper_in)

   implicit none
 
   type(elsi_handle), intent(inout) :: elsi_h      !< Handle
   real(kind=r8),     intent(in)    :: mu_lower_in !< Lower bound of chemical potential
   real(kind=r8),     intent(in)    :: mu_upper_in !< Upper bound of chemical potential

   real(kind=r8)    :: mu_left    ! Left bound of chemical potential interval
   real(kind=r8)    :: mu_right   ! Right bound of chemical potential interval
   real(kind=r8)    :: mu_mid     ! Middle point of chemical potential interval
   real(kind=r8)    :: diff_left  ! Difference in number of electrons on left bound
   real(kind=r8)    :: diff_right ! Difference in number of electrons on right bound
   real(kind=r8)    :: diff_mid   ! Difference in number of electrons on middle point
   logical          :: found_mu
   integer(kind=i4) :: n_steps

   character*40, parameter :: caller = "elsi_find_mu"

   n_steps = 0
   found_mu = .false.

   mu_left = mu_lower_in
   mu_right = mu_upper_in

   do while((.not.found_mu) .and. (n_steps < elsi_h%max_mu_steps))
      call elsi_check_electrons(elsi_h,mu_left,diff_left)
      call elsi_check_electrons(elsi_h,mu_right,diff_right)

      if(abs(diff_left) < elsi_h%occ_tolerance) then
         elsi_h%mu = mu_left
         found_mu = .true.
      elseif(abs(diff_right) < elsi_h%occ_tolerance) then
         elsi_h%mu = mu_right
         found_mu = .true.
      else
         n_steps = n_steps+1

         mu_mid = 0.5_r8*(mu_left+mu_right)

         call elsi_check_electrons(elsi_h,mu_mid,diff_mid)

         if(abs(diff_mid) < elsi_h%occ_tolerance) then
            elsi_h%mu = mu_mid
            found_mu = .true.
         elseif(diff_mid < 0) then
            mu_left = mu_mid
         elseif(diff_mid > 0) then
            mu_right = mu_mid
         endif
      endif
   enddo

   ! Special treatment if mu cannot reach the required accuracy
   if(.not. found_mu) then
      ! Use the chemical potential of the right bound...
      call elsi_check_electrons(elsi_h,mu_right,diff_right)

      elsi_h%mu = mu_right

      ! ...with adjusted occupation numbers
      call elsi_statement_print("  Chemical potential cannot reach the"//&
              " required accuracy by bisection method. The error will be"//&
              " arbitrarily removed from the highest occupied states.",elsi_h)

      call elsi_adjust_occ(elsi_h,diff_right)
   endif

end subroutine

!>
!! This routine cancels the small error in number of electrons.
!!
subroutine elsi_adjust_occ(elsi_h,diff_ne)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h  !< Handle
   real(kind=r8),     intent(inout) :: diff_ne !< Error in electron count

   real(kind=r8),    allocatable :: eval_aux(:)

   real(kind=r8)    :: min_eval
   integer(kind=i4) :: max_id
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: i_val
   integer(kind=i4) :: n_total

   character*40, parameter :: caller = "elsi_adjust_occ"

   n_total = elsi_h%n_states*elsi_h%n_spins*elsi_h%n_kpts

   call elsi_allocate(elsi_h,eval_aux,n_total,"eval_aux",caller)

   ! Put eigenvalues into a 1D array
   i_val = 0

   do i_kpt = 1,elsi_h%n_kpts
      do i_spin = 1,elsi_h%n_spins
         do i_state = 1,elsi_h%n_states
            i_val = i_val+1
            eval_aux(i_val) = elsi_h%eval_all(i_state,i_spin,i_kpt)
         enddo
      enddo
   enddo

   min_eval = minval(eval_aux,1)

   ! Remove error
   do i_val = 1,n_total
      max_id           = maxloc(eval_aux,1)
      eval_aux(max_id) = min_eval-1.0_r8

      i_kpt   = (i_val-1)/(elsi_h%n_spins*elsi_h%n_states)+1
      i_spin  = mod((i_val-1)/elsi_h%n_states,elsi_h%n_spins)+1
      i_state = mod(i_val-1,elsi_h%n_states)+1

      if(elsi_h%k_weight(i_kpt)*elsi_h%occ_num(i_state,i_spin,i_kpt) > diff_ne) then
         elsi_h%occ_num(i_state,i_spin,i_kpt) = elsi_h%occ_num(i_state,i_spin,i_kpt)-&
                                                   diff_ne/elsi_h%k_weight(i_kpt)
         diff_ne = 0.0_r8
      else
         diff_ne = diff_ne-elsi_h%k_weight(i_kpt)*elsi_h%occ_num(i_state,i_spin,i_kpt)
         elsi_h%occ_num(i_state,i_spin,i_kpt) = 0.0_r8
      endif

      if(diff_ne .le. elsi_h%occ_tolerance) exit
   enddo

   call elsi_deallocate(elsi_h,eval_aux,"eval_aux")

end subroutine

!>
!! This routine collects the eigenvalues of all spins and k-points.
!!
subroutine elsi_get_eval_all(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h  !< Handle

   real(kind=r8), allocatable :: tmp_real1(:,:,:)
   real(kind=r8), allocatable :: tmp_real2(:)
   integer(kind=i4)           :: mpierr

   character*40, parameter :: caller = "elsi_get_eval_all"

   if(.not. allocated(elsi_h%eval_all)) then
      call elsi_allocate(elsi_h,elsi_h%eval_all,elsi_h%n_states,&
              elsi_h%n_spins,elsi_h%n_kpts,"eval_all",caller)
   endif

   elsi_h%eval_all = 0.0_r8
   elsi_h%eval_all(1:elsi_h%n_states,elsi_h%i_spin,elsi_h%i_kpt) = &
      elsi_h%eval_elpa(1:elsi_h%n_states)

   if(.not. allocated(elsi_h%occ_num)) then
      call elsi_allocate(elsi_h,elsi_h%occ_num,elsi_h%n_states,&
              elsi_h%n_spins,elsi_h%n_kpts,"occ_num",caller)
   endif

   elsi_h%occ_num = 0.0_r8

   if(.not. allocated(elsi_h%k_weight)) then
      call elsi_allocate(elsi_h,elsi_h%k_weight,&
              elsi_h%n_kpts,"k_weight",caller)
   endif

   elsi_h%k_weight = 0.0_r8
   elsi_h%k_weight(elsi_h%i_kpt) = elsi_h%i_weight

   if(elsi_h%n_spins*elsi_h%n_kpts > 1) then
      call elsi_allocate(elsi_h,tmp_real1,elsi_h%n_states,&
              elsi_h%n_spins,elsi_h%n_kpts,"tmp_real",caller)

      if(elsi_h%myid /= 0) then
         elsi_h%eval_all = 0.0_r8
         elsi_h%k_weight = 0.0_r8
      endif

      call MPI_Allreduce(elsi_h%eval_all,tmp_real1,&
              elsi_h%n_states*elsi_h%n_spins*elsi_h%n_kpts,&
              mpi_real8,mpi_sum,elsi_h%mpi_comm_all,mpierr)

      elsi_h%eval_all = tmp_real1

      call elsi_deallocate(elsi_h,tmp_real1,"tmp_real")

      call elsi_allocate(elsi_h,tmp_real2,elsi_h%n_kpts,"tmp_real",caller)

      call MPI_Allreduce(elsi_h%k_weight,tmp_real2,elsi_h%n_kpts,&
              mpi_real8,mpi_sum,elsi_h%mpi_comm_all,mpierr)

      elsi_h%k_weight = tmp_real2

      call elsi_deallocate(elsi_h,tmp_real2,"tmp_real")
   endif

end subroutine

end module ELSI_ELPA
