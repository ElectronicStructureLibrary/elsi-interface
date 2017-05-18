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

   use iso_c_binding
   use ELSI_PRECISION, only: r8,i4
   use ELSI_CONSTANTS, only: REAL_VALUES,COMPLEX_VALUES
   use ELSI_DIMENSIONS, only: elsi_handle
   use ELSI_TIMERS
   use ELSI_UTILS
   use ELSI_MU
   use ELPA1
   use ELPA2
   use CHECK_SINGULARITY, only: elpa_check_singularity_real_double, elpa_check_singularity_complex_double

   implicit none
   private

   public :: elsi_get_elpa_comms
   public :: elsi_compute_occ_elpa
   public :: elsi_compute_dm_elpa
   public :: elsi_solve_evp_elpa
   public :: elsi_solve_evp_elpa_sp

contains

!========================
! ELSI routines for ELPA
!========================

!>
!! This routine gets the row and column communicators for ELPA.
!!
subroutine elsi_get_elpa_comms(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   integer(kind=i4) :: success

   character*40, parameter :: caller = "elsi_get_elpa_comms"

   success = elpa_get_communicators(elsi_h%mpi_comm,elsi_h%my_p_row,elsi_h%my_p_col,&
                                    elsi_h%mpi_comm_row,elsi_h%mpi_comm_col)

end subroutine

!>
!! This routine computes the chemical potential and occupation numbers
!! from eigenvalues.
!!
subroutine elsi_compute_occ_elpa(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   real(kind=r8) :: mu           !< Chemical potential
   real(kind=r8) :: k_weights(1) !< Weights of k-points

   real(kind=r8), allocatable :: eval_aux(:,:,:)
   real(kind=r8), allocatable :: occ_aux(:,:,:)

   !< Currently this subroutine only supports 1 spin channel and 1 k-point
   integer(kind=i4), parameter :: n_spin   = 1
   integer(kind=i4), parameter :: n_kpoint = 1

   character*40, parameter :: caller = "elsi_compute_occ_elpa"

   k_weights(1) = 1.0_r8

   if(.not. allocated(elsi_h%occ_elpa)) then
       call elsi_allocate(elsi_h,elsi_h%occ_elpa,elsi_h%n_states,"occ_elpa",caller)
   endif

   call elsi_allocate(elsi_h,eval_aux,elsi_h%n_states,1,1,"eval_aux",caller)
   call elsi_allocate(elsi_h,occ_aux,elsi_h%n_states,1,1,"occ_aux",caller)

   eval_aux(1:elsi_h%n_states,1,1) = elsi_h%eval(1:elsi_h%n_states)
   occ_aux = 0.0_r8

   call elsi_compute_mu_and_occ(elsi_h,elsi_h%n_electrons,elsi_h%n_states,n_spin,&
                                n_kpoint,k_weights,eval_aux,occ_aux,mu)

   elsi_h%occ_elpa(:) = occ_aux(:,1,1)

   deallocate(eval_aux)
   deallocate(occ_aux)

end subroutine

!>
!! This routine constructs the density matrix using eigenvectors from ELPA.
!!
subroutine elsi_compute_dm_elpa(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   real(kind=r8),    allocatable :: tmp_real(:,:)    !< Real eigenvectors, temporary
   complex(kind=r8), allocatable :: tmp_complex(:,:) !< Complex eigenvectors, temporary
   real(kind=r8),    allocatable :: factor(:)        !< Factor to construct density matrix
   integer(kind=i4)              :: i,i_col,i_row

   character*40, parameter :: caller = "elsi_compute_dm_elpa"

   call elsi_start_density_matrix_time(elsi_h)

   select case (elsi_h%matrix_data_type)
      case (REAL_VALUES)
         ! Get eigenvectors into tmp_real
         call elsi_allocate(elsi_h,tmp_real,elsi_h%n_l_rows,elsi_h%n_l_cols,"tmp_real",caller)
         tmp_real = elsi_h%evec_real

         ! Compute the factors used to construct density matrix
         call elsi_allocate(elsi_h,factor,elsi_h%n_states,"factor",caller)
         factor = 0.0_r8

         do i = 1,elsi_h%n_states
            if(elsi_h%occ_elpa(i) > 0.0_r8) then
               factor(i) = sqrt(elsi_h%occ_elpa(i))
            endif
         enddo

         do i = 1,elsi_h%n_states
            if(factor(i) > 0.0_r8) then
               if(elsi_h%local_col(i) > 0) then
                  tmp_real(:,elsi_h%local_col(i)) = tmp_real(:,elsi_h%local_col(i))*factor(i)
               endif
            elseif(elsi_h%local_col(i) .ne. 0) then
               tmp_real(:,elsi_h%local_col(i)) = 0.0_r8
            endif
         enddo

         ! Compute density matrix
         elsi_h%den_mat = 0.0_r8

         ! D_elpa = tmp_real*tmp_real^T
         call pdsyrk('U','N',elsi_h%n_g_size,elsi_h%n_states,1.0_r8,tmp_real,&
                     1,1,elsi_h%sc_desc,0.0_r8,elsi_h%den_mat,1,1,elsi_h%sc_desc)

      case (COMPLEX_VALUES)
         ! Get eigenvectors into tmp_complex
         call elsi_allocate(elsi_h,tmp_complex,elsi_h%n_l_rows,elsi_h%n_l_cols,"tmp_complex",caller)
         tmp_complex = elsi_h%evec_complex

         ! Compute the factors used to construct density matrix
         call elsi_allocate(elsi_h,factor,elsi_h%n_states,"factor",caller)
         factor = 0.0_r8

         do i = 1,elsi_h%n_states
            if(elsi_h%occ_elpa(i) > 0.0_r8) then
               factor(i) = sqrt(elsi_h%occ_elpa(i))
            endif
         enddo

         do i = 1,elsi_h%n_states
            if(factor(i) > 0.0_r8) then
               if(elsi_h%local_col(i) > 0) then
                  tmp_complex(:,elsi_h%local_col(i)) = tmp_complex(:,elsi_h%local_col(i))*factor(i)
               endif
            elseif(elsi_h%local_col(i) .ne. 0) then
               tmp_complex(:,elsi_h%local_col(i)) = (0.0_r8,0.0_r8)
            endif
         enddo

         ! Compute density matrix
         elsi_h%den_mat = 0.0_r8
         call elsi_allocate(elsi_h,tmp_real,elsi_h%n_l_rows,elsi_h%n_l_cols,"tmp_real",caller)

         call pdsyrk('U','N',elsi_h%n_g_size,elsi_h%n_states,1.0_r8,real(tmp_complex),&
                     1,1,elsi_h%sc_desc,0.0_r8,elsi_h%den_mat,1,1,elsi_h%sc_desc)
         call pdsyrk('U','N',elsi_h%n_g_size,elsi_h%n_states,1.0_r8,aimag(tmp_complex),&
                     1,1,elsi_h%sc_desc,0.0_r8,tmp_real,1,1,elsi_h%sc_desc)

         elsi_h%den_mat = elsi_h%den_mat+tmp_real
   end select

   deallocate(factor)
   if(allocated(tmp_real))    deallocate(tmp_real)
   if(allocated(tmp_complex)) deallocate(tmp_complex)

   ! Set upper triangle matrix den_mat to full form
   call elsi_allocate(elsi_h,tmp_real,elsi_h%n_l_rows,elsi_h%n_l_cols,"tmp_real",caller)

   call pdtran(elsi_h%n_g_size,elsi_h%n_g_size,1.0_r8,elsi_h%den_mat,1,1,&
               elsi_h%sc_desc,0.0_r8,tmp_real,1,1,elsi_h%sc_desc)

   do i_col = 1,elsi_h%n_g_size-1
      if(elsi_h%local_col(i_col) == 0) cycle
      do i_row = i_col+1,elsi_h%n_g_size
         if(elsi_h%local_row(i_row) > 0) then
            elsi_h%den_mat(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
               tmp_real(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
         endif
      enddo
   enddo

   deallocate(tmp_real)

   call elsi_stop_density_matrix_time(elsi_h)

end subroutine

!> 
!! This routine transforms a generalized eigenvalue problem (Ac = Bcv)
!! to standard form (A'c' = c'v)
!!
!! Starting from Hv = eSv, we first perform a Cholesky decomposition of S
!! S = (U^T)U, resulting in Hv = e(U^T)Uv
!!
!! Using 1=U^-1U we define a new standard eigenvalue problem by
!! H(U^-1)(Uv) = e(U^T)(Uv) => ((U^-1)^T)H(U^-1)(Uv) = e(Uv)
!!
!! On exit, (U^-1) is stored in S, to be used for back-transformation
!!
subroutine elsi_to_standard_evp(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   integer(kind=i4)              :: i_row,i_col
   real(kind=r8),    allocatable :: buffer_real(:,:)
   complex(kind=r8), allocatable :: buffer_complex(:,:)
   logical                       :: success

   character*40, parameter :: caller = "elsi_to_standard_evp"

   select case (elsi_h%matrix_data_type)
      case (COMPLEX_VALUES)
         if(elsi_h%n_elsi_calls == 1) then
            if(.not. elsi_h%no_singularity_check) then
               call elsi_check_singularity(elsi_h)
            endif

            if(elsi_h%n_nonsingular == elsi_h%n_g_size) then ! Not singular
               call elsi_start_cholesky_time(elsi_h)

               elsi_h%overlap_is_singular = .false.

               ! Compute S = (U^T)U, U -> S
               success = elpa_cholesky_complex_double(elsi_h%n_g_size,&
                            elsi_h%ovlp_complex,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                            elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                            .false.)
               if(.not. success) then
                  call elsi_stop(" Cholesky decomposition failed.",elsi_h,caller)
               endif

               ! compute U^-1 -> S
               success = elpa_invert_trm_complex_double(elsi_h%n_g_size,&
                            elsi_h%ovlp_complex,elsi_h%n_l_rows,elsi_h%n_b_rows,&
                            elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                            .false.)
               if(.not. success) then
                  call elsi_stop(" Matrix inversion failed.",elsi_h,caller)
               endif

               call elsi_stop_cholesky_time(elsi_h)
            endif
         endif ! n_elsi_calls == 1

         call elsi_start_transform_evp_time(elsi_h)

         call elsi_allocate(elsi_h,buffer_complex,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)

         if(elsi_h%overlap_is_singular) then ! Use scaled eigenvectors
            ! buffer_complex = H_complex * S_complex
            call pzgemm('N','N',elsi_h%n_g_size,elsi_h%n_nonsingular,elsi_h%n_g_size,&
                        (1.0_r8,0.0_r8),elsi_h%ham_complex,1,1,elsi_h%sc_desc,&
                        elsi_h%ovlp_complex,1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),&
                        buffer_complex,1,1,elsi_h%sc_desc)

            ! H_complex = (S_complex)^* * buffer_complex
            call pzgemm('C','N',elsi_h%n_nonsingular,elsi_h%n_nonsingular,elsi_h%n_g_size,&
                        (1.0_r8,0.0_r8),elsi_h%ovlp_complex,1,1,elsi_h%sc_desc,&
                        buffer_complex,1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),elsi_h%ham_complex,&
                        1,1,elsi_h%sc_desc)

         else ! Use cholesky
            success = elpa_mult_ah_b_complex_double('U','L',elsi_h%n_g_size,elsi_h%n_g_size,&
                         elsi_h%ovlp_complex,elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%ham_complex,&
                         elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%n_b_rows,elsi_h%mpi_comm_row,&
                         elsi_h%mpi_comm_col,buffer_complex,elsi_h%n_l_rows,elsi_h%n_l_cols)

            call pztranc(elsi_h%n_g_size,elsi_h%n_g_size,(1.0_r8,0.0_r8),buffer_complex,1,1,&
                         elsi_h%sc_desc,(0.0_r8,0.0_r8),elsi_h%ham_complex,1,1,elsi_h%sc_desc)

            buffer_complex = elsi_h%ham_complex

            success = elpa_mult_ah_b_complex_double('U','U',elsi_h%n_g_size,elsi_h%n_g_size,&
                         elsi_h%ovlp_complex,elsi_h%n_l_rows,elsi_h%n_l_cols,buffer_complex,&
                         elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%n_b_rows,elsi_h%mpi_comm_row,&
                         elsi_h%mpi_comm_col,elsi_h%ham_complex,elsi_h%n_l_rows,elsi_h%n_l_cols)

            call pztranc(elsi_h%n_g_size,elsi_h%n_g_size,(1.0_r8,0.0_r8),elsi_h%ham_complex,&
                         1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),buffer_complex,1,1,elsi_h%sc_desc)

            ! Set the lower part from the upper
            do i_col = 1,elsi_h%n_g_size-1
               if(elsi_h%local_col(i_col) == 0) cycle
               do i_row = i_col+1,elsi_h%n_g_size
                  if(elsi_h%local_row(i_row) > 0) then
                     elsi_h%ham_complex(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                        buffer_complex(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
                  endif
               enddo
            enddo

            do i_col=1,elsi_h%n_g_size
               if(elsi_h%local_col(i_col) == 0 .or. elsi_h%local_row(i_col) == 0) cycle
               elsi_h%ham_complex(elsi_h%local_row(i_col),elsi_h%local_col(i_col)) = &
                  dble(elsi_h%ham_complex(elsi_h%local_row(i_col),elsi_h%local_col(i_col)))
            enddo
         endif

         call elsi_stop_transform_evp_time(elsi_h)

      case (REAL_VALUES)
         if(elsi_h%n_elsi_calls == 1) then
            if(.not. elsi_h%no_singularity_check) then
               call elsi_check_singularity(elsi_h)
            endif

            if(elsi_h%n_nonsingular == elsi_h%n_g_size) then ! Not singular
               call elsi_start_cholesky_time(elsi_h)

               elsi_h%overlap_is_singular = .false.

               ! Compute S = (U^T)U, U -> S
               success = elpa_cholesky_real_double(elsi_h%n_g_size,elsi_h%ovlp_real,&
                            elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,&
                            elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,.false.)
               if(.not. success) then
                  call elsi_stop(" Cholesky decomposition failed.",elsi_h,caller)
               endif

               ! compute U^-1 -> S
               success = elpa_invert_trm_real_double(elsi_h%n_g_size,elsi_h%ovlp_real,&
                            elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,&
                            elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,.false.)
               if(.not. success) then
                  call elsi_stop(" Matrix inversion failed.",elsi_h,caller)
               endif

               call elsi_stop_cholesky_time(elsi_h)
            endif
         endif ! n_elsi_calls == 1

         call elsi_start_transform_evp_time(elsi_h)

         call elsi_allocate(elsi_h,buffer_real,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)

         if(elsi_h%overlap_is_singular) then ! Use scaled eigenvectors
            ! buffer_real = H_real * S_real
            call pdgemm('N','N',elsi_h%n_g_size,elsi_h%n_nonsingular,elsi_h%n_g_size,&
                        1.0_r8,elsi_h%ham_real,1,1,elsi_h%sc_desc,elsi_h%ovlp_real,1,1,&
                        elsi_h%sc_desc,0.0_r8,buffer_real,1,1,elsi_h%sc_desc)

            ! H_real = (S_real)^T * buffer_real
            call pdgemm('T','N',elsi_h%n_nonsingular,elsi_h%n_nonsingular,elsi_h%n_g_size,&
                        1.0_r8,elsi_h%ovlp_real,1,1,elsi_h%sc_desc,buffer_real,1,1,elsi_h%sc_desc,&
                        0.0_r8,elsi_h%ham_real,1,1,elsi_h%sc_desc)

         else ! Use Cholesky
            success = elpa_mult_at_b_real_double('U','L',elsi_h%n_g_size,elsi_h%n_g_size,&
                         elsi_h%ovlp_real,elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%ham_real,&
                         elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%n_b_rows,elsi_h%mpi_comm_row,&
                         elsi_h%mpi_comm_col,buffer_real,elsi_h%n_l_rows,elsi_h%n_l_cols)

            call pdtran(elsi_h%n_g_size,elsi_h%n_g_size,1.0_r8,buffer_real,1,1,elsi_h%sc_desc,&
                        0.0_r8,elsi_h%ham_real,1,1,elsi_h%sc_desc)

            buffer_real = elsi_h%ham_real

            success = elpa_mult_at_b_real_double('U','U',elsi_h%n_g_size,elsi_h%n_g_size,&
                         elsi_h%ovlp_real,elsi_h%n_l_rows,elsi_h%n_l_cols,buffer_real,&
                         elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%n_b_rows,elsi_h%mpi_comm_row,&
                         elsi_h%mpi_comm_col,elsi_h%ham_real,elsi_h%n_l_rows,elsi_h%n_l_cols)

            call pdtran(elsi_h%n_g_size,elsi_h%n_g_size,1.0_r8,elsi_h%ham_real,1,1,&
                        elsi_h%sc_desc,0.0_r8,buffer_real,1,1,elsi_h%sc_desc)

            ! Set the lower part from the upper
            do i_col = 1,elsi_h%n_g_size-1
               if(elsi_h%local_col(i_col) == 0) cycle
               do i_row = i_col+1,elsi_h%n_g_size
                  if(elsi_h%local_row(i_row) > 0) then
                     elsi_h%ham_real(elsi_h%local_row(i_row),elsi_h%local_col(i_col)) = &
                        buffer_real(elsi_h%local_row(i_row),elsi_h%local_col(i_col))
                  endif
               enddo
            enddo
         endif

         call elsi_stop_transform_evp_time(elsi_h)

   end select

   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

end subroutine

!> 
!! This routine checks the singularity of overlap matrix by computing all
!! its eigenvalues.
!!
!! On exit, S is not modified if not singular, while overwritten by scaled
!! eigenvectors if singular, which is used to transform the generalized
!! eigenvalue problem to standard form without using Cholesky.
!!
subroutine elsi_check_singularity(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   real(kind=r8)                 :: ev_sqrt
   real(kind=r8),    allocatable :: ev_overlap(:)
   real(kind=r8),    allocatable :: buffer_real(:,:)
   complex(kind=r8), allocatable :: buffer_complex(:,:)
   integer(kind=i4)              :: i
   logical                       :: success

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_check_singularity"

   call elsi_start_singularity_check_time(elsi_h)

   select case (elsi_h%matrix_data_type)
      case (COMPLEX_VALUES)
         call elsi_allocate(elsi_h,buffer_complex,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)

         ! Use buffer_complex to store overlap matrix, otherwise it will
         ! be destroyed by eigenvalue calculation
         ! The nonsingular eigenvalues must be the first ones, so find
         ! eigenvalues of negative overlap matrix
         buffer_complex = -elsi_h%ovlp_complex

         call elsi_allocate(elsi_h,ev_overlap,elsi_h%n_g_size,"ev_overlap",caller)

         ! Use customized ELPA 2-stage solver to check overlap singularity
         ! Eigenvectors computed only for singular overlap matrix
         success = elpa_check_singularity_complex_double(elsi_h%n_g_size,elsi_h%n_g_size,&
                      buffer_complex,elsi_h%n_l_rows,ev_overlap,elsi_h%evec_complex,&
                      elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                      elsi_h%mpi_comm_col,elsi_h%mpi_comm,elsi_h%singularity_tolerance,&
                      elsi_h%n_nonsingular)
         if(.not. success) then
            call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                           " Exiting...",elsi_h,caller)
         endif

         ! Stop if n_states is larger that n_nonsingular
         if(elsi_h%n_nonsingular < elsi_h%n_states) then ! Too singular to continue
            call elsi_stop(" Overlap matrix is singular. The number of"//&
                           " basis functions after removing singularity"//&
                           " is smaller than the number of states. Try to"//&
                           " a) decrease the size of basis set, or b)"//&
                           " decrease the number of states, or c) increase"//&
                           " the tolerance of basis singularity."//&
                           " Exiting...",elsi_h,caller)
         elseif(elsi_h%n_nonsingular < elsi_h%n_g_size) then ! Singular
            elsi_h%overlap_is_singular = .true.

            if(elsi_h%stop_singularity) then
               call elsi_stop(" Overlap matrix is singular. This may mean"//&
                              " that a very large basis set is in use."//&
                              " Running with a near-singular basis set"//&
                              " may lead to completely wrong numerical"//&
                              " resutls. The calculation stops here,"//&
                              " because 'stop_singularity' is"//&
                              " set to .true. in elsi_customize."//&
                              " Exiting...",elsi_h,caller)
            endif

            call elsi_statement_print("  Overlap matrix is singular. This"//&
                                      " may mean that a very large basis"//&
                                      " set is in use. The calculation"//&
                                      " will continue. However, please"//&
                                      " note that running with a near-"//&
                                      "singular basis set may lead to"//&
                                      " completely wrong numerical results.",elsi_h)

            write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
               elsi_h%n_nonsingular
            call elsi_statement_print(info_str,elsi_h)

            call elsi_statement_print("  Using scaled eigenvectors of"//&
                                      " overlap matrix for transformation",elsi_h)

            ! Overlap matrix is overwritten with scaled eigenvectors
            do i = 1,elsi_h%n_nonsingular
               ev_sqrt = sqrt(ev_overlap(i))
               if(elsi_h%local_col(i) == 0) cycle
               elsi_h%ovlp_complex(:,elsi_h%local_col(i)) = &
                  elsi_h%evec_complex(:,elsi_h%local_col(i))/ev_sqrt
            enddo

         else ! Nonsingular
            elsi_h%overlap_is_singular = .false.
            call elsi_statement_print("  Overlap matrix is nonsingular",elsi_h)
         endif ! Singular overlap?

      case (REAL_VALUES)
         call elsi_allocate(elsi_h,buffer_real,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)

         ! Use buffer_real to store overlap matrix, otherwise it will be
         ! destroyed by eigenvalue calculation
         ! The nonsingular eigenvalues must be the first ones, so find
         ! eigenvalues of negative overlap matrix
         buffer_real = -elsi_h%ovlp_real

         call elsi_allocate(elsi_h,ev_overlap,elsi_h%n_g_size,"ev_overlap",caller)

         ! Use customized ELPA 2-stage solver to check overlap singularity
         ! Eigenvectors computed only for singular overlap matrix
         success = elpa_check_singularity_real_double(elsi_h%n_g_size,elsi_h%n_g_size,&
                      buffer_real,elsi_h%n_l_rows,ev_overlap,elsi_h%evec_real,elsi_h%n_l_rows,&
                      elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                      elsi_h%mpi_comm,elsi_h%singularity_tolerance,elsi_h%n_nonsingular)
         if(.not. success) then
            call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                           " Exiting...",elsi_h,caller)
         endif

         ! Stop if n_states is larger that n_nonsingular
         if(elsi_h%n_nonsingular < elsi_h%n_states) then ! Too singular to continue
            call elsi_stop(" Overlap matrix is singular. The number of"//&
                           " basis functions after removing singularity"//&
                           " is smaller than the number of states. Try to"//&
                           " a) decrease the size of basis set, or b)"//&
                           " decrease the number of states, or c) increase"//&
                           " the tolerance of basis singularity."//&
                           " Exiting...",elsi_h,caller)
         elseif(elsi_h%n_nonsingular < elsi_h%n_g_size) then ! Singular
            elsi_h%overlap_is_singular = .true.

            if(elsi_h%stop_singularity) then
               call elsi_stop(" Overlap matrix is singular. This may mean"//&
                              " that a very large basis set is in use."//&
                              " Running with a near-singular basis set"//&
                              " may lead to completely wrong numerical"//&
                              " resutls. The calculation stops here,"//&
                              " because 'stop_singularity' is"//&
                              " set to .true. in elsi_customize."//&
                              " Exiting...",elsi_h,caller)
            endif

            call elsi_statement_print("  Overlap matrix is singular. This"//&
                                      " may mean that a very large basis"//&
                                      " set is in use. The calculation"//&
                                      " will continue. However, please"//&
                                      " note that running with a near-"//&
                                      "singular basis set may lead to"//&
                                      " completely wrong numerical results.",elsi_h)

            write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
               elsi_h%n_nonsingular
            call elsi_statement_print(info_str,elsi_h)

            call elsi_statement_print("  Using scaled eigenvectors of"//&
                                      " overlap matrix for transformation",elsi_h)

            ! Overlap matrix is overwritten with scaled eigenvectors
            do i = 1,elsi_h%n_nonsingular
               ev_sqrt = sqrt(ev_overlap(i))
               if(elsi_h%local_col(i) == 0) cycle
               elsi_h%ovlp_real(:,elsi_h%local_col(i)) = &
                  elsi_h%evec_real(:,elsi_h%local_col(i))/ev_sqrt
            enddo

         else ! Nonsingular
            elsi_h%overlap_is_singular = .false.
            call elsi_statement_print("  Overlap matrix is nonsingular",elsi_h)
         endif ! Singular overlap?

   end select ! select matrix_data_type

   if(allocated(ev_overlap))     deallocate(ev_overlap)
   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

   call elsi_stop_singularity_check_time(elsi_h)

end subroutine

!> 
!! This routine does the back-transformation of the eigenvectors in standard
!! form (A'c' = c'v) to the original generalized form (Ac = Bcv)
!!
!! v = (U^-1)v'
!!
subroutine elsi_to_original_ev(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   logical                       :: success
   real(kind=r8),    allocatable :: buffer_real(:,:)
   complex(kind=r8), allocatable :: buffer_complex(:,:)

   character*40, parameter :: caller = "elsi_to_original_ev"

   call elsi_start_back_transform_ev_time(elsi_h)

   select case (elsi_h%matrix_data_type)
      case (COMPLEX_VALUES)
         call elsi_allocate(elsi_h,buffer_complex,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)
         buffer_complex = elsi_h%evec_complex

         if(elsi_h%overlap_is_singular) then
            ! Transform matrix is stored in S_complex after elsi_to_standard_evp
            call pzgemm('N','N',elsi_h%n_g_size,elsi_h%n_states,elsi_h%n_nonsingular,&
                        (1.0_r8,0.0_r8),elsi_h%ovlp_complex,1,1,elsi_h%sc_desc,buffer_complex,&
                        1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),elsi_h%evec_complex,1,1,elsi_h%sc_desc)
         else ! Nonsingular, use Cholesky
            ! (U^-1) is stored in S_complex after elsi_to_standard_evp
            ! C_complex = S_complex * C_complex = S_complex * buffer_complex
            call pztranc(elsi_h%n_g_size,elsi_h%n_g_size,(1.0_r8,0.0_r8),elsi_h%ovlp_complex,&
                         1,1,elsi_h%sc_desc,(0.0_r8,0.0_r8),elsi_h%ham_complex,1,1,elsi_h%sc_desc)

            success = elpa_mult_ah_b_complex_double('L','N',elsi_h%n_g_size,elsi_h%n_states,&
                         elsi_h%ham_complex,elsi_h%n_l_rows,elsi_h%n_l_cols,buffer_complex,&
                         elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%n_b_rows,elsi_h%mpi_comm_row,&
                         elsi_h%mpi_comm_col,elsi_h%evec_complex,elsi_h%n_l_rows,elsi_h%n_l_cols)
         endif

      case (REAL_VALUES)
         call elsi_allocate(elsi_h,buffer_real,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)
         buffer_real = elsi_h%evec_real

         if(elsi_h%overlap_is_singular) then
            ! Transform matrix is stored in S_real after elsi_to_standard_evp
            call pdgemm('N','N',elsi_h%n_g_size,elsi_h%n_states,elsi_h%n_nonsingular,1.0_r8,&
                        elsi_h%ovlp_real,1,1,elsi_h%sc_desc,buffer_real,1,1,elsi_h%sc_desc,0.0_r8,&
                        elsi_h%evec_real,1,1,elsi_h%sc_desc)
         else ! Nonsingular, use Cholesky
            ! (U^-1) is stored in S_real after elsi_to_standard_evp
            ! C_real = S_real * C_real = S_real * buffer_real
            call pdtran(elsi_h%n_g_size,elsi_h%n_g_size,1.0_r8,elsi_h%ovlp_real,1,1,elsi_h%sc_desc,&
                        0.0_r8,elsi_h%ham_real,1,1,elsi_h%sc_desc)

            success = elpa_mult_at_b_real_double('L','N',elsi_h%n_g_size,elsi_h%n_states,&
                         elsi_h%ham_real,elsi_h%n_l_rows,elsi_h%n_l_cols,buffer_real,elsi_h%n_l_rows,&
                         elsi_h%n_l_cols,elsi_h%n_b_rows,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                         elsi_h%evec_real,elsi_h%n_l_rows,elsi_h%n_l_cols)
         endif

   end select

   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

   call elsi_stop_back_transform_ev_time(elsi_h)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_evp_elpa(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   integer(kind=i4) :: mpierr
   logical :: success
   logical :: two_step_solver

   character*40, parameter :: caller = "elsi_solve_evp_elpa"

   elpa_print_times = .false.

   ! Choose 1-stage or 2-stage solver
   if(elsi_h%elpa_one_always) then
      two_step_solver = .false.
   elseif(elsi_h%elpa_two_always) then
      two_step_solver = .true.
   elseif(elsi_h%n_g_size < 256) then
      two_step_solver = .false.
   else
      two_step_solver = .true.
   endif

   ! Transform to standard form
   if(.not. elsi_h%overlap_is_unit) then
      call elsi_to_standard_evp(elsi_h)
   endif

   call elsi_start_standard_evp_time(elsi_h)

   ! Solve evp, return eigenvalues and eigenvectors
   if(two_step_solver) then ! 2-stage solver
      call elsi_statement_print("  Starting ELPA 2-stage eigensolver",elsi_h)
      select case (elsi_h%matrix_data_type)
         case (COMPLEX_VALUES)
            success = elpa_solve_evp_complex_2stage_double(elsi_h%n_nonsingular,elsi_h%n_states,&
                         elsi_h%ham_complex,elsi_h%n_l_rows,elsi_h%eval,elsi_h%evec_complex,&
                         elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                         elsi_h%mpi_comm_col,elsi_h%mpi_comm)
         case (REAL_VALUES)
            success = elpa_solve_evp_real_2stage_double(elsi_h%n_nonsingular,elsi_h%n_states,&
                         elsi_h%ham_real,elsi_h%n_l_rows,elsi_h%eval,elsi_h%evec_real,&
                         elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                         elsi_h%mpi_comm_col,elsi_h%mpi_comm)
      end select
   else ! 1-stage solver
      call elsi_statement_print("  Starting ELPA 1-stage eigensolver",elsi_h)
      select case (elsi_h%matrix_data_type)
         case (COMPLEX_VALUES)
            success = elpa_solve_evp_complex_1stage_double(elsi_h%n_nonsingular,elsi_h%n_states,&
                         elsi_h%ham_complex,elsi_h%n_l_rows,elsi_h%eval,elsi_h%evec_complex,&
                         elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,&
                         elsi_h%mpi_comm_col,elsi_h%mpi_comm)
         case (REAL_VALUES)
            success = elpa_solve_evp_real_1stage_double(elsi_h%n_nonsingular,elsi_h%n_states,&
                         elsi_h%ham_real,elsi_h%n_l_rows,elsi_h%eval,elsi_h%evec_real,elsi_h%n_l_rows,&
                         elsi_h%n_b_rows,elsi_h%n_l_cols,elsi_h%mpi_comm_row,elsi_h%mpi_comm_col,&
                         elsi_h%mpi_comm)
      end select
   endif

   call elsi_stop_standard_evp_time(elsi_h)

   if(.not. success) then
      call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                     " Exiting...",elsi_h,caller)
   endif

   ! Back-transform eigenvectors
   if(.not. elsi_h%overlap_is_unit) then
      call elsi_to_original_ev(elsi_h)
   endif

   call MPI_Barrier(elsi_h%mpi_comm,mpierr)

end subroutine

!> 
!! This routine transforms a generalized eigenvalue problem (Ac = Bcv)
!! to standard form (A'c' = c'v)
!!
!! Starting from Hv = eSv, we first perform a Cholesky decomposition of S
!! S = (U^T)U, resulting in Hv = e(U^T)Uv
!!
!! Using 1=U^-1U we define a new standard eigenvalue problem by
!! H(U^-1)(Uv) = e(U^T)(Uv) => ((U^-1)^T)H(U^-1)(Uv) = e(Uv)
!!
!! On exit, (U^-1) is stored in S, to be used for back-transformation
!!
subroutine elsi_to_standard_evp_sp(elsi_h)

   implicit none
   include 'mpif.h'

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   real(kind=r8),    allocatable :: buffer_real(:,:)
   complex(kind=r8), allocatable :: buffer_complex(:,:)
   logical                       :: success
   integer(kind=i4)              :: nblk=128
   integer(kind=i4)              :: n,nwork

   character*40, parameter :: caller = "elsi_to_standard_evp_sp"

   select case (elsi_h%matrix_data_type)
      case (COMPLEX_VALUES)
         if(elsi_h%n_elsi_calls == 1) then
            if(.not. elsi_h%no_singularity_check) then
               call elsi_check_singularity_sp(elsi_h)
            endif
         endif ! n_elsi_calls == 1

         if(elsi_h%n_nonsingular == elsi_h%n_g_size) then ! Not singular
            call elsi_start_cholesky_time(elsi_h)

            elsi_h%overlap_is_singular = .false.

            ! Compute S = (U^T)U, U -> S
            success = elpa_cholesky_complex_double(elsi_h%n_g_size,elsi_h%ovlp_complex,&
                         elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,mpi_comm_self,&
                         mpi_comm_self,.false.)
            if(.not. success) then
               call elsi_stop(" Cholesky decomposition failed.",elsi_h,caller)
            endif

            ! compute U^-1 -> S
            success = elpa_invert_trm_complex_double(elsi_h%n_g_size,elsi_h%ovlp_complex,&
                         elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,mpi_comm_self,&
                         mpi_comm_self,.false.)
            if(.not. success) then
               call elsi_stop(" Matrix inversion failed.",elsi_h,caller)
            endif

            call elsi_stop_cholesky_time(elsi_h)
         endif

         call elsi_start_transform_evp_time(elsi_h)

         call elsi_allocate(elsi_h,buffer_complex,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)

         if(elsi_h%overlap_is_singular) then ! Use scaled eigenvectors
            ! buffer_complex = H_complex * S_complex
            call zgemm('N','N',elsi_h%n_g_size,elsi_h%n_nonsingular,elsi_h%n_g_size,(1.0_r8,0.0_r8),&
                       elsi_h%ham_complex(1,1),elsi_h%n_g_size,elsi_h%ovlp_complex(1,1),elsi_h%n_g_size,&
                       (0.0_r8,0.0_r8),buffer_complex(1,1),elsi_h%n_g_size)

            ! H_complex = (S_complex)^* * buffer_complex
            call zgemm('C','N',elsi_h%n_nonsingular,elsi_h%n_nonsingular,elsi_h%n_g_size,(1.0_r8,0.0_r8),&
                       elsi_h%ovlp_complex(1,1),elsi_h%n_g_size,buffer_complex(1,1),elsi_h%n_g_size,&
                       (0.0_r8,0.0_r8),elsi_h%ham_complex(1,1),elsi_h%n_g_size)

         else ! Use cholesky
            ! buffer_complex = H_complex * S_complex
            call zgemm('N','N',elsi_h%n_g_size,elsi_h%n_g_size,elsi_h%n_g_size,(1.0_r8,0.0_r8),&
                       elsi_h%ham_complex(1,1),elsi_h%n_g_size,elsi_h%ovlp_complex(1,1),elsi_h%n_g_size,&
                       (0.0_r8,0.0_r8),buffer_complex(1,1),elsi_h%n_g_size)

            ! H_complex = (buffer_complex)^* * S_complex
            call zgemm('C','N',elsi_h%n_g_size,elsi_h%n_g_size,elsi_h%n_g_size,(1.0_r8,0.0_r8),&
                       buffer_complex(1,1),elsi_h%n_g_size,elsi_h%ovlp_complex(1,1),&
                       elsi_h%n_g_size,(0.0_r8,0.0_r8),elsi_h%ham_complex(1,1),elsi_h%n_g_size)
         endif

         call elsi_stop_transform_evp_time(elsi_h)

      case (REAL_VALUES)
         if(elsi_h%n_elsi_calls == 1) then
            if(.not. elsi_h%no_singularity_check) then
               call elsi_check_singularity_sp(elsi_h)
            endif
         endif ! n_elsi_calls == 1

         if(elsi_h%n_nonsingular == elsi_h%n_g_size) then ! Not singular
            call elsi_start_cholesky_time(elsi_h)

            elsi_h%overlap_is_singular = .false.

            ! Compute S = (U^T)U, U -> S
            success = elpa_cholesky_real_double(elsi_h%n_g_size,elsi_h%ovlp_real,elsi_h%n_l_rows,&
                         elsi_h%n_b_rows,elsi_h%n_l_cols,mpi_comm_self,mpi_comm_self,.false.)
            if(.not. success) then
               call elsi_stop(" Cholesky decomposition failed.",elsi_h,caller)
            endif

            ! compute U^-1 -> S
            success = elpa_invert_trm_real_double(elsi_h%n_g_size,elsi_h%ovlp_real,elsi_h%n_l_rows,&
                         elsi_h%n_b_rows,elsi_h%n_l_cols,mpi_comm_self,mpi_comm_self,.false.)
            if(.not. success) then
               call elsi_stop(" Matrix inversion failed.",elsi_h,caller)
            endif

            call elsi_stop_cholesky_time(elsi_h)
         endif

         call elsi_start_transform_evp_time(elsi_h)

         call elsi_allocate(elsi_h,buffer_real,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)

         if(elsi_h%overlap_is_singular) then ! Use scaled eigenvectors
            ! buffer_real = H_real * S_real
            call dgemm('N','N',elsi_h%n_g_size,elsi_h%n_nonsingular,elsi_h%n_g_size,1.0_r8,&
                       elsi_h%ham_real(1,1),elsi_h%n_g_size,elsi_h%ovlp_real(1,1),elsi_h%n_g_size,&
                       0.0_r8,buffer_real(1,1),elsi_h%n_g_size)

            ! H_real = (S_real)^T * buffer_real
            call dgemm('T','N',elsi_h%n_nonsingular,elsi_h%n_nonsingular,elsi_h%n_g_size,1.0_r8,&
                       elsi_h%ovlp_real(1,1),elsi_h%n_g_size,buffer_real(1,1),elsi_h%n_g_size,0.0_r8,&
                       elsi_h%ham_real(1,1),elsi_h%n_g_size)

         else ! Use Cholesky
           ! buffer_real = H_real * S_real
           do n = 1,elsi_h%n_g_size,nblk
              nwork = nblk

              if(n+nwork-1 > elsi_h%n_g_size) nwork = elsi_h%n_g_size-n+1

              call dgemm('N','N',n+nwork-1,nwork,n+nwork-1,1.0_r8,elsi_h%ham_real(1,1),&
                         elsi_h%n_g_size,elsi_h%ovlp_real(1,n),elsi_h%n_g_size,0.0_r8,&
                         buffer_real(1,n),elsi_h%n_g_size) 
           enddo

           ! H_real = (buffer_real)*T * S_real
           do n = 1,elsi_h%n_g_size,nblk

              nwork = nblk

              if(n+nwork-1 > elsi_h%n_g_size) nwork = elsi_h%n_g_size-n+1

              call dgemm('T','N',nwork,elsi_h%n_g_size-n+1,n+nwork-1,1.0_r8,elsi_h%ovlp_real(1,n),&
                         elsi_h%n_g_size,buffer_real(1,n),elsi_h%n_g_size,0.0_r8,elsi_h%ham_real(n,n),&
                         elsi_h%n_g_size)
           enddo
        endif

        call elsi_stop_transform_evp_time(elsi_h)

   end select

   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

end subroutine

subroutine elsi_to_original_ev_sp(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   real(kind=r8),    allocatable :: buffer_real(:,:)
   complex(kind=r8), allocatable :: buffer_complex(:,:)

   character*40, parameter :: caller = "elsi_to_original_ev_sp"

   call elsi_start_back_transform_ev_time(elsi_h)

   select case (elsi_h%matrix_data_type)
      case (COMPLEX_VALUES)
         call elsi_allocate(elsi_h,buffer_complex,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)
         buffer_complex = elsi_h%evec_complex

         if(elsi_h%overlap_is_singular) then
            ! Transform matrix is stored in S_complex after elsi_to_standard_evp
            call zgemm('N','N',elsi_h%n_g_size,elsi_h%n_states,elsi_h%n_nonsingular,(1.0_r8,0.0_r8),&
                       elsi_h%ovlp_complex(1,1),elsi_h%n_g_size,buffer_complex(1,1),elsi_h%n_g_size,&
                       (0.0_r8,0.0_r8),elsi_h%evec_complex(1,1),elsi_h%n_g_size)
         else ! Nonsingular, use Cholesky
            ! (U^-1) is stored in S_complex after elsi_to_standard_evp
            ! C_complex = S_complex * C_complex = S_complex * buffer_complex
            call zgemm('N','N',elsi_h%n_g_size,elsi_h%n_states,elsi_h%n_g_size,(1.0_r8,0.0_r8),&
                       elsi_h%ovlp_complex(1,1),elsi_h%n_g_size,buffer_complex(1,1),elsi_h%n_g_size,&
                       (0.0_r8,0.0_r8),elsi_h%evec_complex(1,1),elsi_h%n_g_size)
         endif

      case (REAL_VALUES)
         call elsi_allocate(elsi_h,buffer_real,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)
         buffer_real = elsi_h%evec_real

         if(elsi_h%overlap_is_singular) then
            ! Transform matrix is stored in S_real after elsi_to_standard_evp
            call dgemm('N','N',elsi_h%n_g_size,elsi_h%n_states,elsi_h%n_nonsingular,1.0_r8,&
                       elsi_h%ovlp_real(1,1),elsi_h%n_g_size,buffer_real(1,1),elsi_h%n_g_size,&
                       0.0_r8,elsi_h%evec_real(1,1),elsi_h%n_g_size)
         else ! Nonsingular, use Cholesky
            ! (U^-1) is stored in S_real after elsi_to_standard_evp
            ! C_real = S_real * C_real = S_real * buffer_real
            call dgemm('N','N',elsi_h%n_g_size,elsi_h%n_states,elsi_h%n_g_size,1.0_r8,&
                       elsi_h%ovlp_real(1,1),elsi_h%n_g_size,buffer_real(1,1),elsi_h%n_g_size,&
                       0.0_r8,elsi_h%evec_real(1,1),elsi_h%n_g_size)
         endif

   end select

   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

   call elsi_stop_back_transform_ev_time(elsi_h)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_evp_elpa_sp(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   real(kind=r8),    allocatable :: d(:),e(:)
   real(kind=r8),    allocatable :: tau_real(:),buffer_real(:,:)
   complex(kind=r8), allocatable :: tau_complex(:),buffer_complex(:,:)
   logical                       :: success
   integer(kind=i4)              :: info

   character*40, parameter :: caller = "elsi_solve_evp_elpa_sp"

   call elsi_allocate(elsi_h,d,elsi_h%n_g_size,"d",caller)
   call elsi_allocate(elsi_h,e,elsi_h%n_g_size,"e",caller)

   ! Transform to standard form
   if(.not. elsi_h%overlap_is_unit) then
      call elsi_to_standard_evp_sp(elsi_h)
   endif

   ! Solve evp, return eigenvalues and eigenvectors
   call elsi_statement_print("  Starting ELPA eigensolver",elsi_h)

   call elsi_start_standard_evp_time(elsi_h)

   select case (elsi_h%matrix_data_type)
      case (COMPLEX_VALUES)
         call elsi_allocate(elsi_h,tau_complex,elsi_h%n_g_size,"tau_complex",caller)
         call elsi_allocate(elsi_h,buffer_real,elsi_h%n_g_size,elsi_h%n_g_size,"buffer_real",caller)
         call elsi_allocate(elsi_h,buffer_complex,elsi_h%n_g_size,elsi_h%n_g_size,"buffer_complex",caller)

         call zhetrd('U',elsi_h%n_g_size,elsi_h%ham_complex,elsi_h%n_g_size,d,e,tau_complex,&
                     buffer_complex,size(buffer_complex),info)
 
         success = elpa_solve_tridi_double(elsi_h%n_g_size,elsi_h%n_states,d,e,buffer_real,&
                                           elsi_h%n_g_size,64,elsi_h%n_g_size,mpi_comm_self,&
                                           mpi_comm_self,.false.)

         elsi_h%eval(1:elsi_h%n_states) = d(1:elsi_h%n_states)
         elsi_h%evec_complex(1:elsi_h%n_g_size,1:elsi_h%n_states) = &
            buffer_real(1:elsi_h%n_g_size,1:elsi_h%n_states)

         call zunmtr('L','U','N',elsi_h%n_g_size,elsi_h%n_states,elsi_h%ham_complex,elsi_h%n_g_size,&
                     tau_complex,elsi_h%evec_complex,elsi_h%n_g_size,buffer_complex,size(buffer_complex),info)

         deallocate(tau_complex)
         deallocate(buffer_real)
         deallocate(buffer_complex)

      case (REAL_VALUES)
         call elsi_allocate(elsi_h,tau_real,elsi_h%n_g_size,"tau_real",caller)
         call elsi_allocate(elsi_h,buffer_real,elsi_h%n_g_size,elsi_h%n_g_size,"buffer_real",caller)

         call dsytrd('U',elsi_h%n_g_size,elsi_h%ham_real,elsi_h%n_g_size,d,e,tau_real,buffer_real,&
                     size(buffer_real),info)

         success = elpa_solve_tridi_double(elsi_h%n_g_size,elsi_h%n_states,d,e,buffer_real,&
                                           elsi_h%n_g_size,64,elsi_h%n_g_size,mpi_comm_self,&
                                           mpi_comm_self,.false.)
 
         elsi_h%eval(1:elsi_h%n_states) = d(1:elsi_h%n_states)
         elsi_h%evec_real(1:elsi_h%n_g_size,1:elsi_h%n_states) = &
            buffer_real(1:elsi_h%n_g_size,1:elsi_h%n_states)
 
         call dormtr('L','U','N',elsi_h%n_g_size,elsi_h%n_states,elsi_h%ham_real,elsi_h%n_g_size,&
                     tau_real,elsi_h%evec_real,elsi_h%n_g_size,buffer_real,size(buffer_real),info)

         deallocate(tau_real)
         deallocate(buffer_real)

   end select

   call elsi_stop_standard_evp_time(elsi_h)

   deallocate(d)
   deallocate(e)

   if(.not. success) then
      call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                     " Exiting...",elsi_h,caller)
   endif

   ! Back-transform eigenvectors
   if(.not. elsi_h%overlap_is_unit) then
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
   include 'mpif.h'   

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   real(kind=r8)                 :: ev_sqrt
   real(kind=r8),    allocatable :: ev_overlap(:)
   real(kind=r8),    allocatable :: buffer_real(:,:)
   complex(kind=r8), allocatable :: buffer_complex(:,:)
   integer(kind=i4)              :: i,i_col
   logical                       :: success

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_check_singularity_sp"

   call elsi_start_singularity_check_time(elsi_h)

   select case (elsi_h%matrix_data_type)
      case (COMPLEX_VALUES)
         call elsi_allocate(elsi_h,buffer_complex,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)

         ! Use buffer_complex to store overlap matrix, otherwise it will
         ! be destroyed by eigenvalue calculation
         ! The nonsingular eigenvalues must be the first ones, so find
         ! eigenvalues of negative overlap matrix
         buffer_complex = -elsi_h%ovlp_complex

         call elsi_allocate(elsi_h,ev_overlap,elsi_h%n_g_size,"ev_overlap",caller)

         ! Use customized ELPA 2-stage solver to check overlap singularity
         ! Eigenvectors computed only for singular overlap matrix
         success = elpa_check_singularity_complex_double(elsi_h%n_g_size,elsi_h%n_g_size,&
                      buffer_complex,elsi_h%n_l_rows,ev_overlap,elsi_h%evec_complex,&
                      elsi_h%n_l_rows,elsi_h%n_b_rows,elsi_h%n_l_cols,mpi_comm_self,mpi_comm_self,&
                      mpi_comm_self,elsi_h%singularity_tolerance,elsi_h%n_nonsingular)
         if(.not. success) then
            call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                           " Exiting...",elsi_h,caller)
         endif

         ! Stop if n_states is larger that n_nonsingular
         if(elsi_h%n_nonsingular < elsi_h%n_states) then ! Too singular to continue
            call elsi_stop(" Overlap matrix is singular. The number of"//&
                           " basis functions after removing singularity"//&
                           " is smaller than the number of states. Try to"//&
                           " a) decrease the size of basis set, or b)"//&
                           " decrease the number of states, or c) increase"//&
                           " the tolerance of basis singularity."//&
                           " Exiting...",elsi_h,caller)
         elseif(elsi_h%n_nonsingular < elsi_h%n_g_size) then ! Singular
            elsi_h%overlap_is_singular = .true.

            if(elsi_h%stop_singularity) then
               call elsi_stop(" Overlap matrix is singular. This may mean"//&
                              " that a very large basis set is in use."//&
                              " Running with a near-singular basis set"//&
                              " may lead to completely wrong numerical"//&
                              " resutls. The calculation stops here,"//&
                              " because 'stop_singularity' is"//&
                              " set to .true. in elsi_customize."//&
                              " Exiting...",elsi_h,caller)
            endif

            call elsi_statement_print("  Overlap matrix is singular. This"//&
                                      " may mean that a very large basis"//&
                                      " set is in use. The calculation"//&
                                      " will continue. However, please"//&
                                      " note that running with a near-"//&
                                      "singular basis set may lead to"//&
                                      " completely wrong numerical results.",elsi_h)

            write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
               elsi_h%n_nonsingular
            call elsi_statement_print(info_str,elsi_h)

            call elsi_statement_print("  Using scaled eigenvectors of"//&
                                      " overlap matrix for transformation",elsi_h)

            ! Overlap matrix is overwritten with scaled eigenvectors
            do i = 1,elsi_h%n_nonsingular
               ev_sqrt = sqrt(ev_overlap(i))
               elsi_h%ovlp_complex(:,i) = elsi_h%evec_complex(:,i)/ev_sqrt
            enddo

         else ! Nonsingular
            elsi_h%overlap_is_singular = .false.
            call elsi_statement_print("  Overlap matrix is nonsingular",elsi_h)
         endif ! Singular overlap?

      case (REAL_VALUES)
         call elsi_allocate(elsi_h,buffer_real,elsi_h%n_l_rows,elsi_h%n_l_cols,"temp",caller)

         ! Use buffer_real to store overlap matrix, otherwise it will be
         ! destroyed by eigenvalue calculation
         ! The nonsingular eigenvalues must be the first ones, so find
         ! eigenvalues of negative overlap matrix
         buffer_real = -elsi_h%ovlp_real

         call elsi_allocate(elsi_h,ev_overlap,elsi_h%n_g_size,"ev_overlap",caller)

         ! Use customized ELPA 2-stage solver to check overlap singularity
         ! Eigenvectors computed only for singular overlap matrix
         success = elpa_check_singularity_real_double(elsi_h%n_g_size,elsi_h%n_g_size,&
                      buffer_real,elsi_h%n_l_rows,ev_overlap,elsi_h%evec_real,elsi_h%n_l_rows,&
                      elsi_h%n_b_rows,elsi_h%n_l_cols,mpi_comm_self,mpi_comm_self,&
                      mpi_comm_self,elsi_h%singularity_tolerance,elsi_h%n_nonsingular)
         if(.not. success) then
            call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                           " Exiting...",elsi_h,caller)
         endif

         ! Stop if n_states is larger that n_nonsingular
         if(elsi_h%n_nonsingular < elsi_h%n_states) then ! Too singular to continue
            call elsi_stop(" Overlap matrix is singular. The number of"//&
                           " basis functions after removing singularity"//&
                           " is smaller than the number of states. Try to"//&
                           " a) decrease the size of basis set, or b)"//&
                           " decrease the number of states, or c) increase"//&
                           " the tolerance of basis singularity."//&
                           " Exiting...",elsi_h,caller)
         elseif(elsi_h%n_nonsingular < elsi_h%n_g_size) then ! Singular
            elsi_h%overlap_is_singular = .true.

            if(elsi_h%stop_singularity) then
               call elsi_stop(" Overlap matrix is singular. This may mean"//&
                              " that a very large basis set is in use."//&
                              " Running with a near-singular basis set"//&
                              " may lead to completely wrong numerical"//&
                              " resutls. The calculation stops here,"//&
                              " because 'stop_singularity' is"//&
                              " set to .true. in elsi_customize."//&
                              " Exiting...",elsi_h,caller)
            endif

            call elsi_statement_print("  Overlap matrix is singular. This"//&
                                      " may mean that a very large basis"//&
                                      " set is in use. The calculation"//&
                                      " will continue. However, please"//&
                                      " note that running with a near-"//&
                                      "singular basis set may lead to"//&
                                      " completely wrong numerical results.",elsi_h)

            write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
               elsi_h%n_nonsingular
            call elsi_statement_print(info_str,elsi_h)

            call elsi_statement_print("  Using scaled eigenvectors of"//&
                                      " overlap matrix for transformation",elsi_h)

            ! Overlap matrix is overwritten with scaled eigenvectors
            do i = 1,elsi_h%n_nonsingular
               ev_sqrt = sqrt(ev_overlap(i))
               elsi_h%ovlp_real(:,i) = elsi_h%evec_real(:,i)/ev_sqrt
            enddo

         else ! Nonsingular
            elsi_h%overlap_is_singular = .false.
            call elsi_statement_print("  Overlap matrix is nonsingular",elsi_h)
         endif ! Singular overlap?

   end select ! select matrix_data_type

   if(allocated(ev_overlap))     deallocate(ev_overlap)
   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

   call elsi_stop_singularity_check_time(elsi_h)

end subroutine

end module ELSI_ELPA
