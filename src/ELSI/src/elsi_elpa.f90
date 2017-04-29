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
   use ELSI_DIMENSIONS
   use ELSI_TIMERS
   use ELSI_UTILS
   use ELSI_MU
   use ELPA1
   use ELPA2

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
subroutine elsi_get_elpa_comms(mpi_comm_row_out,mpi_comm_col_out)

   implicit none

   integer, intent(out) :: mpi_comm_row_out   !< Row MPI communicator
   integer, intent(out) :: mpi_comm_col_out   !< Column MPI communicator

   integer :: success

   success = elpa_get_communicators(mpi_comm_global,my_p_row,my_p_col,&
                                    mpi_comm_row_out,mpi_comm_col_out)

end subroutine

!>
!! This routine computes the chemical potential and occupation numbers
!! from eigenvalues.
!!
subroutine elsi_compute_occ_elpa()

   implicit none

   real*8 :: mu           !< Chemical potential
   real*8 :: k_weights(1) !< Weights of k-points

   real*8, allocatable :: eval_aux(:,:,:)
   real*8, allocatable :: occ_aux(:,:,:)

   !< Currently this subroutine only supports 1 spin channel and 1 k-point
   integer, parameter :: n_spin   = 1
   integer, parameter :: n_kpoint = 1

   character*40, parameter :: caller = "elsi_compute_occ_elpa"

   k_weights(1) = 1.0d0

   if(.not.allocated(occ_elpa)) then
       call elsi_allocate(occ_elpa,n_states,"occ_elpa",caller)
   endif

   call elsi_allocate(eval_aux,n_states,1,1,"eval_aux",caller)
   call elsi_allocate(occ_aux,n_states,1,1,"occ_aux",caller)

   eval_aux(1:n_states,1,1) = eval(1:n_states)
   occ_aux = 0.0d0

   call elsi_compute_mu_and_occ(n_electrons,n_states,n_spin,n_kpoint,&
                                k_weights,eval_aux,occ_aux,mu)

   occ_elpa(:) = occ_aux(:,1,1)

   deallocate(eval_aux)
   deallocate(occ_aux)

end subroutine

!>
!! This routine constructs the density matrix using eigenvectors from ELPA.
!!
subroutine elsi_compute_dm_elpa()

   implicit none

   real*8, allocatable     :: tmp_real(:,:)    !< Real eigenvectors, temporary
   complex*16, allocatable :: tmp_complex(:,:) !< Complex eigenvectors, temporary
   real*8, allocatable :: factor(:) !< Factor to construct density matrix
   integer :: i,i_col,i_row

   character*40, parameter :: caller = "elsi_compute_dm"

   call elsi_start_density_matrix_time()

   select case (mode)
      case (REAL_VALUES)
         ! Get eigenvectors into tmp_real
         call elsi_allocate(tmp_real,n_l_rows,n_l_cols,"tmp_real",caller)
         tmp_real = evec_real

         ! Compute the factors used to construct density matrix
         call elsi_allocate(factor,n_states,"factor",caller)
         factor = 0.0d0

         do i = 1,n_states
            if(occ_elpa(i) > 0.0d0) then
               factor(i) = sqrt(occ_elpa(i))
            endif
         enddo

         do i = 1,n_states
            if(factor(i) > 0.0d0) then
               if(local_col(i) > 0) then
                  tmp_real(:,local_col(i)) = tmp_real(:,local_col(i))*factor(i)
               endif
            elseif(local_col(i) .ne. 0) then
               tmp_real(:,local_col(i)) = 0.0d0
            endif
         enddo

         ! Compute density matrix
         den_mat = 0.0d0

         ! D_elpa = tmp_real*tmp_real^T
         call pdsyrk('U','N',n_g_size,n_states,1.0d0,tmp_real,1,1,sc_desc,&
                     0.0d0,den_mat,1,1,sc_desc)

      case (COMPLEX_VALUES)
         ! Get eigenvectors into tmp_complex
         call elsi_allocate(tmp_complex,n_l_rows,n_l_cols,"tmp_complex",caller)
         tmp_complex = evec_complex

         ! Compute the factors used to construct density matrix
         call elsi_allocate(factor,n_states,"factor",caller)
         factor = 0.0d0

         do i = 1,n_states
            if(occ_elpa(i) > 0.0d0) then
               factor(i) = sqrt(occ_elpa(i))
            endif
         enddo

         do i = 1,n_states
            if(factor(i) > 0.0d0) then
               if(local_col(i) > 0) then
                  tmp_complex(:,local_col(i)) = tmp_complex(:,local_col(i))*factor(i)
               endif
            elseif(local_col(i) .ne. 0) then
               tmp_complex(:,local_col(i)) = (0.0d0,0.0d0)
            endif
         enddo

         ! Compute density matrix
         den_mat = 0.0d0
         call elsi_allocate(tmp_real,n_l_rows,n_l_cols,"tmp_real",caller)

         call pdsyrk('U','N',n_g_size,n_states,1.0d0,real(tmp_complex),1,1,sc_desc,&
                     0.0d0,den_mat,1,1,sc_desc)
         call pdsyrk('U','N',n_g_size,n_states,1.0d0,aimag(tmp_complex),1,1,sc_desc,&
                     0.0d0,tmp_real,1,1,sc_desc)

         den_mat = den_mat+tmp_real
   end select

   deallocate(factor)
   if(allocated(tmp_real))    deallocate(tmp_real)
   if(allocated(tmp_complex)) deallocate(tmp_complex)

   ! Set upper triangle matrix den_mat to full form
   call elsi_allocate(tmp_real,n_l_rows,n_l_cols,"tmp_real",caller)

   call pdtran(n_g_size,n_g_size,1.0d0,den_mat,1,1,sc_desc,0.0d0,tmp_real,1,1,sc_desc)

   do i_col = 1,n_g_size-1
      if(local_col(i_col) == 0) cycle
      do i_row = i_col+1,n_g_size
         if(local_row(i_row) > 0) then
            den_mat(local_row(i_row),local_col(i_col)) = &
               tmp_real(local_row(i_row),local_col(i_col))
         endif
      enddo
   enddo

   deallocate(tmp_real)

   call elsi_stop_density_matrix_time()

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
subroutine elsi_to_standard_evp()

   implicit none

   integer :: i_row,i_col
   real*8, allocatable :: buffer_real(:,:)
   complex*16, allocatable :: buffer_complex(:,:)
   logical :: success

   character*40, parameter :: caller = "elsi_to_standard_evp"

   select case (mode)
      case (COMPLEX_VALUES)
         if(n_elsi_calls == 1) then
            if(.not.no_singularity_check) then
               call elsi_check_singularity()
            endif

            if(n_nonsingular == n_g_size) then ! Not singular
               call elsi_start_cholesky_time()

               overlap_is_singular = .false.

               ! Compute S = (U^T)U, U -> S
               success = elpa_cholesky_complex_double(n_g_size,ovlp_complex,n_l_rows,&
                            n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
               if(.not.success) then
                  call elsi_stop(" Cholesky decomposition failed.",caller)
               endif

               ! compute U^-1 -> S
               success = elpa_invert_trm_complex_double(n_g_size,ovlp_complex,n_l_rows,&
                            n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
               if(.not.success) then
                  call elsi_stop(" Matrix inversion failed.",caller)
               endif

               call elsi_stop_cholesky_time()
            endif
         endif ! n_elsi_calls == 1

         call elsi_start_transform_evp_time()

         call elsi_allocate(buffer_complex,n_l_rows,n_l_cols,"temp",caller)

         if(overlap_is_singular) then ! Use scaled eigenvectors
            ! buffer_complex = H_complex * S_complex
            call pzgemm('N','N',n_g_size,n_nonsingular,n_g_size,(1.0d0,0.0d0),&
                        ham_complex,1,1,sc_desc,ovlp_complex,1,1,sc_desc,&
                        (0.0d0,0.0d0),buffer_complex,1,1,sc_desc)

            ! H_complex = (S_complex)^* * buffer_complex
            call pzgemm('C','N',n_nonsingular,n_nonsingular,n_g_size,(1.0d0,0.0d0),&
                        ovlp_complex,1,1,sc_desc,buffer_complex,1,1,sc_desc,&
                        (0.0d0,0.0d0),ham_complex,1,1,sc_desc)

         else ! Use cholesky
            success = elpa_mult_ah_b_complex_double('U','L',n_g_size,n_g_size,&
                         ovlp_complex,n_l_rows,n_l_cols,ham_complex,n_l_rows,&
                         n_l_cols,n_b_rows,mpi_comm_row,mpi_comm_col,&
                         buffer_complex,n_l_rows,n_l_cols)

            call pztranc(n_g_size,n_g_size,(1.0d0,0.0d0),buffer_complex,1,1,sc_desc,&
                         (0.0d0,0.0d0),ham_complex,1,1,sc_desc)

            buffer_complex = ham_complex

            success = elpa_mult_ah_b_complex_double('U','U',n_g_size,n_g_size,&
                         ovlp_complex,n_l_rows,n_l_cols,buffer_complex,n_l_rows,&
                         n_l_cols,n_b_rows,mpi_comm_row,mpi_comm_col,ham_complex,&
                         n_l_rows,n_l_cols)

            call pztranc(n_g_size,n_g_size,(1.0d0,0.0d0),ham_complex,1,1,sc_desc,&
                         (0.0d0,0.0d0),buffer_complex,1,1,sc_desc)

            ! Set the lower part from the upper
            do i_col = 1,n_g_size-1
               if(local_col(i_col) == 0) cycle
               do i_row = i_col+1,n_g_size
                  if(local_row(i_row) > 0) then
                     ham_complex(local_row(i_row),local_col(i_col)) = &
                        buffer_complex(local_row(i_row),local_col(i_col))
                  endif
               enddo
            enddo

            do i_col=1,n_g_size
               if(local_col(i_col) == 0 .or. local_row(i_col) == 0) cycle
               ham_complex(local_row(i_col),local_col(i_col)) = &
                  dble(ham_complex(local_row(i_col),local_col(i_col)))
            enddo
         endif

         call elsi_stop_transform_evp_time()

      case (REAL_VALUES)
         if(n_elsi_calls == 1) then
            if(.not.no_singularity_check) then
               call elsi_check_singularity()
            endif

            if(n_nonsingular == n_g_size) then ! Not singular
               call elsi_start_cholesky_time()

               overlap_is_singular = .false.

               ! Compute S = (U^T)U, U -> S
               success = elpa_cholesky_real_double(n_g_size,ovlp_real,n_l_rows,&
                            n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
               if(.not.success) then
                  call elsi_stop(" Cholesky decomposition failed.",caller)
               endif

               ! compute U^-1 -> S
               success = elpa_invert_trm_real_double(n_g_size,ovlp_real,n_l_rows,&
                            n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
               if(.not.success) then
                  call elsi_stop(" Matrix inversion failed.",caller)
               endif

               call elsi_stop_cholesky_time()
            endif
         endif ! n_elsi_calls == 1

         call elsi_start_transform_evp_time()

         call elsi_allocate(buffer_real,n_l_rows,n_l_cols,"temp",caller)

         if(overlap_is_singular) then ! Use scaled eigenvectors
            ! buffer_real = H_real * S_real
            call pdgemm('N','N',n_g_size,n_nonsingular,n_g_size,1.0d0,ham_real,1,&
                        1,sc_desc,ovlp_real,1,1,sc_desc,0.0d0,buffer_real,1,1,sc_desc)

            ! H_real = (S_real)^T * buffer_real
            call pdgemm('T','N',n_nonsingular,n_nonsingular,n_g_size,1.0d0,ovlp_real,&
                        1,1,sc_desc,buffer_real,1,1,sc_desc,0.0d0,ham_real,1,1,sc_desc)

         else ! Use Cholesky
            success = elpa_mult_at_b_real_double('U','L',n_g_size,n_g_size,&
                         ovlp_real,n_l_rows,n_l_cols,ham_real,n_l_rows,n_l_cols,&
                         n_b_rows,mpi_comm_row,mpi_comm_col,buffer_real,&
                         n_l_rows,n_l_cols)

            call pdtran(n_g_size,n_g_size,1.0d0,buffer_real,1,1,sc_desc,0.0d0,&
                        ham_real,1,1,sc_desc)

            buffer_real = ham_real

            success = elpa_mult_at_b_real_double('U','U',n_g_size,n_g_size,&
                         ovlp_real,n_l_rows,n_l_cols,buffer_real,n_l_rows,n_l_cols,&
                         n_b_rows,mpi_comm_row,mpi_comm_col,ham_real,n_l_rows,&
                         n_l_cols)

            call pdtran(n_g_size,n_g_size,1.0d0,ham_real,1,1,sc_desc,0.0d0,buffer_real,&
                        1,1,sc_desc)

            ! Set the lower part from the upper
            do i_col = 1,n_g_size-1
               if(local_col(i_col) == 0) cycle
               do i_row = i_col+1,n_g_size
                  if(local_row(i_row) > 0) then
                     ham_real(local_row(i_row),local_col(i_col)) = &
                        buffer_real(local_row(i_row),local_col(i_col))
                  endif
               enddo
            enddo
         endif

         call elsi_stop_transform_evp_time()

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
subroutine elsi_check_singularity()


   real*8 :: ev_sqrt
   real*8, allocatable :: ev_overlap(:)
   real*8, allocatable :: buffer_real(:,:)
   complex*16, allocatable :: buffer_complex(:,:)
   integer :: i
   logical :: success

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_check_singularity"

   call elsi_start_singularity_check_time()

   select case (mode)
      case (COMPLEX_VALUES)
         call elsi_allocate(buffer_complex,n_l_rows,n_l_cols,"temp",caller)

         ! Use buffer_complex to store overlap matrix, otherwise it will
         ! be destroyed by eigenvalue calculation
         ! The nonsingular eigenvalues must be the first ones, so find
         ! eigenvalues of negative overlap matrix
         buffer_complex = -ovlp_complex

         call elsi_allocate(ev_overlap,n_g_size,"ev_overlap",caller)

         ! Use customized ELPA 2-stage solver to check overlap singularity
         ! Eigenvectors computed only for singular overlap matrix
         success = elpa_check_singularity_complex_double(n_g_size,n_g_size,&
                      buffer_complex,n_l_rows,ev_overlap,evec_complex,n_l_rows,&
                      n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,&
                      mpi_comm_global,singularity_tolerance,n_nonsingular)
         if(.not.success) then
            call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                           " Exiting...",caller)
         endif

         ! Stop if n_states is larger that n_nonsingular
         if(n_nonsingular < n_states) then ! Too singular to continue
            call elsi_stop(" Overlap matrix is singular. The number of"//&
                           " basis functions after removing singularity"//&
                           " is smaller than the number of states. Try to"//&
                           " a) decrease the size of basis set, or b)"//&
                           " decrease the number of states, or c) increase"//&
                           " the tolerance of basis singularity."//&
                           " Exiting...",caller)
         elseif(n_nonsingular < n_g_size) then ! Singular
            overlap_is_singular = .true.

            if(stop_singularity) then
               call elsi_stop(" Overlap matrix is singular. This may mean"//&
                              " that a very large basis set is in use."//&
                              " Running with a near-singular basis set"//&
                              " may lead to completely wrong numerical"//&
                              " resutls. The calculation stops here,"//&
                              " because 'force_stop_singularity' is"//&
                              " set to .true. in elsi_customize."//&
                              " Exiting...",caller)
            endif

            call elsi_statement_print("  Overlap matrix is singular. This"//&
                                      " may mean that a very large basis"//&
                                      " set is in use. The calculation"//&
                                      " will continue. However, please"//&
                                      " note that running with a near-"//&
                                      "singular basis set may lead to"//&
                                      " completely wrong numerical results.")

            write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
               n_nonsingular
            call elsi_statement_print(info_str)

            call elsi_statement_print("  Using scaled eigenvectors of"//&
                                      " overlap matrix for transformation")

            ! Overlap matrix is overwritten with scaled eigenvectors
            do i = 1,n_nonsingular
               ev_sqrt = sqrt(ev_overlap(i))
               if(local_col(i) == 0) cycle
               ovlp_complex(:,local_col(i)) = evec_complex(:,local_col(i))/ev_sqrt
            enddo

         else ! Nonsingular
            overlap_is_singular = .false.
            call elsi_statement_print("  Overlap matrix is nonsingular")
         endif ! Singular overlap?

      case (REAL_VALUES)
         call elsi_allocate(buffer_real,n_l_rows,n_l_cols,"temp",caller)

         ! Use buffer_real to store overlap matrix, otherwise it will be
         ! destroyed by eigenvalue calculation
         ! The nonsingular eigenvalues must be the first ones, so find
         ! eigenvalues of negative overlap matrix
         buffer_real = -ovlp_real

         call elsi_allocate(ev_overlap,n_g_size,"ev_overlap",caller)

         ! Use customized ELPA 2-stage solver to check overlap singularity
         ! Eigenvectors computed only for singular overlap matrix
         success = elpa_check_singularity_real_double(n_g_size,n_g_size,&
                      buffer_real,n_l_rows,ev_overlap,evec_real,n_l_rows,&
                      n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,&
                      mpi_comm_global,singularity_tolerance,n_nonsingular)
         if(.not.success) then
            call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                           " Exiting...",caller)
         endif

         ! Stop if n_states is larger that n_nonsingular
         if(n_nonsingular < n_states) then ! Too singular to continue
            call elsi_stop(" Overlap matrix is singular. The number of"//&
                           " basis functions after removing singularity"//&
                           " is smaller than the number of states. Try to"//&
                           " a) decrease the size of basis set, or b)"//&
                           " decrease the number of states, or c) increase"//&
                           " the tolerance of basis singularity."//&
                           " Exiting...",caller)
         elseif(n_nonsingular < n_g_size) then ! Singular
            overlap_is_singular = .true.

            if(stop_singularity) then
               call elsi_stop(" Overlap matrix is singular. This may mean"//&
                              " that a very large basis set is in use."//&
                              " Running with a near-singular basis set"//&
                              " may lead to completely wrong numerical"//&
                              " resutls. The calculation stops here,"//&
                              " because 'force_stop_singularity' is"//&
                              " set to .true. in elsi_customize."//&
                              " Exiting...",caller)
            endif

            call elsi_statement_print("  Overlap matrix is singular. This"//&
                                      " may mean that a very large basis"//&
                                      " set is in use. The calculation"//&
                                      " will continue. However, please"//&
                                      " note that running with a near-"//&
                                      "singular basis set may lead to"//&
                                      " completely wrong numerical results.")

            write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
               n_nonsingular
            call elsi_statement_print(info_str)

            call elsi_statement_print("  Using scaled eigenvectors of"//&
                                      " overlap matrix for transformation")

            ! Overlap matrix is overwritten with scaled eigenvectors
            do i = 1,n_nonsingular
               ev_sqrt = sqrt(ev_overlap(i))
               if(local_col(i) == 0) cycle
               ovlp_real(:,local_col(i)) = evec_real(:,local_col(i))/ev_sqrt
            enddo

         else ! Nonsingular
            overlap_is_singular = .false.
            call elsi_statement_print("  Overlap matrix is nonsingular")
         endif ! Singular overlap?

   end select ! select mode

   if(allocated(ev_overlap))     deallocate(ev_overlap)
   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

   call elsi_stop_singularity_check_time()

end subroutine

!> 
!! This routine does the back-transformation of the eigenvectors in standard
!! form (A'c' = c'v) to the original generalized form (Ac = Bcv)
!!
!! v = (U^-1)v'
!!
subroutine elsi_to_original_ev()

   implicit none

   logical :: success
   real*8, allocatable :: buffer_real(:,:)
   complex*16, allocatable :: buffer_complex(:,:)

   character*40, parameter :: caller = "elsi_to_original_ev"

   call elsi_start_back_transform_ev_time()

   select case (mode)
      case (COMPLEX_VALUES)
         call elsi_allocate(buffer_complex,n_l_rows,n_l_cols,"temp",caller)
         buffer_complex = evec_complex

         if(overlap_is_singular) then
            ! Transform matrix is stored in S_complex after elsi_to_standard_evp
            call pzgemm('N','N',n_g_size,n_states,n_nonsingular,(1.0d0,0.0d0),&
                        ovlp_complex,1,1,sc_desc,buffer_complex,1,1,sc_desc,&
                        (0.0d0,0.0d0),evec_complex,1,1,sc_desc)
         else ! Nonsingular, use Cholesky
            ! (U^-1) is stored in S_complex after elsi_to_standard_evp
            ! C_complex = S_complex * C_complex = S_complex * buffer_complex
            call pztranc(n_g_size,n_g_size,(1.0d0,0.0d0),ovlp_complex,1,1,sc_desc,&
                         (0.0d0,0.0d0),ham_complex,1,1,sc_desc)

            success = elpa_mult_ah_b_complex_double('L','N',n_g_size,n_states,&
                         ham_complex,n_l_rows,n_l_cols,buffer_complex,n_l_rows,&
                         n_l_cols,n_b_rows,mpi_comm_row,mpi_comm_col,evec_complex,&
                         n_l_rows,n_l_cols)
         endif

      case (REAL_VALUES)
         call elsi_allocate(buffer_real,n_l_rows,n_l_cols,"temp",caller)
         buffer_real = evec_real

         if(overlap_is_singular) then
            ! Transform matrix is stored in S_real after elsi_to_standard_evp
            call pdgemm('N','N',n_g_size,n_states,n_nonsingular,1.0d0,ovlp_real,1,&
                        1,sc_desc,buffer_real,1,1,sc_desc,0.0d0,evec_real,1,1,sc_desc)
         else ! Nonsingular, use Cholesky
            ! (U^-1) is stored in S_real after elsi_to_standard_evp
            ! C_real = S_real * C_real = S_real * buffer_real
            call pdtran(n_g_size,n_g_size,1.0d0,ovlp_real,1,1,sc_desc,0.0d0,ham_real,1,1,sc_desc)

            success = elpa_mult_at_b_real_double('L','N',n_g_size,n_states,ham_real,&
                         n_l_rows,n_l_cols,buffer_real,n_l_rows,n_l_cols,n_b_rows,&
                         mpi_comm_row,mpi_comm_col,evec_real,n_l_rows,n_l_cols)
         endif

   end select

   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

   call elsi_stop_back_transform_ev_time()

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_evp_elpa()

   implicit none
   include "mpif.h"

   logical :: success
   logical :: two_step_solver

   character*40, parameter :: caller = "elsi_solve_evp_elpa"

   elpa_print_times = .false.

   ! Choose 1-stage or 2-stage solver
   if(elpa_one_always) then
      two_step_solver = .false.
   elseif(elpa_two_always) then
      two_step_solver = .true.
   elseif(n_g_size < 256) then
      two_step_solver = .false.
   else
      two_step_solver = .true.
   endif

   ! Transform to standard form
   if(.not.overlap_is_unit) then
      call elsi_to_standard_evp()
   endif

   call elsi_start_standard_evp_time()

   ! Solve evp, return eigenvalues and eigenvectors
   if(two_step_solver) then ! 2-stage solver
      call elsi_statement_print("  Starting ELPA 2-stage eigensolver")
      select case (mode)
         case (COMPLEX_VALUES)
            success = elpa_solve_evp_complex_2stage_double(n_nonsingular,n_states,&
                         ham_complex,n_l_rows,eval,evec_complex,n_l_rows,n_b_rows,&
                         n_l_cols,mpi_comm_row,mpi_comm_col,mpi_comm_global)
         case (REAL_VALUES)
            success = elpa_solve_evp_real_2stage_double(n_nonsingular,n_states,&
                         ham_real,n_l_rows,eval,evec_real,n_l_rows,n_b_rows,&
                         n_l_cols,mpi_comm_row,mpi_comm_col,mpi_comm_global)
      end select
   else ! 1-stage solver
      call elsi_statement_print("  Starting ELPA 1-stage eigensolver")
      select case (mode)
         case (COMPLEX_VALUES)
            success = elpa_solve_evp_complex_1stage_double(n_nonsingular,n_states,&
                         ham_complex,n_l_rows,eval,evec_complex,n_l_rows,n_b_rows,&
                         n_l_cols,mpi_comm_row,mpi_comm_col,mpi_comm_global)
         case (REAL_VALUES)
            success = elpa_solve_evp_real_1stage_double(n_nonsingular,n_states,&
                         ham_real,n_l_rows,eval,evec_real,n_l_rows,n_b_rows,&
                         n_l_cols,mpi_comm_row,mpi_comm_col,mpi_comm_global)
      end select
   endif

   call elsi_stop_standard_evp_time()

   if(.not.success) then
      call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                     " Exiting...",caller)
   endif

   ! Back-transform eigenvectors
   if(.not.overlap_is_unit) then
      call elsi_to_original_ev()
   endif

   call MPI_Barrier(mpi_comm_global,mpierr)

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
subroutine elsi_to_standard_evp_sp()

   implicit none
   include 'mpif.h'

   real*8, allocatable :: buffer_real(:,:)
   complex*16, allocatable :: buffer_complex(:,:)
   logical :: success
   integer :: nblk=128
   integer :: n,nwork

   character*40, parameter :: caller = "elsi_to_standard_evp_sp"

   select case (mode)
      case (COMPLEX_VALUES)
         if(n_elsi_calls == 1) then
            if(.not.no_singularity_check) then
               call elsi_check_singularity_sp()
            endif
         endif ! n_elsi_calls == 1

         if(n_nonsingular == n_g_size) then ! Not singular
            call elsi_start_cholesky_time()

            overlap_is_singular = .false.

            ! Compute S = (U^T)U, U -> S
            success = elpa_cholesky_complex_double(n_g_size,ovlp_complex,n_l_rows,&
                         n_b_rows,n_l_cols,mpi_comm_self,mpi_comm_self,.false.)
            if(.not.success) then
               call elsi_stop(" Cholesky decomposition failed.",caller)
            endif

            ! compute U^-1 -> S
            success = elpa_invert_trm_complex_double(n_g_size,ovlp_complex,n_l_rows,&
                         n_b_rows,n_l_cols,mpi_comm_self,mpi_comm_self,.false.)
            if(.not.success) then
               call elsi_stop(" Matrix inversion failed.",caller)
            endif

            call elsi_stop_cholesky_time()
         endif

         call elsi_start_transform_evp_time()

         call elsi_allocate(buffer_complex,n_l_rows,n_l_cols,"temp",caller)

         if(overlap_is_singular) then ! Use scaled eigenvectors
            ! buffer_complex = H_complex * S_complex
            call zgemm('N','N',n_g_size,n_nonsingular,n_g_size,(1.0d0,0.0d0),&
                       ham_complex(1,1),n_g_size,ovlp_complex(1,1),n_g_size,&
                       (0.0d0,0.0d0),buffer_complex(1,1),n_g_size)

            ! H_complex = (S_complex)^* * buffer_complex
            call zgemm('C','N',n_nonsingular,n_nonsingular,n_g_size,(1.0d0,0.0d0),&
                       ovlp_complex(1,1),n_g_size,buffer_complex(1,1),n_g_size,&
                       (0.0d0,0.0d0),ham_complex(1,1),n_g_size)

         else ! Use cholesky
            ! buffer_complex = H_complex * S_complex
            call zgemm('N','N',n_g_size,n_g_size,n_g_size,(1.0d0,0.0d0),&
                       ham_complex(1,1),n_g_size,ovlp_complex(1,1),n_g_size,&
                       (0.0d0,0.0d0),buffer_complex(1,1),n_g_size)

            ! H_complex = (buffer_complex)^* * S_complex
            call zgemm('C','N',n_g_size,n_g_size,n_g_size,(1.0d0,0.0d0),&
                       buffer_complex(1,1),n_g_size,ovlp_complex(1,1),&
                       n_g_size,(0.0d0,0.0d0),ham_complex(1,1),n_g_size)
         endif

         call elsi_stop_transform_evp_time()

      case (REAL_VALUES)
         if(n_elsi_calls == 1) then
            if(.not.no_singularity_check) then
               call elsi_check_singularity_sp()
            endif
         endif ! n_elsi_calls == 1

         if(n_nonsingular == n_g_size) then ! Not singular
            call elsi_start_cholesky_time()

            overlap_is_singular = .false.

            ! Compute S = (U^T)U, U -> S
            success = elpa_cholesky_real_double(n_g_size,ovlp_real,n_l_rows,&
                         n_b_rows,n_l_cols,mpi_comm_self,mpi_comm_self,.false.)
            if(.not.success) then
               call elsi_stop(" Cholesky decomposition failed.",caller)
            endif

            ! compute U^-1 -> S
            success = elpa_invert_trm_real_double(n_g_size,ovlp_real,n_l_rows,&
                         n_b_rows,n_l_cols,mpi_comm_self,mpi_comm_self,.false.)
            if(.not.success) then
               call elsi_stop(" Matrix inversion failed.",caller)
            endif

            call elsi_stop_cholesky_time()
         endif

         call elsi_start_transform_evp_time()

         call elsi_allocate(buffer_real,n_l_rows,n_l_cols,"temp",caller)

         if(overlap_is_singular) then ! Use scaled eigenvectors
            ! buffer_real = H_real * S_real
            call dgemm('N','N',n_g_size,n_nonsingular,n_g_size,1.0d0,ham_real(1,1),&
                       n_g_size,ovlp_real(1,1),n_g_size,0.0d0,buffer_real(1,1),n_g_size)

            ! H_real = (S_real)^T * buffer_real
            call dgemm('T','N',n_nonsingular,n_nonsingular,n_g_size,1.0d0,ovlp_real(1,1),&
                       n_g_size,buffer_real(1,1),n_g_size,0.0d0,ham_real(1,1),n_g_size)

         else ! Use Cholesky
           ! buffer_real = H_real * S_real
           do n=1,n_g_size,nblk
              nwork = nblk
              if(n+nwork-1 > n_g_size) nwork=n_g_size-n+1
              call dgemm('N','N',n+nwork-1,nwork,n+nwork-1,1.0d0,ham_real(1,1),&
                         n_g_size,ovlp_real(1,n),n_g_size,0.0d0,buffer_real(1,n),n_g_size) 
           enddo

           ! H_real = (buffer_real)*T * S_real
           do n=1,n_g_size,nblk
              nwork=nblk
              if(n+nwork-1>n_g_size) nwork=n_g_size-n+1
              call dgemm('T','N',nwork,n_g_size-n+1,n+nwork-1,1.0d0,ovlp_real(1,n),n_g_size,&
                         buffer_real(1,n),n_g_size,0.0d0,ham_real(n,n),n_g_size)
           enddo
        endif

        call elsi_stop_transform_evp_time()

   end select

   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

end subroutine

subroutine elsi_to_original_ev_sp()

   implicit none

   real*8, allocatable :: buffer_real(:,:)
   complex*16, allocatable :: buffer_complex(:,:)

   character*40, parameter :: caller = "elsi_to_original_ev_sp"

   call elsi_start_back_transform_ev_time()

   select case (mode)
      case (COMPLEX_VALUES)
         call elsi_allocate(buffer_complex,n_l_rows,n_l_cols,"temp",caller)
         buffer_complex = evec_complex

         if(overlap_is_singular) then
            ! Transform matrix is stored in S_complex after elsi_to_standard_evp
            call zgemm('N','N',n_g_size,n_states,n_nonsingular,(1.0d0,0.0d0),&
                       ovlp_complex(1,1),n_g_size,buffer_complex(1,1),n_g_size,&
                       (0.0d0,0.0d0),evec_complex(1,1),n_g_size)
         else ! Nonsingular, use Cholesky
            ! (U^-1) is stored in S_complex after elsi_to_standard_evp
            ! C_complex = S_complex * C_complex = S_complex * buffer_complex
            call zgemm('N','N',n_g_size,n_states,n_g_size,(1.0d0,0.0d0),ovlp_complex(1,1),&
                       n_g_size,buffer_complex(1,1),n_g_size,(0.0d0,0.0d0),&
                       evec_complex(1,1),n_g_size)
         endif

      case (REAL_VALUES)
         call elsi_allocate(buffer_real,n_l_rows,n_l_cols,"temp",caller)
         buffer_real = evec_real

         if(overlap_is_singular) then
            ! Transform matrix is stored in S_real after elsi_to_standard_evp
            call dgemm('N','N',n_g_size,n_states,n_nonsingular,1.0d0,ovlp_real(1,1),&
                       n_g_size,buffer_real(1,1),n_g_size,0.0d0,evec_real(1,1),n_g_size)
         else ! Nonsingular, use Cholesky
            ! (U^-1) is stored in S_real after elsi_to_standard_evp
            ! C_real = S_real * C_real = S_real * buffer_real
            call dgemm('N','N',n_g_size,n_states,n_g_size,1.0d0,ovlp_real(1,1),&
                       n_g_size,buffer_real(1,1),n_g_size,0.0d0,evec_real(1,1),n_g_size)
         endif

   end select

   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

   call elsi_stop_back_transform_ev_time()

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_evp_elpa_sp()

   implicit none
   include "mpif.h"

   real*8, allocatable :: d(:),e(:)
   real*8, allocatable :: tau_real(:),buffer_real(:,:)
   complex*16, allocatable :: tau_complex(:),buffer_complex(:,:)
   logical :: success
   integer :: info

   character*40, parameter :: caller = "elsi_solve_evp_elpa_sp"

   call elsi_allocate(d,n_g_size,"d",caller)
   call elsi_allocate(e,n_g_size,"e",caller)

   ! Transform to standard form
   if(.not.overlap_is_unit) then
      call elsi_to_standard_evp_sp()
   endif

   ! Solve evp, return eigenvalues and eigenvectors
   call elsi_statement_print("  Starting ELPA eigensolver")

   call elsi_start_standard_evp_time()

   select case (mode)
      case (COMPLEX_VALUES)
         call elsi_allocate(tau_complex,n_g_size,"tau_complex",caller)
         call elsi_allocate(buffer_real,n_g_size,n_g_size,"buffer_real",caller)
         call elsi_allocate(buffer_complex,n_g_size,n_g_size,"buffer_complex",caller)

         call zhetrd('U',n_g_size,ham_complex,n_g_size,d,e,tau_complex,&
                     buffer_complex,size(buffer_complex),info)
 
         success = elpa_solve_tridi_double(n_g_size,n_states,d,e,buffer_real,&
                                           n_g_size,64,n_g_size,mpi_comm_self,&
                                           mpi_comm_self,.false.)

         eval(1:n_states) = d(1:n_states)
         evec_complex(1:n_g_size,1:n_states) = buffer_real(1:n_g_size,1:n_states)

         call zunmtr('L','U','N',n_g_size,n_states,ham_complex,n_g_size,tau_complex,&
                     evec_complex,n_g_size,buffer_complex,size(buffer_complex),info)

         deallocate(tau_complex)
         deallocate(buffer_real)
         deallocate(buffer_complex)

      case (REAL_VALUES)
         call elsi_allocate(tau_real,n_g_size,"tau_real",caller)
         call elsi_allocate(buffer_real,n_g_size,n_g_size,"buffer_real",caller)

         call dsytrd('U',n_g_size,ham_real,n_g_size,d,e,tau_real,buffer_real,&
                     size(buffer_real),info)

         success = elpa_solve_tridi_double(n_g_size,n_states,d,e,buffer_real,&
                                           n_g_size,64,n_g_size,mpi_comm_self,&
                                           mpi_comm_self,.false.)
 
         eval(1:n_states) = d(1:n_states)
         evec_real(1:n_g_size,1:n_states) = buffer_real(1:n_g_size,1:n_states)
 
         call dormtr('L','U','N',n_g_size,n_states,ham_real,n_g_size,tau_real,&
                     evec_real,n_g_size,buffer_real,size(buffer_real),info)

         deallocate(tau_real)
         deallocate(buffer_real)

   end select

   call elsi_stop_standard_evp_time()

   deallocate(d)
   deallocate(e)

   if(.not.success) then
      call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                     " Exiting...",caller)
   endif

   ! Back-transform eigenvectors
   if(.not.overlap_is_unit) then
      call elsi_to_original_ev_sp()
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
subroutine elsi_check_singularity_sp()

   implicit none
   include 'mpif.h'   

   real*8 :: ev_sqrt
   real*8, allocatable :: ev_overlap(:)
   real*8, allocatable :: buffer_real(:,:)
   complex*16, allocatable :: buffer_complex(:,:)
   integer :: i,i_col
   logical :: success

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_check_singularity_sp"

   call elsi_start_singularity_check_time()

   select case (mode)
      case (COMPLEX_VALUES)
         call elsi_allocate(buffer_complex,n_l_rows,n_l_cols,"temp",caller)

         ! Use buffer_complex to store overlap matrix, otherwise it will
         ! be destroyed by eigenvalue calculation
         ! The nonsingular eigenvalues must be the first ones, so find
         ! eigenvalues of negative overlap matrix
         buffer_complex = -ovlp_complex

         call elsi_allocate(ev_overlap,n_g_size,"ev_overlap",caller)

         ! Use customized ELPA 2-stage solver to check overlap singularity
         ! Eigenvectors computed only for singular overlap matrix
         success = elpa_check_singularity_complex_double(n_g_size,n_g_size,&
                      buffer_complex,n_l_rows,ev_overlap,evec_complex,&
                      n_l_rows,n_b_rows,n_l_cols,mpi_comm_self,mpi_comm_self,&
                      mpi_comm_self,singularity_tolerance,n_nonsingular)
         if(.not.success) then
            call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                           " Exiting...",caller)
         endif

         ! Stop if n_states is larger that n_nonsingular
         if(n_nonsingular < n_states) then ! Too singular to continue
            call elsi_stop(" Overlap matrix is singular. The number of"//&
                           " basis functions after removing singularity"//&
                           " is smaller than the number of states. Try to"//&
                           " a) decrease the size of basis set, or b)"//&
                           " decrease the number of states, or c) increase"//&
                           " the tolerance of basis singularity."//&
                           " Exiting...",caller)
         elseif(n_nonsingular < n_g_size) then ! Singular
            overlap_is_singular = .true.

            if(stop_singularity) then
               call elsi_stop(" Overlap matrix is singular. This may mean"//&
                              " that a very large basis set is in use."//&
                              " Running with a near-singular basis set"//&
                              " may lead to completely wrong numerical"//&
                              " resutls. The calculation stops here,"//&
                              " because 'force_stop_singularity' is"//&
                              " set to .true. in elsi_customize."//&
                              " Exiting...",caller)
            endif

            call elsi_statement_print("  Overlap matrix is singular. This"//&
                                      " may mean that a very large basis"//&
                                      " set is in use. The calculation"//&
                                      " will continue. However, please"//&
                                      " note that running with a near-"//&
                                      "singular basis set may lead to"//&
                                      " completely wrong numerical results.")

            write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
               n_nonsingular
            call elsi_statement_print(info_str)

            call elsi_statement_print("  Using scaled eigenvectors of"//&
                                      " overlap matrix for transformation")

            ! Overlap matrix is overwritten with scaled eigenvectors
            do i = 1,n_nonsingular
               ev_sqrt = sqrt(ev_overlap(i))
               ovlp_complex(:,i) = evec_complex(:,i)/ev_sqrt
            enddo

         else ! Nonsingular
            overlap_is_singular = .false.
            call elsi_statement_print("  Overlap matrix is nonsingular")
         endif ! Singular overlap?

      case (REAL_VALUES)
         call elsi_allocate(buffer_real,n_l_rows,n_l_cols,"temp",caller)

         ! Use buffer_real to store overlap matrix, otherwise it will be
         ! destroyed by eigenvalue calculation
         ! The nonsingular eigenvalues must be the first ones, so find
         ! eigenvalues of negative overlap matrix
         buffer_real = -ovlp_real

         call elsi_allocate(ev_overlap,n_g_size,"ev_overlap",caller)

         ! Use customized ELPA 2-stage solver to check overlap singularity
         ! Eigenvectors computed only for singular overlap matrix
         success = elpa_check_singularity_real_double(n_g_size,n_g_size,&
                      buffer_real,n_l_rows,ev_overlap,evec_real,n_l_rows,&
                      n_b_rows,n_l_cols,mpi_comm_self,mpi_comm_self,&
                      mpi_comm_self,singularity_tolerance,n_nonsingular)
         if(.not.success) then
            call elsi_stop(" ELPA failed when solving eigenvalue problem."//&
                           " Exiting...",caller)
         endif

         ! Stop if n_states is larger that n_nonsingular
         if(n_nonsingular < n_states) then ! Too singular to continue
            call elsi_stop(" Overlap matrix is singular. The number of"//&
                           " basis functions after removing singularity"//&
                           " is smaller than the number of states. Try to"//&
                           " a) decrease the size of basis set, or b)"//&
                           " decrease the number of states, or c) increase"//&
                           " the tolerance of basis singularity."//&
                           " Exiting...",caller)
         elseif(n_nonsingular < n_g_size) then ! Singular
            overlap_is_singular = .true.

            if(stop_singularity) then
               call elsi_stop(" Overlap matrix is singular. This may mean"//&
                              " that a very large basis set is in use."//&
                              " Running with a near-singular basis set"//&
                              " may lead to completely wrong numerical"//&
                              " resutls. The calculation stops here,"//&
                              " because 'force_stop_singularity' is"//&
                              " set to .true. in elsi_customize."//&
                              " Exiting...",caller)
            endif

            call elsi_statement_print("  Overlap matrix is singular. This"//&
                                      " may mean that a very large basis"//&
                                      " set is in use. The calculation"//&
                                      " will continue. However, please"//&
                                      " note that running with a near-"//&
                                      "singular basis set may lead to"//&
                                      " completely wrong numerical results.")

            write(info_str,"(A,I13)") "  | Number of basis functions reduced to: ",&
               n_nonsingular
            call elsi_statement_print(info_str)

            call elsi_statement_print("  Using scaled eigenvectors of"//&
                                      " overlap matrix for transformation")

            ! Overlap matrix is overwritten with scaled eigenvectors
            do i = 1,n_nonsingular
               ev_sqrt = sqrt(ev_overlap(i))
               ovlp_real(:,i) = evec_real(:,i)/ev_sqrt
            enddo

         else ! Nonsingular
            overlap_is_singular = .false.
            call elsi_statement_print("  Overlap matrix is nonsingular")
         endif ! Singular overlap?

   end select ! select mode

   if(allocated(ev_overlap))     deallocate(ev_overlap)
   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

   call elsi_stop_singularity_check_time()

end subroutine

end module 
