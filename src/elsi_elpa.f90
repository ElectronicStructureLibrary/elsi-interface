! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to ELPA.
!!
module ELSI_ELPA

   use ELSI_CONSTANTS, only: UT_MAT,UNSET
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_IO, only: elsi_say,elsi_get_time
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI, only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_integer4
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS, only: elsi_get_nnz,elsi_set_full_mat
   use CHECK_SINGULARITY, only: elpa_check_singularity_real_double,&
       elpa_check_singularity_complex_double
   use ELPA1, only: elpa_print_times,elpa_solve_evp_real_1stage_double,&
       elpa_solve_evp_complex_1stage_double,elpa_cholesky_real_double,&
       elpa_cholesky_complex_double,elpa_invert_trm_real_double,&
       elpa_invert_trm_complex_double,elpa_mult_at_b_real_double,&
       elpa_mult_ah_b_complex_double
   use ELPA2, only: elpa_solve_evp_real_2stage_double,&
       elpa_solve_evp_complex_2stage_double

   implicit none

   private

   public :: elsi_init_elpa
   public :: elsi_cleanup_elpa
   public :: elsi_solve_elpa

   interface elsi_solve_elpa
      module procedure elsi_solve_elpa_real
      module procedure elsi_solve_elpa_cmplx
   end interface

contains

!>
!! This routine initializes ELPA.
!!
subroutine elsi_init_elpa(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh

   integer(kind=i4) :: ierr
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_init_elpa"

   if(.not. ph%elpa_started) then
      call MPI_Comm_split(bh%comm,bh%my_pcol,bh%my_prow,ph%elpa_comm_row,ierr)

      call elsi_check_mpi(bh,"MPI_Comm_split",ierr,caller)

      call MPI_Comm_split(bh%comm,bh%my_prow,bh%my_pcol,ph%elpa_comm_col,ierr)

      call elsi_check_mpi(bh,"MPI_Comm_split",ierr,caller)

      ph%elpa_started = .true.

      if(ph%elpa_gpu) then
         write(msg,"(2X,A)") "No ELPA GPU acceleration available"
         call elsi_say(bh,msg)
      end if

      if(ph%elpa_n_single > 0) then
         write(msg,"(2X,A)") "No single precision ELPA available"
         call elsi_say(bh,msg)
      end if

      if(ph%elpa_autotune) then
         write(msg,"(2X,A)") "No ELPA auto-tuning available"
         call elsi_say(bh,msg)
      end if
   end if

end subroutine

!>
!! This routine transforms a generalized eigenproblem to standard and returns
!! the Cholesky factor for later use.
!!
subroutine elsi_to_standard_evp_real(ph,bh,row_map,col_map,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_to_standard_evp_real"

   if(ph%elpa_first) then
      if(ph%check_sing) then
         call elsi_check_singularity_real(ph,bh,col_map,ovlp,eval,evec)
      end if

      if(ph%n_good == ph%n_basis) then ! Not singular
         call elsi_get_time(t0)

         ph%ovlp_is_sing = .false.

         ! S = U
         success = elpa_cholesky_real_double(ph%n_basis,ovlp,bh%n_lrow,bh%blk,&
                      bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,.false.)

         if(.not. success) then
            call elsi_stop(bh,"Cholesky failed.",caller)
         end if

         ! S = U^(-1)
         success = elpa_invert_trm_real_double(ph%n_basis,ovlp,bh%n_lrow,&
                      bh%blk,bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,&
                      .false.)

         if(.not. success) then
            call elsi_stop(bh,"Matrix inversion failed.",caller)
         end if

         call elsi_get_time(t1)

         write(msg,"(2X,A)") "Finished Cholesky decomposition"
         call elsi_say(bh,msg)
         write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
         call elsi_say(bh,msg)
      end if
   end if

   call elsi_get_time(t0)

   ! H = U^(-T) H U^(-1)
   if(ph%ovlp_is_sing) then
      call pdgemm("N","N",ph%n_basis,ph%n_good,ph%n_basis,1.0_r8,ham,1,1,&
              bh%desc,ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,0.0_r8,evec,1,1,&
              bh%desc)

      call pdgemm("T","N",ph%n_good,ph%n_good,ph%n_basis,1.0_r8,ovlp,1,&
              ph%n_basis-ph%n_good+1,bh%desc,evec,1,1,bh%desc,0.0_r8,ham,1,1,&
              bh%desc)
   else ! Use Cholesky
      success = elpa_mult_at_b_real_double("U","L",ph%n_basis,ph%n_basis,ovlp,&
                   bh%n_lrow,bh%n_lcol,ham,bh%n_lrow,bh%n_lcol,bh%blk,&
                   ph%elpa_comm_row,ph%elpa_comm_col,evec,bh%n_lrow,bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if

      call pdtran(ph%n_basis,ph%n_basis,1.0_r8,evec,1,1,bh%desc,0.0_r8,ham,1,1,&
              bh%desc)

      evec = ham

      success = elpa_mult_at_b_real_double("U","U",ph%n_basis,ph%n_basis,ovlp,&
                   bh%n_lrow,bh%n_lcol,evec,bh%n_lrow,bh%n_lcol,bh%blk,&
                   ph%elpa_comm_row,ph%elpa_comm_col,ham,bh%n_lrow,bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if

      call elsi_set_full_mat(ph,bh,UT_MAT,row_map,col_map,ham)
   end if

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished transformation to standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! This routine checks the singularity of overlap matrix by computing all its
!! eigenvalues. On exit, S is not modified if not singular, or is overwritten by
!! scaled eigenvectors if singular, which can be used to transform the
!! generalized eigenproblem to the standard form.
!!
subroutine elsi_check_singularity_real(ph,bh,col_map,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: ev_sqrt
   integer(kind=i4) :: i
   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   real(kind=r8), allocatable :: copy(:,:)

   character(len=*), parameter :: caller = "elsi_check_singularity_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,copy,bh%n_lrow,bh%n_lcol,"copy",caller)

   ! Overlap will be destroyed by eigenvalue calculation
   copy = ovlp

   ! Use customized ELPA 2-stage solver to check overlap singularity
   ! Eigenvectors computed only for singular overlap matrix
   success = elpa_check_singularity_real_double(ph%n_basis,ph%n_basis,copy,&
                bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,bh%n_lcol,&
                ph%elpa_comm_row,ph%elpa_comm_col,bh%comm,ph%sing_tol,ph%n_good)

   if(.not. success) then
      call elsi_stop(bh,"Singularity check failed.",caller)
   end if

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
      do i = ph%n_basis-ph%n_good+1,ph%n_basis
         ev_sqrt = sqrt(eval(i))

         if(col_map(i) > 0) then
            ovlp(:,col_map(i)) = evec(:,col_map(i))/ev_sqrt
         end if
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
!! This routine back-transforms eigenvectors in the standard form to the
!! original generalized form.
!!
subroutine elsi_to_original_ev_real(ph,bh,ham,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(out) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   real(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_to_original_ev_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)

   tmp = evec

   if(ph%ovlp_is_sing) then
      call pdgemm("N","N",ph%n_basis,ph%n_states_solve,ph%n_good,1.0_r8,ovlp,1,&
              ph%n_basis-ph%n_good+1,bh%desc,tmp,1,1,bh%desc,0.0_r8,evec,1,1,&
              bh%desc)
   else ! Nonsingular, use Cholesky
      call pdtran(ph%n_basis,ph%n_basis,1.0_r8,ovlp,1,1,bh%desc,0.0_r8,ham,1,1,&
              bh%desc)

      success = elpa_mult_at_b_real_double("L","N",ph%n_basis,ph%n_states,ham,&
                   bh%n_lrow,bh%n_lcol,tmp,bh%n_lrow,bh%n_lcol,bh%blk,&
                   ph%elpa_comm_row,ph%elpa_comm_col,evec,bh%n_lrow,bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if
   end if

   call elsi_deallocate(bh,tmp,"tmp")

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished back-transformation of eigenvectors"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_elpa_real(ph,bh,row_map,col_map,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: ierr
   logical :: success
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_solve_elpa_real"

   elpa_print_times = ph%elpa_output

   ! Compute sparsity
   if(bh%nnz_g == UNSET) then
      if(bh%nnz_l == UNSET) then
         call elsi_get_nnz(bh%def0,ham,bh%n_lrow,bh%n_lcol,bh%nnz_l)
      end if

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   end if

   ! Transform to standard form
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_standard_evp_real(ph,bh,row_map,col_map,ham,ovlp,eval,evec)
   end if

   call elsi_get_time(t0)

   write(msg,"(2X,A)") "Starting ELPA eigensolver"
   call elsi_say(bh,msg)

   ! Solve evp, return eigenvalues and eigenvectors
   if(ph%elpa_solver == 2) then
      success = elpa_solve_evp_real_2stage_double(ph%n_good,ph%n_states_solve,&
                   ham,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,bh%n_lcol,&
                   ph%elpa_comm_row,ph%elpa_comm_col,bh%comm)
   else
      success = elpa_solve_evp_real_1stage_double(ph%n_good,ph%n_states_solve,&
                   ham,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,bh%n_lcol,&
                   ph%elpa_comm_row,ph%elpa_comm_col,bh%comm)
   end if

   if(.not. success) then
      call elsi_stop(bh,"ELPA solver failed.",caller)
   end if

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
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
      call elsi_to_original_ev_real(ph,bh,ham,ovlp,evec)
   end if

   ph%elpa_first = .false.

end subroutine

!>
!! This routine transforms a generalized eigenproblem to standard and returns
!! the Cholesky factor for later use.
!!
subroutine elsi_to_standard_evp_cmplx(ph,bh,row_map,col_map,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   complex(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_to_standard_evp_cmplx"

   if(ph%elpa_first) then
      if(ph%check_sing) then
         call elsi_check_singularity_cmplx(ph,bh,col_map,ovlp,eval,evec)
      end if

      if(ph%n_good == ph%n_basis) then ! Not singular
         call elsi_get_time(t0)

         ph%ovlp_is_sing = .false.

         ! S = U
         success = elpa_cholesky_complex_double(ph%n_basis,ovlp,bh%n_lrow,&
                      bh%blk,bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,&
                      .false.)

         if(.not. success) then
            call elsi_stop(bh,"Cholesky failed.",caller)
         end if

         ! S = U^(-1)
         success = elpa_invert_trm_complex_double(ph%n_basis,ovlp,bh%n_lrow,&
                      bh%blk,bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,&
                      .false.)

         if(.not. success) then
            call elsi_stop(bh,"Matrix inversion failed.",caller)
         end if

         call elsi_get_time(t1)

         write(msg,"(2X,A)") "Finished Cholesky decomposition"
         call elsi_say(bh,msg)
         write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
         call elsi_say(bh,msg)
      end if
   end if

   call elsi_get_time(t0)

   ! H = U^(-T) H U^(-1)
   if(ph%ovlp_is_sing) then
      call pzgemm("N","N",ph%n_basis,ph%n_good,ph%n_basis,(1.0_r8,0.0_r8),ham,&
              1,1,bh%desc,ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,&
              (0.0_r8,0.0_r8),evec,1,1,bh%desc)

      call pzgemm("C","N",ph%n_good,ph%n_good,ph%n_basis,(1.0_r8,0.0_r8),ovlp,&
              1,ph%n_basis-ph%n_good+1,bh%desc,evec,1,1,bh%desc,&
              (0.0_r8,0.0_r8),ham,1,1,bh%desc)
   else ! Use cholesky
      success = elpa_mult_ah_b_complex_double("U","L",ph%n_basis,ph%n_basis,&
                   ovlp,bh%n_lrow,bh%n_lcol,ham,bh%n_lrow,bh%n_lcol,bh%blk,&
                   ph%elpa_comm_row,ph%elpa_comm_col,evec,bh%n_lrow,bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if

      call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),evec,1,1,bh%desc,&
              (0.0_r8,0.0_r8),ham,1,1,bh%desc)

      evec = ham

      success = elpa_mult_ah_b_complex_double("U","U",ph%n_basis,ph%n_basis,&
                   ovlp,bh%n_lrow,bh%n_lcol,evec,bh%n_lrow,bh%n_lcol,bh%blk,&
                   ph%elpa_comm_row,ph%elpa_comm_col,ham,bh%n_lrow,bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if

      call elsi_set_full_mat(ph,bh,UT_MAT,row_map,col_map,ham)
   end if

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished transformation to standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! This routine checks the singularity of overlap matrix by computing all its
!! eigenvalues. On exit, S is not modified if not singular, or is overwritten by
!! scaled eigenvectors if singular, which can be used to transform the
!! generalized eigenproblem to the standard form.
!!
subroutine elsi_check_singularity_cmplx(ph,bh,col_map,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: ev_sqrt
   integer(kind=i4) :: i
   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   complex(kind=r8), allocatable :: copy(:,:)

   character(len=*), parameter :: caller = "elsi_check_singularity_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,copy,bh%n_lrow,bh%n_lcol,"copy",caller)

   ! Overlap will be destroyed by eigenvalue calculation
   copy = ovlp

   ! Use customized ELPA 2-stage solver to check overlap singularity
   ! Eigenvectors computed only for singular overlap matrix
   success = elpa_check_singularity_complex_double(ph%n_basis,ph%n_basis,copy,&
                bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,bh%n_lcol,&
                ph%elpa_comm_row,ph%elpa_comm_col,bh%comm,ph%sing_tol,ph%n_good)

   if(.not. success) then
      call elsi_stop(bh,"Singularity check failed.",caller)
   end if

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
      do i = ph%n_basis-ph%n_good+1,ph%n_basis
         ev_sqrt = sqrt(eval(i))

         if(col_map(i) > 0) then
            ovlp(:,col_map(i)) = evec(:,col_map(i))/ev_sqrt
         end if
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
!! This routine back-transforms eigenvectors in the standard form to the
!! original generalized form.
!!
subroutine elsi_to_original_ev_cmplx(ph,bh,ham,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(out) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   logical :: success
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   complex(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_to_original_ev_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)

   tmp = evec

   if(ph%ovlp_is_sing) then
      call pzgemm("N","N",ph%n_basis,ph%n_states_solve,ph%n_good,&
              (1.0_r8,0.0_r8),ovlp,1,ph%n_basis-ph%n_good+1,bh%desc,tmp,1,1,&
              bh%desc,(0.0_r8,0.0_r8),evec,1,1,bh%desc)
   else ! Nonsingular, use Cholesky
      call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),ovlp,1,1,bh%desc,&
              (0.0_r8,0.0_r8),ham,1,1,bh%desc)

      success = elpa_mult_ah_b_complex_double("L","N",ph%n_basis,ph%n_states,&
                   ham,bh%n_lrow,bh%n_lcol,tmp,bh%n_lrow,bh%n_lcol,bh%blk,&
                   ph%elpa_comm_row,ph%elpa_comm_col,evec,bh%n_lrow,bh%n_lcol)

      if(.not. success) then
         call elsi_stop(bh,"Matrix multiplication failed.",caller)
      end if
   end if

   call elsi_deallocate(bh,tmp,"tmp")

   call elsi_get_time(t1)

   write(msg,"(2X,A)") "Finished back-transformation of eigenvectors"
   call elsi_say(bh,msg)
   write(msg,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! This routine interfaces to ELPA.
!!
subroutine elsi_solve_elpa_cmplx(ph,bh,row_map,col_map,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4), intent(in) :: row_map(ph%n_basis)
   integer(kind=i4), intent(in) :: col_map(ph%n_basis)
   complex(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: ierr
   logical :: success
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_solve_elpa_cmplx"

   elpa_print_times = ph%elpa_output

   ! Compute sparsity
   if(bh%nnz_g == UNSET) then
      if(bh%nnz_l == UNSET) then
         call elsi_get_nnz(bh%def0,ham,bh%n_lrow,bh%n_lcol,bh%nnz_l)
      end if

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   end if

   ! Transform to standard form
   if(.not. ph%ovlp_is_unit) then
      call elsi_to_standard_evp_cmplx(ph,bh,row_map,col_map,ham,ovlp,eval,evec)
   end if

   call elsi_get_time(t0)

   write(msg,"(2X,A)") "Starting ELPA eigensolver"
   call elsi_say(bh,msg)

   ! Solve evp, return eigenvalues and eigenvectors
   if(ph%elpa_solver == 2) then
      success = elpa_solve_evp_complex_2stage_double(ph%n_good,&
                   ph%n_states_solve,ham,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,&
                   bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,bh%comm)
   else
      success = elpa_solve_evp_complex_1stage_double(ph%n_good,&
                   ph%n_states_solve,ham,bh%n_lrow,eval,evec,bh%n_lrow,bh%blk,&
                   bh%n_lcol,ph%elpa_comm_row,ph%elpa_comm_col,bh%comm)
   end if

   if(.not. success) then
      call elsi_stop(bh,"ELPA solver failed.",caller)
   end if

   ! Dummy eigenvalues for correct chemical potential, no physical meaning!
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
      call elsi_to_original_ev_cmplx(ph,bh,ham,ovlp,evec)
   end if

   ph%elpa_first = .false.

end subroutine

!>
!! This routine cleans up ELPA.
!!
subroutine elsi_cleanup_elpa(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_cleanup_elpa"

   if(ph%elpa_started) then
      call MPI_Comm_free(ph%elpa_comm_row,ierr)
      call MPI_Comm_free(ph%elpa_comm_col,ierr)
   end if

   ph%elpa_started = .false.

end subroutine

end module ELSI_ELPA
