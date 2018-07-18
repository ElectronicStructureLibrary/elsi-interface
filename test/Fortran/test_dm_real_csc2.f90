! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests real density matrix solver, PEXSI_CSC format.
!!
subroutine test_dm_real_csc2(mpi_comm,solver,h_file,s_file)

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   integer(kind=i4), intent(in) :: mpi_comm
   integer(kind=i4), intent(in) :: solver
   character(len=*), intent(in) :: h_file
   character(len=*), intent(in) :: s_file

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: myid
   integer(kind=i4) :: ierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: n_states
   integer(kind=i4) :: matrix_size
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: n_l_cols
   integer(kind=i4) :: header(8) = 0

   real(kind=r8) :: n_electrons
   real(kind=r8) :: n_test
   real(kind=r8) :: tmp
   real(kind=r8) :: e_test = 0.0_r8
   real(kind=r8) :: e_ref  = 0.0_r8
   real(kind=r8) :: tol    = 0.0_r8
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   real(kind=r8),    allocatable :: ham(:)
   real(kind=r8),    allocatable :: ovlp(:)
   real(kind=r8),    allocatable :: dm(:)
   real(kind=r8),    allocatable :: edm(:)
   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)

   real(kind=r8), external :: ddot

   type(elsi_handle)    :: e_h
   type(elsi_rw_handle) :: rw_h

   ! Reference values
   real(kind=r8), parameter :: e_elpa   = -2622.88214509316_r8
   real(kind=r8), parameter :: e_omm    = -2622.88214509316_r8
   real(kind=r8), parameter :: e_pexsi  = -2622.88194292325_r8
   real(kind=r8), parameter :: e_sips   = -2622.88214509316_r8
   real(kind=r8), parameter :: e_ntpoly = -2622.88214509311_r8

   call MPI_Comm_size(mpi_comm,n_proc,ierr)
   call MPI_Comm_rank(mpi_comm,myid,ierr)

   if(myid == 0) then
      tol = 1.0e-8_r8
      write(*,"(2X,A)") "################################"
      write(*,"(2X,A)") "##     ELSI TEST PROGRAMS     ##"
      write(*,"(2X,A)") "################################"
      write(*,*)
      if(solver == 1) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_real_sparse + ELPA"
         e_ref = e_elpa
      elseif(solver == 2) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_real_sparse + libOMM"
         e_ref = e_omm
      elseif(solver == 3) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_real_sparse + PEXSI"
         e_ref = e_pexsi
         tol   = 1.0e-3_r8
      elseif(solver == 5) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_real_sparse + SLEPc-SIPs"
         e_ref = e_sips
      elseif(solver == 6) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_real_sparse + NTPoly"
         e_ref = e_ntpoly
      endif
      write(*,*)
   endif

   ! Read H and S matrices
   call elsi_init_rw(rw_h,0,1,0,0.0_r8)
   call elsi_set_rw_mpi(rw_h,mpi_comm)

   call elsi_read_mat_dim_sparse(rw_h,h_file,n_electrons,matrix_size,nnz_g,&
           nnz_l,n_l_cols)
   call elsi_get_rw_header(rw_h,header)

   ! This is a limitation of this program, not a limitation of ELSI
   if(mod(matrix_size,n_proc) /= 0) then
      blk = 0

      write(*,"(2X,A)") "Internal error"
      call MPI_Abort(mpi_comm,0,ierr)
      stop
   else
      blk = matrix_size/n_proc
   endif

   allocate(ham(nnz_l))
   allocate(ovlp(nnz_l))
   allocate(dm(nnz_l))
   allocate(edm(nnz_l))
   allocate(row_ind(nnz_l))
   allocate(col_ptr(n_l_cols+1))

   t1 = MPI_Wtime()

   call elsi_read_mat_real_sparse(rw_h,h_file,row_ind,col_ptr,ham)
   call elsi_read_mat_real_sparse(rw_h,s_file,row_ind,col_ptr,ovlp)

   call elsi_finalize_rw(rw_h)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished reading H and S matrices"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   endif

   ! Initialize ELSI
   n_states = int(n_electrons,kind=i4)

   call elsi_init(e_h,solver,1,2,matrix_size,n_electrons,n_states)
   call elsi_set_mpi(e_h,mpi_comm)
   call elsi_set_csc(e_h,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)
   call elsi_set_csc_blk(e_h,blk)

   ! Customize ELSI
   call elsi_set_output(e_h,2)
   call elsi_set_output_log(e_h,1)
   call elsi_set_sing_check(e_h,0)
   call elsi_set_mu_broaden_width(e_h,1.0e-6_r8)
   call elsi_set_omm_n_elpa(e_h,1)
   call elsi_set_pexsi_delta_e(e_h,80.0_r8)
   call elsi_set_pexsi_np_per_pole(e_h,2)
   call elsi_set_sips_n_elpa(e_h,1)

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 1)
   call elsi_dm_real_sparse(e_h,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #1"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   endif

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 2, with the same H)
   call elsi_dm_real_sparse(e_h,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   ! Compute energy density matrix
   if(solver == 1 .or. solver == 2 .or. solver == 3) then
      call elsi_get_edm_real_sparse(e_h,edm)
   endif

   ! Compute electron count
   tmp = ddot(nnz_l,ovlp,1,dm,1)

   call MPI_Reduce(tmp,n_test,1,mpi_real8,mpi_sum,0,mpi_comm,ierr)

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #2"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
      write(*,"(2X,A)") "Finished test program"
      write(*,*)
      if(header(8) == 1111) then
         if(abs(e_test-e_ref) < tol .and. abs(n_test-n_electrons) < tol) then
            write(*,"(2X,A)") "Passed."
         else
            write(*,"(2X,A)") "Failed."
         endif
      endif
      write(*,*)
   endif

   ! Finalize ELSI
   call elsi_finalize(e_h)

   deallocate(ham)
   deallocate(ovlp)
   deallocate(dm)
   deallocate(edm)
   deallocate(row_ind)
   deallocate(col_ptr)

end subroutine
