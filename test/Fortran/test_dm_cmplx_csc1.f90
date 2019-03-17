! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests complex density matrix solver, PEXSI_CSC format.
!!
subroutine test_dm_cmplx_csc1(mpi_comm,solver,h_file,s_file)

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
   integer(kind=i4) :: comm_in_task
   integer(kind=i4) :: comm_among_task
   integer(kind=i4) :: ierr
   integer(kind=i4) :: n_states
   integer(kind=i4) :: n_basis
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: n_l_cols
   integer(kind=i4) :: id_in_task
   integer(kind=i4) :: task_id
   integer(kind=i4) :: buffer(5)
   integer(kind=i4) :: header(8) = 0

   real(kind=r8) :: n_electrons
   real(kind=r8) :: n_test
   real(kind=r8) :: tmp
   real(kind=r8) :: e_test = 0.0_r8
   real(kind=r8) :: e_ref = 0.0_r8
   real(kind=r8) :: tol = 0.0_r8
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   logical :: file_exist

   complex(kind=r8), allocatable :: ham(:)
   complex(kind=r8), allocatable :: ovlp(:)
   complex(kind=r8), allocatable :: dm(:)
   complex(kind=r8), allocatable :: edm(:)

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)

   type(elsi_handle) :: eh
   type(elsi_rw_handle) :: rwh

   ! Reference values
   real(kind=r8), parameter :: e_elpa = -2622.88214509316_r8
   real(kind=r8), parameter :: e_omm = -2622.88214509316_r8
   real(kind=r8), parameter :: e_pexsi = -2622.88194292325_r8
   real(kind=r8), parameter :: e_ntpoly = -2622.88214509311_r8

   character(len=*), parameter :: file_name = "elsi.in"

   complex(kind=r8), external :: zdotc

   call MPI_Comm_size(mpi_comm,n_proc,ierr)
   call MPI_Comm_rank(mpi_comm,myid,ierr)

   if(myid == 0) then
      tol = 1.0e-8_r8
      write(*,"(2X,A)") "################################"
      write(*,"(2X,A)") "##     ELSI TEST PROGRAMS     ##"
      write(*,"(2X,A)") "################################"
      write(*,*)
      if(solver == 1) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex_sparse + ELPA"
         e_ref = e_elpa
      else if(solver == 2) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex_sparse + libOMM"
         e_ref = e_omm
      else if(solver == 3) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex_sparse + PEXSI"
         e_ref = e_pexsi
         tol = 1.0e-3_r8
      else if(solver == 6) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex_sparse + NTPoly"
         e_ref = e_ntpoly
         tol = 1.0e-7_r8
      end if
      write(*,*)
   end if

   if(solver == 3) then
      task_id = myid/2
      id_in_task = mod(myid,2)
   else
      task_id = 0
      id_in_task = myid
   end if

   call MPI_Comm_split(mpi_comm,task_id,id_in_task,comm_in_task,ierr)
   call MPI_Comm_split(mpi_comm,id_in_task,task_id,comm_among_task,ierr)

   if(task_id == 0) then
      ! Read H and S matrices
      call elsi_init_rw(rwh,0,1,0,0.0_r8)
      call elsi_set_rw_mpi(rwh,comm_in_task)

      call elsi_read_mat_dim_sparse(rwh,h_file,n_electrons,n_basis,nnz_g,nnz_l,&
              n_l_cols)
      call elsi_get_rw_header(rwh,header)
   end if

   if(solver == 3) then
      if(task_id == 0) then
         buffer(1) = int(n_electrons,kind=i4)
         buffer(2) = n_basis
         buffer(3) = nnz_g
         buffer(4) = nnz_l
         buffer(5) = n_l_cols
      end if

      call MPI_Bcast(buffer,5,mpi_integer4,0,comm_among_task,ierr)

      n_electrons = real(buffer(1),kind=r8)
      n_basis = buffer(2)
      nnz_g = buffer(3)
      nnz_l = buffer(4)
      n_l_cols = buffer(5)
   end if

   allocate(ham(nnz_l))
   allocate(ovlp(nnz_l))
   allocate(dm(nnz_l))
   allocate(edm(nnz_l))
   allocate(row_ind(nnz_l))
   allocate(col_ptr(n_l_cols+1))

   t1 = MPI_Wtime()

   if(task_id == 0) then
      call elsi_read_mat_complex_sparse(rwh,h_file,row_ind,col_ptr,ham)
      call elsi_read_mat_complex_sparse(rwh,s_file,row_ind,col_ptr,ovlp)

      call elsi_finalize_rw(rwh)
   end if

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished reading H and S matrices"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Initialize ELSI
   n_states = int(n_electrons,kind=i4)

   call elsi_init(eh,solver,1,1,n_basis,n_electrons,n_states)
   call elsi_set_mpi(eh,mpi_comm)
   call elsi_set_csc(eh,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)

   ! Customize ELSI
   call elsi_set_output(eh,2)
   call elsi_set_output_log(eh,1)
   call elsi_set_illcond_check(eh,0)
   call elsi_set_mu_broaden_width(eh,1.0e-6_r8)
   call elsi_set_omm_n_elpa(eh,1)
   call elsi_set_pexsi_delta_e(eh,80.0_r8)
   call elsi_set_pexsi_np_per_pole(eh,2)

   inquire(file=file_name,exist=file_exist)

   if(file_exist) then
      call elsi_set_input_file(eh,file_name)
   end if

   t1 = MPI_Wtime()

   ! Solve
   call elsi_dm_complex_sparse(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #1"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex_sparse(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #2"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex_sparse(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #3"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Reinit for a new geometry
   call elsi_reinit(eh)
   call elsi_set_csc(eh,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)
   call elsi_set_elpa_solver(eh,1)
   call elsi_set_omm_flavor(eh,2)
   call elsi_set_ntpoly_method(eh,1)

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex_sparse(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #4"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex_sparse(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   ! Compute energy density matrix
   call elsi_get_edm_complex_sparse(eh,edm)

   ! Compute electron count
   if(task_id == 0) then
      tmp = real(zdotc(nnz_l,ovlp,1,dm,1),kind=r8)

      call MPI_Reduce(tmp,n_test,1,mpi_real8,mpi_sum,0,comm_in_task,ierr)
   end if

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #5"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
      write(*,"(2X,A)") "Finished test program"
      write(*,*)
      if(header(8) == 1111) then
         write(*,"(2X,A)") "Band energy"
         write(*,"(2X,A,F15.8)") "| This test :",e_test
         write(*,"(2X,A,F15.8)") "| Reference :",e_ref
         write(*,*)
         write(*,"(2X,A)") "Electron count"
         write(*,"(2X,A,F15.8)") "| This test :",n_test
         write(*,"(2X,A,F15.8)") "| Reference :",n_electrons
         write(*,*)
         if(abs(e_test-e_ref) < tol .and. abs(n_test-n_electrons) < tol) then
            write(*,"(2X,A)") "Passed."
         else
            write(*,"(2X,A)") "Failed."
         end if
      end if
      write(*,*)
   end if

   ! Finalize ELSI
   call elsi_finalize(eh)

   deallocate(ham)
   deallocate(ovlp)
   deallocate(dm)
   deallocate(edm)
   deallocate(row_ind)
   deallocate(col_ptr)

   call MPI_Comm_free(comm_in_task,ierr)
   call MPI_Comm_free(comm_among_task,ierr)

end subroutine
