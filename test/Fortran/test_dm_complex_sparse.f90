! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests elsi_dm_complex_sparse.
!!
subroutine test_dm_complex_sparse(mpi_comm,solver,h_file,s_file)

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   integer(kind=i4), intent(in) :: mpi_comm
   integer(kind=i4), intent(in) :: solver
   character(*),     intent(in) :: h_file
   character(*),     intent(in) :: s_file

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: myid
   integer(kind=i4) :: comm_in_task
   integer(kind=i4) :: comm_among_task
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: n_states
   integer(kind=i4) :: matrix_size
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: n_l_cols
   integer(kind=i4) :: id_in_task
   integer(kind=i4) :: task_id
   integer(kind=i4) :: buffer(5)
   integer(kind=i4) :: header(8) = 0

   real(kind=r8) :: n_electrons
   real(kind=r8) :: e_test = 0.0_r8
   real(kind=r8) :: e_ref  = 0.0_r8
   real(kind=r8) :: e_tol  = 0.0_r8
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   complex(kind=r8), allocatable :: ham(:)
   complex(kind=r8), allocatable :: ham_save(:)
   complex(kind=r8), allocatable :: ovlp(:)
   complex(kind=r8), allocatable :: dm(:)
   complex(kind=r8), allocatable :: edm(:)
   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)

   type(elsi_handle)    :: e_h
   type(elsi_rw_handle) :: rw_h

   ! Reference values from calculations on November 20, 2017.
   real(kind=r8), parameter :: e_elpa  = -2622.88214509316_r8
   real(kind=r8), parameter :: e_omm   = -2622.88214509316_r8
   real(kind=r8), parameter :: e_pexsi = -2622.88194292325_r8

   call MPI_Comm_size(mpi_comm,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm,myid,mpierr)

   if(myid == 0) then
      e_tol = 1.0e-8_r8
      write(*,'("  ################################")')
      write(*,'("  ##     ELSI TEST PROGRAMS     ##")')
      write(*,'("  ################################")')
      write(*,*)
      if(solver == 1) then
         write(*,'("  Now start testing  elsi_dm_complex_sparse + ELPA")')
         e_ref = e_elpa
      elseif(solver == 2) then
         write(*,'("  Now start testing  elsi_dm_complex_sparse + libOMM")')
         e_ref = e_omm
      elseif(solver == 3) then
         write(*,'("  Now start testing  elsi_dm_complex_sparse + PEXSI")')
         e_ref = e_pexsi
         e_tol = 1.0e-4_r8
      endif
      write(*,*)
   endif

   if(solver == 3) then
      task_id    = myid/2
      id_in_task = mod(myid,2)
   else
      task_id    = 0
      id_in_task = myid
   endif

   call MPI_Comm_split(mpi_comm,task_id,id_in_task,comm_in_task,mpierr)
   call MPI_Comm_split(mpi_comm,id_in_task,task_id,comm_among_task,mpierr)

   if(task_id == 0) then
      ! Read H and S matrices
      call elsi_init_rw(rw_h,0,1,0,0.0_r8)
      call elsi_set_rw_mpi(rw_h,comm_in_task)

      call elsi_read_mat_dim_sparse(rw_h,h_file,n_electrons,matrix_size,nnz_g,&
              nnz_l,n_l_cols)
      call elsi_get_rw_header(rw_h,header)
   endif

   if(solver == 3) then
      if(task_id == 0) then
         buffer(1) = int(n_electrons,kind=i4)
         buffer(2) = matrix_size
         buffer(3) = nnz_g
         buffer(4) = nnz_l
         buffer(5) = n_l_cols
      endif

      call MPI_Bcast(buffer,5,mpi_integer4,0,comm_among_task,mpierr)

      n_electrons = real(buffer(1),kind=r8)
      matrix_size = buffer(2)
      nnz_g       = buffer(3)
      nnz_l       = buffer(4)
      n_l_cols    = buffer(5)
   endif

   allocate(ham(nnz_l))
   allocate(ham_save(nnz_l))
   allocate(ovlp(nnz_l))
   allocate(dm(nnz_l))
   allocate(edm(nnz_l))
   allocate(row_ind(nnz_l))
   allocate(col_ptr(n_l_cols+1))

   t1 = MPI_Wtime()

   if(task_id == 0) then
      call elsi_read_mat_complex_sparse(rw_h,h_file,row_ind,col_ptr,ham)
      call elsi_read_mat_complex_sparse(rw_h,s_file,row_ind,col_ptr,ovlp)

      call elsi_finalize_rw(rw_h)

      ham_save = ham
   endif

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished reading H and S matrices")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ! Initialize ELSI
   n_states = int(n_electrons,kind=i4)

   call elsi_init(e_h,solver,1,1,matrix_size,n_electrons,n_states)
   call elsi_set_mpi(e_h,mpi_comm)
   call elsi_set_csc(e_h,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)

   ! Customize ELSI
   call elsi_set_output(e_h,2)
   call elsi_set_output_log(e_h,1)
   call elsi_set_sing_check(e_h,0)
   call elsi_set_mu_broaden_width(e_h,1.0e-6_r8)
   call elsi_set_omm_n_elpa(e_h,1)
   call elsi_set_pexsi_delta_e(e_h,80.0_r8)
   call elsi_set_pexsi_np_per_pole(e_h,2)

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 1)
   call elsi_dm_complex_sparse(e_h,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished SCF #1")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   if(task_id == 0) then
      ham = ham_save
   endif

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 2, with the same H)
   call elsi_dm_complex_sparse(e_h,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   ! Compute energy density matrix
   if(solver == 1 .or. solver == 2 .or. solver == 3) then
      call elsi_get_edm_complex_sparse(e_h,edm)
   endif

   if(myid == 0) then
      write(*,'("  Finished SCF #2")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
      write(*,'("  Finished test program")')
      write(*,*)
      if(header(8) == 1111) then
         if(abs(e_test-e_ref) < e_tol) then
            write(*,'("  Passed.")')
         else
            write(*,'("  Failed.")')
         endif
      endif
      write(*,*)
   endif

   ! Finalize ELSI
   call elsi_finalize(e_h)

   deallocate(ham)
   deallocate(ham_save)
   deallocate(ovlp)
   deallocate(dm)
   deallocate(edm)
   deallocate(row_ind)
   deallocate(col_ptr)

   call MPI_Comm_free(comm_in_task,mpierr)
   call MPI_Comm_free(comm_among_task,mpierr)

end subroutine