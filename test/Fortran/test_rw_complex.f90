! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests reading and writing matrices.
!!
subroutine test_rw_complex(mpi_comm,h_file,s_file)

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   integer(kind=i4), intent(in) :: mpi_comm
   character(*),     intent(in) :: h_file
   character(*),     intent(in) :: s_file

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: myid
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: matrix_size
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l

   real(kind=r8) :: n_electrons
   real(kind=r8) :: err
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)

   complex(kind=r8), allocatable :: ham(:,:)
   complex(kind=r8), allocatable :: ham_save(:,:)
   complex(kind=r8), allocatable :: ovlp(:,:)
   complex(kind=r8), allocatable :: ovlp_save(:,:)
   complex(kind=r8), allocatable :: ham_csc(:)
   complex(kind=r8), allocatable :: ham_csc_save(:)
   complex(kind=r8), allocatable :: ovlp_csc(:)
   complex(kind=r8), allocatable :: ovlp_csc_save(:)

   type(elsi_rw_handle) :: rw_h

   real(kind=r8), parameter :: tol = 1.0e-20_r8

   call MPI_Comm_size(mpi_comm,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm,myid,mpierr)

   if(myid == 0) then
      write(*,'("  ################################")')
      write(*,'("  ##     ELSI TEST PROGRAMS     ##")')
      write(*,'("  ################################")')
      write(*,*)
   endif

   ! Set up square-like processor grid
   do npcol = nint(sqrt(real(n_proc))),2,-1
      if(mod(n_proc,npcol) == 0) exit
   enddo
   nprow = n_proc/npcol

   ! Set block size
   blk = 32

   ! Set up BLACS
   blacs_ctxt = mpi_comm
   call BLACS_Gridinit(blacs_ctxt,'r',nprow,npcol)

   ! Read H and S matrices
   if(n_proc == 1) then
      ! Test SINGLE_PROC mode
      call elsi_init_rw(rw_h,0,0,0,0.0_r8)
   else
      ! Test MULTI_PROC mode
      call elsi_init_rw(rw_h,0,1,0,0.0_r8)
      call elsi_set_rw_mpi(rw_h,mpi_comm)
      call elsi_set_rw_blacs(rw_h,blacs_ctxt,blk)
   endif

   call elsi_read_mat_dim(rw_h,h_file,n_electrons,matrix_size,l_rows,l_cols)

   allocate(ham(l_rows,l_cols))
   allocate(ham_save(l_rows,l_cols))
   allocate(ovlp(l_rows,l_cols))
   allocate(ovlp_save(l_rows,l_cols))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex(rw_h,h_file,ham)
   call elsi_read_mat_complex(rw_h,s_file,ovlp)

   call elsi_finalize_rw(rw_h)

   ham_save  = ham
   ovlp_save = ovlp

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished reading H and S matrices")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   if(n_proc == 1) then
      ! Test SINGLE_PROC mode
      call elsi_init_rw(rw_h,1,0,matrix_size,n_electrons)
   else
      ! Test MULTI_PROC mode
      call elsi_init_rw(rw_h,1,1,matrix_size,n_electrons)
      call elsi_set_rw_mpi(rw_h,mpi_comm)
      call elsi_set_rw_blacs(rw_h,blacs_ctxt,blk)
   endif

   call elsi_set_rw_zero_def(rw_h,tol)

   call elsi_write_mat_complex(rw_h,"H_complex.tmp",ham)
   call elsi_write_mat_complex(rw_h,"S_complex.tmp",ovlp)

   call elsi_finalize_rw(rw_h)

   t1 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished writing H and S matrices")')
      write(*,'("  | Time :",F10.3,"s")') t1-t2
      write(*,*)
   endif

   ! Read H and S matrices
   if(n_proc == 1) then
      ! Test SINGLE_PROC mode
      call elsi_init_rw(rw_h,0,0,0,0.0_r8)
   else
      ! Test MULTI_PROC mode
      call elsi_init_rw(rw_h,0,1,0,0.0_r8)
      call elsi_set_rw_mpi(rw_h,mpi_comm)
      call elsi_set_rw_blacs(rw_h,blacs_ctxt,blk)
   endif

   call elsi_read_mat_dim(rw_h,"H_complex.tmp",n_electrons,matrix_size,l_rows,&
           l_cols)

   call elsi_read_mat_complex(rw_h,"H_complex.tmp",ham)
   call elsi_read_mat_complex(rw_h,"S_complex.tmp",ovlp)

   call elsi_finalize_rw(rw_h)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished reading H and S matrices")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   err = max(maxval(abs(ham-ham_save)),maxval(abs(ovlp-ovlp_save)))

   if(myid == 0) then
      if(err <= tol) then
         write(*,'("  Passed.")')
      endif
      write(*,*)
   endif

   deallocate(ham)
   deallocate(ham_save)
   deallocate(ovlp)
   deallocate(ovlp_save)

   ! Read H and S matrices
   call elsi_init_rw(rw_h,0,1,0,0.0_r8)
   call elsi_set_rw_mpi(rw_h,mpi_comm)

   call elsi_read_mat_dim_sparse(rw_h,h_file,n_electrons,matrix_size,nnz_g,&
           nnz_l,l_cols)

   allocate(ham_csc(nnz_l))
   allocate(ham_csc_save(nnz_l))
   allocate(ovlp_csc(nnz_l))
   allocate(ovlp_csc_save(nnz_l))
   allocate(row_ind(nnz_l))
   allocate(col_ptr(l_cols+1))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex_sparse(rw_h,h_file,row_ind,col_ptr,ham_csc)
   call elsi_read_mat_complex_sparse(rw_h,s_file,row_ind,col_ptr,ovlp_csc)

   call elsi_finalize_rw(rw_h)

   ham_csc_save  = ham_csc
   ovlp_csc_save = ovlp_csc

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished reading H and S matrices")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ! Test MULTI_PROC mode
   call elsi_init_rw(rw_h,1,1,matrix_size,n_electrons)
   call elsi_set_rw_mpi(rw_h,mpi_comm)
   call elsi_set_rw_csc(rw_h,nnz_g,nnz_l,l_cols)

   call elsi_write_mat_complex_sparse(rw_h,"H_complex.tmp",row_ind,col_ptr,ham_csc)
   call elsi_write_mat_complex_sparse(rw_h,"S_complex.tmp",row_ind,col_ptr,ovlp_csc)

   call elsi_finalize_rw(rw_h)

   t1 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished writing H and S matrices")')
      write(*,'("  | Time :",F10.3,"s")') t1-t2
      write(*,*)
   endif

   ! Read H and S matrices
   call elsi_init_rw(rw_h,0,1,0,0.0_r8)
   call elsi_set_rw_mpi(rw_h,mpi_comm)

   call elsi_read_mat_dim_sparse(rw_h,"H_complex.tmp",n_electrons,matrix_size,&
           nnz_g,nnz_l,l_cols)

   call elsi_read_mat_complex_sparse(rw_h,"H_complex.tmp",row_ind,col_ptr,ham_csc)
   call elsi_read_mat_complex_sparse(rw_h,"S_complex.tmp",row_ind,col_ptr,ovlp_csc)

   call elsi_finalize_rw(rw_h)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished reading H and S matrices")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   err = max(maxval(abs(ham_csc-ham_csc_save)),&
            maxval(abs(ovlp_csc-ovlp_csc_save)))

   if(myid == 0) then
      if(err <= tol) then
         write(*,'("  Passed.")')
      endif
      write(*,*)
   endif

   deallocate(ham_csc)
   deallocate(ham_csc_save)
   deallocate(ovlp_csc)
   deallocate(ovlp_csc_save)
   deallocate(row_ind)
   deallocate(col_ptr)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)

end subroutine
