! Copyright (c) 2015-2018, the ELSI team. All rights reserved.
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
!! This program tests reading and writing matrices.
!!
program test_rw_complex

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   character(128) :: arg1 ! H file
   character(128) :: arg2 ! S file

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: myid
   integer(kind=i4) :: mpi_comm_global
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

   complex(kind=r8), allocatable :: ham(:,:)
   complex(kind=r8), allocatable :: ham_save(:,:)
   complex(kind=r8), allocatable :: ovlp(:,:)
   complex(kind=r8), allocatable :: ovlp_save(:,:)
   complex(kind=r8), allocatable :: ham_csc(:)
   complex(kind=r8), allocatable :: ham_csc_save(:)
   complex(kind=r8), allocatable :: ovlp_csc(:)
   complex(kind=r8), allocatable :: ovlp_csc_save(:)
   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)

   type(elsi_rw_handle) :: rw_h

   real(kind=r8), parameter :: tol = 1.0e-20_r8

   ! Initialize MPI
   call MPI_Init(mpierr)
   mpi_comm_global = MPI_COMM_WORLD
   call MPI_Comm_size(mpi_comm_global,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm_global,myid,mpierr)

   ! Read command line arguments
   if(COMMAND_ARGUMENT_COUNT() == 2) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
   else
      if(myid == 0) then
         write(*,'("  ################################################")')
         write(*,'("  ##  Wrong number of command line arguments!!  ##")')
         write(*,'("  ##  Arg#1: H matrix file.                     ##")')
         write(*,'("  ##  Arg#2: S matrix file.                     ##")')
         write(*,'("  ################################################")')
         call MPI_Abort(mpi_comm_global,0,mpierr)
         stop
      endif
   endif

   if(myid == 0) then
      write(*,'("  ################################")')
      write(*,'("  ##     ELSI TEST PROGRAMS     ##")')
      write(*,'("  ################################")')
      write(*,*)
      write(*,'("  This test program repeats the following steps:")')
      write(*,*)
      write(*,'("  1) Reads Hamiltonian and overlap matrices;")')
      write(*,'("  2) Writes Hamiltonian and overlap matrices.")')
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
   blacs_ctxt = mpi_comm_global
   call BLACS_Gridinit(blacs_ctxt,'r',nprow,npcol)

   ! Read H and S matrices
   if(n_proc == 1) then
      ! Test SINGLE_PROC mode
      call elsi_init_rw(rw_h,0,0,1,0,0.0_r8)
   else
      ! Test MULTI_PROC mode
      call elsi_init_rw(rw_h,0,1,1,0,0.0_r8)
      call elsi_set_rw_mpi(rw_h,mpi_comm_global)
      call elsi_set_rw_blacs(rw_h,blacs_ctxt,blk)
   endif

   call elsi_read_mat_dim(rw_h,arg1,n_electrons,matrix_size,l_rows,l_cols)

   allocate(ham(l_rows,l_cols))
   allocate(ham_save(l_rows,l_cols))
   allocate(ovlp(l_rows,l_cols))
   allocate(ovlp_save(l_rows,l_cols))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex(rw_h,arg1,ham)
   call elsi_read_mat_complex(rw_h,arg2,ovlp)

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
      call elsi_init_rw(rw_h,1,0,1,matrix_size,n_electrons)
   else
      ! Test MULTI_PROC mode
      call elsi_init_rw(rw_h,1,1,1,matrix_size,n_electrons)
      call elsi_set_rw_mpi(rw_h,mpi_comm_global)
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
      call elsi_init_rw(rw_h,0,0,1,0,0.0_r8)
   else
      ! Test MULTI_PROC mode
      call elsi_init_rw(rw_h,0,1,1,0,0.0_r8)
      call elsi_set_rw_mpi(rw_h,mpi_comm_global)
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
      else
         write(*,'("  Failed!!")')
      endif
      write(*,*)
   endif

   deallocate(ham)
   deallocate(ham_save)
   deallocate(ovlp)
   deallocate(ovlp_save)

   ! Read H and S matrices
   call elsi_init_rw(rw_h,0,1,1,0,0.0_r8)
   call elsi_set_rw_mpi(rw_h,mpi_comm_global)

   call elsi_read_mat_dim_sparse(rw_h,arg1,n_electrons,matrix_size,nnz_g,&
           nnz_l,l_cols)

   allocate(ham_csc(nnz_l))
   allocate(ham_csc_save(nnz_l))
   allocate(ovlp_csc(nnz_l))
   allocate(ovlp_csc_save(nnz_l))
   allocate(row_ind(nnz_l))
   allocate(col_ptr(l_cols+1))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex_sparse(rw_h,arg1,row_ind,col_ptr,ham_csc)
   call elsi_read_mat_complex_sparse(rw_h,arg2,row_ind,col_ptr,ovlp_csc)

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
   call elsi_init_rw(rw_h,1,1,1,matrix_size,n_electrons)
   call elsi_set_rw_mpi(rw_h,mpi_comm_global)
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
   call elsi_init_rw(rw_h,0,1,1,0,0.0_r8)
   call elsi_set_rw_mpi(rw_h,mpi_comm_global)

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
      else
         write(*,'("  Failed!!")')
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
   call MPI_Finalize(mpierr)

end program
