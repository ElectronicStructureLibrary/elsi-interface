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
!! This program tests elsi_dm_real_sparse.
!!
program test_dm_real_sparse

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   character(128) :: arg1 ! solver
   character(128) :: arg2 ! H file
   character(128) :: arg3 ! S file
   character(128) :: arg4 ! make check?

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: myid
   integer(kind=i4) :: mpi_comm_global
   integer(kind=i4) :: comm_in_task
   integer(kind=i4) :: comm_among_task
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: n_states
   integer(kind=i4) :: matrix_size
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: n_l_cols
   integer(kind=i4) :: solver
   integer(kind=i4) :: id_in_task
   integer(kind=i4) :: task_id
   integer(kind=i4) :: buffer(5)

   real(kind=r8) :: n_electrons
   real(kind=r8) :: e_test
   real(kind=r8) :: e_ref
   real(kind=r8) :: e_tol
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   logical :: make_check ! Are we running "make check"?

   real(kind=r8),    allocatable :: ham(:)
   real(kind=r8),    allocatable :: ham_save(:)
   real(kind=r8),    allocatable :: ovlp(:)
   real(kind=r8),    allocatable :: dm(:)
   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)

   type(elsi_handle) :: elsi_h

   ! VY: Reference values from calculations on August 31, 2017.
   real(kind=r8), parameter :: e_elpa  = -1833.07932666530_r8
   real(kind=r8), parameter :: e_omm   = -1833.07932666692_r8
   real(kind=r8), parameter :: e_pexsi = -1833.07836578201_r8

   ! Initialize MPI
   call MPI_Init(mpierr)
   mpi_comm_global = MPI_COMM_WORLD
   call MPI_Comm_size(mpi_comm_global,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm_global,myid,mpierr)

   ! Read command line arguments
   if(COMMAND_ARGUMENT_COUNT() == 4) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
      call GET_COMMAND_ARGUMENT(3,arg3)
      call GET_COMMAND_ARGUMENT(4,arg4)
      read(arg1,*) solver
      make_check = .true.
   elseif(COMMAND_ARGUMENT_COUNT() == 3) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
      call GET_COMMAND_ARGUMENT(3,arg3)
      read(arg1,*) solver
      make_check = .false.
   else
      if(myid == 0) then
         write(*,'("  ################################################")')
         write(*,'("  ##  Wrong number of command line arguments!!  ##")')
         write(*,'("  ##  Arg#1: Choice of solver.                  ##")')
         write(*,'("  ##         (ELPA = 1; libOMM = 2; PEXSI = 3)  ##")')
         write(*,'("  ##  Arg#2: H matrix file.                     ##")')
         write(*,'("  ##  Arg#3: S matrix file.                     ##")')
         write(*,'("  ################################################")')
         call MPI_Abort(mpi_comm_global,0,mpierr)
         stop
      endif
   endif

   if(myid == 0) then
      e_tol = 1.0e-8_r8
      write(*,'("  ################################")')
      write(*,'("  ##     ELSI TEST PROGRAMS     ##")')
      write(*,'("  ################################")')
      write(*,*)
      if(solver == 1) then
         write(*,'("  This test program performs the following computational steps:")')
         write(*,*)
         write(*,'("  1) Reads Hamiltonian and overlap matrices;")')
         write(*,'("  2) Converts the matrices to 2D block-cyclic dense format;")')
         write(*,'("  3) Transforms the generalized eigenproblem to the standard")')
         write(*,'("     form by using Cholesky factorization;")')
         write(*,'("  4) Solves the standard eigenproblem;")')
         write(*,'("  5) Back-transforms the eigenvectors to the generalized problem;")')
         write(*,'("  6) Constructs the density matrix from the eigen-solutions.")')
         write(*,*)
         write(*,'("  Now start testing  elsi_dm_real_sparse + ELPA")')
         e_ref = e_elpa
      elseif(solver == 2) then
         write(*,'("  This test program performs the following computational steps:")')
         write(*,*)
         write(*,'("  1) Reads Hamiltonian and overlap matrices;")')
         write(*,'("  2) Converts the matrices to 2D block-cyclic dense format;")')
         write(*,'("  3) Computes the Cholesky factorization of the overlap matrix;")')
         write(*,'("  4) Computes the density matrix with orbital minimization method.")')
         write(*,*)
         write(*,'("  Now start testing  elsi_dm_real_sparse + libOMM")')
         e_ref = e_omm
      else
         write(*,'("  This test program performs the following computational steps:")')
         write(*,*)
         write(*,'("  1) Reads Hamiltonian and overlap matrices;")')
         write(*,'("  2) Computes the density matrix with pole expansion and selected")')
         write(*,'("     inversion method.")')
         write(*,*)
         write(*,'("  Now start testing  elsi_dm_real_sparse + PEXSI")')
         e_ref = e_pexsi
         e_tol = 1.0e-4_r8
      endif
      write(*,*)
   endif

   t1 = MPI_Wtime()

   if((solver == 1) .or. (solver == 2)) then
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

      task_id = 0
   elseif(solver == 3) then
      task_id    = myid/2
      id_in_task = mod(myid,2)
   endif

   call MPI_Comm_split(mpi_comm_global,task_id,id_in_task,&
           comm_in_task,mpierr)
   call MPI_Comm_split(mpi_comm_global,id_in_task,task_id,&
           comm_among_task,mpierr)

   if(task_id == 0) then
      ! Read H and S matrices
      call elsi_read_mat_dim_sparse(arg2,comm_in_task,n_electrons,&
              matrix_size,nnz_g,nnz_l,n_l_cols)
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
   allocate(row_ind(nnz_l))
   allocate(col_ptr(n_l_cols+1))

   if(task_id == 0) then
      call elsi_read_mat_real_sparse(arg2,comm_in_task,matrix_size,nnz_g,&
              nnz_l,n_l_cols,row_ind,col_ptr,ham)

      call elsi_read_mat_real_sparse(arg3,comm_in_task,matrix_size,nnz_g,&
              nnz_l,n_l_cols,row_ind,col_ptr,ovlp)

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

   call elsi_init(elsi_h,solver,1,1,matrix_size,n_electrons,n_states)
   call elsi_set_mpi(elsi_h,mpi_comm_global)
   call elsi_set_csc(elsi_h,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)

   if((solver == 1) .or. (solver == 2)) then
      call elsi_set_blacs(elsi_h,blacs_ctxt,blk)
   endif

   ! Customize ELSI
   call elsi_set_output(elsi_h,2)
   call elsi_set_sing_check(elsi_h,0)
   call elsi_set_omm_n_elpa(elsi_h,1)
   call elsi_set_pexsi_np_per_pole(elsi_h,2)

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 1)
   call elsi_dm_real_sparse(elsi_h,ham,ovlp,dm,e_test)

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
   call elsi_dm_real_sparse(elsi_h,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished SCF #2")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
      write(*,'("  Finished test program")')
      write(*,*)
      if(make_check) then
         if(abs(e_test-e_ref) < e_tol) then
            write(*,'("  Passed.")')
         else
            write(*,'("  Failed!!")')
         endif
         write(*,*)
      endif
   endif

   ! Finalize ELSI
   call elsi_finalize(elsi_h)

   deallocate(ham)
   deallocate(ham_save)
   deallocate(ovlp)
   deallocate(dm)
   deallocate(row_ind)
   deallocate(col_ptr)

   call MPI_Finalize(mpierr)

end program
