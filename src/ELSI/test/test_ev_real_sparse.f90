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
!! This program tests elsi_ev_real_sparse.
!!
program test_ev_real_sparse

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
   integer(kind=i4) :: myprow
   integer(kind=i4) :: mypcol
   integer(kind=i4) :: mpi_comm_global
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: n_states
   integer(kind=i4) :: matrix_size
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: n_l_cols
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols
   integer(kind=i4) :: solver
   integer(kind=i4) :: i

   real(kind=r8) :: n_electrons
   real(kind=r8) :: mu
   real(kind=r8) :: weight(1)
   real(kind=r8) :: e_test = 0.0_r8
   real(kind=r8) :: e_ref = 0.0_r8
   real(kind=r8) :: e_tol = 0.0_r8
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   logical :: make_check = .false. ! Are we running "make check"?

   real(kind=r8),    allocatable :: ham(:)
   real(kind=r8),    allocatable :: ham_save(:)
   real(kind=r8),    allocatable :: ovlp(:)
   real(kind=r8),    allocatable :: evec(:,:)
   real(kind=r8),    allocatable :: eval(:)
   real(kind=r8),    allocatable :: occ(:)
   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)

   type(elsi_handle)    :: e_h
   type(elsi_rw_handle) :: rw_h

   ! VY: Reference value from calculations on August 31, 2017.
   real(kind=r8), parameter :: e_elpa = -1833.07932666530_r8

   integer(kind=i4), external :: numroc

   ! Initialize MPI
   call MPI_Init(mpierr)
   mpi_comm_global = MPI_COMM_WORLD
   call MPI_Comm_size(mpi_comm_global,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm_global,myid,mpierr)

   ! Read command line arguments
   if(COMMAND_ARGUMENT_COUNT() == 4) then
      make_check = .true.
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
      call GET_COMMAND_ARGUMENT(3,arg3)
      call GET_COMMAND_ARGUMENT(4,arg4)
      read(arg1,*) solver
   elseif(COMMAND_ARGUMENT_COUNT() == 3) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
      call GET_COMMAND_ARGUMENT(3,arg3)
      read(arg1,*) solver
   else
      if(myid == 0) then
         write(*,'("  ################################################")')
         write(*,'("  ##  Wrong number of command line arguments!!  ##")')
         write(*,'("  ##  Arg#1: Choice of solver.                  ##")')
         write(*,'("  ##         (ELPA = 1; SIPs = 5)               ##")')
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
         write(*,'("  3) Checks the singularity of the overlap matrix by computing")')
         write(*,'("     all its eigenvalues;")')
         write(*,'("  4) Transforms the generalized eigenproblem to the standard")')
         write(*,'("     form by using Cholesky factorization;")')
         write(*,'("  5) Solves the standard eigenproblem;")')
         write(*,'("  6) Back-transforms the eigenvectors to the generalized problem.")')
         write(*,*)
         write(*,'("  Now start testing  elsi_ev_real_sparse + ELPA")')
      elseif(solver == 5) then
         write(*,'("  This test program performs the following computational steps:")')
         write(*,*)
         write(*,'("  1) Reads Hamiltonian and overlap matrices;")')
         write(*,'("  2) Solves the generalized eigenproblem with shift-and-invert")')
         write(*,'("     parallel spectral transformation.")')
         write(*,*)
         write(*,'("  Now start testing  elsi_ev_real_sparse + SIPs")')
      endif
      write(*,*)
   endif

   e_ref = e_elpa

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
   call BLACS_Gridinfo(blacs_ctxt,nprow,npcol,myprow,mypcol)

   ! Read H and S matrices
   call elsi_init_rw(rw_h,0,1,1,0,0.0_r8)
   call elsi_set_rw_mpi(rw_h,mpi_comm_global)

   call elsi_read_mat_dim_sparse(rw_h,arg2,n_electrons,matrix_size,nnz_g,nnz_l,&
           n_l_cols)

   l_rows = numroc(matrix_size,blk,myprow,0,nprow)
   l_cols = numroc(matrix_size,blk,mypcol,0,npcol)

   allocate(ham(nnz_l))
   allocate(ham_save(nnz_l))
   allocate(ovlp(nnz_l))
   allocate(row_ind(nnz_l))
   allocate(col_ptr(n_l_cols+1))
   allocate(evec(l_rows,l_cols))
   allocate(eval(matrix_size))
   allocate(occ(matrix_size))

   t1 = MPI_Wtime()

   call elsi_read_mat_real_sparse(rw_h,arg2,row_ind,col_ptr,ham)
   call elsi_read_mat_real_sparse(rw_h,arg3,row_ind,col_ptr,ovlp)

   call elsi_finalize_rw(rw_h)

   ham_save = ham

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished reading H and S matrices")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ! Initialize ELSI
   n_states = int(n_electrons,kind=i4)
   weight(1) = 1.0_r8

   call elsi_init(e_h,solver,1,1,matrix_size,n_electrons,n_states)
   call elsi_set_mpi(e_h,mpi_comm_global)
   call elsi_set_csc(e_h,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)
   call elsi_set_blacs(e_h,blacs_ctxt,blk)

   ! Customize ELSI
   call elsi_set_output(e_h,2)

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 1)
   call elsi_ev_real_sparse(e_h,ham,ovlp,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished SCF #1")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ham = ham_save

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 2, with the same H)
   call elsi_ev_real_sparse(e_h,ham,ovlp,eval,evec)

   t2 = MPI_Wtime()

   call elsi_compute_mu_and_occ(e_h,n_electrons,n_states,1,1,weight,eval,occ,mu)

   e_test = 0.0_r8

   do i = 1,n_states
      e_test = e_test+eval(i)*occ(i)
   enddo

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
   call elsi_finalize(e_h)

   deallocate(ham)
   deallocate(ham_save)
   deallocate(ovlp)
   deallocate(evec)
   deallocate(eval)
   deallocate(occ)
   deallocate(row_ind)
   deallocate(col_ptr)

   call MPI_Finalize(mpierr)

end program
