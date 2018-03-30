! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This program tests elsi_ev_real.
!!
program test_standard_ev_real

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   character(128) :: arg1
   character(128) :: arg2
   character(128) :: arg3

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: myid
   integer(kind=i4) :: myprow
   integer(kind=i4) :: mypcol
   integer(kind=i4) :: mpi_comm
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: sc_desc(9)
   integer(kind=i4) :: info
   integer(kind=i4) :: matrix_size
   integer(kind=i4) :: n_states
   integer(kind=i4) :: solver
   integer(kind=i4) :: local_row
   integer(kind=i4) :: local_col
   integer(kind=i4) :: ldm
   integer(kind=i4) :: n
   integer(kind=i4), allocatable :: seed(:)
   integer(kind=i4), external :: numroc

   real(kind=r8) :: t1
   real(kind=r8) :: t2
   real(kind=r8), allocatable :: mat_a(:,:)
   real(kind=r8), allocatable :: mat_b(:,:)
   real(kind=r8), allocatable :: mat_tmp(:,:)
   real(kind=r8), allocatable :: evec(:,:)
   real(kind=r8), allocatable :: eval(:)

   type(elsi_handle) :: e_h

   ! Initialize MPI
   call MPI_Init(mpierr)
   mpi_comm = MPI_COMM_WORLD
   call MPI_Comm_size(mpi_comm,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm,myid,mpierr)

   ! Read command line arguments
   if(COMMAND_ARGUMENT_COUNT() == 3) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
      call GET_COMMAND_ARGUMENT(3,arg3)

      read(arg1,*) matrix_size
      if(matrix_size <= 0) then
         matrix_size = 1000
      endif

      read(arg2,*) n_states
      if(n_states < 0 .or. n_states > matrix_size) then
         n_states = matrix_size
      endif

      read(arg3,*) solver
   else
      if(myid == 0) then
         write(*,'("  ################################################")')
         write(*,'("  ##  Wrong number of command line arguments!!  ##")')
         write(*,'("  ##  Arg#1: Size of test matrix.               ##")')
         write(*,'("  ##  Arg#2: Number of eigenvectors to compute. ##")')
         write(*,'("  ##  Arg#3: Choice of solver.                  ##")')
         write(*,'("  ##         (ELPA = 1; SIPS = 5)               ##")')
         write(*,'("  ################################################")')
         call MPI_Abort(mpi_comm,0,mpierr)
         stop
      endif
   endif

   if(myid == 0) then
      write(*,'("  ################################")')
      write(*,'("  ##     ELSI TEST PROGRAMS     ##")')
      write(*,'("  ################################")')
      write(*,*)
      if(solver == 1) then
         write(*,'("  Now start testing  elsi_ev_real + ELPA")')
      elseif(solver == 5) then
         write(*,'("  Now start testing  elsi_ev_real + SIPS")')
      endif
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
   call BLACS_Gridinfo(blacs_ctxt,nprow,npcol,myprow,mypcol)

   local_row = numroc(matrix_size,blk,myprow,0,nprow)
   local_col = numroc(matrix_size,blk,mypcol,0,npcol)

   ldm = max(local_row,1)

   call descinit(sc_desc,matrix_size,matrix_size,blk,blk,&
                 0,0,blacs_ctxt,ldm,info)

   ! Generate a random matrix
   call random_seed(size=n)

   allocate(seed(n))
   allocate(mat_a(local_row,local_col))
   allocate(mat_tmp(local_row,local_col))

   seed = myid

   call random_seed(put=seed)
   call random_number(mat_a)

   ! Symmetrize test matrix
   mat_tmp = mat_a

   call pdtran(matrix_size,matrix_size,1.0_r8,mat_tmp,1,1,&
               sc_desc,1.0_r8,mat_a,1,1,sc_desc)

   deallocate(mat_tmp)

   if(myid == 0) then
      write(*,'("  Finished test matrices generation")')
      write(*,*)
   endif

   ! Initialize ELSI
   call elsi_init(e_h,solver,1,0,matrix_size,0.0_r8,n_states)
   call elsi_set_mpi(e_h,mpi_comm)
   call elsi_set_blacs(e_h,blacs_ctxt,blk)

   allocate(mat_b(1,1)) ! Dummy allocation
   allocate(evec(local_row,local_col))
   allocate(eval(matrix_size))

   ! Customize ELSI
   call elsi_set_output(e_h,2)
   call elsi_set_output_log(e_h,1)
   call elsi_set_unit_ovlp(e_h,1)

   t1 = MPI_Wtime()

   ! Solve problem
   call elsi_ev_real(e_h,mat_a,mat_b,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished test program")')
      write(*,'("  | Total computation time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ! Finalize ELSI
   call elsi_finalize(e_h)

   deallocate(mat_a)
   deallocate(mat_b)
   deallocate(eval)
   deallocate(evec)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)
   call MPI_Finalize(mpierr)

end program
