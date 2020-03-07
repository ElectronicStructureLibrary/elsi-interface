! Copyright (c) 2015-2020, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This program tests real eigensolver for a standard eigenvalue problem,
!! BLACS_DENSE format.
!!
program test_standard_ev_real

   use ELSI
   use ELSI_MPI
   use ELSI_PRECISION, only: r8,i4

   implicit none

   character(len=100) :: arg1
   character(len=100) :: arg2
   character(len=100) :: arg3

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: myid
   integer(kind=i4) :: myprow
   integer(kind=i4) :: mypcol
   integer(kind=i4) :: comm
   integer(kind=i4) :: ierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: desc(9)
   integer(kind=i4) :: n_basis
   integer(kind=i4) :: n_states
   integer(kind=i4) :: solver
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols
   integer(kind=i4) :: ldm
   integer(kind=i4) :: n

   real(kind=r8) :: t1
   real(kind=r8) :: t2

   integer(kind=i4), allocatable :: seed(:)

   real(kind=r8), allocatable :: mat_a(:,:)
   real(kind=r8), allocatable :: mat_b(:,:)
   real(kind=r8), allocatable :: mat_tmp(:,:)
   real(kind=r8), allocatable :: evec(:,:)
   real(kind=r8), allocatable :: eval(:)

   integer(kind=i4), external :: numroc

   type(elsi_handle) :: eh

   ! Initialize MPI
   call MPI_Init(ierr)
   comm = MPI_COMM_WORLD
   call MPI_Comm_size(comm,n_proc,ierr)
   call MPI_Comm_rank(comm,myid,ierr)

   ! Read command line arguments
   if(command_argument_count() == 3) then
      call get_command_argument(1,arg1)
      call get_command_argument(2,arg2)
      call get_command_argument(3,arg3)

      read(arg1,*) n_basis
      if(n_basis <= 0) then
         n_basis = 1000
      end if

      read(arg2,*) n_states
      if(n_states < 0 .or. n_states > n_basis) then
         n_states = n_basis
      end if

      read(arg3,*) solver
   else
      if(myid == 0) then
         write(*,"(2X,A)") "################################################"
         write(*,"(2X,A)") "##  Wrong number of command line arguments!!  ##"
         write(*,"(2X,A)") "##  Arg#1: Size of test matrix.               ##"
         write(*,"(2X,A)") "##  Arg#2: Number of eigenvectors to compute. ##"
         write(*,"(2X,A)") "##  Arg#3: 1 = ELPA                           ##"
         write(*,"(2X,A)") "##         4 = EigenExa                       ##"
         write(*,"(2X,A)") "##         5 = SLEPc-SIPs                     ##"
         write(*,"(2X,A)") "################################################"
         call MPI_Abort(comm,0,ierr)
         stop
      end if
   end if

   if(myid == 0) then
      write(*,"(2X,A)") "################################"
      write(*,"(2X,A)") "##     ELSI TEST PROGRAMS     ##"
      write(*,"(2X,A)") "################################"
      write(*,*)
      if(solver == 1) then
         write(*,"(2X,A)") "Now start testing  elsi_ev_real + ELPA"
      else if(solver == 4) then
         write(*,"(2X,A)") "Now start testing  elsi_ev_real + EigenExa"
      else if(solver == 5) then
         write(*,"(2X,A)") "Now start testing  elsi_ev_real + SLEPc-SIPs"
      end if
      write(*,*)
   end if

   ! Set up square-like processor grid
   do npcol = nint(sqrt(real(n_proc))),2,-1
      if(mod(n_proc,npcol) == 0) exit
   end do
   nprow = n_proc/npcol

   ! Set block size
   blk = 32

   ! Set up BLACS
   blacs_ctxt = comm

   call BLACS_Gridinit(blacs_ctxt,'r',nprow,npcol)
   call BLACS_Gridinfo(blacs_ctxt,nprow,npcol,myprow,mypcol)

   l_rows = numroc(n_basis,blk,myprow,0,nprow)
   l_cols = numroc(n_basis,blk,mypcol,0,npcol)

   ldm = max(l_rows,1)

   call descinit(desc,n_basis,n_basis,blk,blk,0,0,blacs_ctxt,ldm,ierr)

   ! Generate a random matrix
   call random_seed(size=n)

   allocate(seed(n))
   allocate(mat_a(l_rows,l_cols))
   allocate(mat_tmp(l_rows,l_cols))

   seed(:) = myid

   call random_seed(put=seed)
   call random_number(mat_a)

   ! Symmetrize test matrix
   mat_tmp(:,:) = mat_a

   call pdtran(n_basis,n_basis,1.0_r8,mat_tmp,1,1,desc,1.0_r8,mat_a,1,1,desc)

   deallocate(seed)
   deallocate(mat_tmp)

   if(myid == 0) then
      write(*,"(2X,A)") "Finished generating test matrices"
      write(*,*)
   end if

   ! Initialize ELSI
   call elsi_init(eh,solver,1,0,n_basis,0.0_r8,n_states)
   call elsi_set_mpi(eh,comm)
   call elsi_set_blacs(eh,blacs_ctxt,blk)

   allocate(mat_b(1,1)) ! Dummy allocation
   allocate(evec(l_rows,l_cols))
   allocate(eval(n_basis))

   ! Customize ELSI
   call elsi_set_output(eh,2)
   call elsi_set_output_log(eh,1)
   call elsi_set_unit_ovlp(eh,1)

   t1 = MPI_Wtime()

   ! Solve problem
   call elsi_ev_real(eh,mat_a,mat_b,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished test program"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
      write(*,*)
   end if

   ! Finalize ELSI
   call elsi_finalize(eh)

   deallocate(mat_a)
   deallocate(mat_b)
   deallocate(eval)
   deallocate(evec)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)
   call MPI_Finalize(ierr)

end program
