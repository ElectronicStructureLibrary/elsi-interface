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
!! This program tests ELSI eigensolver.
!!
program test_standard_ev_real

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   character(128) :: arg1
   character(128) :: arg2
   character(128) :: arg3

   integer(kind=i4) :: n_proc,nprow,npcol,myid,myprow,mypcol
   integer(kind=i4) :: mpi_comm_global,mpierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: BLACS_CTXT
   integer(kind=i4) :: sc_desc(9)
   integer(kind=i4) :: info
   integer(kind=i4) :: matrix_size
   integer(kind=i4) :: n_states
   integer(kind=i4) :: solver
   integer(kind=i4) :: local_row,local_col,ldm
   integer(kind=i4) :: n
   integer(kind=i4), allocatable :: seed(:)
   integer(kind=i4), external :: numroc

   real(kind=r8) :: t1,t2
   real(kind=r8), allocatable :: mat_a(:,:),mat_b(:,:),mat_tmp(:,:),e_vec(:,:)
   real(kind=r8), allocatable :: e_val(:)

   type(elsi_handle) :: elsi_h

   ! Initialize MPI
   call MPI_Init(mpierr)
   mpi_comm_global = MPI_COMM_WORLD
   call MPI_Comm_size(mpi_comm_global,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm_global,myid,mpierr)

   ! Read command line arguments
   if(COMMAND_ARGUMENT_COUNT() == 3) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
      call GET_COMMAND_ARGUMENT(3,arg3)

      read(arg1,*) matrix_size
      if(matrix_size .le. 0) then
         matrix_size = 1000
      endif

      read(arg2,*) n_states
      if((n_states < 0) .or. n_states > matrix_size) then
         n_states = matrix_size
      endif

      read(arg3,*) solver
      if((solver .ne. 1) .and. (solver .ne. 5)) then
         if(myid == 0) then
            write(*,'("  ################################################")')
            write(*,'("  ##  Wrong choice of solver!!                  ##")')
            write(*,'("  ##  Please choose:                            ##")')
            write(*,'("  ##  ELPA = 1; SIPs = 5                        ##")')
            write(*,'("  ################################################")')
            call MPI_Abort(mpi_comm_global,0,mpierr)
            stop
         endif
      endif
   else
      if(myid == 0) then
         write(*,'("  ################################################")')
         write(*,'("  ##  Wrong number of command line arguments!!  ##")')
         write(*,'("  ##  Arg#1: Size of test matrix.               ##")')
         write(*,'("  ##  Arg#2: Number of eigenvectors to compute. ##")')
         write(*,'("  ##  Arg#3: Choice of solver.                  ##")')
         write(*,'("  ##         (ELPA = 1; SIPs = 5)               ##")')
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
      write(*,'("  This test program performs the following computational steps:")')
      write(*,*)
      write(*,'("  1) Generates a random square matrix;")')
      write(*,'("  2) Symmetrize the random matrix;")')
      write(*,'("  3) Computes its eigenvalues and eigenvectors.")')
      write(*,*)
      if(solver == 1) then
         write(*,'("  Now start testing  elsi_ev_real + ELPA")')
      elseif(solver == 5) then
         write(*,'("  Now start testing  elsi_ev_real + SIPs")')
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
   BLACS_CTXT = mpi_comm_global

   call BLACS_Gridinit(BLACS_CTXT,'r',nprow,npcol)
   call BLACS_Gridinfo(BLACS_CTXT,nprow,npcol,myprow,mypcol)

   local_row = numroc(matrix_size,blk,myprow,0,nprow)
   local_col = numroc(matrix_size,blk,mypcol,0,npcol)

   ldm = max(local_row,1)

   call descinit(sc_desc,matrix_size,matrix_size,blk,blk,&
                 0,0,BLACS_CTXT,ldm,info)

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
   endif

   ! Initialize ELSI
   call elsi_init(elsi_h,solver,1,0,matrix_size,0.0_r8,1,1,n_states)
   call elsi_set_mpi(elsi_h,mpi_comm_global,mpi_comm_global)
   call elsi_set_blacs(elsi_h,BLACS_CTXT,blk)

   allocate(mat_b(1,1)) ! Dummy allocation
   allocate(e_vec(local_row,local_col))
   allocate(e_val(matrix_size))

   ! Customize ELSI
   call elsi_set_output(elsi_h,2)
   call elsi_set_unit_ovlp(elsi_h,1)
   
   t1 = MPI_Wtime()

   ! Solve problem
   call elsi_ev_real(elsi_h,mat_a,mat_b,e_val,e_vec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished test program")')
      write(*,'("  | Total computation time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ! Finalize ELSI
   call elsi_finalize(elsi_h)

   deallocate(mat_a)
   deallocate(mat_b)
   deallocate(e_val)
   deallocate(e_vec)

   call MPI_Finalize(mpierr)

end program
