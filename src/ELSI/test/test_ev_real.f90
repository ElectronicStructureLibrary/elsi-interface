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
program test_ev_real

   use ELSI_PRECISION, only: r8,i4
   use ELSI
   use MatrixSwitch ! Only for test matrices generation

   implicit none

   include "mpif.h"

   character(5) :: m_storage
   character(3) :: m_operation
   character(128) :: arg1
   character(128) :: arg2

   integer(kind=i4) :: n_proc,nprow,npcol,myid,myprow,mypcol
   integer(kind=i4) :: mpi_comm_global,mpierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: BLACS_CTXT
   integer(kind=i4) :: n_basis,n_states
   integer(kind=i4) :: matrix_size,l_rows,l_cols,supercell(3)
   integer(kind=i4) :: solver

   real(kind=r8) :: n_electrons,frac_occ,sparsity,orb_r_cut
   real(kind=r8) :: k_point(3)
   real(kind=r8) :: e_test,e_ref
   real(kind=r8) :: t1,t2
   real(kind=r8), allocatable :: ham(:,:),ovlp(:,:),e_val(:),e_vec(:,:)

   type(matrix)      :: H,S
   type(elsi_handle) :: elsi_h

   ! VY: Reference value from calculations on May 7, 2017.
   real(kind=r8), parameter :: e_elpa = -126.817462901838_r8
   real(kind=r8), parameter :: e_tol  = 1e-10_r8

   integer(kind=i4), external :: numroc

   ! Initialize MPI
   call MPI_Init(mpierr)
   mpi_comm_global = MPI_COMM_WORLD
   call MPI_Comm_size(mpi_comm_global,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm_global,myid,mpierr)

   ! Read command line arguments
   if(COMMAND_ARGUMENT_COUNT() == 2) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
      read(arg2,*) solver
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
         write(*,'("  ##  Arg#1: Path to Tomato seed folder.        ##")')
         write(*,'("  ##  Arg#2: Choice of solver.                  ##")')
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
      if(solver == 1) then
         write(*,'("  This test program performs the following computational steps:")')
         write(*,*)
         write(*,'("  1) Generates Hamiltonian and overlap matrices;")')
         write(*,'("  2) Checks the singularity of the overlap matrix by computing")')
         write(*,'("     all its eigenvalues;")')
         write(*,'("  3) Transforms the generalized eigenproblem to the standard")')
         write(*,'("     form by using Cholesky factorization;")')
         write(*,'("  4) Solves the standard eigenproblem;")')
         write(*,'("  5) Back-transforms the eigenvectors to the generalized problem;")')
         write(*,*)
         write(*,'("  Now start testing  elsi_ev_real + ELPA")')
      elseif(solver == 5) then
         write(*,'("  This test program performs the following computational steps:")')
         write(*,*)
         write(*,'("  1) Generates Hamiltonian and overlap matrices;")')
         write(*,'("  2) Converts the matrices to 1D block distributed CSC format;")')
         write(*,'("  3) Solves the generalized eigenproblem with shift-and-invert")')
         write(*,'("     parallel spectral transformation.")')
         write(*,*)
         write(*,'("  Now start testing  elsi_ev_real + SIPs")')
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
   BLACS_CTXT = mpi_comm_global
   call BLACS_Gridinit(BLACS_CTXT,'r',nprow,npcol)
   call BLACS_Gridinfo(BLACS_CTXT,nprow,npcol,myprow,mypcol)

   call ms_scalapack_setup(mpi_comm_global,nprow,'r',blk,icontxt=BLACS_CTXT)

   ! Set parameters
   m_storage = 'pddbc'
   m_operation = 'lap'
   n_basis = 22
   supercell = (/3,3,3/)
   orb_r_cut = 0.5_r8
   k_point(1:3) = (/0.0_r8,0.0_r8,0.0_r8/)

   t1 = MPI_Wtime()

   ! Generate test matrices
   call tomato_TB(arg1,'silicon',.false.,frac_occ,n_basis,.false.,matrix_size,&
                  supercell,.false.,sparsity,orb_r_cut,n_states,.true.,k_point,&
                  .true.,0.0_r8,H,S,m_storage,.true.)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished test matrices generation")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   l_rows = numroc(matrix_size,blk,myprow,0,nprow)
   l_cols = numroc(matrix_size,blk,mypcol,0,npcol)

   allocate(ham(l_rows,l_cols))
   allocate(ovlp(l_rows,l_cols))
   allocate(e_vec(l_rows,l_cols))
   allocate(e_val(matrix_size))

   ! Initialize ELSI
   n_electrons = 2.0_r8*n_states

   call elsi_init(elsi_h,solver,1,0,matrix_size,n_electrons,n_states)
   call elsi_set_mpi(elsi_h,mpi_comm_global)
   call elsi_set_blacs(elsi_h,BLACS_CTXT,blk)

   ! Customize ELSI
   call elsi_set_output(elsi_h,2)

   ham = H%dval
   ovlp = S%dval

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 1)
   call elsi_ev_real(elsi_h,ham,ovlp,e_val,e_vec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished SCF #1")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ham = H%dval

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 2, with the same H)
   call elsi_ev_real(elsi_h,ham,ovlp,e_val,e_vec)

   t2 = MPI_Wtime()

   e_test = 2.0_r8*sum(e_val(1:n_states))

   if(myid == 0) then
      write(*,'("  Finished SCF #2")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
      write(*,'("  Finished test program")')
      write(*,*)
      if(abs(e_test-e_ref) < e_tol) then
         write(*,'("  Passed.")')
      else
         write(*,'("  Failed!!")')
      endif
      write(*,*)
   endif

   ! Finalize ELSI
   call elsi_finalize(elsi_h)

   call m_deallocate(H)
   call m_deallocate(S)
   deallocate(ham)
   deallocate(ovlp)
   deallocate(e_vec)
   deallocate(e_val)

   call MPI_Finalize(mpierr)

end program
