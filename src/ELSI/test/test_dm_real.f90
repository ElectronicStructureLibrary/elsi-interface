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
!! This program tests ELSI density matrix solver elsi_dm_real (multi-processor).
!!
program test_dm_real

   use ELSI
   use MatrixSwitch

   implicit none

   include "mpif.h"

   character(5) :: m_storage
   character(3) :: m_operation
   character(128) :: arg1
   character(128) :: arg2

   integer :: n_proc,nprow,npcol,myid
   integer :: mpi_comm_global,mpi_comm_row,mpi_comm_col,mpierr
   integer :: blk
   integer :: BLACS_CTXT
   integer :: n_basis,n_states
   integer :: matrix_size,supercell(3)
   integer :: solver

   real*8 :: n_electrons,frac_occ,sparsity,orb_r_cut
   real*8 :: k_point(3)
   real*8 :: e_test,e_ref,e_tol
   real*8 :: t1,t2

   ! VY: Reference values from calculations on Feb 24, 2016
   !     Note that PEXSI result is incorrect, since only 2 PEXSI
   !     poles are used for this quick test. Accurate result can
   !     be expected with at least 40 poles.
   real*8, parameter :: e_elpa  = -126.817462901838d0
   real*8, parameter :: e_omm   = -126.817462499630d0
   real*8, parameter :: e_pexsi = -1885.46836305848d0

   type(matrix) :: H,S,D

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
      if(solver < 1 .or. solver > 3) then
         if(myid == 0) then
            write(*,'("  ################################################")')
            write(*,'("  ##  Wrong choice of solver!!                  ##")')
            write(*,'("  ##  Please choose:                            ##")')
            write(*,'("  ##  ELPA = 1; libOMM = 2; PEXSI = 3           ##")')
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
         write(*,'("  ##         (ELPA = 1; libOMM = 2; PEXSI = 3)  ##")')
         write(*,'("  ################################################")')
         call MPI_Abort(mpi_comm_global,0,mpierr)
         stop
      endif
   endif

   ! Set up square-like grid
   do npcol = NINT(SQRT(REAL(n_proc))),2,-1
      if(MOD(n_proc,npcol) == 0) exit
   enddo
   nprow = n_proc/npcol

   ! Set block size
   blk = 128

   if(myid == 0) then
      write(*,'("  ################################")')
      write(*,'("  ##     ELSI TEST PROGRAMS     ##")')
      write(*,'("  ################################")')
      write(*,*)
      write(*,'("  This test program computes the density matrix")')
      write(*,'("  from a Hamiltonian matrix and an overlap matrix")')
      write(*,'("  generated on-the-fly.")')
      write(*,*)
      if(solver == 1) then
         write(*,'("  Now testing  elsi_dm_real + ELPA")')
         e_ref = e_elpa
         e_tol = 1d-10
      elseif(solver == 2) then
         ! Note that the energy tolerance for libOMM is larger, due to the
         ! random initial guess used in the code. In production calculations,
         ! a much higher accuracy can be obtained by using ELPA eigenvectors
         ! as the initial guess, or by completing the SCF cycle to make the
         ! influence of the random initial guess fade away.
         write(*,'("  Now testing  elsi_dm_real + libOMM")')
         e_ref = e_omm
         e_tol = 1d-6
      else
         write(*,'("  Now testing  elsi_dm_real + PEXSI")')
         e_ref = e_pexsi
         e_tol = 1d-10
      endif
      write(*,*)
   endif

   ! Set up BLACS
   BLACS_CTXT = mpi_comm_global
   call BLACS_Gridinit(BLACS_CTXT,'r',nprow,npcol)
   call ms_scalapack_setup(mpi_comm_global,nprow,'r',blk,icontxt=BLACS_CTXT)

   ! Set parameters
   m_storage = 'pddbc'
   m_operation = 'lap'
   n_basis = 22
   supercell = (/3,3,3/)
   orb_r_cut = 0.5d0
   k_point(1:3) = (/0d0,0d0,0d0/)

   t1 = MPI_Wtime()

   ! Generate test matrices
   call tomato_TB(arg1,'silicon',.false.,frac_occ,n_basis,.false.,matrix_size,&
                  supercell,.false.,sparsity,orb_r_cut,n_states,.true.,k_point,&
                  .true.,0d0,H,S,m_storage,.true.)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished test matrices generation")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
   endif

   ! Initialize ELSI
   n_electrons = 2d0*n_states

   call elsi_init(solver,1,0,matrix_size,n_electrons,n_states)
   call elsi_set_mpi(mpi_comm_global)
   call elsi_set_blacs(BLACS_CTXT,blk)

   ! Solve problem
   call m_allocate(D,matrix_size,matrix_size,m_storage)

   ! Disable ELPA if using libOMM
   if(solver == 2) call elsi_customize_omm(n_elpa_steps_omm=0)

   ! Only 2 PEXSI poles for quick test
   if(solver == 3) call elsi_customize_pexsi(n_poles=2)

   ! Get more output
   call elsi_customize(print_detail=.true.)

   t1 = MPI_Wtime()

   call elsi_dm_real(H%dval,S%dval,D%dval,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished test program")')
      write(*,'("  | Total computation time :",F10.3,"s")') t2-t1
      write(*,*)
      if(ABS(e_test-e_ref) < e_tol) then
         write(*,'("  Passed.")')
      else
         write(*,'("  Failed!!")')
      endif
      write(*,*)
   endif

   ! Finalize ELSI
   call elsi_finalize()

   call m_deallocate(H)
   call m_deallocate(S)
   call m_deallocate(D)

   call MPI_Finalize(mpierr)

end program
