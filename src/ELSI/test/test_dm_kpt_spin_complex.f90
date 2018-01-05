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
!! This program tests elsi_dm_complex with 2 k-points and 2 spin channels.
!!
program test_dm_kpt_spin_complex

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
   integer(kind=i4) :: my_spin
   integer(kind=i4) :: my_kpt
   integer(kind=i4) :: my_group
   integer(kind=i4) :: myid_in_group
   integer(kind=i4) :: n_group
   integer(kind=i4) :: group_size
   integer(kind=i4) :: mpi_comm_global
   integer(kind=i4) :: mpi_comm_group
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: n_states
   integer(kind=i4) :: matrix_size
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols
   integer(kind=i4) :: solver

   real(kind=r8) :: n_electrons
   real(kind=r8) :: e_test = 0.0_r8
   real(kind=r8) :: e_ref = 0.0_r8
   real(kind=r8) :: e_tol = 0.0_r8
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   logical :: make_check = .false. ! Are we running "make check"?

   complex(kind=r8), allocatable :: ham(:,:)
   complex(kind=r8), allocatable :: ham_save(:,:)
   complex(kind=r8), allocatable :: ovlp(:,:)
   complex(kind=r8), allocatable :: dm(:,:)
   complex(kind=r8), allocatable :: edm(:,:)

   type(elsi_handle)    :: e_h
   type(elsi_rw_handle) :: rw_h

   integer(kind=i4), parameter :: n_spin = 2
   integer(kind=i4), parameter :: n_kpt  = 2

   real(kind=r8), parameter :: k_weights(2) = 0.5_r8

   ! VY: Reference values from calculations on November 20, 2017.
   real(kind=r8), parameter :: e_elpa  = -2622.88214509316_r8
   real(kind=r8), parameter :: e_omm   = -2622.88214509316_r8
   real(kind=r8), parameter :: e_pexsi = -2622.88194292325_r8

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
         write(*,'("  ##         (ELPA = 1; libOMM = 2; PEXSI = 3)  ##")')
         write(*,'("  ##  Arg#2: H matrix file.                     ##")')
         write(*,'("  ##  Arg#3: S matrix file.                     ##")')
         write(*,'("  ################################################")')
         call MPI_Abort(mpi_comm_global,0,mpierr)
         stop
      endif
   endif

   ! Require at least 4 MPI tasks
   if(mod(n_proc,4) /= 0) then
      if(myid == 0) then
         write(*,'("  #########################################################")')
         write(*,'("  ##  Number of MPI tasks needs to be a multiple of 4!!  ##")')
         write(*,'("  #########################################################")')
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
         write(*,'("  2) Transforms the generalized eigenproblem to the standard")')
         write(*,'("     form by using Cholesky factorization;")')
         write(*,'("  3) Solves the standard eigenproblem;")')
         write(*,'("  4) Back-transforms the eigenvectors to the generalized problem;")')
         write(*,'("  5) Constructs the density matrix from the eigen-solutions.")')
         write(*,*)
         write(*,'("  Now start testing  elsi_dm_complex + ELPA")')
         e_ref = e_elpa
      elseif(solver == 2) then
         write(*,'("  This test program performs the following computational steps:")')
         write(*,*)
         write(*,'("  1) Reads Hamiltonian and overlap matrices;")')
         write(*,'("  2) Computes the density matrix with orbital minimization method.")')
         write(*,*)
         write(*,'("  Now start testing  elsi_dm_complex + libOMM")')
         e_ref = e_omm
      elseif(solver == 3) then
         write(*,'("  This test program performs the following computational steps:")')
         write(*,*)
         write(*,'("  1) Reads Hamiltonian and overlap matrices;")')
         write(*,'("  2) Converts the matrices to 1D block CSC format;")')
         write(*,'("  3) Computes the density matrix with pole expansion and selected")')
         write(*,'("     inversion method.")')
         write(*,*)
         write(*,'("  Now start testing  elsi_dm_complex + PEXSI")')
         e_ref = e_pexsi
         e_tol = 1.0e-4_r8
      endif
      write(*,*)
   endif

   ! Create groups of MPI tasks for k-points and spin channels. Note: In this
   ! example, the number of MPI tasks is forced to be a multiple of 4 for
   ! simplicity. Therefore, the 2 k-points and 2 spin channels can be evenly
   ! distributed among tasks. In practice, the number of MPI tasks assigned to
   ! different k-points and spin channels can be different.
   n_group       = n_spin*n_kpt
   group_size    = n_proc/(n_spin*n_kpt)
   my_group      = myid/group_size
   myid_in_group = mod(myid,group_size)
   my_spin       = 1+my_group/n_spin
   my_kpt        = 1+mod(my_group,n_kpt)

   call MPI_Comm_split(mpi_comm_global,my_group,myid_in_group,mpi_comm_group,&
           mpierr)

   write(*,'("  Task ",I4," solving spin channel ",I2," and k-point",I2)') &
      myid,my_spin,my_kpt

   ! Set up square-like processor grid within each group
   do npcol = nint(sqrt(real(group_size))),2,-1
      if(mod(group_size,npcol) == 0) exit
   enddo
   nprow = group_size/npcol

   ! Set block size
   blk = 16

   ! Set up BLACS
   blacs_ctxt = mpi_comm_group
   call BLACS_Gridinit(blacs_ctxt,'r',nprow,npcol)

   ! Read H and S matrices
   ! H and S are the same for all k-points and spin channels in this example
   call elsi_init_rw(rw_h,0,1,1,0,0.0_r8)
   call elsi_set_rw_mpi(rw_h,mpi_comm_group)
   call elsi_set_rw_blacs(rw_h,blacs_ctxt,blk)

   call elsi_read_mat_dim(rw_h,arg2,n_electrons,matrix_size,l_rows,l_cols)

   allocate(ham(l_rows,l_cols))
   allocate(ham_save(l_rows,l_cols))
   allocate(ovlp(l_rows,l_cols))
   allocate(dm(l_rows,l_cols))
   allocate(edm(l_rows,l_cols))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex(rw_h,arg2,ham)
   call elsi_read_mat_complex(rw_h,arg3,ovlp)

   call elsi_finalize_rw(rw_h)

   ham_save = ham

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,*)
      write(*,'("  Finished reading H and S matrices")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ! Initialize ELSI
   n_states = int(n_electrons,kind=i4)

   call elsi_init(e_h,solver,1,0,matrix_size,n_electrons,n_states)
   call elsi_set_mpi(e_h,mpi_comm_group)
   call elsi_set_blacs(e_h,blacs_ctxt,blk)
   ! Required for spin/kpt calculations
   call elsi_set_mpi_global(e_h,mpi_comm_global)
   call elsi_set_kpoint(e_h,n_kpt,my_kpt,k_weights(my_kpt))
   call elsi_set_spin(e_h,n_spin,my_spin)

   ! Customize ELSI
   if(my_group == 0) then ! Let one group talk
      call elsi_set_output(e_h,2)
   endif
   call elsi_set_sing_check(e_h,0)
   call elsi_set_mu_broaden_width(e_h,1.0e-6_r8)
   call elsi_set_omm_n_elpa(e_h,1)
   call elsi_set_pexsi_delta_e(e_h,80.0_r8)
   call elsi_set_pexsi_np_per_pole(e_h,2)
   call elsi_set_solver_timings_unit(e_h,67)
   call elsi_set_solver_timings_file(e_h,"dm_kpt_spin_complex_timings.json")

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 1)
   call elsi_set_solver_timing_tag(e_h,"TEST")
   call elsi_dm_complex(e_h,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished SCF #1")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ham = ham_save

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 2, with the same H)
   call elsi_set_solver_timing_tag(e_h,"TEST")
   call elsi_dm_complex(e_h,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   ! Compute energy density matrix
   if(solver == 1 .or. solver == 2 .or. solver == 3) then
      call elsi_get_edm_complex(e_h,edm)
   endif

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
   deallocate(dm)
   deallocate(edm)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)
   call MPI_Comm_free(mpi_comm_group,mpierr)
   call MPI_Finalize(mpierr)

end program
