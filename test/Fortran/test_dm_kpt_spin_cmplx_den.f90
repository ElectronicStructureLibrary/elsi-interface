! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This program tests complex density matrix solver, BLACS_DENSE format, with 2
!! k-points and 2 spin channels.
!!
program test_dm_kpt_spin_cmplx_den

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   character(len=100) :: arg1 ! solver
   character(len=999) :: arg2 ! H file
   character(len=999) :: arg3 ! S file

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
   integer(kind=i4) :: mpi_comm
   integer(kind=i4) :: mpi_comm_group
   integer(kind=i4) :: ierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: n_states
   integer(kind=i4) :: n_basis
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols
   integer(kind=i4) :: solver
   integer(kind=i4) :: header(8)

   real(kind=r8) :: n_electrons
   real(kind=r8) :: e_test = 0.0_r8
   real(kind=r8) :: e_ref = 0.0_r8
   real(kind=r8) :: tol = 0.0_r8
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   complex(kind=r8), allocatable :: ham(:,:)
   complex(kind=r8), allocatable :: ham_save(:,:)
   complex(kind=r8), allocatable :: ovlp(:,:)
   complex(kind=r8), allocatable :: dm(:,:)
   complex(kind=r8), allocatable :: edm(:,:)

   type(elsi_handle) :: eh
   type(elsi_rw_handle) :: rwh

   integer(kind=i4), parameter :: n_spin = 2
   integer(kind=i4), parameter :: n_kpt = 2

   real(kind=r8), parameter :: k_weights(2) = 0.5_r8

   ! Reference values
   real(kind=r8), parameter :: e_elpa = -2622.88214509316_r8
   real(kind=r8), parameter :: e_omm = -2622.88214509316_r8
   real(kind=r8), parameter :: e_pexsi = -2622.88194292325_r8
   real(kind=r8), parameter :: e_ntpoly = -2622.88214509311_r8

   ! Initialize MPI
   call MPI_Init(ierr)
   mpi_comm = MPI_COMM_WORLD
   call MPI_Comm_size(mpi_comm,n_proc,ierr)
   call MPI_Comm_rank(mpi_comm,myid,ierr)

   ! Read command line arguments
   if(COMMAND_ARGUMENT_COUNT() == 3) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
      call GET_COMMAND_ARGUMENT(3,arg3)
      read(arg1,*) solver
   else
      if(myid == 0) then
         write(*,"(2X,A)") "################################################"
         write(*,"(2X,A)") "##  Wrong number of command line arguments!!  ##"
         write(*,"(2X,A)") "##  Arg#1: 1 = ELPA                           ##"
         write(*,"(2X,A)") "##         2 = libOMM                         ##"
         write(*,"(2X,A)") "##         3 = PEXSI                          ##"
         write(*,"(2X,A)") "##         6 = NTPoly                         ##"
         write(*,"(2X,A)") "##  Arg#2: H matrix file                      ##"
         write(*,"(2X,A)") "##  Arg#3: S matrix file                      ##"
         write(*,"(2X,A)") "################################################"
         call MPI_Abort(mpi_comm,0,ierr)
         stop
      end if
   end if

   ! Require at least 4 MPI tasks
   if(mod(n_proc,4) /= 0) then
      if(myid == 0) then
         write(*,"(2X,A)") "#########################################################"
         write(*,"(2X,A)") "##  Number of MPI tasks needs to be a multiple of 4!!  ##"
         write(*,"(2X,A)") "#########################################################"
         call MPI_Abort(mpi_comm,0,ierr)
         stop
      end if
   end if

   if(myid == 0) then
      tol = 1.0e-8_r8
      write(*,"(2X,A)") "################################"
      write(*,"(2X,A)") "##     ELSI TEST PROGRAMS     ##"
      write(*,"(2X,A)") "################################"
      write(*,*)
      if(solver == 1) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + ELPA"
         e_ref = e_elpa
      else if(solver == 2) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + libOMM"
         e_ref = e_omm
      else if(solver == 3) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + PEXSI"
         e_ref = e_pexsi
         tol = 1.0e-3_r8
      else if(solver == 6) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + NTPoly"
         e_ref = e_ntpoly
         tol = 1.0e-7_r8
      end if
      write(*,*)
   end if

   ! Create groups of MPI tasks for k-points and spin channels. Note: In this
   ! example, the number of MPI tasks is forced to be a multiple of 4 for
   ! simplicity. Therefore, the 2 k-points and 2 spin channels can be evenly
   ! distributed among tasks. In practice, the number of MPI tasks assigned to
   ! different k-points and spin channels can be different.
   n_group = n_spin*n_kpt
   group_size = n_proc/(n_spin*n_kpt)
   my_group = myid/group_size
   myid_in_group = mod(myid,group_size)
   my_spin = 1+my_group/n_spin
   my_kpt = 1+mod(my_group,n_kpt)

   call MPI_Comm_split(mpi_comm,my_group,myid_in_group,mpi_comm_group,ierr)

   write(*,"(2X,A,I4,A,I2,A,I2)") "| Task",myid," solving spin channel",&
      my_spin," and k-point",my_kpt

   ! Set up square-like processor grid within each group
   do npcol = nint(sqrt(real(group_size))),2,-1
      if(mod(group_size,npcol) == 0) exit
   end do
   nprow = group_size/npcol

   ! Set block size
   blk = 16

   ! Set up BLACS
   blacs_ctxt = mpi_comm_group
   call BLACS_Gridinit(blacs_ctxt,'r',nprow,npcol)

   ! Read H and S matrices
   ! H and S are the same for all k-points and spin channels in this example
   call elsi_init_rw(rwh,0,1,0,0.0_r8)
   call elsi_set_rw_mpi(rwh,mpi_comm_group)
   call elsi_set_rw_blacs(rwh,blacs_ctxt,blk)

   call elsi_read_mat_dim(rwh,arg2,n_electrons,n_basis,l_rows,l_cols)
   call elsi_get_rw_header(rwh,header)

   allocate(ham(l_rows,l_cols))
   allocate(ham_save(l_rows,l_cols))
   allocate(ovlp(l_rows,l_cols))
   allocate(dm(l_rows,l_cols))
   allocate(edm(l_rows,l_cols))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex(rwh,arg2,ham)
   call elsi_read_mat_complex(rwh,arg3,ovlp)

   call elsi_finalize_rw(rwh)

   ham_save = ham

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,*)
      write(*,"(2X,A)") "Finished reading H and S matrices"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Initialize ELSI
   n_states = int(n_electrons,kind=i4)

   call elsi_init(eh,solver,1,0,n_basis,n_electrons,n_states)
   call elsi_set_mpi(eh,mpi_comm_group)
   call elsi_set_blacs(eh,blacs_ctxt,blk)
   ! Required for spin/kpt calculations
   call elsi_set_mpi_global(eh,mpi_comm)
   call elsi_set_kpoint(eh,n_kpt,my_kpt,k_weights(my_kpt))
   call elsi_set_spin(eh,n_spin,my_spin)

   ! Customize ELSI
   if(my_group == 0) then ! Let one group talk
      call elsi_set_output(eh,2)
   end if
   call elsi_set_sing_check(eh,0)
   call elsi_set_mu_broaden_width(eh,1.0e-6_r8)
   call elsi_set_omm_n_elpa(eh,1)
   call elsi_set_pexsi_delta_e(eh,80.0_r8)
   call elsi_set_pexsi_np_per_pole(eh,2)

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 1)
   call elsi_dm_complex(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #1"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ham = ham_save

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 2, with the same H)
   call elsi_dm_complex(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   ! Compute energy density matrix
   call elsi_get_edm_complex(eh,edm)

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #2"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
      write(*,"(2X,A)") "Finished test program"
      write(*,*)
      if(header(8) == 1111) then
         if(abs(e_test-e_ref) < tol) then
            write(*,"(2X,A)") "Passed."
         else
            write(*,"(2X,A)") "Failed."
         end if
      end if
      write(*,*)
   end if

   ! Finalize ELSI
   call elsi_finalize(eh)

   deallocate(ham)
   deallocate(ham_save)
   deallocate(ovlp)
   deallocate(dm)
   deallocate(edm)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)
   call MPI_Comm_free(mpi_comm_group,ierr)
   call MPI_Finalize(ierr)

end program
