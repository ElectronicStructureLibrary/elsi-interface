! Copyright (c) 2015-2020, the ELSI team.
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
   real(kind=r8) :: n_sp
   real(kind=r8) :: tmp
   real(kind=r8) :: e_hp
   real(kind=r8) :: e_sq
   real(kind=r8) :: e_ref
   real(kind=r8) :: tol
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   complex(kind=r8), allocatable :: ham(:,:)
   complex(kind=r8), allocatable :: ham_save(:,:)
   complex(kind=r8), allocatable :: ovlp(:,:)
   complex(kind=r8), allocatable :: ovlp_save(:,:)
   complex(kind=r8), allocatable :: dm(:,:)
   complex(kind=r8), allocatable :: edm(:,:)

   complex(kind=r8), external :: zdotc

   type(elsi_handle) :: eh
   type(elsi_rw_handle) :: rwh

   integer(kind=i4), parameter :: n_spin = 2
   integer(kind=i4), parameter :: n_kpt = 2

   real(kind=r8), parameter :: k_weights(2) = 0.5_r8

   ! Reference values
   real(kind=r8), parameter :: e_elpa = -2622.88214509316_r8
   real(kind=r8), parameter :: e_pexsi = -2622.88194292325_r8

   ! Initialize MPI
   call MPI_Init(ierr)
   mpi_comm = MPI_COMM_WORLD
   call MPI_Comm_size(mpi_comm,n_proc,ierr)
   call MPI_Comm_rank(mpi_comm,myid,ierr)

   e_hp = 0.0_r8
   e_sq = 0.0_r8
   e_ref = e_elpa
   tol = 1.0e-8_r8
   header(:) = 0

   ! Read command line arguments
   if(command_argument_count() == 3) then
      call get_command_argument(1,arg1)
      call get_command_argument(2,arg2)
      call get_command_argument(3,arg3)
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

   ! Below is a requirement of this test (not a requirement of ELSI)
   if(solver == 3) then
      if(mod(n_proc,8) /= 0) then
         if(myid == 0) then
            write(*,"(2X,A)") "#############################################"//&
               "#########################"
            write(*,"(2X,A)") "##  This test requires number of MPI tasks to"//&
               " be a multiple of 8!!  ##"
            write(*,"(2X,A)") "#############################################"//&
               "#########################"
            call MPI_Abort(mpi_comm,0,ierr)
            stop
         end if
      end if
   else
      if(mod(n_proc,4) /= 0) then
         if(myid == 0) then
            write(*,"(2X,A)") "#############################################"//&
               "#########################"
            write(*,"(2X,A)") "##  This test requires number of MPI tasks to"//&
               " be a multiple of 4!!  ##"
            write(*,"(2X,A)") "#############################################"//&
               "#########################"
            call MPI_Abort(mpi_comm,0,ierr)
            stop
         end if
      end if
   end if

   if(myid == 0) then
      write(*,"(2X,A)") "################################"
      write(*,"(2X,A)") "##     ELSI TEST PROGRAMS     ##"
      write(*,"(2X,A)") "################################"
      write(*,*)
      if(solver == 1) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + ELPA"
      else if(solver == 2) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + libOMM"
      else if(solver == 3) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + PEXSI"
         e_ref = e_pexsi
         tol = 1.0e-3_r8
      else if(solver == 6) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + NTPoly"
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
   allocate(ovlp_save(l_rows,l_cols))
   allocate(dm(l_rows,l_cols))
   allocate(edm(l_rows,l_cols))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex(rwh,arg2,ham)
   call elsi_read_mat_complex(rwh,arg3,ovlp)

   call elsi_finalize_rw(rwh)

   ham_save(:,:) = ham
   ovlp_save(:,:) = ovlp

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,*)
      write(*,"(2X,A)") "Finished reading H and S matrices"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Initialize ELSI
   n_states = min(int(n_electrons,kind=i4),n_basis)

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
   call elsi_set_illcond_check(eh,0)
   call elsi_set_mu_broaden_width(eh,1.0e-6_r8)
   call elsi_set_omm_n_elpa(eh,1)
   call elsi_set_pexsi_delta_e(eh,80.0_r8)
!   call elsi_set_pexsi_np_per_pole(eh,2)

   t1 = MPI_Wtime()

   ! Solve
   call elsi_dm_complex(eh,ham,ovlp,dm,e_hp)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #1"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ham(:,:) = ham_save

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex(eh,ham,ovlp,dm,e_hp)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #2"
      write(*,*)
   end if

   ham(:,:) = ham_save

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex(eh,ham,ovlp,dm,e_hp)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #3"
      write(*,*)
   end if

   ! Reinit for a new geometry
   call elsi_reinit(eh)
   call elsi_set_elpa_solver(eh,1)
   call elsi_set_omm_flavor(eh,2)
   call elsi_set_ntpoly_method(eh,1)

   ham(:,:) = ham_save
   ovlp(:,:) = ovlp_save

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex(eh,ham,ovlp,dm,e_hp)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #4"
      write(*,*)
   end if

   ham(:,:) = ham_save

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex(eh,ham,ovlp,dm,e_hp)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #5"
      write(*,*)
   end if

   ! Compute energy density matrix
   call elsi_get_edm_complex(eh,edm)

   ! Compute band energy
   tmp = real(zdotc(l_rows*l_cols,ovlp_save,1,edm,1),kind=r8)*k_weights(my_kpt)

   call MPI_Reduce(tmp,e_sq,1,mpi_real8,mpi_sum,0,mpi_comm,ierr)

   ! Compute electron count
   tmp = real(zdotc(l_rows*l_cols,ovlp_save,1,dm,1),kind=r8)*k_weights(my_kpt)

   call MPI_Reduce(tmp,n_sp,1,mpi_real8,mpi_sum,0,mpi_comm,ierr)

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #5"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
      write(*,"(2X,A)") "Finished test program"
      write(*,*)
      if(header(8) == 1111) then
         write(*,"(2X,A)") "Band energy"
         write(*,"(2X,A,F15.8)") "| Tr[H*P]   :",e_hp
         write(*,"(2X,A,F15.8)") "| Tr[S*Q]   :",e_sq
         write(*,"(2X,A,F15.8)") "| Reference :",e_ref
         write(*,*)
         write(*,"(2X,A)") "Electron count"
         write(*,"(2X,A,F15.8)") "| Tr[S*P]   :",n_sp
         write(*,"(2X,A,F15.8)") "| Reference :",n_electrons
         write(*,*)
         if(abs(e_hp-e_ref) < tol .and. abs(e_sq-e_ref) < tol&
            .and. abs(n_sp-n_electrons) < tol) then
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
   deallocate(ovlp_save)
   deallocate(dm)
   deallocate(edm)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)
   call MPI_Comm_free(mpi_comm_group,ierr)
   call MPI_Finalize(ierr)

end program
