! Copyright (c) 2015-2020, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests complex eigensolver, BLACS_DENSE format.
!!
subroutine test_ev_cmplx_den(mpi_comm,solver,h_file,s_file)

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   integer(kind=i4), intent(in) :: mpi_comm
   integer(kind=i4), intent(in) :: solver
   character(len=*), intent(in) :: h_file
   character(len=*), intent(in) :: s_file

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: myid
   integer(kind=i4) :: ierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: n_states
   integer(kind=i4) :: n_basis
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols
   integer(kind=i4) :: i
   integer(kind=i4) :: header(8)

   real(kind=r8) :: n_electrons
   real(kind=r8) :: mu
   real(kind=r8) :: weight(1)
   real(kind=r8) :: e_test = 0.0_r8
   real(kind=r8) :: tol = 0.0_r8
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   logical :: file_exist

   complex(kind=r8), allocatable :: ham(:,:)
   complex(kind=r8), allocatable :: ham_save(:,:)
   complex(kind=r8), allocatable :: ovlp(:,:)
   complex(kind=r8), allocatable :: ovlp_save(:,:)
   complex(kind=r8), allocatable :: evec(:,:)

   real(kind=r8), allocatable :: eval(:)
   real(kind=r8), allocatable :: occ(:)

   type(elsi_handle) :: eh
   type(elsi_rw_handle) :: rwh

   ! Reference values
   real(kind=r8), parameter :: e_ref = -2622.88214509316_r8

   character(len=*), parameter :: file_name = "elsi.in"

   call MPI_Comm_size(mpi_comm,n_proc,ierr)
   call MPI_Comm_rank(mpi_comm,myid,ierr)

   if(myid == 0) then
      tol = 1.0e-8_r8
      write(*,"(2X,A)") "################################"
      write(*,"(2X,A)") "##     ELSI TEST PROGRAMS     ##"
      write(*,"(2X,A)") "################################"
      write(*,*)
      if(solver == 1) then
         write(*,"(2X,A)") "Now start testing  elsi_ev_complex + ELPA"
      else if(solver == 7) then
         write(*,"(2X,A)") "Now start testing  elsi_ev_complex + MAGMA"
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
   blacs_ctxt = mpi_comm
   call BLACS_Gridinit(blacs_ctxt,'r',nprow,npcol)

   ! Read H and S matrices
   if(n_proc == 1) then
      ! Test SINGLE_PROC mode
      call elsi_init_rw(rwh,0,0,0,0.0_r8)
   else
      ! Test MULTI_PROC mode
      call elsi_init_rw(rwh,0,1,0,0.0_r8)
      call elsi_set_rw_mpi(rwh,mpi_comm)
      call elsi_set_rw_blacs(rwh,blacs_ctxt,blk)
   end if

   call elsi_read_mat_dim(rwh,h_file,n_electrons,n_basis,l_rows,l_cols)
   call elsi_get_rw_header(rwh,header)

   allocate(ham(l_rows,l_cols))
   allocate(ham_save(l_rows,l_cols))
   allocate(ovlp(l_rows,l_cols))
   allocate(ovlp_save(l_rows,l_cols))
   allocate(evec(l_rows,l_cols))
   allocate(eval(n_basis))
   allocate(occ(n_basis))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex(rwh,h_file,ham)
   call elsi_read_mat_complex(rwh,s_file,ovlp)

   call elsi_finalize_rw(rwh)

   ham_save = ham
   ovlp_save = ovlp

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished reading H and S matrices"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Initialize ELSI
   n_states = min(int(n_electrons,kind=i4),n_basis)
   weight(1) = 1.0_r8

   if(n_proc == 1) then
      ! Test SINGLE_PROC mode
      call elsi_init(eh,solver,0,0,n_basis,n_electrons,n_states)
   else
      ! Test MULTI_PROC mode
      call elsi_init(eh,solver,1,0,n_basis,n_electrons,n_states)
      call elsi_set_mpi(eh,mpi_comm)
      call elsi_set_blacs(eh,blacs_ctxt,blk)
   end if

   ! Customize ELSI
   call elsi_set_output(eh,2)
   call elsi_set_output_log(eh,1)
   call elsi_set_mu_broaden_width(eh,1.0e-6_r8)

   inquire(file=file_name,exist=file_exist)

   if(file_exist) then
      call elsi_set_input_file(eh,file_name)
   end if

   t1 = MPI_Wtime()

   ! Solve
   call elsi_ev_complex(eh,ham,ovlp,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #1"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ham = ham_save

   if(n_proc == 1) then
      ovlp = ovlp_save
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_ev_complex(eh,ham,ovlp,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #2"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ham = ham_save

   if(n_proc == 1) then
      ovlp = ovlp_save
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_ev_complex(eh,ham,ovlp,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #3"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Reinit for a new geometry
   call elsi_reinit(eh)
   call elsi_set_elpa_solver(eh,1)
   call elsi_set_magma_solver(eh,2)

   ham = ham_save
   ovlp = ovlp_save

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_ev_complex(eh,ham,ovlp,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #4"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ham = ham_save

   if(n_proc == 1) then
      ovlp = ovlp_save
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_ev_complex(eh,ham,ovlp,eval,evec)

   t2 = MPI_Wtime()

   call elsi_compute_mu_and_occ(eh,n_electrons,n_states,1,1,weight,eval,occ,mu)

   e_test = 0.0_r8

   do i = 1,n_states
      e_test = e_test+eval(i)*occ(i)
   end do

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #5"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
      write(*,"(2X,A)") "Finished test program"
      write(*,*)
      if(header(8) == 1111) then
         write(*,"(2X,A)") "Band energy"
         write(*,"(2X,A,F15.8)") "| This test :",e_test
         write(*,"(2X,A,F15.8)") "| Reference :",e_ref
         write(*,*)
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
   deallocate(ovlp_save)
   deallocate(evec)
   deallocate(eval)
   deallocate(occ)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)

end subroutine
