! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests complex eigensolver, PEXSI_CSC format.
!!
subroutine test_ev_cmplx_csc1(mpi_comm,solver,h_file,s_file)

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
   integer(kind=i4) :: myprow
   integer(kind=i4) :: mypcol
   integer(kind=i4) :: ierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: n_states
   integer(kind=i4) :: n_basis
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: n_l_cols
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols
   integer(kind=i4) :: i
   integer(kind=i4) :: header(8)

   real(kind=r8) :: n_electrons
   real(kind=r8) :: mu
   real(kind=r8) :: weight(1)
   real(kind=r8) :: e_test = 0.0_r8
   real(kind=r8) :: e_ref = 0.0_r8
   real(kind=r8) :: tol = 0.0_r8
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   complex(kind=r8), allocatable :: ham(:)
   complex(kind=r8), allocatable :: ovlp(:)
   complex(kind=r8), allocatable :: evec(:,:)

   real(kind=r8), allocatable :: eval(:)
   real(kind=r8), allocatable :: occ(:)

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)

   type(elsi_handle) :: eh
   type(elsi_rw_handle) :: rwh

   ! Reference values
   real(kind=r8), parameter :: e_elpa = -2622.88214509316_r8

   integer(kind=i4), external :: numroc

   call MPI_Comm_size(mpi_comm,n_proc,ierr)
   call MPI_Comm_rank(mpi_comm,myid,ierr)

   if(myid == 0) then
      tol = 1.0e-8_r8
      write(*,"(2X,A)") "################################"
      write(*,"(2X,A)") "##     ELSI TEST PROGRAMS     ##"
      write(*,"(2X,A)") "################################"
      write(*,*)
      if(solver == 1) then
         write(*,"(2X,A)") "Now start testing  elsi_ev_complex_sparse + ELPA"
      end if
      write(*,*)
   end if

   e_ref = e_elpa

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
   call BLACS_Gridinfo(blacs_ctxt,nprow,npcol,myprow,mypcol)

   ! Read H and S matrices
   call elsi_init_rw(rwh,0,1,0,0.0_r8)
   call elsi_set_rw_mpi(rwh,mpi_comm)

   call elsi_read_mat_dim_sparse(rwh,h_file,n_electrons,n_basis,nnz_g,nnz_l,&
           n_l_cols)
   call elsi_get_rw_header(rwh,header)

   l_rows = numroc(n_basis,blk,myprow,0,nprow)
   l_cols = numroc(n_basis,blk,mypcol,0,npcol)

   allocate(ham(nnz_l))
   allocate(ovlp(nnz_l))
   allocate(row_ind(nnz_l))
   allocate(col_ptr(n_l_cols+1))
   allocate(evec(l_rows,l_cols))
   allocate(eval(n_basis))
   allocate(occ(n_basis))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex_sparse(rwh,h_file,row_ind,col_ptr,ham)
   call elsi_read_mat_complex_sparse(rwh,s_file,row_ind,col_ptr,ovlp)

   call elsi_finalize_rw(rwh)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished reading H and S matrices"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Initialize ELSI
   n_states = int(n_electrons,kind=i4)
   weight(1) = 1.0_r8

   call elsi_init(eh,solver,1,1,n_basis,n_electrons,n_states)
   call elsi_set_mpi(eh,mpi_comm)
   call elsi_set_csc(eh,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)
   call elsi_set_blacs(eh,blacs_ctxt,blk)

   ! Customize ELSI
   call elsi_set_output(eh,2)
   call elsi_set_output_log(eh,1)
   call elsi_set_mu_broaden_width(eh,1.0e-6_r8)

   t1 = MPI_Wtime()

   ! Solve
   call elsi_ev_complex_sparse(eh,ham,ovlp,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #1"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_ev_complex_sparse(eh,ham,ovlp,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #2"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_ev_complex_sparse(eh,ham,ovlp,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #3"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Reinit for a new geometry
   call elsi_reinit(eh)
   call elsi_set_csc(eh,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)
   call elsi_set_elpa_solver(eh,1)

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_ev_complex_sparse(eh,ham,ovlp,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #4"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_ev_complex_sparse(eh,ham,ovlp,eval,evec)

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
   end if

   ! Finalize ELSI
   call elsi_finalize(eh)

   deallocate(ham)
   deallocate(ovlp)
   deallocate(evec)
   deallocate(eval)
   deallocate(occ)
   deallocate(row_ind)
   deallocate(col_ptr)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)

end subroutine
