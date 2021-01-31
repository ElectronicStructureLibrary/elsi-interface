! Copyright (c) 2015-2021, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests real BSE solver, BLACS_DENSE format.
!!
subroutine test_bse_real_den(comm,solver,a_file,b_file)

   use ELSI
   use ELSI_MPI
   use ELSI_PRECISION, only: r8,i4

   implicit none

   integer(kind=i4), intent(in) :: comm
   integer(kind=i4), intent(in) :: solver
   character(len=*), intent(in) :: a_file
   character(len=*), intent(in) :: b_file

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: myid
   integer(kind=i4) :: myprow
   integer(kind=i4) :: mypcol
   integer(kind=i4) :: ierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: n_basis
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols
   integer(kind=i4) :: l_rows2
   integer(kind=i4) :: l_cols2
   integer(kind=i4) :: header(8)

   real(kind=r8) :: n_electrons
   real(kind=r8) :: tol
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   logical :: file_exist

   real(kind=r8), allocatable :: mat_a(:,:)
   real(kind=r8), allocatable :: mat_b(:,:)
   real(kind=r8), allocatable :: evec(:,:)
   real(kind=r8), allocatable :: eval(:)

   type(elsi_handle) :: eh
   type(elsi_rw_handle) :: rwh

   integer(kind=i4), external :: numroc

   ! Reference values
   real(kind=r8), parameter :: e_ref = 0.31442451613_r8

   character(len=*), parameter :: file_name = "elsi.in"

   call MPI_Comm_size(comm,n_proc,ierr)
   call MPI_Comm_rank(comm,myid,ierr)

   tol = 1.0e-8_r8
   header(:) = 0

   if(myid == 0) then
      write(*,"(2X,A)") "################################"
      write(*,"(2X,A)") "##     ELSI TEST PROGRAMS     ##"
      write(*,"(2X,A)") "################################"
      write(*,*)
      if(solver == 8) then
         write(*,"(2X,A)") "Now start testing  elsi_bse_real + BSEPACK"
      end if
      write(*,*)
   end if

   ! Set up square-like processor grid
   do npcol = nint(sqrt(real(n_proc))),2,-1
      if(mod(n_proc,npcol) == 0) exit
   end do
   nprow = n_proc/npcol

   ! Set block size
   blk = 8

   ! Set up BLACS
   blacs_ctxt = comm
   call BLACS_Gridinit(blacs_ctxt,'r',nprow,npcol)
   call BLACS_Gridinfo(blacs_ctxt,nprow,npcol,myprow,mypcol)

   ! Read A and B matrices
   call elsi_init_rw(rwh,0,1,0,0.0_r8)
   call elsi_set_rw_mpi(rwh,comm)
   call elsi_set_rw_blacs(rwh,blacs_ctxt,blk)

   call elsi_read_mat_dim(rwh,a_file,n_electrons,n_basis,l_rows,l_cols)
   call elsi_get_rw_header(rwh,header)

   l_rows2 = numroc(2*n_basis,blk,myprow,0,nprow)
   l_cols2 = numroc(2*n_basis,blk,mypcol,0,npcol)

   allocate(mat_a(l_rows,l_cols))
   allocate(mat_b(l_rows,l_cols))
   allocate(evec(l_rows2,l_cols2))
   allocate(eval(n_basis))

   t1 = MPI_Wtime()

   call elsi_read_mat_real(rwh,a_file,mat_a)
   call elsi_read_mat_real(rwh,b_file,mat_b)

   call elsi_finalize_rw(rwh)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished reading A and B matrices"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Initialize ELSI
   call elsi_init(eh,solver,1,0,n_basis,n_electrons,n_basis)
   call elsi_set_mpi(eh,comm)
   call elsi_set_blacs(eh,blacs_ctxt,blk)

   ! Customize ELSI
   call elsi_set_output(eh,2)
   call elsi_set_output_log(eh,1)

   inquire(file=file_name,exist=file_exist)

   if(file_exist) then
      call elsi_set_input_file(eh,file_name)
   end if

   t1 = MPI_Wtime()

   ! Solve
   call elsi_bse_real(eh,mat_a,mat_b,eval,evec)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished BSE"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
      write(*,"(2X,A)") "Finished test program"
      write(*,*)
      if(header(8) == 2222) then
         write(*,"(2X,A)") "Lowest eigenvalue"
         write(*,"(2X,A,F15.8)") "| This test :",eval(1)
         write(*,"(2X,A,F15.8)") "| Reference :",e_ref
         write(*,*)
         if(abs(eval(1)-e_ref) < tol) then
            write(*,"(2X,A)") "Passed."
         else
            write(*,"(2X,A)") "Failed."
         end if
      end if
      write(*,*)
   end if

   ! Finalize ELSI
   call elsi_finalize(eh)

   deallocate(mat_a)
   deallocate(mat_b)
   deallocate(evec)
   deallocate(eval)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)

end subroutine
