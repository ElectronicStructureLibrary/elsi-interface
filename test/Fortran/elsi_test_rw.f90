! Copyright (c) 2015-2020, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This program tests ELSI matrix I/O.
!!
program elsi_test_rw

   use ELSI_MPI
   use ELSI_PRECISION, only: i4

   implicit none

   character(len=100) :: arg1 ! real or cmplx
   character(len=999) :: arg2 ! H file
   character(len=999) :: arg3 ! S file

   integer(kind=i4) :: myid
   integer(kind=i4) :: ierr

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierr)

   ! Read command line arguments
   if(command_argument_count() == 3) then
      call get_command_argument(1,arg1)
      call get_command_argument(2,arg2)
      call get_command_argument(3,arg3)
   else
      call test_die()
   end if

   select case(arg1(1:1))
   case("r") ! real
      call test_rw_real(MPI_COMM_WORLD,arg2,arg3)
   case("c") ! complex
      call test_rw_cmplx(MPI_COMM_WORLD,arg2,arg3)
   case default
      call test_die()
   end select

   call MPI_Finalize(ierr)

contains

subroutine test_die()

   implicit none

   if(myid == 0) then
      write(*,"(A)") "  ########################################"
      write(*,"(A)") "  ##                                    ##"
      write(*,"(A)") "  ##  Wrong command line argument(s)!!  ##"
      write(*,"(A)") "  ##                                    ##"
      write(*,"(A)") "  ##  Arg #1: 'real' or 'complex'       ##"
      write(*,"(A)") "  ##  Arg #2: H matrix file             ##"
      write(*,"(A)") "  ##  Arg #3: S matrix file             ##"
      write(*,"(A)") "  ##                                    ##"
      write(*,"(A)") "  ########################################"
      call MPI_Abort(MPI_COMM_WORLD,0,ierr)
      stop
   end if

end subroutine

end program
