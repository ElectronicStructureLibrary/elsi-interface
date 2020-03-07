! Copyright (c) 2015-2020, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This program tests ELSI.
!!
program elsi_test

   use ELSI_MPI
   use ELSI_PRECISION, only: i4

   implicit none

   character(len=100) :: arg1 ! ev, dm, or bse
   character(len=100) :: arg2 ! dense or sparse
   character(len=100) :: arg3 ! real or cmplx
   character(len=100) :: arg4 ! solver
   character(len=999) :: arg5 ! H file
   character(len=999) :: arg6 ! S file

   integer(kind=i4) :: solver
   integer(kind=i4) :: myid
   integer(kind=i4) :: ierr

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierr)

   ! Read command line arguments
   if(command_argument_count() == 6) then
      call get_command_argument(1,arg1)
      call get_command_argument(2,arg2)
      call get_command_argument(3,arg3)
      call get_command_argument(4,arg4)
      call get_command_argument(5,arg5)
      call get_command_argument(6,arg6)
      read(arg4,*) solver
   else
      call test_die()
   end if

   select case(arg1(1:1))
   case("e") ! ev
      select case(arg2(1:1))
      case("0") ! BLACS_DENSE
         select case(arg3(1:1))
         case("r") ! real
            call test_ev_real_den(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_ev_cmplx_den(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case("1") ! PEXSI_CSC
         select case(arg3(1:1))
         case("r") ! real
            call test_ev_real_csc1(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_ev_cmplx_csc1(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case("2") ! SIESTA_CSC
         select case(arg3(1:1))
         case("r") ! real
            call test_ev_real_csc2(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_ev_cmplx_csc2(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case("3") ! GENERIC_COO
         select case(arg3(1:1))
         case("r") ! real
            call test_ev_real_coo(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_ev_cmplx_coo(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case default
         call test_die()
      end select
   case("d") ! dm
      select case(arg2(1:1))
      case("0") ! BLACS_DENSE
         select case(arg3(1:1))
         case("r") ! real
            call test_dm_real_den(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_dm_cmplx_den(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case("1") ! PEXSI_CSC
         select case(arg3(1:1))
         case("r") ! real
            call test_dm_real_csc1(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_dm_cmplx_csc1(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case("2") ! SIESTA_CSC
         select case(arg3(1:1))
         case("r") ! real
            call test_dm_real_csc2(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_dm_cmplx_csc2(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case("3") ! GENERIC_COO
         select case(arg3(1:1))
         case("r") ! real
            call test_dm_real_coo(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_dm_cmplx_coo(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case default
         call test_die()
      end select
   case("b") ! bse
      select case(arg3(1:1))
      case("r") ! real
         call test_bse_real_den(MPI_COMM_WORLD,solver,arg5,arg6)
      case("c") ! complex
         call test_bse_cmplx_den(MPI_COMM_WORLD,solver,arg5,arg6)
      case default
         call test_die()
      end select
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
      write(*,"(A)") "  ##  Arg #1: 'ev', 'dm', or 'bse'      ##"
      write(*,"(A)") "  ##  Arg #2: 0 = BLACS_DENSE           ##"
      write(*,"(A)") "  ##          1 = PEXSI_CSC             ##"
      write(*,"(A)") "  ##          2 = SIESTA_CSC            ##"
      write(*,"(A)") "  ##          3 = GENERIC_COO           ##"
      write(*,"(A)") "  ##  Arg #3: 'real' or 'complex'       ##"
      write(*,"(A)") "  ##  Arg #4: 1 = ELPA                  ##"
      write(*,"(A)") "  ##          2 = libOMM                ##"
      write(*,"(A)") "  ##          3 = PEXSI                 ##"
      write(*,"(A)") "  ##          4 = EigenExa              ##"
      write(*,"(A)") "  ##          5 = SLEPc-SIPs            ##"
      write(*,"(A)") "  ##          6 = NTPoly                ##"
      write(*,"(A)") "  ##          7 = MAGMA                 ##"
      write(*,"(A)") "  ##          8 = BSEPACK               ##"
      write(*,"(A)") "  ##  Arg #5: H matrix file             ##"
      write(*,"(A)") "  ##  Arg #6: S matrix file             ##"
      write(*,"(A)") "  ##                                    ##"
      write(*,"(A)") "  ########################################"
      call MPI_Abort(MPI_COMM_WORLD,0,ierr)
      stop
   end if

end subroutine

end program
