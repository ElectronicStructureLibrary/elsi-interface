! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This program tests ELSI.
!!
program elsi_test

   use ELSI_PRECISION, only: r8,i4

   implicit none

   include "mpif.h"

   character(128) :: arg1 ! ev or dm
   character(128) :: arg2 ! dense or sparse
   character(128) :: arg3 ! real or cmplx
   character(128) :: arg4 ! solver
   character(128) :: arg5 ! H file
   character(128) :: arg6 ! S file

   integer(kind=i4) :: solver
   integer(kind=i4) :: myid
   integer(kind=i4) :: mpierr

   call MPI_Init(mpierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,myid,mpierr)

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
   endif

   select case(arg1(1:1))
   case("e") ! ev
      select case(arg2(1:1))
      case("d") ! dense
         select case(arg3(1:1))
         case("r") ! real
            call test_ev_real(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_ev_complex(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case("s") ! sparse
         select case(arg3(1:1))
         case("r") ! real
            call test_ev_real_sparse(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_ev_complex_sparse(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case default
         call test_die()
      end select
   case("d") ! dm
      select case(arg2(1:1))
      case("d") ! dense
         select case(arg3(1:1))
         case("r") ! real
            call test_dm_real(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_dm_complex(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case("s") ! sparse
         select case(arg3(1:1))
         case("r") ! real
            call test_dm_real_sparse(MPI_COMM_WORLD,solver,arg5,arg6)
         case("c") ! complex
            call test_dm_complex_sparse(MPI_COMM_WORLD,solver,arg5,arg6)
         case default
            call test_die()
         end select
      case default
         call test_die()
      end select
   case default
      call test_die()
   end select

   call MPI_Finalize(mpierr)

contains

subroutine test_die()

   implicit none

   if(myid == 0) then
      write(*,'(A)') "  ########################################"
      write(*,'(A)') "  ##                                    ##"
      write(*,'(A)') "  ##  Wrong command line argument(s)!!  ##"
      write(*,'(A)') "  ##                                    ##"
      write(*,'(A)') "  ##  Arg #1: 'ev' or 'dm'              ##"
      write(*,'(A)') "  ##  Arg #2: 'dense' or 'sparse'       ##"
      write(*,'(A)') "  ##  Arg #3: 'real' or 'complex'       ##"
      write(*,'(A)') "  ##  Arg #4: 1 = ELPA                  ##"
      write(*,'(A)') "  ##          2 = libOMM                ##"
      write(*,'(A)') "  ##          3 = PEXSI                 ##"
      write(*,'(A)') "  ##          5 = SIPs                  ##"
      write(*,'(A)') "  ##  Arg #5: H matrix file             ##"
      write(*,'(A)') "  ##  Arg #6: S matrix file             ##"
      write(*,'(A)') "  ##                                    ##"
      write(*,'(A)') "  ########################################"
      call MPI_Abort(MPI_COMM_WORLD,0,mpierr)
      stop
   endif

end subroutine

end program
