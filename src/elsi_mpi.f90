! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains a collection of basic utility routines related to MPI.
!!
module ELSI_MPI

   use ELSI_DATATYPE, only: elsi_basic_t
   use ELSI_PRECISION, only: i4

   implicit none

   include "mpif.h"

   public :: elsi_stop
   public :: elsi_check_mpi

contains

!>
!! Clean shutdown in case of errors.
!!
subroutine elsi_stop(bh,info,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   character(len=*), intent(in) :: info
   character(len=*), intent(in) :: caller

   character(len=200) :: info_str
   integer(kind=i4) :: ierr

   if(bh%mpi_all_ready) then
      write(info_str,"(A,I7,4A)") "**Error! MPI task ",bh%myid_all," in ",&
         trim(caller),": ",trim(info)
      write(*,"(A)") trim(info_str)

      call MPI_Abort(bh%comm_all,0,ierr)
   else if(bh%mpi_ready) then
      write(info_str,"(A,I7,4A)") "**Error! MPI task ",bh%myid," in ",&
         trim(caller),": ",trim(info)
      write(*,"(A)") trim(info_str)

      call MPI_Abort(bh%comm,0,ierr)
   else
      write(info_str,"(4A)") "**Error! ",trim(caller),": ",trim(info)
      write(*,"(A)") trim(info_str)
   end if

   stop

end subroutine

!>
!! Checks if an MPI call is successful.
!!
subroutine elsi_check_mpi(bh,routine,ierr,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   character(len=*), intent(in) :: routine
   integer(kind=i4), intent(in) :: ierr
   character(len=*), intent(in) :: caller

   character(len=200) :: info_str

   if(ierr /= MPI_SUCCESS) then
      write(info_str,"(2A)") routine," failed."

      call elsi_stop(bh,info_str,caller)
   end if

end subroutine

end module ELSI_MPI
