! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains a collection of basic utility routines related to MPI.
!!
module ELSI_MPI

   use ELSI_DATATYPE,  only: elsi_handle
   use ELSI_PRECISION, only: i4

   implicit none

   include "mpif.h"

   public :: elsi_stop
   public :: elsi_check_mpi

contains

!>
!! Clean shutdown in case of errors.
!!
subroutine elsi_stop(e_h,info,caller)

   implicit none

   type(elsi_handle), intent(in) :: e_h
   character(len=*),  intent(in) :: info
   character(len=*),  intent(in) :: caller

   character(len=200) :: info_str
   integer(kind=i4)   :: ierr

   if(e_h%global_mpi_ready) then
      write(info_str,"(A,I7,5A)") "**Error! MPI task ",e_h%myid_all," in ",&
         trim(caller),": ",trim(info)," Exiting..."
      write(e_h%stdio%print_unit,"(A)") trim(info_str)

      if(e_h%n_procs_all > 1) then
         call MPI_Abort(e_h%mpi_comm_all,0,ierr)
      endif
   elseif(e_h%mpi_ready) then
      write(info_str,"(A,I7,5A)") "**Error! MPI task ",e_h%myid," in ",&
         trim(caller),": ",trim(info)," Exiting..."
      write(e_h%stdio%print_unit,"(A)") trim(info_str)

      if(e_h%n_procs > 1) then
         call MPI_Abort(e_h%mpi_comm,0,ierr)
      endif
   else
      write(info_str,"(5A)") "**Error! ",trim(caller),": ",trim(info),&
         " Exiting..."
      write(e_h%stdio%print_unit,"(A)") trim(info_str)
   endif

   stop

end subroutine

!>
!! Checks if an MPI call is successful.
!!
subroutine elsi_check_mpi(e_h,routine,ierr,caller)

   implicit none

   type(elsi_handle), intent(in) :: e_h
   character(len=*),  intent(in) :: routine
   integer(kind=i4),  intent(in) :: ierr
   character(len=*),  intent(in) :: caller

   if(ierr /= MPI_SUCCESS) then
      call elsi_stop(e_h,routine,caller)
   endif

end subroutine

end module ELSI_MPI
