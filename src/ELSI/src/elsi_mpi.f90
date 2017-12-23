! Copyright (c) 2015-2017, the ELSI team. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!  * Neither the name of the "ELectronic Structure Infrastructure" project nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!>
!! This module contains a collection of basic utility routines.
!!
module ELSI_MPI

   use ELSI_DATATYPE
   use ELSI_PRECISION, only: i4

   implicit none

   include "mpif.h"

   public :: elsi_stop

contains

!>
!! Clean shutdown in case of errors.
!!
subroutine elsi_stop(info,e_h,caller)

   implicit none

   character(len=*),  intent(in) :: info
   type(elsi_handle), intent(in) :: e_h
   character(len=*),  intent(in) :: caller

   character*800    :: info_str
   integer(kind=i4) :: mpierr

   if(e_h%global_mpi_ready) then
      write(info_str,"(A,I7,5A)") "**Error! MPI task ",e_h%myid_all," in ",&
         trim(caller),": ",trim(info)," Exiting..."
      write(e_h%print_unit,"(A)") trim(info_str)

      if(e_h%n_procs_all > 1) then
         call MPI_Abort(e_h%mpi_comm_all,0,mpierr)
      endif
   elseif(e_h%mpi_ready) then
      write(info_str,"(A,I7,5A)") "**Error! MPI task ",e_h%myid," in ",&
         trim(caller),": ",trim(info)," Exiting..."
      write(e_h%print_unit,"(A)") trim(info_str)

      if(e_h%n_procs > 1) then
         call MPI_Abort(e_h%mpi_comm,0,mpierr)
      endif
   else
      write(info_str,"(5A)") "**Error! ",trim(caller),": ",trim(info),&
         " Exiting..."
      write(e_h%print_unit,"(A)") trim(info_str)
   endif

   stop

end subroutine

end module ELSI_MPI
