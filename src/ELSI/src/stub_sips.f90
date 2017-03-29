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

! This module is only compiled when the actual SIPs is not available, to make
! the SIPs part of ELSI compile.

module m_qetsc

   implicit none

   public :: initialize_qetsc
   public :: finalize_qetsc
   public :: load_elsi_ham
   public :: load_elsi_ovlp
   public :: set_eps
   public :: compute_subintervals
   public :: set_eps_subintervals
   public :: run_eps_inertias_check
   public :: inertias_to_eigenvalues
   public :: solve_eps_check
   public :: get_eps_eigenvalues
   public :: get_eps_interval

contains

subroutine initialize_qetsc()

   implicit none

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine  finalize_qetsc()

   implicit none

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine load_elsi_ham(global_size,local_size,local_nnz,col_idx,row_ptr,ham_val)

   implicit none

   integer :: global_size
   integer :: local_size
   integer :: local_nnz
   integer :: col_idx(local_nnz)
   integer :: row_ptr(local_size+1)
   real*8  :: ham_val(local_size+1)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine load_elsi_ovlp(global_size,local_size,local_nnz,col_idx,row_ptr,ovlp_val)

   implicit none

   integer :: global_size
   integer :: local_size
   integer :: local_nnz
   integer :: col_idx(local_nnz)
   integer :: row_ptr(local_size+1)
   real*8  :: ovlp_val(local_size+1)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine set_eps(eps_type)

   implicit none

   integer :: eps_type

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine compute_subintervals(nsub,subtype,lefttype,interval,buffer,&
                                subbuffer,subs,evals)

   implicit none

   integer          :: nsub
   real*8           :: subs(nsub+1)
   integer          :: subtype
   integer          :: lefttype
   real*8           :: interval(2)
   real*8           :: buffer
   real*8           :: subbuffer
   real*8, optional :: evals(:)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine set_eps_subintervals(nslice,subs)

   implicit none

   integer :: nslice
   real*8  :: subs(nslice+1)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine run_eps_inertias_check(unbound,nreq,nsub,subs,shifts,inertias,nsolve)

   implicit none

   integer :: unbound
   integer :: nreq
   integer :: nsub
   real*8  :: subs(:)
   real*8  :: shifts(:)
   integer :: inertias(:)
   integer :: nsolve

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine inertias_to_eigenvalues(nshift,nreq,buffer,shifts,inertias,evals)

   implicit none

   integer :: nshift
   integer :: nreq
   real*8  :: buffer
   real*8  :: shifts(nshift)
   integer :: inertias(nshift)
   real*8  :: evals(nreq)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine solve_eps_check(nreq,nsub,subs,nsolve)

   implicit none

   integer :: nreq
   integer :: nsub
   real*8  :: subs(nsub)
   integer :: nsolve

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

function get_eps_eigenvalues(nreq) result(eigs)

   implicit none

   integer :: nreq
   real*8  :: eigs(nreq)

   eigs = 0.0d0

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end function

function get_eps_interval() result(interval)

   implicit none

   real*8 :: interval(2)

   interval = 0.0d0

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end function

end module
