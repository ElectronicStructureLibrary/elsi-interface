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

   use elsi_precision, only : dp

   implicit none

   public :: init_sips
   public :: finalize_sips
   public :: sips_load_ham
   public :: sips_load_ovlp
   public :: sips_set_evp
   public :: sips_get_interval
   public :: sips_get_slicing
   public :: sips_solve_evp
   public :: sips_get_eigenvalues

contains

subroutine init_sips()

   implicit none

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine  finalize_sips()

   implicit none

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_load_ham(global_size,local_size,local_nnz,col_idx,row_ptr,ham_val)

   implicit none

   integer :: global_size
   integer :: local_size
   integer :: local_nnz
   integer :: col_idx(local_nnz)
   integer :: row_ptr(local_size+1)
   real(kind=dp)  :: ham_val(local_size+1)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_load_ovlp(global_size,local_size,local_nnz,col_idx,row_ptr,ovlp_val)

   implicit none

   integer :: global_size
   integer :: local_size
   integer :: local_nnz
   integer :: col_idx(local_nnz)
   integer :: row_ptr(local_size+1)
   real(kind=dp)  :: ovlp_val(local_size+1)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_set_evp(eps_type)

   implicit none

   integer :: eps_type

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_get_slicing(nsub,subtype,lefttype,interval,buffer,&
                            subbuffer,subs,n,evals)

   implicit none

   integer          :: nsub
   real(kind=dp)           :: subs(nsub+1)
   integer          :: subtype
   integer          :: lefttype
   real(kind=dp)           :: interval(2)
   real(kind=dp)           :: buffer
   real(kind=dp)           :: subbuffer
   integer          :: n
   real(kind=dp), optional :: evals(:)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_solve_evp(nreq,nsub,subs,nsolve)

   implicit none

   integer :: nreq
   integer :: nsub
   real(kind=dp)  :: subs(nsub)
   integer :: nsolve

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_get_eigenvalues(nreq,eigs)

   implicit none

   integer :: nreq
   real(kind=dp)  :: eigs(nreq)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_get_eigenvectors(n_basis,idx,evec)

   implicit none

   integer :: n_basis
   integer :: idx
   real(kind=dp)  :: evec(n_basis)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_get_interval(interval)

   implicit none

   real(kind=dp) :: interval(2)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

end module
