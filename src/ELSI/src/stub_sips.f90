! Copyright (c) 2015-2018, the ELSI team. All rights reserved.
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
!! This module is only compiled when the actual SIPs is not available, to make
!! the SIPs part of ELSI compile.
!!
module M_QETSC

   use ELSI_PRECISION, only: r8,i4

   implicit none

   private

   public :: sips_initialize
   public :: sips_finalize
   public :: sips_load_ham_ovlp
   public :: sips_load_ham
   public :: sips_update_ham
   public :: sips_set_eps
   public :: sips_update_eps
   public :: sips_set_slices
   public :: sips_solve_eps
   public :: sips_get_eigenvalues
   public :: sips_get_eigenvectors

contains

subroutine sips_initialize()

   implicit none

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_finalize()

   implicit none

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_load_ham_ovlp(global_size,local_size,local_nnz,col_idx,row_ptr,&
              ham_val,ovlp_val)

   implicit none

   integer(kind=i4) :: global_size
   integer(kind=i4) :: local_size
   integer(kind=i4) :: local_nnz
   integer(kind=i4) :: col_idx(local_nnz)
   integer(kind=i4) :: row_ptr(local_size+1)
   real(kind=r8)    :: ham_val(local_size+1)
   real(kind=r8)    :: ovlp_val(local_size+1)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_load_ham(global_size,local_size,local_nnz,col_idx,row_ptr,&
              ham_val)

   implicit none

   integer(kind=i4) :: global_size
   integer(kind=i4) :: local_size
   integer(kind=i4) :: local_nnz
   integer(kind=i4) :: col_idx(local_nnz)
   integer(kind=i4) :: row_ptr(local_size+1)
   real(kind=r8)    :: ham_val(local_size+1)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_update_ham(global_size,local_size,local_nnz,col_idx,row_ptr,&
              ham_val)

   implicit none

   integer(kind=i4) :: global_size
   integer(kind=i4) :: local_size
   integer(kind=i4) :: local_nnz
   integer(kind=i4) :: col_idx(local_nnz)
   integer(kind=i4) :: row_ptr(local_size+1)
   real(kind=r8)    :: ham_val(local_size+1)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_set_eps(stdevp)

   implicit none

   integer(kind=i4) :: stdevp

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_update_eps(nsub)

   implicit none

   integer(kind=i4) :: nsub

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_set_slices(tol,buf,nev,eigs,nsub,subs)

   implicit none

   real(kind=r8)    :: tol
   real(kind=r8)    :: buf
   integer(kind=i4) :: nev
   real(kind=r8)    :: eigs(nev)
   integer(kind=i4) :: nsub
   real(kind=r8)    :: subs(nsub+1)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_solve_eps(nreq)

   implicit none

   integer(kind=i4) :: nreq

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_get_eigenvalues(nev,eigs)

   implicit none

   integer(kind=i4) :: nev
   real(kind=r8)    :: eigs(nev)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sips_get_eigenvectors(n_basis,idx,evec)

   implicit none

   integer(kind=i4) :: n_basis
   integer(kind=i4) :: idx
   real(kind=r8)    :: evec(n_basis)

   write(*,"(A)") " A SIPs stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

end module M_QETSC
