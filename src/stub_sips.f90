! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module is only compiled when the actual SLEPc-SIPs is not available, to
!! make the SLEPc-SIPs part of ELSI compile.
!!
module M_SIPS

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
   public :: sips_get_inertias
   public :: sips_get_slices
   public :: sips_get_slices_from_inertias
   public :: sips_set_slices
   public :: sips_solve_eps
   public :: sips_get_eigenvalues
   public :: sips_get_eigenvectors
   public :: sips_get_dm
   public :: sips_get_edm

contains

subroutine sips_initialize()

   implicit none

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_finalize()

   implicit none

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_load_ham_ovlp(ncol_g,ncol_l,nnz_l,col_idx,row_ptr,ham,ovlp)

   implicit none

   integer(kind=i4) :: ncol_g
   integer(kind=i4) :: ncol_l
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: col_idx(nnz_l)
   integer(kind=i4) :: row_ptr(ncol_l+1)
   real(kind=r8) :: ham(nnz_l)
   real(kind=r8) :: ovlp(nnz_l)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_load_ham(ncol_g,ncol_l,nnz_l,col_idx,row_ptr,ham)

   implicit none

   integer(kind=i4) :: ncol_g
   integer(kind=i4) :: ncol_l
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: col_idx(nnz_l)
   integer(kind=i4) :: row_ptr(ncol_l+1)
   real(kind=r8) :: ham(nnz_l)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_update_ham(ncol_g,ncol_l,nnz_l,col_idx,row_ptr,ham)

   implicit none

   integer(kind=i4) :: ncol_g
   integer(kind=i4) :: ncol_l
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: col_idx(nnz_l)
   integer(kind=i4) :: row_ptr(ncol_l+1)
   real(kind=r8) :: ham(nnz_l)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_set_eps(stdevp)

   implicit none

   integer(kind=i4) :: stdevp

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_update_eps(nsub)

   implicit none

   integer(kind=i4) :: nsub

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_set_slices(nsub,subs)

   implicit none

   integer(kind=i4) :: nsub
   real(kind=r8) :: subs(nsub+1)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_solve_eps(nreq,nconv)

   implicit none

   integer(kind=i4) :: nreq
   integer(kind=i4) :: nconv

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_get_eigenvalues(nev,evals)

   implicit none

   integer(kind=i4) :: nev
   real(kind=r8) :: evals(nev)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_get_eigenvectors(nev,lrow,evec)

   implicit none

   integer(kind=i4) :: nev
   integer(kind=i4) :: lrow
   real(kind=r8) :: evec(lrow,nev)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_get_inertias(nsub,subs,inertias)

   implicit none

   integer(kind=i4) :: nsub
   real(kind=r8) :: subs(nsub+1)
   integer(kind=i4) :: inertias(nsub+1)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_get_slices(algr,nev,nsub,buf,subbuf,evals,subs)

   implicit none

   integer(kind=i4) :: algr
   integer(kind=i4) :: nev
   integer(kind=i4) :: nsub
   real(kind=r8) :: buf
   real(kind=r8) :: subbuf
   real(kind=r8) :: evals(nev)
   real(kind=r8) :: subs(nsub+1)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_get_slices_from_inertias(nev,nsub,inertias,subs)

   implicit none

   integer(kind=i4) :: nev
   integer(kind=i4) :: nsub
   integer(kind=i4) :: inertias(nsub+1)
   real(kind=r8) :: subs(nsub+1)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_get_dm(ncol_l,nnz_l,col_idx,row_ptr,nev,occ,dm)

   implicit none

   integer(kind=i4) :: ncol_l
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: col_idx(nnz_l)
   integer(kind=i4) :: row_ptr(ncol_l+1)
   integer(kind=i4) :: nev
   real(kind=r8) :: occ(nev)
   real(kind=r8) :: dm(nnz_l)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

subroutine sips_get_edm(ncol_l,nnz_l,col_idx,row_ptr,nev,occ,edm)

   implicit none

   integer(kind=i4) :: ncol_l
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: col_idx(nnz_l)
   integer(kind=i4) :: row_ptr(ncol_l+1)
   integer(kind=i4) :: nev
   real(kind=r8) :: occ(nev)
   real(kind=r8) :: edm(nnz_l)

   write(*,"(A)") "**Error! A SLEPc-SIPs stub routine was called."
   stop

end subroutine

end module M_SIPS
