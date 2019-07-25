! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide stub routines which are only compiled when the actual EigenExa is not
!! available.
!!
module EIGEN_LIBS_MOD

   use ELSI_PRECISION, only: r8,i4

   implicit none

   private

   public :: eigen_init
   public :: eigen_get_procs
   public :: eigen_get_id
   public :: eigen_get_matdims
   public :: eigen_s
   public :: eigen_sx
   public :: eigen_free

contains

subroutine eigen_init(comm)

   implicit none

   integer(kind=i4) :: comm

   write(*,"(A)") "**Error! An EigenExa stub routine was called"
   stop

end subroutine

subroutine eigen_get_procs(procs,x_procs,y_procs)

   implicit none

   integer(kind=i4) :: procs
   integer(kind=i4) :: x_procs
   integer(kind=i4) :: y_procs

   write(*,"(A)") "**Error! An EigenExa stub routine was called"
   stop

end subroutine

subroutine eigen_get_id(id,x_id,y_id)

   implicit none

   integer(kind=i4) :: id
   integer(kind=i4) :: x_id
   integer(kind=i4) :: y_id

   write(*,"(A)") "**Error! An EigenExa stub routine was called"
   stop

end subroutine

subroutine eigen_get_matdims(n,nx,ny)

   implicit none

   integer(kind=i4) :: n
   integer(kind=i4) :: nx
   integer(kind=i4) :: ny

   write(*,"(A)") "**Error! An EigenExa stub routine was called"
   stop

end subroutine

subroutine eigen_s(n,nvec,a,lda,w,z,ldz,m_forward,m_backward,mode)

   implicit none

   integer(kind=i4) :: n
   integer(kind=i4) :: nvec
   integer(kind=i4) :: lda
   real(kind=r8) :: a(lda,*)
   real(kind=r8) :: w(*)
   integer(kind=i4) :: ldz
   real(kind=r8) :: z(ldz,*)
   integer(kind=i4) :: m_forward
   integer(kind=i4) :: m_backward
   character :: mode

   write(*,"(A)") "**Error! An EigenExa stub routine was called"
   stop

end subroutine

subroutine eigen_sx(n,nvec,a,lda,w,z,ldz,m_forward,m_backward,mode)

   implicit none

   integer(kind=i4) :: n
   integer(kind=i4) :: nvec
   integer(kind=i4) :: lda
   real(kind=r8) :: a(lda,*)
   real(kind=r8) :: w(*)
   integer(kind=i4) :: ldz
   real(kind=r8) :: z(ldz,*)
   integer(kind=i4) :: m_forward
   integer(kind=i4) :: m_backward
   character :: mode

   write(*,"(A)") "**Error! An EigenExa stub routine was called"
   stop

end subroutine

subroutine eigen_free()

   implicit none

   write(*,"(A)") "**Error! An EigenExa stub routine was called"
   stop

end subroutine

end module EIGEN_LIBS_MOD
