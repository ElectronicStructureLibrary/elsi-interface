! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide stub routines which are only compiled when the actual MAGMA is not
!! available.
!!
module MAGMA

   use ELSI_PRECISION, only: r8,i4

   implicit none

   private

   public :: magmaf_init
   public :: magmaf_finalize
   public :: magmaf_num_gpus
   public :: magmaf_dsyevdx_m
   public :: magmaf_dsyevdx_2stage_m
   public :: magmaf_dsygvdx_m
   public :: magmaf_dsygvdx_2stage_m
   public :: magmaf_zheevdx_m
   public :: magmaf_zheevdx_2stage_m
   public :: magmaf_zhegvdx_m
   public :: magmaf_zhegvdx_2stage_m

contains

subroutine magmaf_init()

   implicit none

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

subroutine magmaf_finalize()

   implicit none

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

function magmaf_num_gpus()

   implicit none

   integer(kind=i4) :: magmaf_num_gpus

   magmaf_num_gpus = 0

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end function

subroutine magmaf_dsyevdx_m(ngpu,jobz,range,uplo,n,A,lda,vl,vu,il,iu,mout,w,&
   work,lwork,iwork,liwork,info)

   implicit none

   integer(kind=i4) :: ngpu
   character :: jobz
   character :: range
   character :: uplo
   integer(kind=i4) :: n
   real(kind=r8) :: A(*)
   integer(kind=i4) :: lda
   real(kind=r8) :: vl
   real(kind=r8) :: vu
   integer(kind=i4) :: il
   integer(kind=i4) :: iu
   integer(kind=i4) :: mout(*)
   real(kind=r8) :: w(*)
   real(kind=r8) :: work(*)
   integer(kind=i4) :: lwork
   integer(kind=i4) :: iwork(*)
   integer(kind=i4) :: liwork
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

subroutine magmaf_dsyevdx_2stage_m(ngpu,jobz,range,uplo,n,A,lda,vl,vu,il,iu,&
   mout,w,work,lwork,iwork,liwork,info)

   implicit none

   integer(kind=i4) :: ngpu
   character :: jobz
   character :: range
   character :: uplo
   integer(kind=i4) :: n
   real(kind=r8) :: A(*)
   integer(kind=i4) :: lda
   real(kind=r8) :: vl
   real(kind=r8) :: vu
   integer(kind=i4) :: il
   integer(kind=i4) :: iu
   integer(kind=i4) :: mout(*)
   real(kind=r8) :: w(*)
   real(kind=r8) :: work(*)
   integer(kind=i4) :: lwork
   integer(kind=i4) :: iwork(*)
   integer(kind=i4) :: liwork
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

subroutine magmaf_dsygvdx_m(ngpu,itype,jobz,range,uplo,n,A,lda,B,ldb,vl,vu,il,&
   iu,mout,w,work,lwork,iwork,liwork,info)

   implicit none

   integer(kind=i4) :: ngpu
   integer(kind=i4) :: itype
   character :: jobz
   character :: range
   character :: uplo
   integer(kind=i4) :: n
   real(kind=r8) :: A(*)
   integer(kind=i4) :: lda
   real(kind=r8) :: B(*)
   integer(kind=i4) :: ldb
   real(kind=r8) :: vl
   real(kind=r8) :: vu
   integer(kind=i4) :: il
   integer(kind=i4) :: iu
   integer(kind=i4) :: mout(*)
   real(kind=r8) :: w(*)
   real(kind=r8) :: work(*)
   integer(kind=i4) :: lwork
   integer(kind=i4) :: iwork(*)
   integer(kind=i4) :: liwork
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

subroutine magmaf_dsygvdx_2stage_m(ngpu,itype,jobz,range,uplo,n,A,lda,B,ldb,vl,&
   vu,il,iu,mout,w,work,lwork,iwork,liwork,info)

   implicit none

   integer(kind=i4) :: ngpu
   integer(kind=i4) :: itype
   character :: jobz
   character :: range
   character :: uplo
   integer(kind=i4) :: n
   real(kind=r8) :: A(*)
   integer(kind=i4) :: lda
   real(kind=r8) :: B(*)
   integer(kind=i4) :: ldb
   real(kind=r8) :: vl
   real(kind=r8) :: vu
   integer(kind=i4) :: il
   integer(kind=i4) :: iu
   integer(kind=i4) :: mout(*)
   real(kind=r8) :: w(*)
   real(kind=r8) :: work(*)
   integer(kind=i4) :: lwork
   integer(kind=i4) :: iwork(*)
   integer(kind=i4) :: liwork
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

subroutine magmaf_zheevdx_m(ngpu,jobz,range,uplo,n,A,lda,vl,vu,il,iu,mout,w,&
   work,lwork,rwork,lrwork,iwork,liwork,info)

   implicit none

   integer(kind=i4) :: ngpu
   character :: jobz
   character :: range
   character :: uplo
   integer(kind=i4) :: n
   complex(kind=r8) :: A(*)
   integer(kind=i4) :: lda
   real(kind=r8) :: vl
   real(kind=r8) :: vu
   integer(kind=i4) :: il
   integer(kind=i4) :: iu
   integer(kind=i4) :: mout(*)
   real(kind=r8) :: w(*)
   complex(kind=r8) :: work(*)
   integer(kind=i4) :: lwork
   real(kind=r8) :: rwork(*)
   integer(kind=i4) :: lrwork
   integer(kind=i4) :: iwork(*)
   integer(kind=i4) :: liwork
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

subroutine magmaf_zheevdx_2stage_m(ngpu,jobz,range,uplo,n,A,lda,vl,vu,il,iu,&
   mout,w,work,lwork,rwork,lrwork,iwork,liwork,info)

   implicit none

   integer(kind=i4) :: ngpu
   character :: jobz
   character :: range
   character :: uplo
   integer(kind=i4) :: n
   complex(kind=r8) :: A(*)
   integer(kind=i4) :: lda
   real(kind=r8) :: vl
   real(kind=r8) :: vu
   integer(kind=i4) :: il
   integer(kind=i4) :: iu
   integer(kind=i4) :: mout(*)
   real(kind=r8) :: w(*)
   complex(kind=r8) :: work(*)
   integer(kind=i4) :: lwork
   real(kind=r8) :: rwork(*)
   integer(kind=i4) :: lrwork
   integer(kind=i4) :: iwork(*)
   integer(kind=i4) :: liwork
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

subroutine magmaf_zhegvdx_m(ngpu,itype,jobz,range,uplo,n,A,lda,B,ldb,vl,vu,il,&
   iu,mout,w,work,lwork,rwork,lrwork,iwork,liwork,info)

   implicit none

   integer(kind=i4) :: ngpu
   integer(kind=i4) :: itype
   character :: jobz
   character :: range
   character :: uplo
   integer(kind=i4) :: n
   complex(kind=r8) :: A(*)
   integer(kind=i4) :: lda
   complex(kind=r8) :: B(*)
   integer(kind=i4) :: ldb
   real(kind=r8) :: vl
   real(kind=r8) :: vu
   integer(kind=i4) :: il
   integer(kind=i4) :: iu
   integer(kind=i4) :: mout(*)
   real(kind=r8) :: w(*)
   complex(kind=r8) :: work(*)
   integer(kind=i4) :: lwork
   real(kind=r8) :: rwork(*)
   integer(kind=i4) :: lrwork
   integer(kind=i4) :: iwork(*)
   integer(kind=i4) :: liwork
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

subroutine magmaf_zhegvdx_2stage_m(ngpu,itype,jobz,range,uplo,n,A,lda,B,ldb,vl,&
   vu,il,iu,mout,w,work,lwork,rwork,lrwork,iwork,liwork,info)

   implicit none

   integer(kind=i4) :: ngpu
   integer(kind=i4) :: itype
   character :: jobz
   character :: range
   character :: uplo
   integer(kind=i4) :: n
   complex(kind=r8) :: A(*)
   integer(kind=i4) :: lda
   complex(kind=r8) :: B(*)
   integer(kind=i4) :: ldb
   real(kind=r8) :: vl
   real(kind=r8) :: vu
   integer(kind=i4) :: il
   integer(kind=i4) :: iu
   integer(kind=i4) :: mout(*)
   real(kind=r8) :: w(*)
   complex(kind=r8) :: work(*)
   integer(kind=i4) :: lwork
   real(kind=r8) :: rwork(*)
   integer(kind=i4) :: lrwork
   integer(kind=i4) :: iwork(*)
   integer(kind=i4) :: liwork
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A MAGMA stub routine was called"
   stop

end subroutine

end module MAGMA
