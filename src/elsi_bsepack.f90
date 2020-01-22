! Copyright (c) 2015-2020, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Interface to BSEPACK.
!!
module ELSI_BSEPACK

   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI, only: elsi_stop
   use ELSI_OUTPUT, only: elsi_say,elsi_get_time
   use ELSI_PRECISION, only: r8,i4

   implicit none

   private

   public :: elsi_solve_bsepack

   interface elsi_solve_bsepack
      module procedure elsi_solve_bsepack_real
      module procedure elsi_solve_bsepack_cmplx
   end interface

contains

!>
!! Interface to BSEPACK.
!!
subroutine elsi_solve_bsepack_real(ph,bh,mat_a,mat_b,eval,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: mat_a(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: mat_b(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(ph%bse_n_lrow,ph%bse_n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: lwork
   integer(kind=i4) :: liwork
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   real(kind=r8), allocatable :: work(:)
   integer(kind=i4), allocatable :: iwork(:)

   character(len=*), parameter :: caller = "elsi_solve_bsepack_real"

   write(msg,"(A)") "Starting BSEPACK eigensolver"
   call elsi_say(bh,msg)

   call elsi_get_time(t0)

   lwork = -1
   liwork = -1

   call elsi_allocate(bh,work,1,"work",caller)
   call elsi_allocate(bh,iwork,1,"iwork",caller)

   call pdbseig(0,ph%n_basis,mat_a,1,1,bh%desc,mat_b,1,1,bh%desc,eval,evec,1,1,&
        ph%bse_desc,work,lwork,iwork,liwork,ierr)

   lwork = ceiling(work(1),kind=i4)
   liwork = iwork(1)

   call elsi_deallocate(bh,work,"work")
   call elsi_deallocate(bh,iwork,"iwork")
   call elsi_allocate(bh,work,lwork,"work",caller)
   call elsi_allocate(bh,iwork,liwork,"iwork",caller)

   call pdbseig(0,ph%n_basis,mat_a,1,1,bh%desc,mat_b,1,1,bh%desc,eval,evec,1,1,&
        ph%bse_desc,work,lwork,iwork,liwork,ierr)

   call elsi_deallocate(bh,work,"work")
   call elsi_deallocate(bh,iwork,"iwork")

   if(ierr /= 0) then
      write(msg,"(A)") "BSEPACK eigensolver failed"
      call elsi_stop(bh,msg,caller)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished solving BSE eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Interface to BSEPACK.
!!
subroutine elsi_solve_bsepack_cmplx(ph,bh,mat_a,mat_b,eval,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: mat_a(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: mat_b(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(ph%bse_n_lrow,ph%bse_n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: lwork
   integer(kind=i4) :: lrwork
   integer(kind=i4) :: liwork
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   complex(kind=r8), allocatable :: work(:)
   real(kind=r8), allocatable :: rwork(:)
   integer(kind=i4), allocatable :: iwork(:)

   character(len=*), parameter :: caller = "elsi_solve_bsepack_cmplx"

   write(msg,"(A)") "Starting BSEPACK eigensolver"
   call elsi_say(bh,msg)

   call elsi_get_time(t0)

   lwork = -1
   lrwork = -1
   liwork = -1

   call elsi_allocate(bh,work,1,"work",caller)
   call elsi_allocate(bh,rwork,1,"rwork",caller)
   call elsi_allocate(bh,iwork,1,"iwork",caller)

   call pzbseig(0,ph%n_basis,mat_a,1,1,bh%desc,mat_b,1,1,bh%desc,eval,evec,1,1,&
        ph%bse_desc,work,lwork,rwork,lrwork,iwork,liwork,ierr)

   lwork = ceiling(real(work(1),kind=r8),kind=i4)
   lrwork = ceiling(rwork(1),kind=i4)
   liwork = iwork(1)

   call elsi_deallocate(bh,work,"work")
   call elsi_deallocate(bh,rwork,"rwork")
   call elsi_deallocate(bh,iwork,"iwork")
   call elsi_allocate(bh,work,lwork,"work",caller)
   call elsi_allocate(bh,rwork,lrwork,"rwork",caller)
   call elsi_allocate(bh,iwork,liwork,"iwork",caller)

   call pzbseig(0,ph%n_basis,mat_a,1,1,bh%desc,mat_b,1,1,bh%desc,eval,evec,1,1,&
        ph%bse_desc,work,lwork,rwork,lrwork,iwork,liwork,ierr)

   call elsi_deallocate(bh,work,"work")
   call elsi_deallocate(bh,rwork,"rwork")
   call elsi_deallocate(bh,iwork,"iwork")

   if(ierr /= 0) then
      write(msg,"(A)") "BSEPACK eigensolver failed"
      call elsi_stop(bh,msg,caller)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished solving BSE eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

end module ELSI_BSEPACK
