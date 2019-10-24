! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Interface to MAGMA.
!!
module ELSI_MAGMA

   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI, only: elsi_stop
   use ELSI_OUTPUT, only: elsi_say,elsi_get_time
   use ELSI_PRECISION, only: r8,i4
   use MAGMA, only: magmaf_init,magmaf_finalize,magmaf_num_gpus,&
       magmaf_dsyevdx_m,magmaf_dsyevdx_2stage_m,magmaf_dsygvdx_m,&
       magmaf_dsygvdx_2stage_m,magmaf_zheevdx_m,magmaf_zheevdx_2stage_m,&
       magmaf_zhegvdx_m,magmaf_zhegvdx_2stage_m

   implicit none

   private

   public :: elsi_init_magma
   public :: elsi_cleanup_magma
   public :: elsi_solve_magma

   interface elsi_solve_magma
      module procedure elsi_solve_magma_real
      module procedure elsi_solve_magma_cmplx
   end interface

contains

!>
!! Initialize MAGMA.
!!
subroutine elsi_init_magma(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   character(len=*), parameter :: caller = "elsi_init_magma"

   if(.not. ph%magma_started) then
      ph%magma_started = .true.

      call magmaf_init()

      ph%magma_n_gpus = magmaf_num_gpus()
   end if

end subroutine

!>
!! Interface to MAGMA.
!!
subroutine elsi_solve_magma_real(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: n_solved(1)
   integer(kind=i4) :: lwork
   integer(kind=i4) :: liwork
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   real(kind=r8), allocatable :: work(:)
   integer(kind=i4), allocatable :: iwork(:)

   character(len=*), parameter :: caller = "elsi_solve_magma_real"

   write(msg,"(A)") "Starting MAGMA eigensolver"
   call elsi_say(bh,msg)

   call elsi_get_time(t0)

   call elsi_allocate(bh,work,1,"work",caller)
   call elsi_allocate(bh,iwork,1,"iwork",caller)

   ! Solve
   if(ph%unit_ovlp) then
      select case(ph%magma_solver)
      case(1)
         call magmaf_dsyevdx_m(ph%magma_n_gpus,"V","I","L",ph%n_basis,ham,&
              ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,eval,work,-1,&
              iwork,-1,ierr)
      case(2)
         call magmaf_dsyevdx_2stage_m(ph%magma_n_gpus,"V","I","L",ph%n_basis,&
              ham,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,eval,work,-1,&
              iwork,-1,ierr)
      end select
   else
      select case(ph%magma_solver)
      case(1)
         call magmaf_dsygvdx_m(ph%magma_n_gpus,1,"V","I","L",ph%n_basis,ham,&
              ph%n_basis,ovlp,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,&
              eval,work,-1,iwork,-1,ierr)
      case(2)
         call magmaf_dsygvdx_2stage_m(ph%magma_n_gpus,1,"V","I","L",ph%n_basis,&
              ham,ph%n_basis,ovlp,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,&
              n_solved,eval,work,-1,iwork,-1,ierr)
      end select
   end if

   lwork = floor(work(1),kind=i4)
   liwork = iwork(1)

   call elsi_deallocate(bh,work,"work")
   call elsi_deallocate(bh,iwork,"iwork")
   call elsi_allocate(bh,work,lwork,"work",caller)
   call elsi_allocate(bh,iwork,liwork,"iwork",caller)

   if(ph%unit_ovlp) then
      select case(ph%magma_solver)
      case(1)
         call magmaf_dsyevdx_m(ph%magma_n_gpus,"V","I","L",ph%n_basis,ham,&
              ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,eval,work,lwork,&
              iwork,liwork,ierr)
      case(2)
         call magmaf_dsyevdx_2stage_m(ph%magma_n_gpus,"V","I","L",ph%n_basis,&
              ham,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,eval,work,&
              lwork,iwork,liwork,ierr)
      end select
   else
      select case(ph%magma_solver)
      case(1)
         call magmaf_dsygvdx_m(ph%magma_n_gpus,1,"V","I","L",ph%n_basis,ham,&
              ph%n_basis,ovlp,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,&
              eval,work,lwork,iwork,liwork,ierr)
      case(2)
         call magmaf_dsygvdx_2stage_m(ph%magma_n_gpus,1,"V","I","L",ph%n_basis,&
              ham,ph%n_basis,ovlp,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,&
              n_solved,eval,work,lwork,iwork,liwork,ierr)
      end select
   end if

   if(ierr /= 0) then
      write(msg,"(A)") "MAGMA eigensolver failed"
      call elsi_stop(bh,msg,caller)
   end if

   call elsi_deallocate(bh,work,"work")
   call elsi_deallocate(bh,iwork,"iwork")

   evec = ham

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished solving eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Interface to MAGMA.
!!
subroutine elsi_solve_magma_cmplx(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   complex(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: n_solved(1)
   integer(kind=i4) :: lwork
   integer(kind=i4) :: lrwork
   integer(kind=i4) :: liwork
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   complex(kind=r8), allocatable :: work(:)
   real(kind=r8), allocatable :: rwork(:)
   integer(kind=i4), allocatable :: iwork(:)

   character(len=*), parameter :: caller = "elsi_solve_magma_cmplx"

   write(msg,"(A)") "Starting MAGMA eigensolver"
   call elsi_say(bh,msg)

   call elsi_get_time(t0)

   call elsi_allocate(bh,work,1,"work",caller)
   call elsi_allocate(bh,rwork,1,"rwork",caller)
   call elsi_allocate(bh,iwork,1,"iwork",caller)

   ! Solve
   if(ph%unit_ovlp) then
      select case(ph%magma_solver)
      case(1)
         call magmaf_zheevdx_m(ph%magma_n_gpus,"V","I","L",ph%n_basis,ham,&
              ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,eval,work,-1,&
              rwork,-1,iwork,-1,ierr)
      case(2)
         call magmaf_zheevdx_2stage_m(ph%magma_n_gpus,"V","I","L",ph%n_basis,&
              ham,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,eval,work,-1,&
              rwork,-1,iwork,-1,ierr)
      end select
   else
      select case(ph%magma_solver)
      case(1)
         call magmaf_zhegvdx_m(ph%magma_n_gpus,1,"V","I","L",ph%n_basis,ham,&
              ph%n_basis,ovlp,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,&
              eval,work,-1,rwork,-1,iwork,-1,ierr)
      case(2)
         call magmaf_zhegvdx_2stage_m(ph%magma_n_gpus,1,"V","I","L",ph%n_basis,&
              ham,ph%n_basis,ovlp,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,&
              n_solved,eval,work,-1,rwork,-1,iwork,-1,ierr)
      end select
   end if

   lwork = floor(real(work(1),kind=r8),kind=i4)
   lrwork = floor(rwork(1),kind=i4)
   liwork = iwork(1)

   call elsi_deallocate(bh,work,"work")
   call elsi_deallocate(bh,rwork,"rwork")
   call elsi_deallocate(bh,iwork,"iwork")
   call elsi_allocate(bh,work,lwork,"work",caller)
   call elsi_allocate(bh,rwork,lrwork,"rwork",caller)
   call elsi_allocate(bh,iwork,liwork,"iwork",caller)

   if(ph%unit_ovlp) then
      select case(ph%magma_solver)
      case(1)
         call magmaf_zheevdx_m(ph%magma_n_gpus,"V","I","L",ph%n_basis,ham,&
              ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,eval,work,lwork,&
              rwork,lrwork,iwork,liwork,ierr)
      case(2)
         call magmaf_zheevdx_2stage_m(ph%magma_n_gpus,"V","I","L",ph%n_basis,&
              ham,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,eval,work,&
              lwork,rwork,lrwork,iwork,liwork,ierr)
      end select
   else
      select case(ph%magma_solver)
      case(1)
         call magmaf_zhegvdx_m(ph%magma_n_gpus,1,"V","I","L",ph%n_basis,ham,&
              ph%n_basis,ovlp,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,n_solved,&
              eval,work,lwork,rwork,lrwork,iwork,liwork,ierr)
      case(2)
         call magmaf_zhegvdx_2stage_m(ph%magma_n_gpus,1,"V","I","L",ph%n_basis,&
              ham,ph%n_basis,ovlp,ph%n_basis,0.0_r8,0.0_r8,1,ph%n_states,&
              n_solved,eval,work,lwork,rwork,lrwork,iwork,liwork,ierr)
      end select
   end if

   if(ierr /= 0 .or. n_solved(1) < ph%n_states) then
      write(msg,"(A)") "MAGMA eigensolver failed"
      call elsi_stop(bh,msg,caller)
   end if

   call elsi_deallocate(bh,work,"work")
   call elsi_deallocate(bh,rwork,"rwork")
   call elsi_deallocate(bh,iwork,"iwork")

   evec = ham

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished solving eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Clean up MAGMA.
!!
subroutine elsi_cleanup_magma(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   character(len=*), parameter :: caller = "elsi_cleanup_magma"

   if(ph%magma_started) then
      call magmaf_finalize()

      ph%magma_started = .false.
   end if

end subroutine

end module ELSI_MAGMA
