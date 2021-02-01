! Copyright (c) 2015-2021, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests occ number normal occupation.
!!
subroutine test_occ_non_aufbau(comm,mu_width)

   use ELSI
   use ELSI_MPI
   use ELSI_PRECISION, only: r8,i4

   implicit none

   integer(kind=i4), intent(in) :: comm
   real(kind=r8), intent(in) :: mu_width

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: myid
   integer(kind=i4) :: myprow
   integer(kind=i4) :: mypcol
   integer(kind=i4) :: ierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: n_basis
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols
   integer(kind=i4) :: l_rows2
   integer(kind=i4) :: l_cols2
   integer(kind=i4) :: header(8)
   integer(kind=i4) :: i_state

   integer(kind=i4) :: n_kpt
   real(kind=r8) :: k_wt(1)
   real(kind=r8)  :: occ(100,1,1) !< Occupation members
   real(kind=r8)  :: mu !< Chemical potential

   real(kind=r8) :: n_electrons
   real(kind=r8) :: tol
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   logical :: file_exist

   real(kind=r8), allocatable :: eval(:,:,:)

   type(elsi_handle) :: eh

   integer(kind=i4), external :: numroc

   ! Reference values
   real(kind=r8), parameter :: e_ref = 0.31442451613_r8

   character(len=*), parameter :: file_name = "elsi.in"

   tol = 1.0e-8_r8
   header(:) = 0

   if(myid == 0) then
      write(*,"(2X,A)") "################################"
      write(*,"(2X,A)") "##     ELSI TEST PROGRAMS     ##"
      write(*,"(2X,A)") "################################"
      write(*,*)
      write(*,*)
   end if

   n_kpt = 1
   n_basis = 100
   n_electrons = 100
   allocate(eval(n_basis,1,1))
   do i_state = 1, 100
     eval(i_state,1,1) = dble(i_state)
   end do
   k_wt = 1.0


   ! Initialize ELSI
   call elsi_init(eh,1,1,0,n_basis,n_electrons,n_basis)
   call elsi_set_mpi(eh,comm)
   eh%ph%occ_non_aufbau = .true.
   eh%ph%mu_width = mu_width

   ! Customize ELSI
   call elsi_set_output(eh,2)
   call elsi_set_output_log(eh,1)

   call elsi_compute_mu_and_occ(eh,n_electrons,n_basis,1,n_kpt,k_wt,eval,occ,mu)


   if(myid == 0) then
      write(*,"(2X,A)") "Finished occ normal"
      write(*,*) "#chemical potential mu", mu
      write(*,*) "# state eigenvalue occ"
      do i_state = 1, 100
        write(*,*) i_state, eval(i_state,1,1), occ(i_state,1,1)
      end do
      write(*,*)
   end if

   ! Finalize ELSI
   call elsi_finalize(eh)

   deallocate(eval)

end subroutine
