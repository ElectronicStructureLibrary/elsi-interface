! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests real density matrix solver, BLACS_DENSE format.
!!
subroutine test_dm_real_den(mpi_comm,solver,h_file,s_file)

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   integer(kind=i4), intent(in) :: mpi_comm
   integer(kind=i4), intent(in) :: solver
   character(len=*), intent(in) :: h_file
   character(len=*), intent(in) :: s_file

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: myid
   integer(kind=i4) :: ierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: n_states
   integer(kind=i4) :: matrix_size
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols
   integer(kind=i4) :: header(8)

   real(kind=r8) :: n_electrons
   real(kind=r8) :: n_test
   real(kind=r8) :: tmp
   real(kind=r8) :: e_test = 0.0_r8
   real(kind=r8) :: e_ref  = 0.0_r8
   real(kind=r8) :: tol    = 0.0_r8
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   real(kind=r8), allocatable :: ham(:,:)
   real(kind=r8), allocatable :: ham_save(:,:)
   real(kind=r8), allocatable :: ovlp(:,:)
   real(kind=r8), allocatable :: ovlp_save(:,:)
   real(kind=r8), allocatable :: dm(:,:)
   real(kind=r8), allocatable :: edm(:,:)

   real(kind=r8), external :: ddot

   type(elsi_handle)    :: e_h
   type(elsi_rw_handle) :: rw_h

   ! Reference values
   real(kind=r8), parameter :: e_elpa   = -2622.88214509316_r8
   real(kind=r8), parameter :: e_omm    = -2622.88214509316_r8
   real(kind=r8), parameter :: e_pexsi  = -2622.88194292325_r8
   real(kind=r8), parameter :: e_sips   = -2622.88214509316_r8
   real(kind=r8), parameter :: e_ntpoly = -2622.88214509311_r8

   call MPI_Comm_size(mpi_comm,n_proc,ierr)
   call MPI_Comm_rank(mpi_comm,myid,ierr)

   if(myid == 0) then
      tol = 1.0e-8_r8
      write(*,"(2X,A)") "################################"
      write(*,"(2X,A)") "##     ELSI TEST PROGRAMS     ##"
      write(*,"(2X,A)") "################################"
      write(*,*)
      if(solver == 1) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_real + ELPA"
         e_ref = e_elpa
      elseif(solver == 2) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_real + libOMM"
         e_ref = e_omm
      elseif(solver == 3) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_real + PEXSI"
         e_ref = e_pexsi
         tol   = 1.0e-3_r8
      elseif(solver == 5) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_real + SLEPc-SIPs"
         e_ref = e_sips
      elseif(solver == 6) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_real + NTPoly"
         e_ref = e_ntpoly
      endif
      write(*,*)
   endif

   ! Set up square-like processor grid
   do npcol = nint(sqrt(real(n_proc))),2,-1
      if(mod(n_proc,npcol) == 0) exit
   enddo
   nprow = n_proc/npcol

   ! Set block size
   blk = 32

   ! Set up BLACS
   blacs_ctxt = mpi_comm
   call BLACS_Gridinit(blacs_ctxt,'r',nprow,npcol)

   ! Read H and S matrices
   call elsi_init_rw(rw_h,0,1,0,0.0_r8)
   call elsi_set_rw_mpi(rw_h,mpi_comm)
   call elsi_set_rw_blacs(rw_h,blacs_ctxt,blk)

   call elsi_read_mat_dim(rw_h,h_file,n_electrons,matrix_size,l_rows,l_cols)
   call elsi_get_rw_header(rw_h,header)

   allocate(ham(l_rows,l_cols))
   allocate(ham_save(l_rows,l_cols))
   allocate(ovlp(l_rows,l_cols))
   allocate(ovlp_save(l_rows,l_cols))
   allocate(dm(l_rows,l_cols))
   allocate(edm(l_rows,l_cols))

   t1 = MPI_Wtime()

   call elsi_read_mat_real(rw_h,h_file,ham)
   call elsi_read_mat_real(rw_h,s_file,ovlp)

   call elsi_finalize_rw(rw_h)

   ham_save  = ham
   ovlp_save = ovlp

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished reading H and S matrices"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   endif

   ! Initialize ELSI
   n_states = int(n_electrons,kind=i4)

   call elsi_init(e_h,solver,1,0,matrix_size,n_electrons,n_states)
   call elsi_set_mpi(e_h,mpi_comm)
   call elsi_set_blacs(e_h,blacs_ctxt,blk)

   ! Customize ELSI
   call elsi_set_output(e_h,2)
   call elsi_set_output_log(e_h,1)
   call elsi_set_sing_check(e_h,0)
   call elsi_set_mu_broaden_width(e_h,1.0e-6_r8)
   call elsi_set_omm_n_elpa(e_h,1)
   call elsi_set_pexsi_delta_e(e_h,80.0_r8)
   call elsi_set_pexsi_np_per_pole(e_h,2)
   call elsi_set_sips_n_elpa(e_h,1)

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 1)
   call elsi_dm_real(e_h,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #1"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   endif

   ham = ham_save

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 2, with the same H)
   call elsi_dm_real(e_h,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   ! Compute energy density matrix
   call elsi_get_edm_real(e_h,edm)

   ! Compute electron count
   tmp = ddot(l_rows*l_cols,ovlp_save,1,dm,1)

   call MPI_Reduce(tmp,n_test,1,mpi_real8,mpi_sum,0,mpi_comm,ierr)

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #2"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
      write(*,"(2X,A)") "Finished test program"
      write(*,*)
      if(header(8) == 1111) then
         write(*,"(2X,A)") "Band energy"
         write(*,"(2X,A,F15.8)") "| This test :",e_test
         write(*,"(2X,A,F15.8)") "| Reference :",e_ref
         write(*,*)
         write(*,"(2X,A)") "Electron count"
         write(*,"(2X,A,F15.8)") "| This test :",n_test
         write(*,"(2X,A,F15.8)") "| Reference :",n_electrons
         write(*,*)
         if(abs(e_test-e_ref) < tol .and. abs(n_test-n_electrons) < tol) then
            write(*,"(2X,A)") "Passed."
         else
            write(*,"(2X,A)") "Failed."
         endif
      endif
      write(*,*)
   endif

   ! Finalize ELSI
   call elsi_finalize(e_h)

   deallocate(ham)
   deallocate(ham_save)
   deallocate(ovlp)
   deallocate(ovlp_save)
   deallocate(dm)
   deallocate(edm)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)

end subroutine
