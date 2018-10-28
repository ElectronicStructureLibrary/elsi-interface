! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests complex density matrix solver, BLACS_DENSE format.
!!
subroutine test_dm_cmplx_den(mpi_comm,solver,h_file,s_file)

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
   integer(kind=i4) :: n_basis
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols
   integer(kind=i4) :: header(8)

   real(kind=r8) :: n_electrons
   real(kind=r8) :: n_test
   real(kind=r8) :: tmp
   real(kind=r8) :: e_test = 0.0_r8
   real(kind=r8) :: e_ref = 0.0_r8
   real(kind=r8) :: tol = 0.0_r8
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   complex(kind=r8), allocatable :: ham(:,:)
   complex(kind=r8), allocatable :: ham_save(:,:)
   complex(kind=r8), allocatable :: ovlp(:,:)
   complex(kind=r8), allocatable :: ovlp_save(:,:)
   complex(kind=r8), allocatable :: dm(:,:)
   complex(kind=r8), allocatable :: edm(:,:)

   complex(kind=r8), external :: zdotc

   type(elsi_handle) :: eh
   type(elsi_rw_handle) :: rwh

   ! Reference values
   real(kind=r8), parameter :: e_elpa = -2622.88214509316_r8
   real(kind=r8), parameter :: e_omm = -2622.88214509316_r8
   real(kind=r8), parameter :: e_pexsi = -2622.88194292325_r8
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
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + ELPA"
         e_ref = e_elpa
      else if(solver == 2) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + libOMM"
         e_ref = e_omm
      else if(solver == 3) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + PEXSI"
         e_ref = e_pexsi
         tol = 1.0e-3_r8
      else if(solver == 6) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex + NTPoly"
         e_ref = e_ntpoly
         tol = 1.0e-7_r8
      end if
      write(*,*)
   end if

   ! Set up square-like processor grid
   do npcol = nint(sqrt(real(n_proc))),2,-1
      if(mod(n_proc,npcol) == 0) exit
   end do
   nprow = n_proc/npcol

   ! Set block size
   blk = 32

   ! Set up BLACS
   blacs_ctxt = mpi_comm
   call BLACS_Gridinit(blacs_ctxt,'r',nprow,npcol)

   ! Read H and S matrices
   call elsi_init_rw(rwh,0,1,0,0.0_r8)
   call elsi_set_rw_mpi(rwh,mpi_comm)
   call elsi_set_rw_blacs(rwh,blacs_ctxt,blk)

   call elsi_read_mat_dim(rwh,h_file,n_electrons,n_basis,l_rows,l_cols)
   call elsi_get_rw_header(rwh,header)

   allocate(ham(l_rows,l_cols))
   allocate(ham_save(l_rows,l_cols))
   allocate(ovlp(l_rows,l_cols))
   allocate(ovlp_save(l_rows,l_cols))
   allocate(dm(l_rows,l_cols))
   allocate(edm(l_rows,l_cols))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex(rwh,h_file,ham)
   call elsi_read_mat_complex(rwh,s_file,ovlp)

   call elsi_finalize_rw(rwh)

   ham_save = ham
   ovlp_save = ovlp

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished reading H and S matrices"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Initialize ELSI
   n_states = int(n_electrons,kind=i4)

   call elsi_init(eh,solver,1,0,n_basis,n_electrons,n_states)
   call elsi_set_mpi(eh,mpi_comm)
   call elsi_set_blacs(eh,blacs_ctxt,blk)

   ! Customize ELSI
   call elsi_set_output(eh,2)
   call elsi_set_output_log(eh,1)
   call elsi_set_sing_check(eh,0)
   call elsi_set_mu_broaden_width(eh,1.0e-6_r8)
   call elsi_set_omm_n_elpa(eh,1)
   call elsi_set_pexsi_delta_e(eh,80.0_r8)
   call elsi_set_pexsi_np_per_pole(eh,2)

   t1 = MPI_Wtime()

   ! Solve
   call elsi_dm_complex(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #1"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ham = ham_save

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #2"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ham = ham_save

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #3"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Reinit for a new geometry
   call elsi_reinit(eh)
   call elsi_set_elpa_solver(eh,1)
   call elsi_set_omm_flavor(eh,2)
   call elsi_set_ntpoly_method(eh,1)

   ham = ham_save
   ovlp = ovlp_save

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #4"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ham = ham_save

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex(eh,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   ! Compute energy density matrix
   call elsi_get_edm_complex(eh,edm)

   ! Compute electron count
   tmp = real(zdotc(l_rows*l_cols,ovlp_save,1,dm,1),kind=r8)

   call MPI_Reduce(tmp,n_test,1,mpi_real8,mpi_sum,0,mpi_comm,ierr)

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #5"
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
         end if
      end if
      write(*,*)
   end if

   ! Finalize ELSI
   call elsi_finalize(eh)

   deallocate(ham)
   deallocate(ham_save)
   deallocate(ovlp)
   deallocate(ovlp_save)
   deallocate(dm)
   deallocate(edm)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)

end subroutine
