! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests elsi_dm_complex.
!!
subroutine test_dm_complex(mpi_comm,solver,h_file,s_file)

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   integer(kind=i4), intent(in) :: mpi_comm
   integer(kind=i4), intent(in) :: solver
   character(*),     intent(in) :: h_file
   character(*),     intent(in) :: s_file

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: myid
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: n_states
   integer(kind=i4) :: matrix_size
   integer(kind=i4) :: l_rows
   integer(kind=i4) :: l_cols

   real(kind=r8) :: n_electrons
   real(kind=r8) :: e_test = 0.0_r8
   real(kind=r8) :: e_ref  = 0.0_r8
   real(kind=r8) :: e_tol  = 0.0_r8
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   complex(kind=r8), allocatable :: ham(:,:)
   complex(kind=r8), allocatable :: ham_save(:,:)
   complex(kind=r8), allocatable :: ovlp(:,:)
   complex(kind=r8), allocatable :: dm(:,:)
   complex(kind=r8), allocatable :: edm(:,:)

   type(elsi_handle)    :: e_h
   type(elsi_rw_handle) :: rw_h

   ! Reference values from calculations on November 20, 2017.
   real(kind=r8), parameter :: e_elpa  = -2622.88214509316_r8
   real(kind=r8), parameter :: e_omm   = -2622.88214509316_r8
   real(kind=r8), parameter :: e_pexsi = -2622.88194292325_r8

   call MPI_Comm_size(mpi_comm,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm,myid,mpierr)

   if(myid == 0) then
      e_tol = 1.0e-8_r8
      write(*,'("  ################################")')
      write(*,'("  ##     ELSI TEST PROGRAMS     ##")')
      write(*,'("  ################################")')
      write(*,*)
      if(solver == 1) then
         write(*,'("  Now start testing  elsi_dm_complex + ELPA")')
         e_ref = e_elpa
      elseif(solver == 2) then
         write(*,'("  Now start testing  elsi_dm_complex + libOMM")')
         e_ref = e_omm
      elseif(solver == 3) then
         write(*,'("  Now start testing  elsi_dm_complex + PEXSI")')
         e_ref = e_pexsi
         e_tol = 1.0e-4_r8
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

   allocate(ham(l_rows,l_cols))
   allocate(ham_save(l_rows,l_cols))
   allocate(ovlp(l_rows,l_cols))
   allocate(dm(l_rows,l_cols))
   allocate(edm(l_rows,l_cols))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex(rw_h,h_file,ham)
   call elsi_read_mat_complex(rw_h,s_file,ovlp)

   call elsi_finalize_rw(rw_h)

   ham_save = ham

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished reading H and S matrices")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
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

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 1)
   call elsi_dm_complex(e_h,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished SCF #1")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ham = ham_save

   t1 = MPI_Wtime()

   ! Solve (pseudo SCF 2, with the same H)
   call elsi_dm_complex(e_h,ham,ovlp,dm,e_test)

   t2 = MPI_Wtime()

   ! Compute energy density matrix
   if(solver == 1 .or. solver == 2 .or. solver == 3) then
      call elsi_get_edm_complex(e_h,edm)
   endif

   if(myid == 0) then
      write(*,'("  Finished SCF #2")')
      write(*,'("  | Time :",F10.3,"s")') t2-t1
      write(*,*)
      write(*,'("  Finished test program")')
      write(*,*)
      if(abs(e_test-e_ref) < e_tol) then
         write(*,'("  Passed.")')
      endif
      write(*,*)
   endif

   ! Finalize ELSI
   call elsi_finalize(e_h)

   deallocate(ham)
   deallocate(ham_save)
   deallocate(ovlp)
   deallocate(dm)
   deallocate(edm)

   call BLACS_Gridexit(blacs_ctxt)
   call BLACS_Exit(1)

end subroutine
