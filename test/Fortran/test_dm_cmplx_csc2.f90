! Copyright (c) 2015-2020, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This subroutine tests complex density matrix solver, SIESTA_CSC format.
!!
subroutine test_dm_cmplx_csc2(mpi_comm,solver,h_file,s_file)

   use ELSI_PRECISION, only: r8,i4
   use ELSI

   implicit none

   include "mpif.h"

   integer(kind=i4), intent(in) :: mpi_comm
   integer(kind=i4), intent(in) :: solver
   character(len=*), intent(in) :: h_file
   character(len=*), intent(in) :: s_file

   integer(kind=i4) :: n_proc
   integer(kind=i4) :: myid
   integer(kind=i4) :: ierr
   integer(kind=i4) :: blk
   integer(kind=i4) :: n_states
   integer(kind=i4) :: n_basis
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: n_l_cols
   integer(kind=i4) :: header(8)

   real(kind=r8) :: n_electrons
   real(kind=r8) :: n_sp
   real(kind=r8) :: tmp
   real(kind=r8) :: e_hp
   real(kind=r8) :: e_sq
   real(kind=r8) :: e_ref
   real(kind=r8) :: tol
   real(kind=r8) :: t1
   real(kind=r8) :: t2

   logical :: file_exist

   complex(kind=r8), allocatable :: ham(:)
   complex(kind=r8), allocatable :: ovlp(:)
   complex(kind=r8), allocatable :: dm(:)
   complex(kind=r8), allocatable :: edm(:)

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)

   type(elsi_handle) :: eh
   type(elsi_rw_handle) :: rwh

   ! Reference values
   real(kind=r8), parameter :: e_elpa = -2622.88214509316_r8
   real(kind=r8), parameter :: e_pexsi = -2622.88194292325_r8

   character(len=*), parameter :: file_name = "elsi.in"

   complex(kind=r8), external :: zdotc

   call MPI_Comm_size(mpi_comm,n_proc,ierr)
   call MPI_Comm_rank(mpi_comm,myid,ierr)

   e_hp = 0.0_r8
   e_sq = 0.0_r8
   e_ref = e_elpa
   tol = 1.0e-8_r8
   header(:) = 0

   if(myid == 0) then
      write(*,"(2X,A)") "################################"
      write(*,"(2X,A)") "##     ELSI TEST PROGRAMS     ##"
      write(*,"(2X,A)") "################################"
      write(*,*)
      if(solver == 1) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex_sparse + ELPA"
      else if(solver == 2) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex_sparse + libOMM"
      else if(solver == 3) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex_sparse + PEXSI"
         e_ref = e_pexsi
         tol = 1.0e-3_r8
      else if(solver == 6) then
         write(*,"(2X,A)") "Now start testing  elsi_dm_complex_sparse + NTPoly"
         tol = 1.0e-7_r8
      end if
      write(*,*)
   end if

   ! Read H and S matrices
   call elsi_init_rw(rwh,0,1,0,0.0_r8)
   call elsi_set_rw_mpi(rwh,mpi_comm)

   call elsi_read_mat_dim_sparse(rwh,h_file,n_electrons,n_basis,nnz_g,nnz_l,&
        n_l_cols)
   call elsi_get_rw_header(rwh,header)

   ! Below is a requirement of this test (not a requirement of ELSI)
   if(mod(n_basis,n_proc) /= 0) then
      blk = 0

      if(myid == 0) then
         write(*,"(2X,A)") "################################################"//&
            "######################"
         write(*,"(2X,A)") "##  This test requires number of MPI tasks to be"//&
         " a divisor of 288!! ##"
         write(*,"(2X,A)") "################################################"//&
            "######################"
      end if
      call MPI_Abort(mpi_comm,0,ierr)
      stop
   else
      blk = n_basis/n_proc
   end if

   allocate(ham(nnz_l))
   allocate(ovlp(nnz_l))
   allocate(dm(nnz_l))
   allocate(edm(nnz_l))
   allocate(row_ind(nnz_l))
   allocate(col_ptr(n_l_cols+1))

   t1 = MPI_Wtime()

   call elsi_read_mat_complex_sparse(rwh,h_file,row_ind,col_ptr,ham)
   call elsi_read_mat_complex_sparse(rwh,s_file,row_ind,col_ptr,ovlp)

   call elsi_finalize_rw(rwh)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished reading H and S matrices"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Initialize ELSI
   n_states = min(int(n_electrons,kind=i4),n_basis)

   call elsi_init(eh,solver,1,2,n_basis,n_electrons,n_states)
   call elsi_set_mpi(eh,mpi_comm)
   call elsi_set_csc(eh,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)
   call elsi_set_csc_blk(eh,blk)

   ! Customize ELSI
   call elsi_set_output(eh,2)
   call elsi_set_output_log(eh,1)
   call elsi_set_illcond_check(eh,0)
   call elsi_set_mu_broaden_scheme(eh,2)
   call elsi_set_mu_mp_order(eh,1)
   call elsi_set_mu_broaden_width(eh,1.0e-6_r8)
   call elsi_set_omm_n_elpa(eh,1)
   call elsi_set_pexsi_delta_e(eh,80.0_r8)
   call elsi_set_pexsi_np_per_pole(eh,2)

   inquire(file=file_name,exist=file_exist)

   if(file_exist) then
      call elsi_set_input_file(eh,file_name)
   end if

   t1 = MPI_Wtime()

   ! Solve
   call elsi_dm_complex_sparse(eh,ham,ovlp,dm,e_hp)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #1"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex_sparse(eh,ham,ovlp,dm,e_hp)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #2"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex_sparse(eh,ham,ovlp,dm,e_hp)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #3"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   ! Reinit for a new geometry
   call elsi_reinit(eh)
   call elsi_set_csc(eh,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)
   call elsi_set_elpa_solver(eh,1)
   call elsi_set_omm_flavor(eh,2)
   call elsi_set_ntpoly_method(eh,1)

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex_sparse(eh,ham,ovlp,dm,e_hp)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #4"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
   end if

   t1 = MPI_Wtime()

   ! Solve again
   call elsi_dm_complex_sparse(eh,ham,ovlp,dm,e_hp)

   t2 = MPI_Wtime()

   ! Compute energy density matrix
   call elsi_get_edm_complex_sparse(eh,edm)

   ! Compute band energy
   tmp = real(zdotc(nnz_l,ovlp,1,edm,1),kind=r8)

   call MPI_Reduce(tmp,e_sq,1,mpi_real8,mpi_sum,0,mpi_comm,ierr)

   ! Compute electron count
   tmp = real(zdotc(nnz_l,ovlp,1,dm,1),kind=r8)

   call MPI_Reduce(tmp,n_sp,1,mpi_real8,mpi_sum,0,mpi_comm,ierr)

   if(myid == 0) then
      write(*,"(2X,A)") "Finished SCF #5"
      write(*,"(2X,A,F10.3,A)") "| Time :",t2-t1,"s"
      write(*,*)
      write(*,"(2X,A)") "Finished test program"
      write(*,*)
      if(header(8) == 1111) then
         write(*,"(2X,A)") "Band energy"
         write(*,"(2X,A,F15.8)") "| Tr[H*P]   :",e_hp
         write(*,"(2X,A,F15.8)") "| Tr[S*Q]   :",e_sq
         write(*,"(2X,A,F15.8)") "| Reference :",e_ref
         write(*,*)
         write(*,"(2X,A)") "Electron count"
         write(*,"(2X,A,F15.8)") "| Tr[S*P]   :",n_sp
         write(*,"(2X,A,F15.8)") "| Reference :",n_electrons
         write(*,*)
         if(abs(e_hp-e_ref) < tol .and. abs(e_sq-e_ref) < tol&
            .and. abs(n_sp-n_electrons) < tol) then
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
   deallocate(ovlp)
   deallocate(dm)
   deallocate(edm)
   deallocate(row_ind)
   deallocate(col_ptr)

end subroutine
