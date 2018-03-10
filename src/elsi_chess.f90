! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to CheSS.
!!
module ELSI_CHESS

   use ELSI_DATATYPE,          only: elsi_handle
   use ELSI_IO,                only: elsi_say
   use ELSI_MPI,               only: elsi_check_mpi
   use ELSI_PRECISION,         only: r8,i4
   use ELSI_TIMINGS,           only: elsi_get_time
   use FOE_BASE,               only: foe_data_get_real
   use FOE_COMMON,             only: init_foe
   use SPARSEMATRIX_BASE,      only: sparsematrix_init_errors,&
                                     sparsematrix_initialize_timing_categories,&
                                     SPARSE_TASKGROUP
   use SPARSEMATRIX_HIGHLEVEL, only: matrices_init,&
                                     matrix_fermi_operator_expansion,&
                                     sparse_matrix_init_from_data_ccs
   use SPARSEMATRIX_INIT,      only: init_matrix_taskgroups_wrapper

   implicit none

   private

   public :: elsi_set_chess_default
   public :: elsi_init_chess
   public :: elsi_solve_evp_chess_real
   public :: elsi_compute_edm_chess_real

contains

!>
!! This routine initializes CheSS.
!!
subroutine elsi_init_chess(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   real(kind=r8) :: n_electron(1)

   character(len=40), parameter :: caller = "elsi_init_chess"

   if(e_h%handle_ready) then
      e_h%handle_changed = .true.
   endif

   if(e_h%n_elsi_calls == 1) then
      ! Initialize f_lib
      call f_lib_initialize_stub()

      ! Initialize sparsematrix error handling and timing
      call sparsematrix_init_errors()
      call sparsematrix_initialize_timing_categories()

      ! Initialize sparse matrices
      call sparse_matrix_init_from_data_ccs(e_h%myid,e_h%n_procs,e_h%mpi_comm,&
              e_h%n_basis,e_h%nnz_g,e_h%row_ind_chess,e_h%col_ptr_chess,&
              e_h%sparse_mat_chess(1),init_matmul=.false.)

      ! TODO: Sparsity buffer
      call sparse_matrix_init_from_data_ccs(e_h%myid,e_h%n_procs,e_h%mpi_comm,&
              e_h%n_basis,e_h%nnz_g,e_h%row_ind_chess,e_h%col_ptr_chess,&
              e_h%sparse_mat_chess(2),init_matmul=.true.,nvctr_mult=e_h%nnz_g,&
              row_ind_mult=e_h%row_ind_chess,&
              col_ptr_mult=e_h%col_ptr_chess(1:e_h%n_basis))

      ! Initialize task groups
      call init_matrix_taskgroups_wrapper(e_h%myid,e_h%n_procs,e_h%mpi_comm,&
              .false.,2,e_h%sparse_mat_chess)

      n_electron = e_h%n_electrons

      ! Initialize FOE objects
      call init_foe(e_h%myid,e_h%n_procs,1,n_electron,e_h%chess_foe,&
              fscale=e_h%chess_erf_decay,fscale_lowerbound=e_h%chess_erf_min,&
              fscale_upperbound=e_h%chess_erf_max,evlow=e_h%chess_ev_ham_min,&
              evhigh=e_h%chess_ev_ham_max,betax=e_h%chess_beta)

      call init_foe(e_h%myid,e_h%n_procs,1,n_electron,e_h%chess_ice,&
              evlow=e_h%chess_ev_ovlp_min,evhigh=e_h%chess_ev_ovlp_max,&
              betax=e_h%chess_beta)

      ! Allocate CheSS matrices
      call matrices_init(e_h%sparse_mat_chess(1),e_h%ham_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(e_h%sparse_mat_chess(1),e_h%ovlp_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(e_h%sparse_mat_chess(2),e_h%dm_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(e_h%sparse_mat_chess(2),e_h%edm_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(e_h%sparse_mat_chess(2),e_h%ovlp_inv_sqrt_chess(1),&
              matsize=SPARSE_TASKGROUP)

      e_h%ovlp_chess%matrix_compr = e_h%ovlp_real_chess

      e_h%chess_started = .true.
   endif

   e_h%ham_chess%matrix_compr = e_h%ham_real_chess

end subroutine

!>
!! This routine interfaces to CheSS.
!!
subroutine elsi_solve_evp_chess_real(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   logical            :: calc_ovlp_inv_sqrt
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   integer(kind=i4)   :: ierr
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_solve_evp_chess_real"

   call elsi_get_time(t0)

   if(e_h%n_elsi_calls == 1) then
      calc_ovlp_inv_sqrt = .true.
   else
      calc_ovlp_inv_sqrt = .false.
   endif

   call elsi_say(e_h,"  Starting CheSS density matrix solver")

   call matrix_fermi_operator_expansion(e_h%myid,e_h%n_procs,e_h%mpi_comm,&
           e_h%chess_foe,e_h%chess_ice,e_h%sparse_mat_chess(1),&
           e_h%sparse_mat_chess(1),e_h%sparse_mat_chess(2),e_h%ovlp_chess,&
           e_h%ham_chess,e_h%ovlp_inv_sqrt_chess,e_h%dm_chess,e_h%energy_hdm,&
           calculate_minusonehalf=calc_ovlp_inv_sqrt,foe_verbosity=0,&
           symmetrize_kernel=.true.,calculate_energy_density_kernel=.false.,&
           energy_kernel=e_h%edm_chess)

   e_h%mu = foe_data_get_real(e_h%chess_foe,"ef",1)

   call MPI_Barrier(e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Barrier",ierr,caller)

   call elsi_get_time(t1)

   write(info_str,"('  Finished density matrix calculation')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine computes the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_chess_real(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   character(len=40), parameter :: caller = "elsi_compute_edm_chess_real"

   ! TODO: Compute edm

end subroutine

!>
!! This routine sets default CheSS parameters.
!!
subroutine elsi_set_chess_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   character(len=40), parameter :: caller = "elsi_set_chess_default"

   if(e_h%handle_ready) then
      e_h%handle_changed = .true.
   endif

   ! Initial guess of the decay length of the error function
   e_h%chess_erf_decay = 5.0e-2_r8

   ! Lower bound of the decay length
   e_h%chess_erf_min = 5.0e-3_r8

   ! Upper bound of the decay length
   e_h%chess_erf_max = 5.0e-2_r8

   ! Lower bound of the eigenvalues of H
   e_h%chess_ev_ham_min = -2.0_r8

   ! Upper bound of the eigenvalues of H
   e_h%chess_ev_ham_max = 2.0_r8

   ! Lower bound of the eigenvalues of S
   e_h%chess_ev_ovlp_min = 1.0e-4_r8

   ! Upper bound of the eigenvalues of S
   e_h%chess_ev_ovlp_max = 2.0_r8

   ! A patameter used to estimate eigenspectrum
   e_h%chess_beta = -1.0e3_r8

end subroutine

end module ELSI_CHESS
