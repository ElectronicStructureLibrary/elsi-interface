! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides routines for solving the Kohn-Sham electronic structure
!! using ELPA, libOMM, PEXSI, SIPS, DMP.
!!
module ELSI

   use ELSI_DATATYPE, only: elsi_handle,elsi_io_handle,elsi_rw_handle
   use ELSI_IO
   use ELSI_MAT_IO
   use ELSI_OCC
   use ELSI_MUTATOR
   use ELSI_SETUP
   use ELSI_SOLVER

   implicit none

   private

   ! Data type
   public :: elsi_handle
   public :: elsi_io_handle
   public :: elsi_rw_handle

   ! Setup
   public :: elsi_init
   public :: elsi_finalize
   public :: elsi_set_mpi
   public :: elsi_set_mpi_global
   public :: elsi_set_spin
   public :: elsi_set_kpoint
   public :: elsi_set_blacs
   public :: elsi_set_csc

   ! Mutator
   public :: elsi_set_output
   public :: elsi_set_write_unit
   public :: elsi_set_unit_ovlp
   public :: elsi_set_zero_def
   public :: elsi_set_sing_check
   public :: elsi_set_sing_tol
   public :: elsi_set_sing_stop
   public :: elsi_set_uplo
   public :: elsi_set_csc_blk
   public :: elsi_set_elpa_solver
   public :: elsi_set_omm_flavor
   public :: elsi_set_omm_n_elpa
   public :: elsi_set_omm_tol
   public :: elsi_set_omm_ev_shift
   public :: elsi_set_pexsi_n_mu
   public :: elsi_set_pexsi_n_pole
   public :: elsi_set_pexsi_np_per_pole
   public :: elsi_set_pexsi_np_symbo
   public :: elsi_set_pexsi_ordering
   public :: elsi_set_pexsi_temp
   public :: elsi_set_pexsi_gap
   public :: elsi_set_pexsi_delta_e
   public :: elsi_set_pexsi_mu_min
   public :: elsi_set_pexsi_mu_max
   public :: elsi_set_pexsi_inertia_tol
   public :: elsi_set_sips_n_elpa
   public :: elsi_set_sips_n_slice
   public :: elsi_set_sips_slice_type
   public :: elsi_set_sips_buffer
   public :: elsi_set_sips_inertia_tol
   public :: elsi_set_sips_interval
   public :: elsi_set_sips_first_ev
   public :: elsi_set_dmp_method
   public :: elsi_set_dmp_max_step
   public :: elsi_set_dmp_tol
   public :: elsi_set_mu_broaden_scheme
   public :: elsi_set_mu_broaden_width
   public :: elsi_set_mu_tol
   public :: elsi_set_mu_spin_degen
   public :: elsi_set_mu_mp_order
   public :: elsi_set_output_log
   public :: elsi_set_log_unit
   public :: elsi_set_log_file
   public :: elsi_set_log_tag
   public :: elsi_set_uuid
   public :: elsi_set_calling_code
   public :: elsi_get_pexsi_mu_min
   public :: elsi_get_pexsi_mu_max
   public :: elsi_get_n_sing
   public :: elsi_get_mu
   public :: elsi_get_entropy
   public :: elsi_get_edm_real
   public :: elsi_get_edm_complex
   public :: elsi_get_edm_real_sparse
   public :: elsi_get_edm_complex_sparse

   ! Solver
   public :: elsi_ev_real
   public :: elsi_ev_complex
   public :: elsi_ev_real_sparse
   public :: elsi_ev_complex_sparse
   public :: elsi_dm_real
   public :: elsi_dm_complex
   public :: elsi_dm_real_sparse
   public :: elsi_dm_complex_sparse
   public :: elsi_compute_mu_and_occ
   public :: elsi_compute_entropy

   ! Reading and writing matrix
   public :: elsi_init_rw
   public :: elsi_finalize_rw
   public :: elsi_set_rw_mpi
   public :: elsi_set_rw_blacs
   public :: elsi_set_rw_csc
   public :: elsi_set_rw_output
   public :: elsi_set_rw_write_unit
   public :: elsi_set_rw_zero_def
   public :: elsi_set_rw_header
   public :: elsi_get_rw_header
   public :: elsi_read_mat_dim
   public :: elsi_read_mat_dim_sparse
   public :: elsi_read_mat_real
   public :: elsi_read_mat_real_sparse
   public :: elsi_read_mat_complex
   public :: elsi_read_mat_complex_sparse
   public :: elsi_write_mat_real
   public :: elsi_write_mat_real_sparse
   public :: elsi_write_mat_complex
   public :: elsi_write_mat_complex_sparse

end module ELSI
