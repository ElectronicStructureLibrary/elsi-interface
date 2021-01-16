! Copyright (c) 2015-2021, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This is the public interface module of ELSI.
!!
module ELSI

   use ELSI_DATATYPE, only: elsi_handle,elsi_rw_handle
   use ELSI_GEO
   use ELSI_GET
   use ELSI_INPUT
   use ELSI_RW
   use ELSI_SET
   use ELSI_SETUP
   use ELSI_SOLVER

   implicit none

   private

   ! Data type
   public :: elsi_handle
   public :: elsi_rw_handle

   ! Setup
   public :: elsi_init
   public :: elsi_set_mpi
   public :: elsi_set_mpi_global
   public :: elsi_set_spin
   public :: elsi_set_kpoint
   public :: elsi_set_blacs
   public :: elsi_set_csc
   public :: elsi_set_csc_blk
   public :: elsi_set_coo
   public :: elsi_reinit
   public :: elsi_finalize

   ! Set
   public :: elsi_set_input_file
   public :: elsi_set_output
   public :: elsi_set_output_unit
   public :: elsi_set_output_log
   public :: elsi_set_save_ovlp
   public :: elsi_set_unit_ovlp
   public :: elsi_set_zero_def
   public :: elsi_set_sparsity_mask
   public :: elsi_set_illcond_check
   public :: elsi_set_illcond_tol
   public :: elsi_set_spin_degeneracy
   public :: elsi_set_energy_gap
   public :: elsi_set_spectrum_width
   public :: elsi_set_dimensionality
   public :: elsi_set_elpa_solver
   public :: elsi_set_elpa_n_single
   public :: elsi_set_elpa_gpu
   public :: elsi_set_elpa_autotune
   public :: elsi_set_omm_flavor
   public :: elsi_set_omm_n_elpa
   public :: elsi_set_omm_tol
   public :: elsi_set_pexsi_method
   public :: elsi_set_pexsi_n_mu
   public :: elsi_set_pexsi_n_pole
   public :: elsi_set_pexsi_np_per_pole
   public :: elsi_set_pexsi_np_symbo
   public :: elsi_set_pexsi_temp
   public :: elsi_set_pexsi_mu_min
   public :: elsi_set_pexsi_mu_max
   public :: elsi_set_pexsi_inertia_tol
   public :: elsi_set_eigenexa_method
   public :: elsi_set_sips_n_elpa
   public :: elsi_set_sips_n_slice
   public :: elsi_set_sips_inertia_tol
   public :: elsi_set_sips_ev_min
   public :: elsi_set_sips_ev_max
   public :: elsi_set_ntpoly_method
   public :: elsi_set_ntpoly_isr
   public :: elsi_set_ntpoly_tol
   public :: elsi_set_ntpoly_filter
   public :: elsi_set_ntpoly_max_iter
   public :: elsi_set_magma_solver
   public :: elsi_set_mu_broaden_scheme
   public :: elsi_set_mu_broaden_width
   public :: elsi_set_mu_tol
   public :: elsi_set_mu_mp_order

   ! Get
   public :: elsi_get_version
   public :: elsi_get_datestamp
   public :: elsi_get_initialized
   public :: elsi_get_n_illcond
   public :: elsi_get_ovlp_ev_min
   public :: elsi_get_ovlp_ev_max
   public :: elsi_get_pexsi_mu_min
   public :: elsi_get_pexsi_mu_max
   public :: elsi_get_mu
   public :: elsi_get_entropy
   public :: elsi_get_edm_real
   public :: elsi_get_edm_complex
   public :: elsi_get_edm_real_sparse
   public :: elsi_get_edm_complex_sparse
   public :: elsi_get_eval
   public :: elsi_get_evec_real
   public :: elsi_get_evec_complex
   public :: elsi_get_occ

   ! Deprecated
   public :: elsi_set_write_unit
   public :: elsi_set_sing_check
   public :: elsi_set_sing_tol
   public :: elsi_set_elpa_gpu_kernels
   public :: elsi_set_pexsi_gap
   public :: elsi_set_pexsi_delta_e
   public :: elsi_set_sips_interval
   public :: elsi_set_mu_spin_degen
   public :: elsi_set_uuid
   public :: elsi_get_n_sing

   ! Solve
   public :: elsi_ev_real
   public :: elsi_ev_complex
   public :: elsi_ev_real_sparse
   public :: elsi_ev_complex_sparse
   public :: elsi_dm_real
   public :: elsi_dm_complex
   public :: elsi_dm_real_sparse
   public :: elsi_dm_complex_sparse
   public :: elsi_bse_real
   public :: elsi_bse_complex

   ! Tool
   public :: elsi_orthonormalize_ev_real
   public :: elsi_orthonormalize_ev_complex
   public :: elsi_orthonormalize_ev_real_sparse
   public :: elsi_orthonormalize_ev_complex_sparse
   public :: elsi_extrapolate_dm_real
   public :: elsi_extrapolate_dm_complex
   public :: elsi_extrapolate_dm_real_sparse
   public :: elsi_extrapolate_dm_complex_sparse
   public :: elsi_compute_dm_real
   public :: elsi_compute_dm_complex
   public :: elsi_compute_edm_real
   public :: elsi_compute_edm_complex
   public :: elsi_compute_mu_and_occ
   public :: elsi_compute_entropy

   ! Read and write matrix
   public :: elsi_init_rw
   public :: elsi_finalize_rw
   public :: elsi_set_rw_mpi
   public :: elsi_set_rw_blacs
   public :: elsi_set_rw_csc
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
