include "mpi4py/mpi.pxi"

ctypedef elsi_handle__* elsi_t

ctypedef elsi_rw_handle__* elsi_rw_t

cdef import from "elsi.h":
    cdef struct elsi_handle__:
        pass
    cdef struct elsi_rw_handle__:
        pass

    void c_elsi_init(elsi_t *handle_c,
         int solver,
         int parallel_mode,
         int matrix_format,
         int n_basis,
         double n_electron,
         int n_state);


    void c_elsi_init(elsi_t *handle_c,
         int solver,
         int parallel_mode,
         int matrix_format,
         int n_basis,
         double n_electron,
         int n_state);
    
    void c_elsi_set_mpi(elsi_t handle_c,
         MPI_Fint comm);
    
    void c_elsi_set_mpi_global(elsi_t handle_c,
         MPI_Fint comm_global);
    
    void c_elsi_set_spin(elsi_t handle_c,
         int n_spin,
         int i_spin);
    
    void c_elsi_set_kpoint(elsi_t handle_c,
         int n_kpt,
         int i_kpt,
         double weight);
    
    void c_elsi_set_blacs(elsi_t handle_c,
         int blacs_ctxt,
         int block_size);
    
    void c_elsi_set_csc(elsi_t handle_c,
         int nnz,
         int nnz_l,
         int n_lcol,
         int *row_ind,
         int *col_ptr);
    
    void c_elsi_set_csc_blk(elsi_t handle_c,
         int blk);
    
    void c_elsi_set_coo(elsi_t handle_c,
         int nnz,
         int nnz_l,
         int *row_ind,
         int *col_ind);
    
    void c_elsi_reinit(elsi_t handle_c);
    
    void c_elsi_finalize(elsi_t handle_c);
    
    void c_elsi_ev_real(elsi_t handle_c,
         double *ham,
         double *ovlp,
         double *eval,
         double *evec);
    
    void c_elsi_ev_complex(elsi_t handle_c,
         double complex *ham,
         double complex *ovlp,
         double *eval,
         double complex *evec);
    
    void c_elsi_ev_real_sparse(elsi_t handle_c,
         double *ham,
         double *ovlp,
         double *eval,
         double *evec);
    
    void c_elsi_ev_complex_sparse(elsi_t handle_c,
         double complex *ham,
         double complex *ovlp,
         double *eval,
         double complex *evec);
    
    void c_elsi_dm_real(elsi_t handle_c,
         double *ham,
         double *ovlp,
         double *dm,
         double *energy);
    
    void c_elsi_dm_complex(elsi_t handle_c,
         double complex *ham,
         double complex *ovlp,
         double complex *dm,
         double *energy);
    
    void c_elsi_dm_real_sparse(elsi_t handle_c,
         double *ham,
         double *ovlp,
         double *dm,
         double *energy);
    
    void c_elsi_dm_complex_sparse(elsi_t handle_c,
         double complex *ham,
         double complex *ovlp,
         double complex *dm,
         double *energy);
    
    void c_elsi_bse_real(elsi_t handle_c,
         double *mat_a,
         double *mat_b,
         double *eval,
         double *evec);
    
    void c_elsi_bse_complex(elsi_t handle_c,
         double complex *mat_a,
         double complex *mat_b,
         double *eval,
         double complex *evec);
    
    void c_elsi_set_input_file(elsi_t handle_c,
         char *name_c);
    
    void c_elsi_set_output(elsi_t handle_c,
         int output);
    
    void c_elsi_set_output_log(elsi_t handle_c,
         int output_log);
    
    void c_elsi_set_save_ovlp(elsi_t handle_c,
         int save_ovlp);
    
    void c_elsi_set_unit_ovlp(elsi_t handle_c,
         int unit_ovlp);
    
    void c_elsi_set_zero_def(elsi_t handle_c,
         double zero_def);
    
    void c_elsi_set_sparsity_mask(elsi_t handle_c,
         int sparsity_mask);
    
    void c_elsi_set_illcond_check(elsi_t handle_c,
         int illcond_check);
    
    void c_elsi_set_illcond_tol(elsi_t handle_c,
         double illcond_tol);
    
    void c_elsi_set_sing_check(elsi_t handle_c,
         int sing_check);
    
    void c_elsi_set_sing_tol(elsi_t handle_c,
         double sing_tol);
    
    void c_elsi_set_spin_degeneracy(elsi_t handle_c,
         double spin_degeneracy);
    
    void c_elsi_set_energy_gap(elsi_t handle_c,
         double gap);
    
    void c_elsi_set_spectrum_width(elsi_t handle_c,
         double width);
    
    void c_elsi_set_dimensionality(elsi_t handle_c,
         int dimensionality);
    
    void c_elsi_set_elpa_solver(elsi_t handle_c,
         int solver);
    
    void c_elsi_set_elpa_n_single(elsi_t handle_c,
         int n_single);
    
    void c_elsi_set_elpa_gpu(elsi_t handle_c,
         int gpu);
    
    void c_elsi_set_elpa_gpu_kernels(elsi_t handle_c,
         int gpu_kernels);
    
    void c_elsi_set_elpa_autotune(elsi_t handle_c,
         int autotune);
    
    void c_elsi_set_omm_flavor(elsi_t handle_c,
         int flavor);
    
    void c_elsi_set_omm_n_elpa(elsi_t handle_c,
         int n_elpa);
    
    void c_elsi_set_omm_tol(elsi_t handle_c,
         double tol);
    
    void c_elsi_set_pexsi_method(elsi_t handle_c,
         int method);
    
    void c_elsi_set_pexsi_n_mu(elsi_t handle_c,
         int n_mu);
    
    void c_elsi_set_pexsi_n_pole(elsi_t handle_c,
         int n_pole);
    
    void c_elsi_set_pexsi_np_per_pole(elsi_t handle_c,
         int np_per_pole);
    
    void c_elsi_set_pexsi_np_symbo(elsi_t handle_c,
         int np_symbo);
    
    void c_elsi_set_pexsi_temp(elsi_t handle_c,
         double temp);
    
    void c_elsi_set_pexsi_gap(elsi_t handle_c,
         double gap);
    
    void c_elsi_set_pexsi_delta_e(elsi_t handle_c,
         double delta_e);
    
    void c_elsi_set_pexsi_mu_min(elsi_t handle_c,
         double mu_min);
    
    void c_elsi_set_pexsi_mu_max(elsi_t handle_c,
         double mu_max);
    
    void c_elsi_set_pexsi_inertia_tol(elsi_t handle_c,
         double inertia_tol);
    
    void c_elsi_set_eigenexa_method(elsi_t handle_c,
         int method);
    
    void c_elsi_set_sips_n_elpa(elsi_t handle_c,
         int n_elpa);
    
    void c_elsi_set_sips_n_slice(elsi_t handle_c,
         int n_slice);
    
    void c_elsi_set_sips_inertia_tol(elsi_t handle_c,
         double inertia_tol);
    
    void c_elsi_set_sips_ev_min(elsi_t handle_c,
         double ev_min);
    
    void c_elsi_set_sips_ev_max(elsi_t handle_c,
         double ev_max);
    
    void c_elsi_set_sips_interval(elsi_t handle_c,
         double lower,
         double upper);
    
    void c_elsi_set_ntpoly_method(elsi_t handle_c,
         int method);
    
    void c_elsi_set_ntpoly_isr(elsi_t handle_c,
         int isr);
    
    void c_elsi_set_ntpoly_tol(elsi_t handle_c,
         double tol);
    
    void c_elsi_set_ntpoly_filter(elsi_t handle_c,
         double filter);
    
    void c_elsi_set_ntpoly_max_iter(elsi_t handle_c,
         int max_iter);
    
    void c_elsi_set_magma_solver(elsi_t handle_c,
         int solver);
    
    void c_elsi_set_mu_broaden_scheme(elsi_t handle_c,
         int broaden_scheme);
    
    void c_elsi_set_mu_broaden_width(elsi_t handle_c,
         double broaden_width);
    
    void c_elsi_set_mu_tol(elsi_t handle_c,
         double tol);
    
    void c_elsi_set_mu_mp_order(elsi_t handle_c,
         int mp_order);
    
    void c_elsi_set_n_frozen(elsi_t handle_c,
         int n_frozen);
    
    void c_elsi_set_frozen_method(elsi_t handle_c,
         int frozen_method);
    
    void c_elsi_set_frozen_idx(elsi_t handle_c,
         int *frozen_idx);
    
    void c_elsi_get_version(int *major,
         int *minor,
         int *patch);
    
    void c_elsi_get_datestamp(int *datestamp);
    
    void c_elsi_get_initialized(elsi_t handle_c,
         int *initialized);
    
    void c_elsi_get_n_illcond(elsi_t handle_c,
         int *n_illcond);
    
    void c_elsi_get_n_sing(elsi_t handle_c,
         int *n_sing);
    
    void c_elsi_get_ovlp_ev_min(elsi_t handle_c,
         double *ev_min);
    
    void c_elsi_get_ovlp_ev_max(elsi_t handle_c,
         double *ev_max);
    
    void c_elsi_get_pexsi_mu_min(elsi_t handle_c,
         double *mu_min);
    
    void c_elsi_get_pexsi_mu_max(elsi_t handle_c,
         double *mu_max);
    
    void c_elsi_get_mu(elsi_t handle_c,
         double *mu);
    
    void c_elsi_get_entropy(elsi_t handle_c,
         double *entropy);
    
    void c_elsi_get_edm_real(elsi_t handle_c,
         double *edm);
    
    void c_elsi_get_edm_complex(elsi_t handle_c,
         double complex *edm);
    
    void c_elsi_get_edm_real_sparse(elsi_t handle_c,
         double *edm);
    
    void c_elsi_get_edm_complex_sparse(elsi_t handle_c,
         double complex *edm);
    
    void c_elsi_get_eval(elsi_t handle_c,
         double *eval);
    
    void c_elsi_get_evec_real(elsi_t handle_c,
         double *evec);
    
    void c_elsi_get_evec_complex(elsi_t handle_c,
         double complex *evec);
    
    void c_elsi_get_occ(elsi_t handle_c,
         double *occ);
    
    void c_elsi_orthonormalize_ev_real(elsi_t handle_c,
         double *ovlp,
         double *evec);
    
    void c_elsi_orthonormalize_ev_complex(elsi_t handle_c,
         double complex *ovlp,
         double complex *evec);
    
    void c_elsi_orthonormalize_ev_real_sparse(elsi_t handle_c,
         double *ovlp,
         double *evec);
    
    void c_elsi_orthonormalize_ev_complex_sparse(elsi_t handle_c,
         double complex *ovlp,
         double complex *evec);
    
    void c_elsi_extrapolate_dm_real(elsi_t handle_c,
         double *ovlp,
         double *dm);
    
    void c_elsi_extrapolate_dm_complex(elsi_t handle_c,
         double complex *ovlp,
         double complex *dm);
    
    void c_elsi_extrapolate_dm_real_sparse(elsi_t handle_c,
         double *ovlp,
         double *dm);
    
    void c_elsi_extrapolate_dm_complex_sparse(elsi_t handle_c,
         double complex *ovlp,
         double complex *dm);
    
    void c_elsi_init_rw(elsi_rw_t *handle_c,
         int rw_task,
         int parallel_mode,
         int n_basis,
         double n_electron);
    
    void c_elsi_set_rw_mpi(elsi_rw_t handle_c,
         MPI_Fint comm);
    
    void c_elsi_set_rw_blacs(elsi_rw_t handle_c,
         int blacs_ctxt,
         int block_size);
    
    void c_elsi_set_rw_csc(elsi_rw_t handle_c,
         int nnz,
         int nnz_l,
         int n_lcol);
    
    void c_elsi_finalize_rw(elsi_rw_t handle_c);
    
    void c_elsi_set_rw_zero_def(elsi_rw_t handle_c,
         double zero_def);
    
    void c_elsi_read_mat_dim(elsi_rw_t handle_c,
         char *name_c,
         double *n_electrons,
         int *n_basis,
         int *n_lrow,
         int *n_lcol);
    
    void c_elsi_read_mat_dim_sparse(elsi_rw_t handle_c,
         char *name_c,
         double *n_electrons,
         int *n_basis,
         int *nnz_g,
         int *nnz_l,
         int *n_lcol);
    
    void c_elsi_read_mat_real(elsi_rw_t handle_c,
         char *name_c,
         double *mat);
    
    void c_elsi_read_mat_real_sparse(elsi_rw_t handle_c,
         char *name_c,
         int *row_ind,
         int *col_ptr,
         double *mat);
    
    void c_elsi_write_mat_real(elsi_rw_t handle_c,
         char *name_c,
         double *mat);
    
    void c_elsi_write_mat_real_sparse(elsi_rw_t handle_c,
         char *name_c,
         int *row_ind,
         int *col_ptr,
         double *mat);
    
    void c_elsi_read_mat_complex(elsi_rw_t handle_c,
         char *name_c,
         double complex *mat);
    
    void c_elsi_read_mat_complex_sparse(elsi_rw_t handle_c,
         char *name_c,
         int *row_ind,
         int *col_ptr,
         double complex *mat);
    
    void c_elsi_write_mat_complex(elsi_rw_t handle_c,
         char *name_c,
         double complex *mat);
    
    void c_elsi_write_mat_complex_sparse(elsi_rw_t handle_c,
         char *name_c,
         int *row_ind,
         int *col_ptr,
         double complex *mat);
