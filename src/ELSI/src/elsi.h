#include <complex.h>

void c_elsi_init(int solver,
                 int parallel_mode,
                 int matrix_format,
                 int matrix_size,
                 double n_electrons_in,
                 int n_states_in);

void c_elsi_set_mpi(int mpi_comm_global_in);

void c_elsi_set_blacs(int icontext,
                      int block_size);

void c_elsi_finalize();

void c_elsi_set_sparsity(int nnz_g_in,
                         int nnz_l_in,
                         int nnz_l_cols_in,
                         int* row_ind_in,
                         int* col_ptr_in);

void c_elsi_customize(int print_detail,
                      int unit_overlap,
                      double hartree_to_ev,
                      double numerical_zero,
                      int no_check_singularity,
                      double singularity_threshold,
                      int force_stop_singularity);

void c_elsi_customize_mu(int broadening_scheme,
                         double broadening_width,
                         double mu_accuracy,
                         int mu_max_steps);

void c_elsi_customize_omm(int n_elpa_steps_omm,
                          double eigenspectrum_shift,
                          double omm_tolerance,
                          int use_pspblas);

void c_elsi_customize_pexsi(double temperature,
                            double gap,
                            double delta_E,
                            int n_poles,
                            int max_iteration,
                            double mu_min,
                            double mu_max,
                            double mu0,
                            double mu_inertia_tolerance,
                            double mu_inertia_expansion,
                            double mu_safeguard,
                            double n_electron_accuracy,
                            int matrix_type,
                            int is_symbolic_factorize,
                            int ordering,
                            int np_symbolic_factorize,
                            int verbosity);

void c_elsi_customize_elpa(int elpa_solver);

void c_elsi_ev_real(double *H_in,
                    double *S_in,
                    double *e_val_out,
                    double *e_vec_out);

void c_elsi_ev_complex(double complex *H_in,
                       double complex *S_in,
                       double *e_val_out,
                       double complex *e_vec_out);

double c_elsi_dm_real(double *H_in,
                      double *S_in,
                      double *D_out);

double c_elsi_dm_real_sparse(double *H_in,
                             double *S_in,
                             double *D_out);
