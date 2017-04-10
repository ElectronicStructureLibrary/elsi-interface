/*
 Copyright (c) 2015-2017, the ELSI team. All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  * Neither the name of the "ELectronic Structure Infrastructure" project nor
    the names of its contributors may be used to endorse or promote products
    derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
 INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <complex.h>

#ifdef __cplusplus
extern "C"{
#endif

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
                            int n_procs_per_pole,
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

void c_elsi_ev_complex(double _Complex *H_in,
                       double _Complex *S_in,
                       double *e_val_out,
                       double _Complex *e_vec_out);

double c_elsi_dm_real(double *H_in,
                      double *S_in,
                      double *D_out);

double c_elsi_dm_real_sparse(double *H_in,
                             double *S_in,
                             double *D_out);

#ifdef __cplusplus
}
#endif
