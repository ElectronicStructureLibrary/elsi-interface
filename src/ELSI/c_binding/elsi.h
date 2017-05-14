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

#ifndef ELSI_H_INCLUDED
#define ELSI_H_INCLUDED

#include <complex.h>

#define DECLARE_HANDLE(name) struct name##__ { int unused; }; \
                             typedef struct name##__ *name

DECLARE_HANDLE(elsi_handle);

#ifdef __cplusplus
extern "C"{
#endif

void c_elsi_init(elsi_handle *handle_c,
                 int solver,
                 int parallel_mode,
                 int matrix_storage_format,
                 int matrix_size,
                 double n_electrons,
                 int n_states);

void c_elsi_set_mpi(elsi_handle handle_c,
                    int mpi_comm);

void c_elsi_set_blacs(elsi_handle handle_c,
                      int blacs_ctxt,
                      int block_size);

void c_elsi_finalize(elsi_handle handle_c);

void c_elsi_set_csc(elsi_handle handle_c,
                    int nnz,
                    int nnz_l,
                    int n_l_cols,
                    int* row_ind,
                    int* col_ptr);

void c_elsi_customize(elsi_handle handle_c,
                      int print_detail,
                      int overlap_is_unit,
                      double zero_threshold,
                      int no_singularity_check,
                      double singularity_tolerance,
                      int stop_singularity);

void c_elsi_customize_mu(elsi_handle handle_c,
                         int broadening_scheme,
                         double broadening_width,
                         double occ_accuracy,
                         int mu_max_steps,
                         double spin_degeneracy);

void c_elsi_customize_omm(elsi_handle handle_c,
                          int n_elpa_steps,
                          int omm_flavor,
                          double eigen_shift,
                          double omm_tolerance,
                          int use_pspblas,
                          int omm_output);

void c_elsi_customize_pexsi(elsi_handle handle_c,
                            double temperature,
                            double gap,
                            double delta_e,
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

void c_elsi_customize_elpa(elsi_handle handle_c,
                           int elpa_solver);

void c_elsi_ev_real(elsi_handle handle_c,
                    double *H,
                    double *S,
                    double *e_val,
                    double *e_vec);

void c_elsi_ev_complex(elsi_handle handle_c,
                       double _Complex *H,
                       double _Complex *S,
                       double *e_val,
                       double _Complex *e_vec);

double c_elsi_dm_real(elsi_handle handle_c,
                      double *H,
                      double *S,
                      double *D);

double c_elsi_dm_real_sparse(elsi_handle handle_c,
                             double *H,
                             double *S,
                             double *D);

void c_elsi_collect_pexsi(elsi_handle handle_c,
                          double *mu,
                          double *edm,
                          double *fdm);

#ifdef __cplusplus
}
#endif

#endif
