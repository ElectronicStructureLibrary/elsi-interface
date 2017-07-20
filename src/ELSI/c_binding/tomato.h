#ifndef TOMATO_H_INCLUDED
#define TOMATO_H_INCLUDED

#define DECLARE_HANDLE(name) struct name##__ { int unused; }; \
                             typedef struct name##__ *name

DECLARE_HANDLE(matrix);

#ifdef __cplusplus
extern "C"{
#endif

void c_tomato_tb(char* seed_dir,
                 char* system,
                 int switch1,
                 double* frac_occ,
                 int n_basis_per_atom,
                 int switch2,
                 int* n_basis,
                 int* supercell,
                 int switch3,
                 double* sparsity,
                 double* r_cut,
                 int* n_states,
                 int gamma_only,
                 double* k_point,
                 int defect,
                 double perturbation,
                 matrix *h_c,
                 matrix *s_c,
                 char* ms_storage,
                 int build_matrix);

void c_m_deallocate(matrix handle_c);

void c_ms_scalapack_setup(int mpi_comm,
                          int nprow,
                          char* order,
                          int blk,
                          int ctxt);

void c_get_mat_info(matrix handle_c,
                    int* is_real,
                    int* l_row,
                    int* l_col);

void c_get_mat_real(matrix handle_c,
                    double *mat_real);

#ifdef __cplusplus
}
#endif

#endif
