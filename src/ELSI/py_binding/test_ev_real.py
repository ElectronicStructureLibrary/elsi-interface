#!/opt/local/bin/python

from mpi4py import MPI
import numpy as np
import elsi_python as elsi
import tomato_python as tomato
import ctypes

# Parameterz
e_elpa  = -126.817462901838
e_tol = 1.0e-10
template_basedir = np.array('/Users/wph/opt/ELSI/tomato-seed', dtype="c")
system_label = np.array('silicon', dtype="c")

#m_storage = np.array('pddbc', dtype="c")
m_storage = np.array('sdden', dtype="c")
solver = 1
n_basis = 22
supercell = np.array([3,3,3], dtype=np.int32, order='F')
orb_r_cut = 0.5
k_point = np.array([0.0,0.0,0.0],  dtype=np.float64, order='F')
defect_perturbation = 0.0

# Internal Variables
frac_occ = 0.0
matrix_size = 0
sparsity = 0.0
n_states = 0
n_rows_H = 0
n_cols_H = 0
n_rows_S = 0
n_cols_S = 0

# Initialize matrices
n_states, frac_occ, n_basis, matrix_size, sparsity, orb_r_cut, \
                       n_rows_H, n_cols_H, n_rows_S, n_cols_S = tomato.tomato_tb_get_dims_py_wrapper(template_basedir, system_label,
                                                                False, frac_occ, n_basis,
                                                                False, matrix_size, supercell,
                                                                False, sparsity, orb_r_cut,
                                                                True, k_point,
                                                                True, defect_perturbation,
                                                                m_storage, True)
H = np.zeros( (n_rows_H,n_cols_H), dtype=np.float64, order='F')
S = np.zeros( (n_rows_S,n_cols_S), dtype=np.float64, order='F')

# Generate test matrix
n_states, frac_occ, n_basis, matrix_size, sparsity, orb_r_cut = tomato.tomato_tb_real_py_wrapper(template_basedir, system_label,
                                                                False, frac_occ, n_basis,
                                                                False, matrix_size, supercell,
                                                                False, sparsity, orb_r_cut,
                                                                True, k_point,
                                                                True, defect_perturbation,
                                                                H, S, m_storage, True)

e_val = np.zeros( (matrix_size), dtype=np.float64, order='F')
e_vec = np.zeros( (matrix_size,matrix_size), dtype=np.float64, order='F')

# Initialize ELSI
n_electrons = 2.0*n_states
elsi.elsi_init_py_wrapper(solver,0,0,matrix_size,n_electrons,n_states)
elsi.elsi_customize_py_wrapper(True,False,1.0e-13,False,1.0e-5,False)

# Run ELSI
elsi.elsi_ev_real_py_wrapper(H,S,e_val,e_vec)

# Confirm validity
e_ref = e_elpa
e_test = 2.0*sum(e_val[0:n_states])

print '  Finished test program'
if(abs(e_test-e_ref) < e_tol):
    print '  Passed.'
else:
    print '  Failed!!'
print

elsi.elsi_finalize_py_wrapper()
