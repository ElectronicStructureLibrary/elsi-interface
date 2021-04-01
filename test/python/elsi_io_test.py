###
### read H&S matrices written by ELSI and solve the generalized eigenvalue problem.
### This code is serial for now.
###

from ctypes import *
import numpy as np
import scipy.linalg

#initialize Hamiltonian matrix
libelsi = CDLL("./libelsi.so")
rh = c_void_p()
r_task = c_int(0)
wh = c_void_p()
w_task = c_int(1)
para_mode = c_int(0)
n_basis = c_int(0)
n_elec = c_double(0.0)
n_lrow = c_int()
n_lcol = c_int()

### read Hamiltonian matrix
H_mat_file = create_string_buffer(b"../matrices/H_real.csc")
libelsi.c_elsi_init_rw(byref(rh),r_task,para_mode,n_basis,n_elec)
libelsi.c_elsi_read_mat_dim(rh,H_mat_file,byref(n_elec),byref(n_basis),byref(n_lrow),byref(n_lcol))
print("Number of basis functions:",n_basis.value)
print("Number of electrons:",n_elec.value)
H_mat = np.zeros([n_basis.value,n_basis.value])
libelsi.c_elsi_read_mat_real(rh,H_mat_file,H_mat.ctypes.data_as(POINTER(c_double)))
libelsi.c_elsi_finalize_rw(rh)
print("Finished reading Hamiltonian matrices")
print()

### read overlap matrix
S_mat_file = create_string_buffer(b"../matrices/S_real.csc")
libelsi.c_elsi_init_rw(byref(rh),r_task,para_mode,n_basis,n_elec)
libelsi.c_elsi_read_mat_dim(rh,S_mat_file,byref(n_elec),byref(n_basis),byref(n_lrow),byref(n_lcol))
print("Number of basis functions:",n_basis.value)
print("Number of electrons:",n_elec.value)
S_mat = np.zeros([n_basis.value,n_basis.value])
libelsi.c_elsi_read_mat_real(rh,S_mat_file,S_mat.ctypes.data_as(POINTER(c_double)))
libelsi.c_elsi_finalize_rw(rh)
print("Finished reading overlap matrices")
print()

### solve the generalized eigenvalue problem by scipy
print("Solve the generalized eigenvalue problem")
e_vals, e_vecs = scipy.linalg.eigh(H_mat, S_mat)
print("eigenvalues:")
for e_val in e_vals:
    print(e_val)

### solver eigenvalue problem

### write density matrix
#new_mat_file = create_string_buffer(b"mat_new.csc")
#libelsi.c_elsi_init_rw(byref(wh),w_task,para_mode,n_basis,n_elec)
#libelsi.c_elsi_write_mat_real(wh,new_mat_file,new_mat.ctypes.data_as(POINTER(c_double)))
#libelsi.c_elsi_finalize_rw(wh)
#print("Finished writing matrices")


