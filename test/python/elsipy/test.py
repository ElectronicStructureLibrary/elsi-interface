from elsipy import elsi

elsi_rw_h = elsi.elsirw()
H_matrix, n_electrons, n_basis, n_lrow, n_lcol = elsi_rw_h.elsi_read_real("H_real.csc")
S_matrix, n_electrons, n_basis, n_lrow, n_lcol = elsi_rw_h.elsi_read_real("S_real.csc")

elsi_h = elsi.elsi(n_basis=n_basis,n_electron=n_electrons,n_state=n_basis)
eig_vals, eig_vecs = elsi_h.elsi_ev_real(H_matrix,S_matrix)
print(eig_vals)
