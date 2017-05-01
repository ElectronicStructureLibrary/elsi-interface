#!/opt/local/bin/python

from mpi4py import MPI
import elsi_python
import numpy as np
import random
import sys

#comm = MPI.COMM_WORLD
#commf = comm.py2f()
#print  MPI.COMM_WORLD

matrix_size = int(sys.argv[1])

H=np.zeros( (matrix_size,matrix_size), dtype=np.float64, order='F')
S=np.identity( matrix_size, dtype=np.float64 )
#S=np.zeros( (matrix_size,matrix_size), dtype=np.float64, order='F')
for i in xrange(0, matrix_size):
    for j in xrange(i, matrix_size):
        H[j,i] = random.random()
#        S[j,i] = random.random()
        if i != j:
            H[i,j] = H[j,i]
#            S[i,j] = S[j,i]

eigenval_lesser, eigenvec_lesser = np.linalg.eig(H)
eigenval_lesser = np.sort(eigenval_lesser)

print
print "THE ONE TRUE MATRIX SIZE"
print matrix_size
print "THE ONE TRUE HAMILTONIAN MATRIX"
print H
print "THE ONE TRUE OVERLAP MATRIX"
print S
print
eigenvec=np.zeros( (matrix_size,matrix_size), dtype=np.float64, order='F')
eigenval=np.zeros( (matrix_size), dtype=np.float64, order='F')
elsi_python.elsi_init_py_wrapper(1,0,0,matrix_size,matrix_size,matrix_size)

elsi_python.elsi_ev_real_py_wrapper(H,S,eigenval,eigenvec)
eigenval_diffs = eigenval - eigenval_lesser

print "THE ONE TRUE SET OF EIGENVALUES"
print eigenval
#print "THE ONE TRUE SET OF EIGENVECTORS"
#print eigenvec
print
print "SOME OTHER SET OF EIGENVALUES"
print eigenval_lesser
# Ordering will not be valid
#print "THE LESSER SET OF EIGENVECTORS"
#print eigenvec_lesser
print "THE DIFFERENCES"
print abs(eigenval_diffs)
print "THE MAXIMUM DIFF"
print max(abs(eigenval_diffs))
print
print "UNLIMITED POWER!"
print
