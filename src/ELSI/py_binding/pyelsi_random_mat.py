#!/opt/local/bin/python

from mpi4py import MPI
import elsi_python as elsi
import elsi_scalapack_python as scalapack
import numpy as np
import random
import sys
import math

matrix_size = int(sys.argv[1])

# Initialize MPI
mpi_comm_global_py = MPI.COMM_WORLD
n_proc = mpi_comm_global_py.Get_size()
myid = mpi_comm_global_py.Get_rank()
mpi_comm_global = mpi_comm_global_py.py2f()

# Initialize BLACS
if n_proc > 1:
    npcol = int(math.ceil(math.sqrt(n_proc)))
    for i in range(int(math.ceil(math.sqrt(n_proc))),2,-1):
        npcol = i
        if n_proc % npcol == 0:
            break
    nprow = n_proc/npcol

    blk = 128

    BLACS_CTXT = mpi_comm_global
    BLACS_CTXT = scalapack.blacs_gridinit_py_wrapper(BLACS_CTXT,'r',nprow,npcol)
    scalapack.ms_scalapack_setup_py_wrapper(mpi_comm_global,nprow,'r',blk,False,0,True,BLACS_CTXT)

# Initialize random (global) Hamiltonian matrix
H_global=np.zeros( (matrix_size,matrix_size), dtype=np.float64, order='F')
S_global=np.identity( matrix_size, dtype=np.float64 )
for i in xrange(0, matrix_size):
    for j in xrange(i, matrix_size):
        H_global[j,i] = random.random()
        if i != j:
            H_global[i,j] = H_global[j,i]

# MPI task 0 broadcasts its Hamiltonian to all other ranks

# Diagonalize the global Hamiltonian matrix using numpy
eigenval_lesser, eigenvec_lesser = np.linalg.eig(H_global)
eigenval_lesser = np.sort(eigenval_lesser)

# Use matrix switch to convert dense serial matrix to 2D block-cyclic

print
print "THE ONE TRUE MATRIX SIZE"
print matrix_size
print "THE ONE TRUE HAMILTONIAN MATRIX"
print H_global
print "THE ONE TRUE OVERLAP MATRIX"
print S_global
print
if n_proc > 1:
    elsi.elsi_init_py_wrapper(1,1,0,matrix_size,matrix_size,matrix_size)
    elsi.elsi_set_mpi_py_wrapper(mpi_comm_global)
    elsi.elsi_set_blacs_py_wrapper(BLACS_CTXT,blk)
else:
    elsi.elsi_init_py_wrapper(1,0,0,matrix_size,matrix_size,matrix_size)

eigenvec=np.zeros( (matrix_size,matrix_size), dtype=np.float64, order='F')
eigenval=np.zeros( (matrix_size), dtype=np.float64, order='F')
elsi.elsi_ev_real_py_wrapper(H_global,S_global,eigenval,eigenvec)

eigenval_diffs = eigenval - eigenval_lesser

print "THE ONE TRUE SET OF EIGENVALUES"
print eigenval
print
print "SOME OTHER SET OF EIGENVALUES"
print eigenval_lesser
print
print "THE DIFFERENCES"
print abs(eigenval_diffs)
print "THE MAXIMUM DIFF"
print max(abs(eigenval_diffs))
print
print "UNLIMITED POWER!"
print
