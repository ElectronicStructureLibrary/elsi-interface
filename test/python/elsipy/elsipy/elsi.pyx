from mpi4py import MPI
from mpi4py cimport MPI
cimport numpy as np
import numpy as np
import sys
from elsi cimport *


cdef class elsi:
    cdef elsi_t eh

    cdef int solver 
    cdef int parallel_mode
    cdef int matrix_format
    cdef int n_basis
    cdef double n_electron
    cdef int n_state

    def __init__(self, int solver = 0, int parallel_mode = 0, int matrix_format = 0, 
                       int n_basis = 0, double n_electron = 0.0, int n_state = 0):
        self.solver = solver
        self.parallel_mode = parallel_mode
        self.matrix_format = matrix_format
        self.n_basis = n_basis
        self.n_electron = n_electron
        self.n_state = n_state
        c_elsi_init(&(self.eh), solver, parallel_mode, matrix_format, n_basis, n_electron, n_state)

    def __del__(self):
        c_elsi_finalize(self.eh)

    def elsi_set_mpi(self,comm):
        if not isinstance(comm, MPI.Comm):
            raise Exception("Parallel mode request. But we did not recieve MPI Communicator.")
        c_elsi_set_mpi(self.eh, MPI_Comm_c2f(<MPI_Comm>(<MPI.Comm>comm).ob_mpi))

    def elsi_set_blacs(self,blacs_ctxt,block_size):
        cdef int ictxt = 0
        if blacs_ctxt==-1:
            raise Exception("Parallel mode request. But we did not recieve BLACS ctxt")
        if block_size<=0:
            raise Exception("Parallel mode request. But illegal BLACS block size.")
        ictxt = <int> blacs_ctxt
        c_elsi_set_blacs(self.eh, ictxt, block_size)

    def elsi_ev_real(self,double[:,:] h_matrix,double[:,:] s_matrix):
        cdef Py_ssize_t n_lrow = h_matrix.shape[0]
        cdef Py_ssize_t n_lcol = h_matrix.shape[1]
        cdef np.ndarray[double, ndim=1] eig_vals = np.zeros((self.n_basis), dtype=np.dtype("d"))
        cdef np.ndarray[double, ndim=2] eig_vecs = np.zeros((n_lrow, n_lcol), dtype=np.dtype("d"))
        cdef double *ptr_h_matrix = &h_matrix[0,0]
        cdef double *ptr_s_matrix = &s_matrix[0,0]
        cdef double *ptr_eig_vals = &eig_vals[0]
        cdef double *ptr_eig_vecs = &eig_vecs[0,0]
        c_elsi_ev_real(self.eh, ptr_h_matrix, ptr_s_matrix, ptr_eig_vals, ptr_eig_vecs)
        return eig_vals, eig_vecs

cdef class elsirw:
    cdef elsi_rw_t rwh

    def __init__(self):
        pass

    def elsi_read_real(self, str filename, int parallel_mode=0, int blacs_ctxt=-1, int block_size=0, comm = 0):
        """read
    
        Parameters
        ----------
    
        Returns
        -------
        """
        cdef int rw_task = 0 # rw_task = 0 means read
        cdef int n_basis = 0 # for initialization
        cdef double n_electrons = 0.0 # for initialization
        cdef int n_lrow = 1
        cdef int n_lcol = 1
        cdef int ictxt = 0

        ##### encode string to char *
        ##### https://groups.google.com/g/cython-users/c/E1E96mpS5cg
        filename_byte_string = filename.encode('UTF-8')
        cdef char* filename_c_string = filename_byte_string

        c_elsi_init_rw(&(self.rwh), rw_task, parallel_mode, n_basis, n_electrons)

        if (parallel_mode == 1): # parallel
            if not isinstance(comm, MPI.Comm):
                raise Exception("Parallel mode request. But we did not recieve MPI Communicator.")
            c_elsi_set_rw_mpi(self.rwh, MPI_Comm_c2f(<MPI_Comm>(<MPI.Comm>comm).ob_mpi))
            if blacs_ctxt==-1:
                raise Exception("Parallel mode request. But we did not recieve BLACS ctxt")
            if block_size<=0:
                raise Exception("Parallel mode request. But illegal BLACS block size.")
            ictxt = <int> blacs_ctxt
            c_elsi_set_rw_blacs(self.rwh, ictxt, block_size)
        print(filename)
        print(filename_byte_string)

        c_elsi_read_mat_dim(self.rwh, filename_c_string, &n_electrons, &n_basis, &n_lrow, &n_lcol)
        cdef np.ndarray[double, ndim=2] local_array = np.zeros((n_lrow, n_lcol), dtype=np.dtype("d"))
        cdef double *ptr_local_array = &local_array[0,0]
        c_elsi_read_mat_real(self.rwh, filename_c_string, ptr_local_array)
        c_elsi_finalize_rw(self.rwh)

        return local_array, n_electrons, n_basis, n_lrow, n_lcol

