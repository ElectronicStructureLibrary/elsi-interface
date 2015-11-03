extern"C" {
 void set_up_blacsgrid_from_fortran(int my_mpi_comm_world, int* my_blacs_ctxt,
      int* np_rows, int* np_cols, int* nprow, int* npcol, int* my_prow,
      int* my_pcol);
}

extern"C" {
 int elpa_get_communicators(int my_mpi_comm_world, int my_prow, int my_pcol,
       int* mpi_comm_rows, int* mpi_comm_cols);
}

extern"C" {
 void set_up_blacs_descriptor_from_fortran(int na, int nblk, int my_prow,
      int my_pcol, int np_rows, int np_cols, int* na_rows, int* na_cols,
      int* sc_desc, int my_blacs_ctxt, int* info);
}

extern"C" {
 void prepare_matrix_real_from_fortran(int na, int myid, int na_rows, 
      int na_cols, int* sc_desc, int* iseed, double* a, double* z,
      double* as);
}

extern"C" {
 int elpa_solve_evp_real_1stage(int na, int nev, int na_cols, double* a,
      int lda, double* ev, double* z, int ldz, int nblk, int mpi_comm_rows,
      int mpi_comm_cols);
}

extern"C" {
 int check_correctness_real_from_fortran(int na, int nev, int na_rows, 
       int na_cols, double* as, double* z, double* ev, int* sc_desc, 
       int myid, double* tmp1, double* tmp2);
}
