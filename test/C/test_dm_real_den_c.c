/* Copyright (c) 2015-2018, the ELSI team.
   All rights reserved.

   This file is part of ELSI and is distributed under the BSD 3-clause license,
   which may be found in the LICENSE file in the ELSI root directory. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <elsi.h>

void test_dm_real_den_c(MPI_Comm comm,
                        int solver,
                        char *h_file,
                        char *s_file) {

   int n_proc;
   int n_prow;
   int n_pcol;
   int myid;
   int mpierr;
   int blacs_ctxt;
   int blk;
   int l_row;
   int l_col;
   int l_size;
   int n_basis;
   int n_states;
   int format;
   int parallel;
   int int_one;
   int int_zero;
   int success;
   int tmp;
   int i;

   double n_electrons;
   double *h;
   double *s;
   double *dm;
   double e_elpa;
   double e_omm;
   double e_pexsi;
   double e_ntpoly;
   double e_test;
   double e_tol;
   double e_ref;

   MPI_Fint comm_f;
   elsi_handle eh;
   elsi_rw_handle rwh;

   e_elpa = -2622.88214509316;
   e_omm = -2622.88214509316;
   e_pexsi = -2622.88143358352;
   e_ntpoly = -2622.88214509311;

   MPI_Comm_size(comm,&n_proc);
   MPI_Comm_rank(comm,&myid);

   // Fortran communicator
   comm_f = MPI_Comm_c2f(comm);

   // Parameters
   blk = 16;
   format = 0; // BLACS_DENSE
   parallel = 1; // MULTI_PROC
   int_one = 1;
   int_zero = 0;

   if (solver == 1) {
       e_ref = e_elpa;
       e_tol = 0.00000001;
   }
   if (solver == 2) {
       e_ref = e_omm;
       e_tol = 0.00000001;
   }
   if (solver == 3) {
       e_ref = e_pexsi;
       e_tol = 0.001;
   }
   if (solver == 6) {
       e_ref = e_ntpoly;
       e_tol = 0.0000001;
   }

   tmp = (int) round(sqrt((double) n_proc));
   for (n_pcol=tmp; n_pcol>1; n_pcol--) {
	if (n_proc%n_pcol == 0) {
	break;
	}
   }
   n_prow = n_proc/n_pcol;

   // Set up BLACS
   blacs_ctxt = Csys2blacs_handle(comm);
   Cblacs_gridinit(&blacs_ctxt,"R",n_prow,n_pcol);

   // Read H and S matrices
   c_elsi_init_rw(&rwh,0,1,0,0.0);
   c_elsi_set_rw_mpi(rwh,comm_f);
   c_elsi_set_rw_blacs(rwh,blacs_ctxt,blk);

   c_elsi_read_mat_dim(rwh,h_file,&n_electrons,&n_basis,&l_row,&l_col);

   l_size = l_row * l_col;
   h = malloc(l_size * sizeof(double));
   s = malloc(l_size * sizeof(double));
   dm = malloc(l_size * sizeof(double));

   c_elsi_read_mat_real(rwh,h_file,h);
   c_elsi_read_mat_real(rwh,s_file,s);

   c_elsi_finalize_rw(rwh);

   n_states = n_electrons;

   // Initialize ELSI
   c_elsi_init(&eh,solver,parallel,format,n_basis,n_electrons,n_states);
   c_elsi_set_mpi(eh,comm_f);
   c_elsi_set_blacs(eh,blacs_ctxt,blk);

   // Customize ELSI
   c_elsi_set_output(eh,2);
   c_elsi_set_sing_check(eh,0);
   c_elsi_set_mu_broaden_width(eh,0.000001);
   c_elsi_set_omm_n_elpa(eh,1);
   c_elsi_set_pexsi_delta_e(eh,80.0);
   c_elsi_set_pexsi_np_per_pole(eh,2);

   // Call ELSI density matrix solver
   c_elsi_dm_real(eh,h,s,dm,&e_test);

   // Finalize ELSI
   c_elsi_finalize(eh);

   if (myid == 0) {
       if (fabs(e_test-e_ref) < e_tol) {
           printf("  Passed.\n");
       }
   }

   free(h);
   free(s);
   free(dm);

   Cblacs_gridexit(blacs_ctxt);

   return;
}
