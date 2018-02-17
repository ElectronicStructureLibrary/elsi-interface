/* Copyright (c) 2015-2018, the ELSI team.
   All rights reserved.

   This file is part of ELSI and is distributed under the BSD 3-clause license,
   which may be found in the LICENSE file in the ELSI root directory. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <elsi.h>

void test_dm_complex_c(int mpi_comm,
                       int solver,
                       char *h_file,
                       char *s_file) {

   int n_proc,n_prow,n_pcol,myid;
   int mpierr;
   int blacs_ctxt;
   int blk,l_row,l_col,l_size;
   int n_basis,n_states;
   int format,parallel;
   int int_one,int_zero;
   int success;
   int tmp;
   int i;

   double n_electrons;
   double _Complex *h,*s,*dm;
   double e_elpa,e_omm,e_pexsi,e_test,e_tol,e_ref;

   elsi_handle    e_h;
   elsi_rw_handle rw_h;

   e_elpa  = -2622.88214509316;
   e_omm   = -2622.88214509316;
   e_pexsi = -2622.88143358352;

   MPI_Comm_size(mpi_comm,&n_proc);
   MPI_Comm_rank(mpi_comm,&myid);

   // Parameters
   blk      = 16;
   format   = 0; // BLACS_DENSE
   parallel = 1; // MULTI_PROC
   int_one  = 1;
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
       e_tol = 0.0001;
   }

   tmp = (int) round(sqrt((double) n_proc));
   for (n_pcol=tmp; n_pcol>1; n_pcol--) {
	if (n_proc%n_pcol == 0) {
	break;
	}
   }
   n_prow = n_proc/n_pcol;

   // Set up BLACS
   blacs_ctxt = mpi_comm;
   blacs_gridinit_(&blacs_ctxt,"R",&n_prow,&n_pcol);

   // Read H and S matrices
   c_elsi_init_rw(&rw_h,0,1,0,0.0);
   c_elsi_set_rw_mpi(rw_h,mpi_comm);
   c_elsi_set_rw_blacs(rw_h,blacs_ctxt,blk);
   c_elsi_set_rw_output(rw_h,2);

   c_elsi_read_mat_dim(rw_h,h_file,&n_electrons,&n_basis,&l_row,&l_col);

   l_size = l_row * l_col;
   h      = malloc(l_size * sizeof(double _Complex));
   s      = malloc(l_size * sizeof(double _Complex));
   dm     = malloc(l_size * sizeof(double _Complex));

   c_elsi_read_mat_complex(rw_h,h_file,h);
   c_elsi_read_mat_complex(rw_h,s_file,s);

   c_elsi_finalize_rw(rw_h);

   n_states = n_electrons;

   // Initialize ELSI
   c_elsi_init(&e_h,solver,parallel,format,n_basis,n_electrons,n_states);
   c_elsi_set_mpi(e_h,mpi_comm);
   c_elsi_set_blacs(e_h,blacs_ctxt,blk);

   // Customize ELSI
   c_elsi_set_output(e_h,2);
   c_elsi_set_sing_check(e_h,0);
   c_elsi_set_mu_broaden_width(e_h,0.000001);
   c_elsi_set_omm_n_elpa(e_h,1);
   c_elsi_set_pexsi_delta_e(e_h,80.0);
   c_elsi_set_pexsi_np_per_pole(e_h,2);

   // Call ELSI density matrix solver
   c_elsi_dm_complex(e_h,h,s,dm,&e_test);

   // Finalize ELSI
   c_elsi_finalize(e_h);

   if (myid == 0) {
       if (fabs(e_test-e_ref) < e_tol) {
           printf("  Passed.\n");
       }
   }

   free(h);
   free(s);
   free(dm);

   blacs_gridexit_(&blacs_ctxt);

   return;
}
