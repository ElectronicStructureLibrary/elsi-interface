/* Copyright (c) 2015-2018, the ELSI team.
   All rights reserved.

   This file is part of ELSI and is distributed under the BSD 3-clause license,
   which may be found in the LICENSE file in the ELSI root directory. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <elsi.h>

void test_ev_complex_c(int mpi_comm,
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
   double *eval;
   double _Complex *h,*s,*evec;
   double e_elpa,e_test,e_tol;

   elsi_handle    e_h;
   elsi_rw_handle rw_h;

   e_elpa = -2564.61963724048;
   e_tol  = 0.00000001;

   MPI_Comm_size(mpi_comm,&n_proc);
   MPI_Comm_rank(mpi_comm,&myid);

   // Parameters
   blk      = 16;
   format   = 0; // BLACS_DENSE
   int_one  = 1;
   int_zero = 0;

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
   evec   = malloc(l_size * sizeof(double _Complex));
   eval   = malloc(n_basis * sizeof(double));

   c_elsi_read_mat_complex(rw_h,h_file,h);
   c_elsi_read_mat_complex(rw_h,s_file,s);

   c_elsi_finalize_rw(rw_h);

   n_states = n_electrons;

   // Initialize ELSI
   if (n_proc == 1) {
       parallel = 0; // Test SINGLE_PROC mode

       c_elsi_init(&e_h,solver,parallel,format,n_basis,n_electrons,n_states);
   }
   else {
       parallel = 1; // Test MULTI_PROC mode

       c_elsi_init(&e_h,solver,parallel,format,n_basis,n_electrons,n_states);
       c_elsi_set_mpi(e_h,mpi_comm);
       c_elsi_set_blacs(e_h,blacs_ctxt,blk);
   }

   // Customize ELSI
   c_elsi_set_output(e_h,2);

   // Call ELSI eigensolver
   c_elsi_ev_complex(e_h,h,s,eval,evec);

   // Finalize ELSI
   c_elsi_finalize(e_h);

   e_test = 0.0;
   for (i=0; i<n_states; i++) {
        e_test += 2.0*eval[i];
   }

   if (myid == 0) {
       if (fabs(e_test-e_elpa) < e_tol) {
           printf("  Passed.\n");
       }
   }

   free(h);
   free(s);
   free(evec);
   free(eval);

   blacs_gridexit_(&blacs_ctxt);

   return;
}
