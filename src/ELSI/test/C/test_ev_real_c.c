/* Copyright (c) 2015-2018, the ELSI team. All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  * Neither the name of the "ELectronic Structure Infrastructure" project nor
    the names of its contributors may be used to endorse or promote products
    derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
 INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

// This program tests elsi_ev_real.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <elsi.h>

void main(int argc, char** argv) {

   int n_proc,n_prow,n_pcol,myid;
   int mpierr;
   int blacs_ctxt;
   int blk,l_row,l_col,l_size;
   int n_basis,n_states;
   int solver,format,parallel;
   int int_one,int_zero;
   int success;
   int tmp;
   int i;

   MPI_Comm mpi_comm_global;

   double n_electrons;
   double *h,*s,*eval,*evec;
   double e_elpa,e_test,e_tol;

   e_elpa = -2564.61963724048;
   e_tol  = 0.00000001;

   elsi_handle    e_h;
   elsi_rw_handle rw_h;

   // Set up MPI
   MPI_Init(&argc,&argv);
   mpi_comm_global = MPI_COMM_WORLD;
   MPI_Comm_size(mpi_comm_global,&n_proc);
   MPI_Comm_rank(mpi_comm_global,&myid);

   // Parameters
   blk      = 16;
   solver   = 1; // ELPA
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
   blacs_ctxt = mpi_comm_global;
   blacs_gridinit_(&blacs_ctxt,"R",&n_prow,&n_pcol);

   // Read H and S matrices
   c_elsi_init_rw(&rw_h,0,1,0,0.0);
   c_elsi_set_rw_mpi(rw_h,mpi_comm_global);
   c_elsi_set_rw_blacs(rw_h,blacs_ctxt,blk);
   c_elsi_set_rw_output(rw_h,2);

   c_elsi_read_mat_dim(rw_h,argv[2],&n_electrons,&n_basis,&l_row,&l_col);

   l_size = l_row * l_col;
   h      = malloc(l_size * sizeof(double));
   s      = malloc(l_size * sizeof(double));
   evec   = malloc(l_size * sizeof(double));
   eval   = malloc(n_basis * sizeof(double));

   c_elsi_read_mat_real(rw_h,argv[2],h);
   c_elsi_read_mat_real(rw_h,argv[3],s);

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
       c_elsi_set_mpi(e_h,mpi_comm_global);
       c_elsi_set_blacs(e_h,blacs_ctxt,blk);
   }

   // Customize ELSI
   c_elsi_set_output(e_h,2);

   // Call ELSI eigensolver
   c_elsi_ev_real(e_h,h,s,eval,evec);

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
       else {
           printf("  Failed!!\n");
       }
   }

   free(h);
   free(s);
   free(evec);
   free(eval);

   blacs_gridexit_(&blacs_ctxt);
   blacs_exit_(1);
   MPI_Finalize();

   return;
}
