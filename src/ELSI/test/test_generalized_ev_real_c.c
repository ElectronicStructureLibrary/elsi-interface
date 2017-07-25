/* Copyright (c) 2015-2017, the ELSI team. All rights reserved.

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

// This program demonstrates how to use the ELSI interface in C to solve a
// generalized eigenvalue problem.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <elsi.h>
#include <tomato.h>

void main(int argc, char** argv) {

   int n_proc,n_prow,n_pcol,myid,my_prow,my_pcol;
   int mpi_comm_global;
   int mpierr;
   int blacs_ctxt;
   int blk,l_row,l_col,l_size;
   int n_basis,n_states;
   int solver,format,parallel;
   int int_one,int_zero;
   int success;
   int tmp;
   int i;
   int *supercell;

   double *k_point;
   double frac_occ,sparsity,r_cut;
   double n_electrons;
   double double_one,double_zero;
   double *h,*s,*eval,*evec;
   double e_elpa,e_test,e_tol;

   e_elpa = -126.817462901838;
   e_tol = 0.0000000001;

   elsi_handle elsi_h;
   matrix h_ms,s_ms;

   // Set up MPI
   MPI_Init(&argc,&argv);
   mpi_comm_global = MPI_COMM_WORLD;
   MPI_Comm_size(mpi_comm_global,&n_proc);
   MPI_Comm_rank(mpi_comm_global,&myid);

   // Parameters
   blk         = 16;
   solver      = 1; // ELPA
   format      = 0; // BLACS_DENSE
   parallel    = 1; // MULTI_PROC
   r_cut       = 0.5; // Tomato only
   int_one     = 1;
   int_zero    = 0;
   double_one  = 1.0;
   double_zero = 0.0;

   supercell = malloc(3 * sizeof(int));
   k_point   = malloc(3 * sizeof(double));
   for (i=0; i<3; i++) {
        supercell[i] = 3;
        k_point[i]   = 0.0;
   }

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
   blacs_gridinfo_(&blacs_ctxt,&n_prow,&n_pcol,&my_prow,&my_pcol);

   // MatrixSwitch
   c_ms_scalapack_setup(mpi_comm_global,n_prow,"r",blk,blacs_ctxt);

   // Prepare matrices by Tomato
   c_tomato_tb(argv[1],"silicon",0,&frac_occ,22,0,&n_basis,supercell,0,
               &sparsity,&r_cut,&n_states,1,k_point,1,0.0,&h_ms,&s_ms,"pddbc",1);

   n_electrons = 2.0*n_states;

   l_row = numroc_(&n_basis,&blk,&my_prow,&int_zero,&n_prow);
   l_col = numroc_(&n_basis,&blk,&my_pcol,&int_zero,&n_pcol);

   l_size = l_row * l_col;
   h      = malloc(l_size * sizeof(double));
   s      = malloc(l_size * sizeof(double));
   evec   = malloc(l_size * sizeof(double));
   eval   = malloc(n_basis * sizeof(double));

   c_get_mat_real(h_ms,h);
   c_get_mat_real(s_ms,s);

   c_m_deallocate(h_ms);
   c_m_deallocate(s_ms);
   free(supercell);
   free(k_point);

   // Initialize ELSI
   c_elsi_init(&elsi_h,solver,parallel,format,n_basis,n_electrons,n_states);
   c_elsi_set_mpi(elsi_h,mpi_comm_global);
   c_elsi_set_blacs(elsi_h,blacs_ctxt,blk);

   // Customize ELSI
   c_elsi_set_output(elsi_h,2);
   c_elsi_set_sing_check(elsi_h,0);

   // Call ELSI eigensolver
   c_elsi_ev_real(elsi_h,h,s,eval,evec);

   // Finalize ELSI
   c_elsi_finalize(elsi_h);

   e_test = 0.0;
   for (i=0; i<n_states; i++) {
        e_test += eval[i];
   }

   e_test *= 2.0;

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

   MPI_Finalize();

   return;
}
