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
// standard eigenvalue problem.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <mkl_blacs.h>
#include <mkl_scalapack.h>
#include <elsi.h>

void main(int argc, char** argv) {

   int n_proc,n_prow,n_pcol,myid,my_prow,my_pcol;
   int mpi_comm_elsi,mpi_comm_row,mpi_comm_col;
   int mpierr;
   int blacs_ctxt;
   int info,*sc_desc;
   int blk,l_row,l_col,l_size,ldm;
   int n_basis,n_states;
   int solver,format,parallel;
   int int_one,int_zero;
   int success;
   int tmp;
   int i;

   double n_electrons;
   double double_one,double_zero;
   double *h,*h_tmp,*s,*eval,*evec;

   int print,no_s;

   elsi_handle elsi_h;

   // Set up MPI
   MPI_Init(&argc,&argv);
   mpi_comm_elsi = MPI_COMM_WORLD;
   MPI_Comm_size(mpi_comm_elsi,&n_proc);
   MPI_Comm_rank(mpi_comm_elsi,&myid);

   // Parameters
   n_basis     = 4000;
   n_states    = 1000;
   n_electrons = 0; // not used at all
   blk         = 16;
   solver      = 1; // ELPA
   format      = 0; // BLACS_DENSE
   parallel    = 1; // MULTI_PROC
   int_one     = 1;
   int_zero    = 0;
   double_one  = 1.0;
   double_zero = 0.0;
   print       = 2; // Print details
   no_s        = 1; // No overlap, i.e. standard eigenproblem

   tmp = (int) round(sqrt((double) n_proc));
   for (n_pcol=tmp; n_pcol>1; n_pcol--) {
	if (n_proc%n_pcol == 0) {
	break;
	}
   }
   n_prow = n_proc/n_pcol;

   // Set up BLACS
   sc_desc = malloc(9*sizeof(int));
   blacs_ctxt = mpi_comm_elsi;
   blacs_gridinit_(&blacs_ctxt,"R",&n_prow,&n_pcol);
   blacs_gridinfo_(&blacs_ctxt,&n_prow,&n_pcol,&my_prow,&my_pcol);
   l_row = numroc_(&n_basis,&blk,&my_prow,&int_zero,&n_prow);
   l_col = numroc_(&n_basis,&blk,&my_pcol,&int_zero,&n_pcol);

   if (l_row > 1) {
	ldm = l_row;
   }
   else {
	ldm = 1;
   }

   descinit_(sc_desc,&n_basis,&n_basis,&blk,&blk,&int_zero,&int_zero,&blacs_ctxt,&ldm,&info);

   // Prepare a symmetric matrix by pdtran
   l_size = l_row * l_col;
   h      = malloc(l_size * sizeof(double));
   h_tmp  = malloc(l_size * sizeof(double));
   s      = NULL;
   evec   = malloc(l_size * sizeof(double));
   eval   = malloc(n_basis * sizeof(double));

   for (i=0; i<l_size; i++) {
	h[i]     = rand();
	h[i]    /= RAND_MAX;
	h_tmp[i] = h[i];
   }

   pdtran_(&n_basis,&n_basis,&double_one,h_tmp,&int_one,&int_one,sc_desc,&double_one,h,&int_one,&int_one,sc_desc);

   free(sc_desc);
   free(h_tmp);

   mpierr = MPI_Barrier(mpi_comm_elsi);

   // Initialize ELSI
   c_elsi_init(&elsi_h,solver,parallel,format,n_basis,n_electrons,n_states);
   c_elsi_set_mpi(elsi_h,mpi_comm_elsi);
   c_elsi_set_blacs(elsi_h,blacs_ctxt,blk);

   // Customize ELSI
   c_elsi_set_output(elsi_h,print);
   c_elsi_set_unit_ovlp(elsi_h,no_s);

   // Call ELSI eigensolver
   c_elsi_ev_real(elsi_h,h,s,eval,evec);

   // Finalize ELSI
   c_elsi_finalize(elsi_h);

   free(h);
   free(evec);
   free(eval);

   MPI_Finalize();

   return;
}
