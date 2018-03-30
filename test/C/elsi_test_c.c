/* Copyright (c) 2015-2018, the ELSI team.
   All rights reserved.

   This file is part of ELSI and is distributed under the BSD 3-clause license,
   which may be found in the LICENSE file in the ELSI root directory. */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void test_dm_real_c(MPI_Comm comm, int solver, char *h_file, char *s_file);
void test_dm_complex_c(MPI_Comm comm, int solver, char *h_file, char *s_file);
void test_ev_real_c(MPI_Comm comm, int solver, char *h_file, char *s_file);
void test_ev_complex_c(MPI_Comm comm, int solver, char *h_file, char *s_file);

void main(int argc, char** argv) {

   int myid;
   int mpierr;
   int solver;

   MPI_Comm comm;

   // Set up MPI
   MPI_Init(&argc,&argv);
   comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm,&myid);

   // Parameters
   solver = atoi(argv[4]);

   switch(argv[1][0]) {
       case 'e' :
           switch(argv[3][0]) {
               case 'r' :
                   test_ev_real_c(comm,solver,argv[5],argv[6]);
                   break;
               case 'c' :
                   test_ev_complex_c(comm,solver,argv[5],argv[6]);
                   break;
               default :
                   if (myid == 0) {
                       printf("  Wrong command line argument(s)!!\n");
                   }
           }
           break;
       case 'd' :
           switch(argv[3][0]) {
               case 'r' :
                   test_dm_real_c(comm,solver,argv[5],argv[6]);
                   break;
               case 'c' :
                   test_dm_complex_c(comm,solver,argv[5],argv[6]);
                   break;
               default :
                   if (myid == 0) {
                       printf("  Wrong command line argument(s)!!\n");
                   }
           }
           break;
       default :
           if (myid == 0) {
               printf("  Wrong command line argument(s)!!\n");
           }
   }

   MPI_Finalize();

   return;
}
