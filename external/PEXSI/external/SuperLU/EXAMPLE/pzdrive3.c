
/*! @file 
 * \brief Driver program for PZGSSVX example
 *
 * <pre>
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 * </pre>
 */

#include <math.h>
#include "superlu_zdefs.h"

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * The driver program PZDRIVE3.
 *
 * This example illustrates how to use PZGSSVX to solve
 * systems repeatedly with the same sparsity pattern and similar
 * numerical values of matrix A.
 * In this case, the column permutation vector and symbolic factorization are
 * computed only once. The following data structures will be reused in the
 * subsequent call to PZGSSVX:
 *        ScalePermstruct : DiagScale, R, C, perm_r, perm_c
 *        LUstruct        : etree, Glu_persist, Llu
 *
 * NOTE:
 * The distributed nonzero structures of L and U remain the same,
 * although the numerical values are different. So 'Llu' is set up once
 * in the first call to PZGSSVX, and reused in the subsequent call.
 *
 * On an IBM SP, the program may be run by typing:
 *    poe pzdrive3 -r <proc rows> -c <proc columns> <input_matrix> -procs <p>
 * </pre>
 */

int main(int argc, char *argv[])
{
    superlu_options_t options;
    SuperLUStat_t stat;
    SuperMatrix A;
    NRformat_loc *Astore;
    ScalePermstruct_t ScalePermstruct;
    LUstruct_t LUstruct;
    SOLVEstruct_t SOLVEstruct;
    gridinfo_t grid;
    double   *berr;
    doublecomplex   *b, *b1, *xtrue, *nzval, *nzval1;
    int_t    *colind, *colind1, *rowptr, *rowptr1;
    int_t    i, j, m, n, nnz_loc, m_loc, fst_row;
    int_t    nprow, npcol;
    int      iam, info, ldb, ldx, nrhs;
    char     **cpp, c;
    FILE *fp, *fopen();


    nprow = 1;  /* Default process rows.      */
    npcol = 1;  /* Default process columns.   */
    nrhs = 1;   /* Number of right-hand side. */

    /* ------------------------------------------------------------
       INITIALIZE MPI ENVIRONMENT. 
       ------------------------------------------------------------*/
    MPI_Init( &argc, &argv );

    /* Parse command line argv[]. */
    for (cpp = argv+1; *cpp; ++cpp) {
	if ( **cpp == '-' ) {
	    c = *(*cpp+1);
	    ++cpp;
	    switch (c) {
	      case 'h':
		  printf("Options:\n");
		  printf("\t-r <int>: process rows    (default %d)\n", nprow);
		  printf("\t-c <int>: process columns (default %d)\n", npcol);
		  exit(0);
		  break;
	      case 'r': nprow = atoi(*cpp);
		        break;
	      case 'c': npcol = atoi(*cpp);
		        break;
	    }
	} else { /* Last arg is considered a filename */
	    if ( !(fp = fopen(*cpp, "r")) ) {
                ABORT("File does not exist");
            }
	    break;
	}
    }

    /* ------------------------------------------------------------
       INITIALIZE THE SUPERLU PROCESS GRID. 
       ------------------------------------------------------------*/
    superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);

    /* Bail out if I do not belong in the grid. */
    iam = grid.iam;
    if ( iam >= nprow * npcol )	goto out;
    if ( !iam ) printf("\tProcess grid\t%d X %d\n", grid.nprow, grid.npcol);
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Enter main()");
#endif

    /* ------------------------------------------------------------
       GET THE MATRIX FROM FILE AND SETUP THE RIGHT HAND SIDE. 
       ------------------------------------------------------------*/
    zcreate_matrix(&A, nrhs, &b, &ldb, &xtrue, &ldx, fp, &grid);

    if ( !(b1 = doublecomplexMalloc_dist(ldb * nrhs)) )
        ABORT("Malloc fails for b1[]");
    for (j = 0; j < nrhs; ++j)
        for (i = 0; i < ldb; ++i) b1[i+j*ldb] = b[i+j*ldb];
    if ( !(berr = doubleMalloc_dist(nrhs)) )
	ABORT("Malloc fails for berr[].");
    m = A.nrow;
    n = A.ncol;

    /* Save a copy of the matrix A. */
    Astore = (NRformat_loc *) A.Store;
    nnz_loc = Astore->nnz_loc;
    m_loc = Astore->m_loc;
    fst_row = Astore->fst_row;
    nzval = Astore->nzval;
    colind = Astore->colind;
    rowptr = Astore->rowptr;
    nzval1 = doublecomplexMalloc_dist(nnz_loc);
    colind1 = intMalloc_dist(nnz_loc);
    rowptr1 = intMalloc_dist(m_loc+1);
    for (i = 0; i < nnz_loc; ++i) {
        nzval1[i] = nzval[i];
        colind1[i] = colind[i];
    }
    for (i = 0; i < m_loc+1; ++i) rowptr1[i] = rowptr[i];

    /* ------------------------------------------------------------
       WE SOLVE THE LINEAR SYSTEM FOR THE FIRST TIME.
       ------------------------------------------------------------*/

    /* Set the default input options:
        options.Fact = DOFACT;
        options.Equil = YES;
        options.ColPerm = METIS_AT_PLUS_A;
        options.RowPerm = LargeDiag;
        options.ReplaceTinyPivot = YES;
        options.Trans = NOTRANS;
        options.IterRefine = DOUBLE;
        options.SolveInitialized = NO;
        options.RefineInitialized = NO;
        options.PrintStat = YES;
     */
    set_default_options_dist(&options);

    /* Initialize ScalePermstruct and LUstruct. */
    ScalePermstructInit(m, n, &ScalePermstruct);
    LUstructInit(m, n, &LUstruct);

    /* Initialize the statistics variables. */
    PStatInit(&stat);

    /* Call the linear equation solver: factorize and solve. */
    pzgssvx(&options, &A, &ScalePermstruct, b, ldb, nrhs, &grid,
            &LUstruct, &SOLVEstruct, berr, &stat, &info);

    /* Check the accuracy of the solution. */
    pzinf_norm_error(iam, m_loc, nrhs, b, ldb, xtrue, ldx, &grid);
    
    PStatPrint(&options, &stat, &grid);        /* Print the statistics. */
    PStatFree(&stat);
    Destroy_CompRowLoc_Matrix_dist(&A); /* Deallocate storage of matrix A.  */
    SUPERLU_FREE(b);                 /* Free storage of right-hand side.    */


    /* ------------------------------------------------------------
       NOW WE SOLVE ANOTHER LINEAR SYSTEM.
       THE MATRIX A HAS THE SAME SPARSITY PATTERN AND THE SIMILAR
       NUMERICAL VALUES AS THAT IN A PREVIOUS SYSTEM.
       ------------------------------------------------------------*/
    options.Fact = SamePattern_SameRowPerm;
    PStatInit(&stat); /* Initialize the statistics variables. */

    /* Set up the local A in NR_loc format */
    zCreate_CompRowLoc_Matrix_dist(&A, m, n, nnz_loc, m_loc, fst_row,
				   nzval1, colind1, rowptr1,
				   SLU_NR_loc, SLU_Z, SLU_GE);

    /* Solve the linear system. */
    pzgssvx(&options, &A, &ScalePermstruct, b1, ldb, nrhs, &grid,
            &LUstruct, &SOLVEstruct, berr, &stat, &info);

    /* Check the accuracy of the solution. */
    if ( !iam )
        printf("Solve a system with the same pattern and similar values.\n");
    pzinf_norm_error(iam, m_loc, nrhs, b1, ldb, xtrue, ldx, &grid);

    /* Print the statistics. */
    PStatPrint(&options, &stat, &grid);

    /* ------------------------------------------------------------
       DEALLOCATE STORAGE.
       ------------------------------------------------------------*/
    PStatFree(&stat);
    Destroy_CompRowLoc_Matrix_dist(&A); /* Deallocate storage of matrix A.  */
    Destroy_LU(n, &grid, &LUstruct); /* Deallocate storage associated with    
					the L and U matrices.               */
    ScalePermstructFree(&ScalePermstruct);
    LUstructFree(&LUstruct);         /* Deallocate the structure of L and U.*/
    if ( options.SolveInitialized ) {
        zSolveFinalize(&options, &SOLVEstruct);
    }
    SUPERLU_FREE(b1);	             /* Free storage of right-hand side.    */
    SUPERLU_FREE(xtrue);             /* Free storage of the exact solution. */
    SUPERLU_FREE(berr);


    /* ------------------------------------------------------------
       RELEASE THE SUPERLU PROCESS GRID.
       ------------------------------------------------------------*/
out:
    superlu_gridexit(&grid);

    /* ------------------------------------------------------------
       TERMINATES THE MPI EXECUTION ENVIRONMENT.
       ------------------------------------------------------------*/
    MPI_Finalize();

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit main()");
#endif

}


int cpp_defs()
{
    printf(".. CPP definitions:\n");
#if ( PRNTlevel>=1 )
    printf("\tPRNTlevel = %d\n", PRNTlevel);
#endif
#if ( DEBUGlevel>=1 )
    printf("\tDEBUGlevel = %d\n", DEBUGlevel);
#endif
#if ( PROFlevel>=1 )
    printf("\tPROFlevel = %d\n", PROFlevel);
#endif
#if ( StaticPivot>=1 )
    printf("\tStaticPivot = %d\n", StaticPivot);
#endif
    printf("....\n");
    return 0;
}


