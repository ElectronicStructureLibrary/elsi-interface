
/*! @file 
 * \brief Performs LU factorization in parallel
 *
 * <pre>
 * -- Distributed SuperLU routine (version 3.2) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * October 22, 2012
 *
 * Modified:
 *     September 1, 1999
 *     Feburary 7, 2001  use MPI_Isend/MPI_Irecv
 *     October 15, 2008  latency-reducing panel factorization
 *     July    12, 2011  static scheduling and arbitrary look-ahead
 *     March   13, 2013  change NTAGS to MPI_TAG_UB value
 *
 * Sketch of the algorithm
 * =======================
 *
 * The following relations hold:
 *     * A_kk = L_kk * U_kk
 *     * L_ik = Aik * U_kk^(-1)
 *     * U_kj = L_kk^(-1) * A_kj
 *
 *              ----------------------------------
 *              |   |                            |
 *              ----|-----------------------------
 *              |   | \ U_kk|                    |
 *              |   |   \   |        U_kj        |
 *              |   |L_kk \ |         ||         |
 *              ----|-------|---------||----------
 *              |   |       |         \/         |
 *              |   |       |                    |
 *              |   |       |                    |
 *              |   |       |                    |
 *              |   | L_ik ==>       A_ij        |
 *              |   |       |                    |
 *              |   |       |                    |
 *              |   |       |                    |
 *              ----------------------------------
 *
 * Handle the first block of columns separately.
 *     * Factor diagonal and subdiagonal blocks and test for exact
 *       singularity. ( pzgstrf2(0), one column at a time )
 *     * Compute block row of U
 *     * Update trailing matrix
 * 
 * Loop over the remaining blocks of columns.
 *   mycol = MYCOL( iam, grid );
 *   myrow = MYROW( iam, grid );
 *   N = nsupers;
 *   For (k = 1; k < N; ++k) {
 *       krow = PROW( k, grid );
 *       kcol = PCOL( k, grid );
 *       Pkk = PNUM( krow, kcol, grid );
 *
 *     * Factor diagonal and subdiagonal blocks and test for exact
 *       singularity.
 *       if ( mycol == kcol ) {
 *           pzgstrf2(k), one column at a time 
 *       }
 *
 *     * Parallel triangular solve
 *       if ( iam == Pkk ) multicast L_k,k to this process row;
 *       if ( myrow == krow && mycol != kcol ) {
 *          Recv L_k,k from process Pkk;
 *          for (j = k+1; j < N; ++j) 
 *              if ( PCOL( j, grid ) == mycol && A_k,j != 0 )
 *                 U_k,j = L_k,k \ A_k,j;
 *       }
 *
 *     * Parallel rank-k update
 *       if ( myrow == krow ) multicast U_k,k+1:N to this process column;
 *       if ( mycol == kcol ) multicast L_k+1:N,k to this process row;
 *       if ( myrow != krow ) {
 *          Pkj = PNUM( krow, mycol, grid );
 *          Recv U_k,k+1:N from process Pkj;
 *       }
 *       if ( mycol != kcol ) {
 *          Pik = PNUM( myrow, kcol, grid );
 *          Recv L_k+1:N,k from process Pik;
 *       }
 *       for (j = k+1; k < N; ++k) {
 *          for (i = k+1; i < N; ++i) 
 *              if ( myrow == PROW( i, grid ) && mycol == PCOL( j, grid )
 *                   && L_i,k != 0 && U_k,j != 0 )
 *                 A_i,j = A_i,j - L_i,k * U_k,j;
 *       }
 *  }
 *
 * </pre>
 */

#include <math.h>
#include "superlu_zdefs.h"

/*
 * Internal prototypes
 */
static void pzgstrf2(superlu_options_t *, int_t, int_t, int_t, double, Glu_persist_t *,
                     gridinfo_t *, LocalLU_t *, MPI_Request *, 
                     SuperLUStat_t *, int *);
#ifdef _CRAY
static void pzgstrs2(int_t, int_t, int_t, Glu_persist_t *, gridinfo_t *,
		     LocalLU_t *, SuperLUStat_t *, _fcd, _fcd, _fcd);
#else
static void pzgstrs2(int_t, int_t, int_t, Glu_persist_t *, gridinfo_t *,
		     LocalLU_t *, SuperLUStat_t *);
#endif

#define ISORT  /* Note: qsort() has bug on Mac */

#ifdef ISORT
extern void isort(int_t N, int_t *ARRAY1, int_t *ARRAY2);
extern void isort1(int_t N, int_t *ARRAY);
#else
int superlu_sort_perm(const void* arg1,const void* arg2);
#endif

typedef struct {
  int id, key;
  void *next;
} etree_node;

/* return the mpi_tag assuming 5 pairs of communications and MPI_TAG_UB >= 5 *
 * for each supernodal column, the five communications are:                  *
 * 0,1: for sending L to "right"                                             *
 * 2,3: for sending off-diagonal blocks of U "down"                          *
 * 4  : for sending the diagonal blcok down (in pxgstrf2)                    */
int tag_ub;
#define SLU_MPI_TAG(id,num) ( (5*(num)+id) % tag_ub )


/************************************************************************/

/*! \brief
 * 
 * <pre>
 * Purpose
 * =======
 *
 *  PDGSTRF performs the LU factorization in parallel.
 *
 * Arguments
 * =========
 * 
 * options (input) superlu_options_t*
 *         The structure defines the input parameters to control
 *         how the LU decomposition will be performed.
 *         The following field should be defined:
 *         o ReplaceTinyPivot (yes_no_t)
 *           Specifies whether to replace the tiny diagonals by
 *           sqrt(epsilon)*norm(A) during LU factorization.
 *
 * m      (input) int
 *        Number of rows in the matrix.
 *
 * n      (input) int
 *        Number of columns in the matrix.
 *
 * anorm  (input) double
 *        The norm of the original matrix A, or the scaled A if
 *        equilibration was done.
 *
 * LUstruct (input/output) LUstruct_t*
 *         The data structures to store the distributed L and U factors.
 *         The following fields should be defined:
 *
 *         o Glu_persist (input) Glu_persist_t*
 *           Global data structure (xsup, supno) replicated on all processes,
 *           describing the supernode partition in the factored matrices
 *           L and U:
 *	       xsup[s] is the leading column of the s-th supernode,
 *             supno[i] is the supernode number to which column i belongs.
 *
 *         o Llu (input/output) LocalLU_t*
 *           The distributed data structures to store L and U factors.
 *           See superlu_zdefs.h for the definition of 'LocalLU_t'.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh. It contains the MPI communicator, the number
 *        of process rows (NPROW), the number of process columns (NPCOL),
 *        and my process rank. It is an input argument to all the
 *        parallel routines.
 *        Grid can be initialized by subroutine SUPERLU_GRIDINIT.
 *        See superlu_ddefs.h for the definition of 'gridinfo_t'.
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics on runtime and floating-point operation count.
 *        See util.h for the definition of 'SuperLUStat_t'.
 *
 * info   (output) int*
 *        = 0: successful exit
 *        < 0: if info = -i, the i-th argument had an illegal value
 *        > 0: if info = i, U(i,i) is exactly zero. The factorization has
 *             been completed, but the factor U is exactly singular,
 *             and division by zero will occur if it is used to solve a
 *             system of equations.
 * </pre>
 */
int_t pzgstrf
(
 superlu_options_t *options, int m, int n, double anorm,
 LUstruct_t *LUstruct, gridinfo_t *grid, SuperLUStat_t *stat, int *info
 )
{
#ifdef _CRAY
    _fcd ftcs = _cptofcd("N", strlen("N"));
    _fcd ftcs1 = _cptofcd("L", strlen("L"));
    _fcd ftcs2 = _cptofcd("N", strlen("N"));
    _fcd ftcs3 = _cptofcd("U", strlen("U"));
#endif
    doublecomplex zero = {0.0, 0.0};
    doublecomplex alpha = {1.0, 0.0}, beta = {0.0, 0.0};
    int_t *xsup;
    int_t *lsub, *lsub1, *usub, *Usub_buf;
    int_t  **Lsub_buf_2, **Usub_buf_2; 
    doublecomplex **Lval_buf_2, **Uval_buf_2;        /* pointers to starts of bufs */ 
    doublecomplex *lusup, *lusup1, *uval, *Uval_buf; /* pointer to current buf     */
    int_t fnz, i, ib, ijb, ilst, it, iukp, jb, jj, klst, knsupc,
          lb, lib, ldv, ljb, lptr, lptr0, lptrj, luptr, luptr0, luptrj,
          nlb, nub, nsupc, rel, rukp, il, iu;
    int_t Pc, Pr;
    int   iam, kcol, krow, yourcol, mycol, myrow, pi, pj;
    int   j, k, kk, lk, nsupers;
    int   nsupr, nbrow, segsize;
    int_t  msg0, msg2;
    int_t  **Ufstnz_br_ptr, **Lrowind_bc_ptr;
    doublecomplex **Unzval_br_ptr, **Lnzval_bc_ptr;
    int_t  *index;
    doublecomplex *nzval;
    int_t  *iuip, *ruip;/* Pointers to U index/nzval; size ceil(NSUPERS/Pr). */
    doublecomplex *ucol;
    int_t  *indirect;
    doublecomplex *tempv, *tempv2d;
    int iinfo;
    int_t *ToRecv, *ToSendD, **ToSendR;
    Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
    LocalLU_t *Llu = LUstruct->Llu;
    superlu_scope_t *scp;
    float s_eps;
    double thresh;
    doublecomplex *tempU2d, *tempu;
    int    full, ldt, ldu, lead_zero, ncols, ncb, nrb, p, pr, pc, nblocks;
    int_t *etree_supno_l, *etree_supno, *blocks, *blockr, *Ublock, *Urows, *Lblock, *Lrows,
          *perm_u, *sf_block, *sf_block_l, *nnodes_l, *nnodes_u, *edag_supno_l, *recvbuf,
          **edag_supno;
#ifdef ISORT
    int_t *iperm_u;
#endif
    int *msgcnt;  /* Count the size of the message xfer'd in each buffer:
                   *     0 : transferred in Lsub_buf[]
           	   *     1 : transferred in Lval_buf[]
                   *     2 : transferred in Usub_buf[]
                   *     3 : transferred in Uval_buf[]
                   */
    int **msgcnts, **msgcntsU;
    int k0; /* counter of the next supernode to be factored */
    int kk0, kk1, kk2, jj0, iukp0, rukp0, flag0, flag1;
    int *factored, *factoredU, nnodes, *sendcnts, *sdispls, *recvcnts, *rdispls, *srows, *rrows;
    etree_node *head, *tail, *ptr;
    int *num_child;
    int num_look_aheads, look_id, *look_ahead;
    int_t *perm_c_supno, *iperm_c_supno;
    MPI_Request *recv_req, **recv_reqs, **send_reqs, **send_reqs_u, **recv_reqs_u;
    MPI_Request *send_req, *U_diag_blk_send_req = NULL;
    MPI_Status status;
    void *attr_val;
    int flag;    

#if ( DEBUGlevel>=2 ) 
    int_t num_copy=0, num_update=0;
#endif
#if ( PRNTlevel==3 )
    int  zero_msg = 0, total_msg = 0;
#endif
#if ( PROFlevel>=1 )
    double t1, t2;
    float msg_vol = 0, msg_cnt = 0;
    int_t iword = sizeof(int_t), zword = sizeof(doublecomplex);
#endif

    /* Test the input parameters. */
    *info = 0;
    if ( m < 0 ) *info = -2;
    else if ( n < 0 ) *info = -3;
    if ( *info ) {
	pxerbla("pzgstrf", grid, -*info);
	return (-1);
    }

    /* Quick return if possible. */
    if ( m == 0 || n == 0 ) return 0;

    /*
     * Initialization.
     */
    iam = grid->iam;
    Pc = grid->npcol;
    Pr = grid->nprow;
    myrow = MYROW( iam, grid );
    mycol = MYCOL( iam, grid );
    nsupers = Glu_persist->supno[n-1] + 1;
    xsup = Glu_persist->xsup;
    s_eps = slamch_("Epsilon");
    thresh = s_eps * anorm;

    MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &attr_val, &flag);
    if ( !flag ) {
       fprintf( stderr, "Could not get TAG_UB\n" );
       return (-1);
    }
    tag_ub = *(int*)attr_val;
#if ( PRNTlevel>=1 )
    if ( !iam ) printf("MPI tag upper bound = %d\n", tag_ub);
#endif

#if ( DEBUGlevel>=1 )
    if( s_eps == 0.0 ) printf( " ***** warning s_eps = %e *****\n",s_eps );
    CHECK_MALLOC(iam, "Enter pzgstrf()");
#endif
    stat->ops[FACT] = 0.0;

    /* make sure the range of look-ahead window [0, MAX_LOOKAHEADS-1] */
    num_look_aheads = SUPERLU_MAX( 0, SUPERLU_MIN( options->num_lookaheads, MAX_LOOKAHEADS-1 ) );

    if ( Pr*Pc > 1 ) {
 	if ( !(U_diag_blk_send_req =
 	       (MPI_Request *) SUPERLU_MALLOC(Pr*sizeof(MPI_Request))))
 	    ABORT("Malloc fails for U_diag_blk_send_req[].");
        U_diag_blk_send_req[myrow] = 0; /* flag no outstanding Isend */

	/* allocating buffers for look-ahead */
        i = Llu->bufmax[0];
        if ( i != 0 ) {
            if ( !(Llu->Lsub_buf_2[0] = intMalloc_dist((num_look_aheads+1) * ((size_t)i))) )
                ABORT("Malloc fails for Lsub_buf.");
            for( jj=0; jj<num_look_aheads; jj++ )
                Llu->Lsub_buf_2[jj+1] = Llu->Lsub_buf_2[jj] + i;
        }
        i = Llu->bufmax[1];
        if ( i != 0 ) {
            if ( !(Llu->Lval_buf_2[0] = doublecomplexMalloc_dist((num_look_aheads+1) * ((size_t)i))) )
                ABORT("Malloc fails for Lval_buf[].");
            for( jj=0; jj<num_look_aheads; jj++ )
                Llu->Lval_buf_2[jj+1] = Llu->Lval_buf_2[jj] + i;
        } 
	i = Llu->bufmax[2];
	if ( i != 0 ) {
	    if ( !(Llu->Usub_buf_2[0] = intMalloc_dist((num_look_aheads+1)*i)) )
		ABORT("Malloc fails for Usub_buf_2[].");
	    for( jj=0; jj<num_look_aheads; jj++ )
		Llu->Usub_buf_2[jj+1] = Llu->Usub_buf_2[jj] + i;
	}
	i = Llu->bufmax[3];
	if ( i != 0 ) {
	    if ( !(Llu->Uval_buf_2[0] = doublecomplexMalloc_dist((num_look_aheads+1)*i)) )
		ABORT("Malloc fails for Uval_buf_2[].");
	    for( jj=0; jj<num_look_aheads; jj++ )
		Llu->Uval_buf_2[jj+1] = Llu->Uval_buf_2[jj] + i;
	}
    }
    /* creating pointers to the look-ahead buffers */
    if ( !(Lsub_buf_2 = SUPERLU_MALLOC( (1+num_look_aheads) * sizeof(int_t*) )) )
	ABORT("Malloc fails for Lsub_buf_2[].");
    if ( !(Lval_buf_2 = SUPERLU_MALLOC( (1+num_look_aheads) * sizeof(doublecomplex*) )) )
	ABORT("Malloc fails for Lval_buf_2[].");  
    if( !(Usub_buf_2 = SUPERLU_MALLOC( (1+num_look_aheads) * sizeof(int_t*) )) )
	ABORT("Malloc fails for Uval_buf_2[].");  
    if( !(Uval_buf_2 = SUPERLU_MALLOC( (1+num_look_aheads) * sizeof(doublecomplex*) )) )
	ABORT("Malloc fails for buf_2[].");  
    for( i=0; i<=num_look_aheads; i++ ) {
      Lval_buf_2[i] = Llu->Lval_buf_2[i];
      Lsub_buf_2[i] = Llu->Lsub_buf_2[i];
      Uval_buf_2[i] = Llu->Uval_buf_2[i];
      Usub_buf_2[i] = Llu->Usub_buf_2[i];
    }
    if ( !(msgcnts  = SUPERLU_MALLOC( (1+num_look_aheads) * sizeof(int_t*) )) ) 
	ABORT("Malloc fails for msgcnts[].");
    if ( !(msgcntsU = SUPERLU_MALLOC( (1+num_look_aheads) * sizeof(int_t*) )) )
	ABORT("Malloc fails for msgcntsU[].");
    for( i=0; i<=num_look_aheads; i++ ) {
	if ( !(msgcnts[i]  = SUPERLU_MALLOC( 4 * sizeof(int_t) )) ) 
	ABORT("Malloc fails for msgcnts[].");
	if ( !(msgcntsU[i] = SUPERLU_MALLOC( 4 * sizeof(int_t) )) )
	ABORT("Malloc fails for msgcntsU[].");
    }

    if ( !(recv_reqs_u = SUPERLU_MALLOC( (1+num_look_aheads) * sizeof(MPI_Request*) )) )
	ABORT("Malloc fails for recv_reqs_u[].");
    if ( !(send_reqs_u = SUPERLU_MALLOC( (1+num_look_aheads) * sizeof(MPI_Request*) )) )
	ABORT("Malloc fails for send_reqs_u[].");
    if ( !(send_reqs   = SUPERLU_MALLOC( (1+num_look_aheads) * sizeof(MPI_Request*) )) )
	ABORT("Malloc fails for send_reqs_u[].");
    if ( !(recv_reqs   = SUPERLU_MALLOC( (1+num_look_aheads) * sizeof(MPI_Request*) )) )
	ABORT("Malloc fails for recv_reqs[].");
    for( i=0; i<=num_look_aheads; i++ ) {
      if ( !(recv_reqs_u[i] =
         (MPI_Request *) SUPERLU_MALLOC(2*sizeof(MPI_Request))))
        ABORT("Malloc fails for recv_req_u[i].");
      if ( !(send_reqs_u[i] =
         (MPI_Request *) SUPERLU_MALLOC(2*Pr*sizeof(MPI_Request))))
        ABORT("Malloc fails for send_req_u[i].");
      if ( !(send_reqs[i] =
         (MPI_Request *) SUPERLU_MALLOC(2*Pc*sizeof(MPI_Request))))
        ABORT("Malloc fails for send_reqs[i].");
      if ( !(recv_reqs[i] =
         (MPI_Request *) SUPERLU_MALLOC(4*sizeof(MPI_Request))))
        ABORT("Malloc fails for recv_req[].");
      send_reqs[i][0] = send_reqs[i][1] = MPI_REQUEST_NULL;
      recv_reqs[i][0] = recv_reqs[i][1] = MPI_REQUEST_NULL;
    }
    if ( !(factored  = SUPERLU_MALLOC( nsupers * sizeof(int_t) )) )
	ABORT("Malloc fails for factored[].");
    if ( !(factoredU  = SUPERLU_MALLOC( nsupers * sizeof(int_t) )) )
	ABORT("Malloc fails for factoredU[].");
    for( i=0; i<nsupers; i++ ) factored[i] = factoredU[i] = -1;

#if ( PRNTlevel>=1 )
    double tt1 = SuperLU_timer_();
#endif
    nblocks = 0;
    ncb = nsupers / Pc;
    nrb = nsupers / Pr;
   /* ================================================== *
    * static scheduling of j-th step of LU-factorization *
    * ================================================== */
    if( options->lookahead_etree == YES && /* use e-tree of symmetrized matrix, and      */
       (options->ParSymbFact == NO ||      /* 1) symmetric fact with serial symbolic, or */
       (options->SymPattern == YES &&      /* 2) symmetric pattern, and                  */
        options->RowPerm == NOROWPERM)) )  /*    no rowperm to destroy the symmetry      */
    {
      /* if symmetric pattern or using e-tree of |A^T|+|A|,
         then we can use a simple tree structure for static schduling */

      if( options->ParSymbFact == NO ) {
          /* Use the etree computed from serial symb. fact., and turn it
             into supernodal tree.  */
          int_t *etree = LUstruct->etree;
#if ( PRNTlevel>=1 )
          if( grid->iam == 0 ) printf( " === using column e-tree ===\n" );
#endif

	  /* look for the first off-diagonal blocks */
          etree_supno = SUPERLU_MALLOC( nsupers * sizeof(int_t) );
          for( i=0; i<nsupers; i++ ) etree_supno[i] = nsupers;
          for( j=0, lb=0; lb<nsupers; lb++ ) {
              for( k=0; k<SuperSize(lb); k++ ) {
                  jb = Glu_persist->supno[etree[j+k]];
                  if( jb != lb ) etree_supno[lb] = SUPERLU_MIN( etree_supno[lb], jb );
              }
              j += SuperSize(lb);
          }
      } else { /* ParSymbFACT==YES and SymPattern==YES  and RowPerm == NOROWPERM */
         /* Compute an "etree" based on struct(L), 
            assuming struct(U) = struct(L').   */
#if ( PRNTlevel>=1 )
        if( grid->iam == 0 ) printf( " === using supernodal e-tree ===\n" );
#endif

        /* find the first block in each supernodal-column of local L-factor */
        etree_supno_l = SUPERLU_MALLOC( nsupers * sizeof(int_t) );
        for( i=0; i<nsupers; i++ ) etree_supno_l[i] = nsupers;
        for( lb=0; lb<ncb; lb++ ) {
          jb = lb * grid->npcol + mycol;
          index = Llu->Lrowind_bc_ptr[lb];
          if ( index ) { /* Not an empty column */
            i = index[0];
            k = BC_HEADER; 
	    krow = PROW( jb, grid );
            if( krow == myrow ) { /* skip the diagonal block */
              k += LB_DESCRIPTOR + index[k+1];
              i--;
            }
            if( i > 0 ) {
              etree_supno_l[jb] = index[k];
              k += LB_DESCRIPTOR + index[k+1];
	      i --;
	    }

	    for( j=0; j<i; j++ ) {
              etree_supno_l[jb] = SUPERLU_MIN( etree_supno_l[jb], index[k] );
              k += LB_DESCRIPTOR + index[k+1];
	    }
          }
        }
        if( mycol < nsupers % grid->npcol ) {
          jb = ncb * grid->npcol + mycol;
          index = Llu->Lrowind_bc_ptr[ncb];
          if ( index ) { /* Not an empty column */
            i = index[0];
            k = BC_HEADER; 
            krow = PROW( jb, grid );
            if( krow == myrow ) { /* skip the diagonal block */
              k += LB_DESCRIPTOR + index[k+1];
              i--;
            }
            if( i > 0 ) {
	      etree_supno_l[jb] = index[k];
	      k += LB_DESCRIPTOR + index[k+1];
              i --;
	    }
            for( j=0; j<i; j++ ) {
              etree_supno_l[jb] = SUPERLU_MIN( etree_supno_l[jb], index[k] );
              k += LB_DESCRIPTOR + index[k+1];
            }
          }
        }

        /* form global e-tree */
        etree_supno = SUPERLU_MALLOC( nsupers * sizeof(int_t) );
        MPI_Allreduce( etree_supno_l, etree_supno, nsupers, mpi_int_t, MPI_MIN, grid->comm );
        SUPERLU_FREE(etree_supno_l);
      }

      /* initialize the num of child for each node */
      num_child = SUPERLU_MALLOC( nsupers * sizeof(int_t) );
      for( i=0; i<nsupers; i++ ) num_child[i] = 0;
      for( i=0; i<nsupers; i++ ) if( etree_supno[i] != nsupers ) num_child[etree_supno[i]] ++;

      /* push initial leaves to the fifo queue */
      nnodes = 0;
      for( i=0; i<nsupers; i++ ) {
        if( num_child[i] == 0 ) {
          ptr = SUPERLU_MALLOC( sizeof(etree_node) );
          ptr->id = i;
          ptr->next = NULL;
          /*printf( " == push leaf %d (%d) ==\n",i,nnodes );*/
          nnodes ++;

          if( nnodes == 1 ) {
            head = ptr;
            tail = ptr;
          } else {
            tail->next = ptr;
            tail = ptr;
          }
        }
      }

      /* process fifo queue, and compute the ordering */
      i = 0;
      perm_c_supno  = SUPERLU_MALLOC( nsupers * sizeof(int_t) );
      while( nnodes > 0 ) {
        ptr = head;  j = ptr->id;
        head = ptr->next;
        perm_c_supno[i] = j;
        SUPERLU_FREE(ptr);
        i++; nnodes --;

        if( etree_supno[j] != nsupers ) {
          num_child[etree_supno[j]] --;
          if( num_child[etree_supno[j]] == 0 ) {
            nnodes ++;

            ptr = SUPERLU_MALLOC( sizeof(etree_node) );
            ptr->id = etree_supno[j];
            ptr->next = NULL;

            /*printf( "=== push %d ===\n",ptr->id );*/
            if( nnodes == 1 ) {
              head = ptr;
              tail = ptr;
            } else {
              tail->next = ptr;
              tail = ptr;
            }
          }
        }
        /*printf( "\n" );*/
      }
      SUPERLU_FREE(num_child);
      SUPERLU_FREE(etree_supno);

    } else { /* Unsymmetric pattern */
        /* Need to process both L- and U-factors, use the symmetrically
           pruned graph of L & U instead of tree (very naive implementation) */
      int nrbp1 = nrb + 1;

      /* allocate some workspace */
      if ( !(sendcnts = SUPERLU_MALLOC( (4+2*nrbp1)*Pr*Pc * sizeof(int))) )
	ABORT("Malloc fails for sendcnts[].");
      sdispls  = &sendcnts[Pr*Pc];
      recvcnts = &sdispls [Pr*Pc];
      rdispls  = &recvcnts[Pr*Pc];
      srows    = &rdispls [Pr*Pc];
      rrows    = &srows   [Pr*Pc*nrbp1];

      myrow = MYROW( iam, grid );
#if ( PRNTlevel>=1 )
      if( grid->iam == 0 ) printf( " === using DAG ===\n" );
#endif

     /* send supno block of local U-factor to a processor *
      * who owns the corresponding block of L-factor      */

      /* srows   : # of block to send to a processor from each supno row */
      /* sendcnts: total # of blocks to send to a processor              */
      for (p=0; p<Pr*Pc*nrbp1; p++) srows[p] = 0; 
      for (p=0; p<Pr*Pc; p++ ) sendcnts[p] = 0;

      /* sending blocks of U-factors corresponding to L-factors */
      /* count the number of blocks to send */
      for (lb = 0; lb < nrb; ++lb) {
        jb = lb * Pr + myrow;
        pc = jb % Pc;
        index = Llu->Ufstnz_br_ptr[lb];

        if ( index ) { /* Not an empty row */
            k = BR_HEADER;
            nblocks += index[0];
            for (j = 0; j < index[0]; ++j) {
                ib = index[k];
                pr = ib % Pr;
                p  = pr * Pc + pc;
                sendcnts[p] ++;
                srows[p*nrbp1+lb] ++;

                k += UB_DESCRIPTOR + SuperSize( index[k] );
            }
        }
      }
      if ( myrow < nsupers % grid->nprow ) {
        jb = nrb * Pr + myrow;
        pc = jb % Pc;
        index = Llu->Ufstnz_br_ptr[nrb];

        if ( index ) { /* Not an empty row */
            k = BR_HEADER;
            nblocks += index[0];
            for (j = 0; j < index[0]; ++j) {
                ib = index[k];
                pr = ib % Pr;
                p  = pr * Pc + pc;
                sendcnts[p] ++;
                srows[p*nrbp1+nrb] ++;
                k += UB_DESCRIPTOR + SuperSize( index[k] );
            }
        }
      }

      /* insert blocks to send */
      sdispls[0] = 0;
      for( p=1; p<Pr*Pc; p++ ) sdispls[p] = sdispls[p-1]+sendcnts[p-1];
      if ( !(blocks = intMalloc_dist( nblocks )) ) ABORT("Malloc fails for blocks[].");
      for (lb = 0; lb < nrb; ++lb) {
        jb = lb * Pr + myrow;
        pc = jb % Pc;
        index = Llu->Ufstnz_br_ptr[lb];

        if ( index ) { /* Not an empty row */
            k = BR_HEADER;
            for (j = 0; j < index[0]; ++j) {
                ib = index[k];
                pr = ib % Pr;
                p  = pr * Pc + pc;
                blocks[sdispls[p]] = ib;
                sdispls[p] ++;

                k += UB_DESCRIPTOR + SuperSize( index[k] );
            }
        }
      }
      if ( myrow < nsupers % grid->nprow ) {
        jb = nrb * Pr + myrow;
        pc = jb % Pc;
        index = Llu->Ufstnz_br_ptr[nrb];

        if ( index ) { /* Not an empty row */
            k = BR_HEADER;
            for (j = 0; j < index[0]; ++j) {
                ib = index[k];
                pr = ib % Pr;
                p  = pr * Pc + pc;
                blocks[sdispls[p]] = ib;
                sdispls[p] ++;

                k += UB_DESCRIPTOR + SuperSize( index[k] );
            }
        }
      }
 
      /* communication */
      MPI_Alltoall( sendcnts,  1, MPI_INT, recvcnts,   1, MPI_INT, grid->comm );
      MPI_Alltoall( srows, nrbp1, MPI_INT, rrows,  nrbp1, MPI_INT, grid->comm );

      nblocks = recvcnts[0];
      rdispls[0] = sdispls[0] = 0;
      for( p=1; p<Pr*Pc; p++ ) {
        rdispls[p] = rdispls[p-1]+recvcnts[p-1];
        sdispls[p] = sdispls[p-1]+sendcnts[p-1];
        nblocks += recvcnts[p];
      }

      if ( !(blockr =  intMalloc_dist( nblocks )) )  ABORT("Malloc fails for blockr[].");
      MPI_Alltoallv( blocks, sendcnts, sdispls, mpi_int_t, blockr, recvcnts, rdispls, mpi_int_t, grid->comm );
      SUPERLU_FREE( blocks );

      /* store the received U-blocks by rows */
      nlb = nsupers / Pc; 
      if ( !(Ublock = intMalloc_dist( nblocks )) )  ABORT("Malloc fails for Ublock[].");
      if ( !(Urows  = intMalloc_dist( 1+nlb )) )  ABORT("Malloc fails for Urows[].");
      k = 0;
      for (jb=0; jb<nlb; jb++ ) {
        j = jb*Pc + mycol;
        pr = j % Pr;
        lb = j / Pr;
        Urows[jb] = 0;

        for( pc=0; pc<Pc; pc++ ) {
          p = pr*Pc+pc; /* the processor owning this block of U-factor */

          for( i=rdispls[p]; i<rdispls[p]+rrows[p*nrbp1+lb]; i++) {
            Ublock[k] = blockr[i];
            k++; Urows[jb] ++;
          }
          rdispls[p] += rrows[p*nrbp1+lb];
        }
        /* sort by the column indices to make things easier for later on */

#ifdef ISORT
        isort1( Urows[jb], &(Ublock[k-Urows[jb]]) );
#else
	qsort( &(Ublock[k-Urows[jb]]), (size_t)(Urows[jb]), sizeof(int_t), &superlu_sort_perm );
#endif
      }
      if ( mycol < nsupers % grid->npcol ) {
        j = nlb*Pc + mycol;
        pr = j % Pr;
        lb = j / Pr;
        Urows[nlb] = 0;

        for( pc=0; pc<Pc; pc++ ) {
          p = pr*Pc+pc;
          for( i=rdispls[p]; i<rdispls[p]+rrows[p*nrbp1+lb]; i++) {
            Ublock[k] = blockr[i];
            k++; Urows[nlb] ++;
          }
          rdispls[p] += rrows[p*nrb+lb];
        }
#ifdef ISORT
        isort1( Urows[nlb], &(Ublock[k-Urows[nlb]]) );
#else
	qsort( &(Ublock[k-Urows[nlb]]), (size_t)(Urows[nlb]), sizeof(int_t), &superlu_sort_perm );
#endif
      }
      SUPERLU_FREE( blockr );

      /* sort the block in L-factor */
      nblocks = 0;
      for( lb=0; lb<ncb; lb++ ) {
        jb = lb * Pc + mycol;
        index = Llu->Lrowind_bc_ptr[lb];
        if ( index ) { /* Not an empty column */
            nblocks += index[0];
        }
      }
      if( mycol < nsupers % grid->npcol ) {
        jb = ncb * Pc + mycol;
        index = Llu->Lrowind_bc_ptr[ncb];
        if ( index ) { /* Not an empty column */
          nblocks += index[0];
        }
      }

      if ( !(Lblock = intMalloc_dist( nblocks )) ) ABORT("Malloc fails for Lblock[].");
      if ( !(Lrows  = intMalloc_dist( 1+ncb )) ) ABORT("Malloc fails for Lrows[].");
      for( lb=0; lb<=ncb; lb++ ) Lrows[lb] = 0;
      nblocks = 0;
      for( lb=0; lb<ncb; lb++ ) {
        Lrows[lb] = 0;

        jb = lb * Pc + mycol;
        index = Llu->Lrowind_bc_ptr[lb];
        if ( index ) { /* Not an empty column */
            i = index[0];
            k = BC_HEADER; 
	    krow = PROW( jb, grid );
            if( krow == myrow ) { /* skip the diagonal block */
              k += LB_DESCRIPTOR + index[k+1];
              i--;
            }

	    for( j=0; j<i; j++ ) {
                Lblock[nblocks] = index[k];
                Lrows[lb] ++;
                nblocks++;

		k += LB_DESCRIPTOR + index[k+1];
	    }
        }
#ifdef ISORT
        isort1( Lrows[lb], &(Lblock[nblocks-Lrows[lb]]) );
#else
	qsort( &(Lblock[nblocks-Lrows[lb]]), (size_t)(Lrows[lb]), sizeof(int_t), &superlu_sort_perm );
#endif
      }
      if( mycol < nsupers % grid->npcol ) {
        Lrows[ncb] = 0;
        jb = ncb * Pc + mycol;
        index = Llu->Lrowind_bc_ptr[ncb];
        if ( index ) { /* Not an empty column */
          i = index[0];
          k = BC_HEADER; 
          krow = PROW( jb, grid );
          if( krow == myrow ) { /* skip the diagonal block */
            k += LB_DESCRIPTOR + index[k+1];
            i--;
          }
          for( j=0; j<i; j++ ) {
            Lblock[nblocks] = index[k];
            Lrows[ncb] ++;
            nblocks++;
            k += LB_DESCRIPTOR + index[k+1];
          }
#ifdef ISORT
          isort1( Lrows[ncb], &(Lblock[nblocks-Lrows[ncb]]) );
#else
	  qsort( &(Lblock[nblocks-Lrows[ncb]]), (size_t)(Lrows[ncb]), sizeof(int_t), &superlu_sort_perm );
#endif
        }
      }

      /* look for the first local symmetric nonzero block match */
      if ( !(sf_block   = intMalloc_dist( nsupers )) )
	ABORT("Malloc fails for sf_block[].");
      if ( !(sf_block_l = intMalloc_dist( nsupers )) )
	ABORT("Malloc fails for sf_block_l[].");
      for( lb=0; lb<nsupers; lb++ ) sf_block_l[lb] = nsupers;
      i = 0; j = 0; 
      for( jb=0; jb<nlb; jb++ ) {
        if( Urows[jb] > 0 ) {
          ib = i + Urows[jb];
          lb = jb * Pc + mycol;
          for( k=0; k<Lrows[jb]; k++ ) {
            while( Ublock[i] < Lblock[j] && i+1 < ib ) i++; 

            if( Ublock[i] == Lblock[j] ) {
              sf_block_l[lb] = Lblock[j];
              j += (Lrows[jb]-k);
              k = Lrows[jb];
            } else {
              j++;
            }
          }
          i = ib;
        } else {
          j += Lrows[jb];
        }
      }
      if( mycol < nsupers % grid->npcol ) {
        if( Urows[nlb] > 0 ) {
          ib = i + Urows[nlb];
          lb = nlb * Pc + mycol;
          for( k=0; k<Lrows[nlb]; k++ ) {
            while( Ublock[i] < Lblock[j] && i+1 < ib ) i++; 

            if( Ublock[i] == Lblock[j] ) {
              sf_block_l[lb] = Lblock[j];
              j += (Lrows[nlb]-k);
              k = Lrows[nlb];
            } else {
              j++;
            }
          }
          i = ib;
        } else {
          j += Lrows[nlb];
        }
      }
      /* compute the first global symmetric matchs */
      MPI_Allreduce( sf_block_l, sf_block, nsupers, mpi_int_t, MPI_MIN, grid->comm );
      SUPERLU_FREE( sf_block_l );

      /* count number of nodes in DAG (i.e., the number of blocks on and above the first match) */
      if ( !(nnodes_l = intMalloc_dist( nsupers )) )
	ABORT("Malloc fails for nnodes_l[].");
      if ( !(nnodes_u = intMalloc_dist( nsupers )) )
	ABORT("Malloc fails for nnodes_u[].");
      for( lb=0; lb<nsupers; lb++ ) nnodes_l[lb] = 0;
      for( lb=0; lb<nsupers; lb++ ) nnodes_u[lb] = 0;

      nblocks = 0; 
      /* from U-factor */
      for (i=0, jb=0; jb<nlb; jb++ ) {
        lb = jb*Pc + mycol;
        ib = i + Urows[jb];
        while( i < ib ) {
          if( Ublock[i] <= sf_block[lb] ) {
            nnodes_u[lb] ++;
            i++; nblocks++;
          } else { /* get out*/
            i = ib;
          }
        }
        i = ib;
      }
      if( mycol < nsupers % grid->npcol ) {
        lb = nlb*Pc + mycol;
        ib = i + Urows[nlb];
        while( i < ib ) {
          if( Ublock[i] <= sf_block[lb] ) {
            nnodes_u[lb] ++;
            i++; nblocks++;
          } else { /* get out*/
            i = ib;
          }
        }
        i = ib;
      }

      /* from L-factor */
      for (i=0, jb=0; jb<nlb; jb++ ) {
        lb = jb*Pc + mycol;
        ib = i + Lrows[jb];
        while( i < ib ) {
          if( Lblock[i] < sf_block[lb] ) {
            nnodes_l[lb] ++;
            i++; nblocks++;
          } else {
            i = ib;
          }
        }
        i = ib;
      }
      if( mycol < nsupers % grid->npcol ) {
        lb = nlb*Pc + mycol;
        ib = i + Lrows[nlb];
        while( i < ib ) {
          if( Lblock[i] < sf_block[lb] ) {
            nnodes_l[lb] ++;
            i++; nblocks++;
          } else {
            i = ib;
          }
        }
        i = ib;
      }

#ifdef USE_ALLGATHER
      /* insert local nodes in DAG */
      if ( !(edag_supno_l = intMalloc_dist( nsupers+nblocks )) )
	ABORT("Malloc fails for edag_supno_l[].");
      iu = il = nblocks = 0;
      for( lb=0; lb<nsupers; lb++ ) {
        j  = lb / Pc;
        pc = lb % Pc;

        edag_supno_l[nblocks] = nnodes_l[lb]+nnodes_u[lb]; nblocks ++;
        if( mycol == pc ) {
          /* from U-factor */
          ib = iu + Urows[j];
          for( jb=0; jb<nnodes_u[lb]; jb++ ) {
            edag_supno_l[nblocks] = Ublock[iu];
            iu++; nblocks++;
          }
          iu = ib;

          /* from L-factor */
          ib = il + Lrows[j];
          for( jb=0; jb<nnodes_l[lb]; jb++ ) {
            edag_supno_l[nblocks] = Lblock[il];
            il++; nblocks++;
          }
          il = ib;
        }
      }
      SUPERLU_FREE( nnodes_u );

      /* form global DAG on each processor */
      MPI_Allgather( &nblocks, 1, MPI_INT, recvcnts, 1, MPI_INT, grid->comm );
      nblocks = recvcnts[0];
      rdispls[0] = 0;
      for( lb=1; lb<Pc*Pr; lb++ ) {
        rdispls[lb] = nblocks;
        nblocks += recvcnts[lb];
      }
      if ( !(recvbuf = intMalloc_dist( nblocks )) )
	ABORT("Malloc fails for recvbuf[].");
      MPI_Allgatherv( edag_supno_l, recvcnts[iam], mpi_int_t, 
                      recvbuf, recvcnts, rdispls, mpi_int_t, grid->comm );
      SUPERLU_FREE(edag_supno_l);

      if ( !(edag_supno = SUPERLU_MALLOC( nsupers * sizeof(int_t*) )) )
        ABORT("Malloc fails for edag_supno[].");
      k = 0;
      for( lb=0; lb<nsupers; lb++ ) nnodes_l[lb] = 0;
      for( p=0; p<Pc*Pr; p++ ) {
        for( lb=0; lb<nsupers; lb++ ) {
          nnodes_l[lb] += recvbuf[k];
          k += (1+recvbuf[k]);
        }
      }
      for( lb=0; lb<nsupers; lb++ ) {
        if( nnodes_l[lb] > 0 ) 
          if( !(edag_supno[lb] = intMalloc_dist( nnodes_l[lb] )) )
	    ABORT("Malloc fails for edag_supno[lb].");
        nnodes_l[lb] = 0;
      }
      k = 0;
      for( p=0; p<Pc*Pr; p++ ) {
        for( lb=0; lb<nsupers; lb++ ) {
          jb = k + recvbuf[k] + 1;
          k ++;
          for( ; k<jb; k++ ) {
            edag_supno[lb][nnodes_l[lb]] = recvbuf[k];
            nnodes_l[lb] ++;
          }
        }
      }
      SUPERLU_FREE(recvbuf);
#else
      int nlsupers = nsupers/Pc;
      if( mycol < nsupers%Pc ) nlsupers ++;

      /* insert local nodes in DAG */
      if ( !(edag_supno_l = intMalloc_dist( nlsupers+nblocks )) )
        ABORT("Malloc fails for edag_supno_l[].");
      iu = il = nblocks = 0;
      for( lb=0; lb<nsupers; lb++ ) {
        j  = lb / Pc;
        pc = lb % Pc;
        if( mycol == pc ) {
          edag_supno_l[nblocks] = nnodes_l[lb]+nnodes_u[lb]; nblocks ++;
          /* from U-factor */
          ib = iu + Urows[j];
          for( jb=0; jb<nnodes_u[lb]; jb++ ) {
            edag_supno_l[nblocks] = Ublock[iu];
            iu++; nblocks++;
          }
          iu = ib;

          /* from L-factor */
          ib = il + Lrows[j];
          for( jb=0; jb<nnodes_l[lb]; jb++ ) {
            edag_supno_l[nblocks] = Lblock[il];
            il++; nblocks++;
          }
          il = ib;
        } else if( nnodes_l[lb]+nnodes_u[lb] != 0 ) printf( " # %d: nnodes[%d]=%d+%d\n",grid->iam,lb,nnodes_l[lb],nnodes_u[lb] );
      }
      SUPERLU_FREE( nnodes_u );
      /* form global DAG on each processor */
      MPI_Allgather( &nblocks, 1, MPI_INT, recvcnts, 1, MPI_INT, grid->comm );
      nblocks = recvcnts[0];
      rdispls[0] = 0;
      for( lb=1; lb<Pc*Pr; lb++ ) {
        rdispls[lb] = nblocks;
        nblocks += recvcnts[lb];
      }
      if ( !(recvbuf = intMalloc_dist( nblocks )) )
        ABORT("Malloc fails for recvbuf[].");

      MPI_Allgatherv( edag_supno_l, recvcnts[iam], mpi_int_t,
                      recvbuf, recvcnts, rdispls, mpi_int_t, grid->comm );
      SUPERLU_FREE(edag_supno_l);

      if ( !(edag_supno = SUPERLU_MALLOC( nsupers * sizeof(int_t*) )) )
        ABORT("Malloc fails for edag_supno[].");
      k = 0;
      for( lb=0; lb<nsupers; lb++ ) nnodes_l[lb] = 0;
      for( p=0; p<Pc*Pr; p++ ) {
        yourcol = MYCOL( p, grid );

        for( lb=0; lb<nsupers; lb++ ) {
          j  = lb / Pc;
          pc = lb % Pc;
          if( yourcol == pc ) {
            nnodes_l[lb] += recvbuf[k];
            k += (1+recvbuf[k]);
          }
        }
      }
      for( lb=0; lb<nsupers; lb++ ) {
        if( nnodes_l[lb] > 0 )
          if( !(edag_supno[lb] = intMalloc_dist( nnodes_l[lb] )) )
            ABORT("Malloc fails for edag_supno[lb].");
        nnodes_l[lb] = 0;
      }
      k = 0;
      for( p=0; p<Pc*Pr; p++ ) {
        yourcol = MYCOL( p, grid );

        for( lb=0; lb<nsupers; lb++ ) {
          j  = lb / Pc;
          pc = lb % Pc;
          if( yourcol == pc ) {
            jb = k + recvbuf[k] + 1;
            k ++;
            for( ; k<jb; k++ ) {
              edag_supno[lb][nnodes_l[lb]] = recvbuf[k];
              nnodes_l[lb] ++;
            }
          }
        }
      }
      SUPERLU_FREE(recvbuf);
#endif

      /* initialize the num of child for each node */
      num_child = SUPERLU_MALLOC( nsupers * sizeof(int_t) );
      for( i=0; i<nsupers; i++ ) num_child[i] = 0;
      for( i=0; i<nsupers; i++ ) {
        for( jb=0; jb<nnodes_l[i]; jb++ ) {
          num_child[edag_supno[i][jb]]++;
        }
      }

      /* push initial leaves to the fifo queue */
      nnodes = 0;
      for( i=0; i<nsupers; i++ ) {
        if( num_child[i] == 0 ) {
          ptr = SUPERLU_MALLOC( sizeof(etree_node) );
          ptr->id = i;
          ptr->next = NULL;
          /*printf( " == push leaf %d (%d) ==\n",i,nnodes );*/
          nnodes ++;

          if( nnodes == 1 ) {
            head = ptr;
            tail = ptr;
          } else {
            tail->next = ptr;
            tail = ptr;
          }
        }
      }

      /* process fifo queue, and compute the ordering */
      i = 0;
      perm_c_supno  = SUPERLU_MALLOC( nsupers * sizeof(int_t) );
      while( nnodes > 0 ) {

        /*printf( "=== pop %d (%d) ===\n",head->id,i );*/
        ptr = head;  j = ptr->id;
        head = ptr->next;

        perm_c_supno[i] = j;
        SUPERLU_FREE(ptr);
        i++; nnodes --;

        for( jb=0; jb<nnodes_l[j]; jb++ ) {
          num_child[edag_supno[j][jb]]--;
          if( num_child[edag_supno[j][jb]] == 0 ) {
            nnodes ++;

            ptr = SUPERLU_MALLOC( sizeof(etree_node) );
            ptr->id = edag_supno[j][jb];
            ptr->next = NULL;

            /*printf( "=== push %d ===\n",ptr->id );*/
            if( nnodes == 1 ) {
              head = ptr;
              tail = ptr;
            } else {
              tail->next = ptr;
              tail = ptr;
            }
          }
        }
        /*printf( "\n" );*/
      }
      SUPERLU_FREE(num_child);

      for( lb=0; lb<nsupers; lb++ ) if( nnodes_l[lb] > 0 ) SUPERLU_FREE(edag_supno[lb] );
      SUPERLU_FREE(edag_supno);
      SUPERLU_FREE(nnodes_l);
      SUPERLU_FREE(sendcnts);
      SUPERLU_FREE(sf_block);
      SUPERLU_FREE(Ublock);
      SUPERLU_FREE(Urows);
      SUPERLU_FREE(Lblock);
      SUPERLU_FREE(Lrows);
    }
   /* ======================== *
    * end of static scheduling *
    * ======================== */

    iperm_c_supno = SUPERLU_MALLOC( nsupers * sizeof(int_t) );
    for( lb=0; lb<nsupers; lb++ ) iperm_c_supno[perm_c_supno[lb]] = lb;

    /* constructing look-ahead table to indicate the last dependency */
    int *look_ahead_l; 
    stat->num_look_aheads = num_look_aheads;

    look_ahead_l = SUPERLU_MALLOC( nsupers * sizeof(int) );
    look_ahead   = SUPERLU_MALLOC( nsupers * sizeof(int) );
    for( lb=0; lb<nsupers; lb++ ) look_ahead_l[lb] = -1;

    /* go through U-factor */
    for (lb = 0; lb < nrb; ++lb) {
      ib = lb * Pr + myrow;
      index = Llu->Ufstnz_br_ptr[lb];
      if ( index ) { /* Not an empty row */
        k = BR_HEADER;
        for (j = 0; j < index[0]; ++j) {
            jb = index[k];
            if( jb != ib ) look_ahead_l[jb] = SUPERLU_MAX( iperm_c_supno[ib], look_ahead_l[jb] );
            k += UB_DESCRIPTOR + SuperSize( index[k] );
        }
      }
    }
    if ( myrow < nsupers % grid->nprow ) {
      ib = nrb * Pr + myrow;
      index = Llu->Ufstnz_br_ptr[nrb];
      if ( index ) { /* Not an empty row */
        k = BR_HEADER;
        for (j = 0; j < index[0]; ++j) {
            jb = index[k];
            if( jb != ib ) look_ahead_l[jb] = SUPERLU_MAX( iperm_c_supno[ib], look_ahead_l[jb] );
            k += UB_DESCRIPTOR + SuperSize( index[k] );
        }
      }
    }

    if( options->SymPattern == NO ) {
      /* go through L-factor */
      for( lb=0; lb<ncb; lb++ ) {
        ib = lb * Pc + mycol;
        index = Llu->Lrowind_bc_ptr[lb];
        if ( index ) { 
            k = BC_HEADER;
            for( j=0; j<index[0]; j++ ) {
                jb = index[k];
                if( jb != ib ) look_ahead_l[jb] = SUPERLU_MAX( iperm_c_supno[ib], look_ahead_l[jb] ); 
                k += LB_DESCRIPTOR + index[k+1];
            }
        }
      }
      if( mycol < nsupers % grid->npcol ) {
        ib = ncb * Pc + mycol;
        index = Llu->Lrowind_bc_ptr[ncb];
        if ( index ) { 
          k = BC_HEADER; 
          for( j=0; j<index[0]; j++ ) {
             jb = index[k];
             if( jb != ib ) look_ahead_l[jb] = SUPERLU_MAX( iperm_c_supno[ib], look_ahead_l[jb] ); 
             k += LB_DESCRIPTOR + index[k+1];
          }
        }
      }
    }
    MPI_Allreduce( look_ahead_l , look_ahead , nsupers, MPI_INT, MPI_MAX, grid->comm );
    SUPERLU_FREE(look_ahead_l);
#ifdef ISORT
    iperm_u = SUPERLU_MALLOC( nsupers * sizeof(int_t) );
    perm_u  = SUPERLU_MALLOC( nsupers * sizeof(int_t) );
#else
    perm_u  = SUPERLU_MALLOC( 2*nsupers * sizeof(int_t) );
#endif
#if ( PRNTlevel>=1 )
    if( grid->iam == 0 ) printf( " * init: %e seconds\n",SuperLU_timer_()-tt1 );
#endif

    k = sp_ienv_dist(3); /* max supernode size */
    if ( !(Llu->ujrow = doublecomplexMalloc_dist(k*(k+1)/2)) )
	ABORT("Malloc fails for ujrow[].");

#if ( PRNTlevel>=1 )
    if ( !iam ) {
	printf(".. thresh = s_eps %e * anorm %e = %e\n", s_eps, anorm, thresh);
	printf(".. Buffer size: Lsub %ld\tLval %ld\tUsub %ld\tUval %ld\tLDA %ld\n",
	       (long int)Llu->bufmax[0], (long int)Llu->bufmax[1], 
	       (long int)Llu->bufmax[2], (long int)Llu->bufmax[3], (long int)Llu->bufmax[4]);
    }
#endif

    Lrowind_bc_ptr = Llu->Lrowind_bc_ptr;
    Lnzval_bc_ptr = Llu->Lnzval_bc_ptr;
    Ufstnz_br_ptr = Llu->Ufstnz_br_ptr;
    Unzval_br_ptr = Llu->Unzval_br_ptr;
    ToRecv = Llu->ToRecv;
    ToSendD = Llu->ToSendD;
    ToSendR = Llu->ToSendR;

    ldt = sp_ienv_dist(3); /* Size of maximum supernode */
    k = CEILING( nsupers, Pr ); /* Number of local block rows */

    if ( !(tempv2d = doublecomplexCalloc_dist(2*((size_t)ldt)*ldt)) )
        ABORT("Calloc fails for tempv2d[].");
    tempU2d = tempv2d + ldt*ldt;
    if ( !(indirect = intMalloc_dist(ldt)) )
        ABORT("Malloc fails for indirect[].");
    if ( !(iuip = intMalloc_dist(k)) )
        ABORT("Malloc fails for iuip[].");
    if ( !(ruip = intMalloc_dist(k)) )
        ABORT("Malloc fails for ruip[].");


  /* --------------------------------------------------------------- *
     Handle the first block column separately to start the pipeline.
   * --------------------------------------------------------------- */
    look_id = 0;
    msgcnt = msgcnts[0];
    send_req = send_reqs[0];
    recv_req = recv_reqs[0];

    k0 = 0;
    k  = perm_c_supno[0];
    kcol = PCOL( k, grid );
    krow = PROW( k, grid );
    if ( mycol == kcol ) {
	/* double ttt1 = SuperLU_timer_(); */
	pzgstrf2(options, nsupers, k0, k, thresh, Glu_persist, grid, Llu, 
                  U_diag_blk_send_req, stat, info);
	/* stat->time7 += SuperLU_timer_()-ttt1; */
	scp = &grid->rscp; /* The scope of process row. */

	/* Process column *kcol* multicasts numeric values of L(:,k) 
	   to process rows. */
	lk = LBj( k, grid ); /* Local block number. */
	lsub = Lrowind_bc_ptr[lk];
	lusup = Lnzval_bc_ptr[lk];
	if ( lsub ) {
	    msgcnt[0] = lsub[1] + BC_HEADER + lsub[0]*LB_DESCRIPTOR;
	    msgcnt[1] = lsub[1] * SuperSize( k );
	} else {
	    msgcnt[0] = msgcnt[1] = 0;
	}
	
	for (pj = 0; pj < Pc; ++pj) {
	    if ( ToSendR[lk][pj] != EMPTY ) {
#if ( PROFlevel>=1 )
                TIC(t1);
#endif
		MPI_Isend( lsub, msgcnt[0], mpi_int_t, pj, 
                           SLU_MPI_TAG(0,0) /* 0 */, 
                           scp->comm, &send_req[pj] );
		MPI_Isend( lusup, msgcnt[1], SuperLU_MPI_DOUBLE_COMPLEX, pj, 
                           SLU_MPI_TAG(1,0) /* 1 */, 
                           scp->comm, &send_req[pj+Pc] );
#if ( DEBUGlevel>=2 )
                printf("(%d) Send L(:,%4d): lsub %4d, lusup %4d to Pc %2d\n",
                       iam, 0, msgcnt[0], msgcnt[1], pj);
#endif
#if ( PROFlevel>=1 )
                TOC(t2, t1);
                stat->utime[COMM] += t2;
                msg_cnt += 2;
                msg_vol += msgcnt[0]*iword + msgcnt[1]*zword;
#endif
	    }
	} /* for pj ... */
    } else { /* Post immediate receives. */
	if ( ToRecv[k] >= 1 ) { /* Recv block column L(:,0). */
	    scp = &grid->rscp; /* The scope of process row. */
	    MPI_Irecv( Lsub_buf_2[0], Llu->bufmax[0], mpi_int_t, kcol,
		       SLU_MPI_TAG(0,0) /* 0 */, 
                       scp->comm, &recv_req[0] );
	    MPI_Irecv( Lval_buf_2[0], Llu->bufmax[1], SuperLU_MPI_DOUBLE_COMPLEX, kcol,
		       SLU_MPI_TAG(1,0) /* 1 */, 
                       scp->comm, &recv_req[1] );
	}
    } /* if mycol == 0 */
    factored[k] = 0;

#define IRECV_U
#ifdef  IRECV_U
    /* post receive of first U-row */
    if( myrow != krow ) {
	if ( ToRecv[k] == 2 ) { /* Recv block row U(k,:). */
	    scp = &grid->cscp; /* The scope of process column. */
	    Usub_buf = Llu->Usub_buf_2[0];
	    Uval_buf = Llu->Uval_buf_2[0];
	    MPI_Irecv( Usub_buf, Llu->bufmax[2], mpi_int_t, krow,
	               SLU_MPI_TAG(2,0) /* 2%tag_ub */, 
                       scp->comm, &recv_reqs_u[0][0] );
	    MPI_Irecv( Uval_buf, Llu->bufmax[3], SuperLU_MPI_DOUBLE_COMPLEX, krow,
	               SLU_MPI_TAG(3,0) /* 3%tag_ub */, 
                       scp->comm, &recv_reqs_u[0][1] ); 
	}
    } 
#endif

    /* ------------------------------------------
       MAIN LOOP: Loop through all block columns.
       ------------------------------------------ */
    for (k0 = 0; k0 < nsupers; ++k0) {
	k = perm_c_supno[k0];

        /* ============================================ *
	 * ======== look-ahead the new columns ======== *
         * ============================================ */
	/* tt1 = SuperLU_timer_(); */
        if( k0 == 0 ) { /* look-ahead all the columns in the window */
          kk1 = k0+1;
          kk2 = SUPERLU_MIN( k0+num_look_aheads, nsupers-1 );
        } else { /* look-ahead one new column after the current window */
          kk1 = k0+num_look_aheads;
	  kk2 = SUPERLU_MIN( kk1, nsupers-1 );
        }

	for( kk0 = kk1; kk0<=kk2; kk0++ ) { /* loop through look-ahead window */
	    kk = perm_c_supno[kk0]; /* use the ordering from static scheduling */
	    look_id = kk0 % (1+num_look_aheads); /* which column in window */

	    if( look_ahead[kk] < k0 ) { /* dose not depends on the current column */

	      kcol = PCOL( kk, grid );
	      if ( mycol == kcol ) {

	        /* Factor diagonal and subdiagonal blocks and test for exact
	           singularity.  */
	        factored[kk] = 0;
		/* double ttt1 = SuperLU_timer_(); */
	        pzgstrf2(options, nsupers, kk0, kk, thresh, Glu_persist, grid, Llu,
                         U_diag_blk_send_req, stat, info);
		/* stat->time7 += SuperLU_timer_() - ttt1; */

	        /* Process column *kcol+1* multicasts numeric values of L(:,k+1) 
	           to process rows. */
		/* ttt1 = SuperLU_timer_(); */
	        msgcnt   = msgcnts[look_id]; /* point to the proper count array */
	        send_req = send_reqs[look_id];

	        lk = LBj( kk, grid ); /* Local block number. */
	        lsub1 = Lrowind_bc_ptr[lk];
 	        if ( lsub1 ) {
		  msgcnt[0] = lsub1[1] + BC_HEADER + lsub1[0]*LB_DESCRIPTOR;
		  msgcnt[1] = lsub1[1] * SuperSize( kk );
	        } else {
		  msgcnt[0] = 0;
		  msgcnt[1] = 0;
	        }
	        scp = &grid->rscp; /* The scope of process row. */
	        for (pj = 0; pj < Pc; ++pj) {
		  if ( ToSendR[lk][pj] != EMPTY ) {
		    lusup1 = Lnzval_bc_ptr[lk];
		    MPI_Isend( lsub1, msgcnt[0], mpi_int_t, pj,
			       SLU_MPI_TAG(0,kk0) /* (4*kk0)%tag_ub */, 
                               scp->comm, &send_req[pj] );
		    MPI_Isend( lusup1, msgcnt[1], SuperLU_MPI_DOUBLE_COMPLEX, pj,
			       SLU_MPI_TAG(1,kk0) /* (4*kk0+1)%tag_ub */, 
                               scp->comm, &send_req[pj+Pc] );
		  } 
	        } 
		/* stat->time9 += SuperLU_timer_() - ttt1; */
	      } else { /* Post Recv of block column L(:,k+1). */
		/* double ttt1 = SuperLU_timer_(); */
	        if ( ToRecv[kk] >= 1 ) {
		  scp = &grid->rscp; /* The scope of process row. */
		  recv_req = recv_reqs[look_id];

		  MPI_Irecv( Lsub_buf_2[look_id], Llu->bufmax[0], mpi_int_t, kcol,
			     SLU_MPI_TAG(0,kk0) /* (4*kk0)%tag_ub */, 
                             scp->comm, &recv_req[0] );
		  MPI_Irecv( Lval_buf_2[look_id], Llu->bufmax[1], SuperLU_MPI_DOUBLE_COMPLEX, kcol, 
			     SLU_MPI_TAG(1,kk0) /* (4*kk0+1)%tag_ub */, 
                             scp->comm, &recv_req[1] );
	        } 
		/* stat->time10 += SuperLU_timer_() - ttt1; */
	      } /* if mycol == Pc(k+1) */
	    } /* if look ahead */

#ifdef IRECV_U
	    /* post irecev for U-row look-ahead */
	    krow = PROW( kk, grid );
	    if ( myrow != krow ) {
	      if ( ToRecv[kk] == 2 ) { /* post iRecv block row U(k,:). */
	        scp = &grid->cscp;     /* The scope of process column. */
	        Usub_buf = Llu->Usub_buf_2[look_id];
	        Uval_buf = Llu->Uval_buf_2[look_id];

	        MPI_Irecv( Usub_buf, Llu->bufmax[2], mpi_int_t, krow,
	                   SLU_MPI_TAG(2,kk0) /* (4*kk0+2)%tag_ub */, 
                           scp->comm, &recv_reqs_u[look_id][0] );
	        MPI_Irecv( Uval_buf, Llu->bufmax[3], SuperLU_MPI_DOUBLE_COMPLEX, krow,
	                   SLU_MPI_TAG(3,kk0) /* (4*kk0+3)%tag_ub */, 
                           scp->comm, &recv_reqs_u[look_id][1] ); 
  	      }
	    }
#endif
	} /* for each column in look-ahead window */
	/* stat->time4 += SuperLU_timer_()-tt1; */

	/* ================================= *
	 * == looking-ahead the U columns == *
	 * ================================= */
        kk1 = k0;
        kk2 = SUPERLU_MIN( k0+num_look_aheads, nsupers-1 );
	for( kk0 = kk1; kk0<kk2; kk0++ ) {
	    kk = perm_c_supno[kk0];
	    if( factoredU[kk0] != 1 && look_ahead[kk] < k0 ) {
	      kcol = PCOL( kk, grid );
	      krow = PROW( kk, grid );
	      lk = LBj( kk, grid ); /* Local block number. */

	      look_id  = kk0 % (1+num_look_aheads);
	      msgcnt   = msgcntsU[look_id];
	      recv_req = recv_reqs[look_id];

	     /* ================================================= *
	      * checking if diagonal block has been received      *
              * for panel factorization of U in look-ahead window *
	      * ================================================= */

	      if ( mycol == kcol ) {
	        flag0 = flag1 = 1;
		msgcnt[0] = msgcnt[1] = -1;
	      } else {
	        flag0 = flag1 = 0;
	        if( ToRecv[kk] >= 1 ) {
		  if( recv_req[0] != MPI_REQUEST_NULL ) {
		    MPI_Test( &recv_req[0], &flag0, &status );
		    if( flag0 ) {
		      MPI_Get_count( &status, mpi_int_t, &msgcnt[0] );
		      recv_req[0] = MPI_REQUEST_NULL;
		    } 
		  } else flag0 = 1;
		  if( recv_req[1] != MPI_REQUEST_NULL ) {
		    MPI_Test( &recv_req[1], &flag1, &status );
		    if( flag1 ) {
		      MPI_Get_count( &status, mpi_int_t, &msgcnt[1] );
		      recv_req[1] = MPI_REQUEST_NULL;
		    }
		  } else flag1 = 1;
	        } else msgcnt[0] = 0;
	      }

	      if( flag0 && flag1 ) {
		/* tt1 = SuperLU_timer_(); */
		scp = &grid->cscp; /* The scope of process column. */
		if ( myrow == krow ) {
		    factoredU[kk0] = 1;
		    /* Parallel triangular solve across process row *krow* --
		       U(k,j) = L(k,k) \ A(k,j).  */
	    	    /* double ttt2 = SuperLU_timer_(); */
#ifdef _CRAY
		    pzgstrs2(n, kk0, kk, Glu_persist, grid, Llu, stat, ftcs1, ftcs2, ftcs3);
#else
		    pzgstrs2(n, kk0, kk, Glu_persist, grid, Llu, stat);
#endif
		    /* stat->time8 += SuperLU_timer_()-ttt2; */

		    /* Multicasts U(k,:) to process columns. */
		    lk = LBi( kk, grid );
		    usub = Ufstnz_br_ptr[lk];
		    uval = Unzval_br_ptr[lk];
		    if ( usub )	{
			msgcnt[2] = usub[2];
			msgcnt[3] = usub[1];
		    } else {
			msgcnt[2] = msgcnt[3] = 0;
		    }

		    if ( ToSendD[lk] == YES ) {
			for (pi = 0; pi < Pr; ++pi) {
			    if ( pi != myrow ) {
#if ( PROFlevel>=1 )
				TIC(t1);
#endif
#if ( VAMPIR>=1 )
				VT_begin(3);
#endif
				MPI_Isend( usub, msgcnt[2], mpi_int_t, pi,
					   SLU_MPI_TAG(2,kk0) /* (4*kk0+2)%tag_ub */, 
                                           scp->comm, &send_reqs_u[look_id][pi] );
				MPI_Isend( uval, msgcnt[3], SuperLU_MPI_DOUBLE_COMPLEX, pi,
					   SLU_MPI_TAG(3,kk0) /* (4*kk0+3)%tag_ub */, 
                                           scp->comm, &send_reqs_u[look_id][pi+Pr] );
#if ( VAMPIR>=1 )
				VT_end(3);
#endif
#if ( PROFlevel>=1 )
				TOC(t2, t1);
				stat->utime[COMM] += t2;
				msg_cnt += 2;
				msg_vol += msgcnt[2]*iword + msgcnt[3]*zword;
#endif
#if ( DEBUGlevel>=2 )
				printf("(%d) Send U(%4d,:) to Pr %2d\n", iam, k, pi);
#endif
			    } /* if pi ... */
			} /* for pi ... */
		    } /* if ToSendD ... */

		    /* stat->time2 += SuperLU_timer_()-tt1; */

		} /* if myrow == krow */
	      } /* if flag0 ... */
	    } /* if factoredU[] ... */
	} /* for kk0 ... */

        /* ============================================== *
	 * == start of processing the current row of U == *
	 * ============================================== */
	knsupc = SuperSize( k );
	krow = PROW( k, grid );
	kcol = PCOL( k, grid );

	/* tt1 = SuperLU_timer_(); */
	look_id = k0 % (1+num_look_aheads);
	recv_req = recv_reqs[look_id];
	send_req = send_reqs[look_id];
	msgcnt   = msgcnts[look_id];
	Usub_buf = Llu->Usub_buf_2[look_id];
	Uval_buf = Llu->Uval_buf_2[look_id];

	if ( mycol == kcol ) {
	    lk = LBj( k, grid ); /* Local block number. */

	    for (pj = 0; pj < Pc; ++pj) {
                /* Wait for Isend to complete before using lsub/lusup. */
		if ( ToSendR[lk][pj] != EMPTY ) {
		    MPI_Wait( &send_req[pj], &status );
		    MPI_Wait( &send_req[pj+Pc], &status );
		}
	    }
	    lsub = Lrowind_bc_ptr[lk];
	    lusup = Lnzval_bc_ptr[lk];
	} else {
	    if ( ToRecv[k] >= 1 ) { /* Recv block column L(:,k). */
		scp = &grid->rscp; /* The scope of process row. */
		/* =========================================== *
		 * waiting for L(:,K) for outer-product uptate *
		 * if iam in U(k,:) then                       *
		 *  the diagonal block did not reach in time   *
		 *  for panel factorization of U(k,:)          *
		 * =========================================== */
#if ( PROFlevel>=1 )
		TIC(t1);
#endif
		if( recv_req[0] != MPI_REQUEST_NULL ) {
		  MPI_Wait( &recv_req[0], &status );
		  MPI_Get_count( &status, mpi_int_t, &msgcnt[0] );
		  recv_req[0] = MPI_REQUEST_NULL;
		} else {
		  msgcnt[0] = msgcntsU[look_id][0];
		}

		if( recv_req[1] != MPI_REQUEST_NULL ) {
		  MPI_Wait( &recv_req[1], &status );
		  MPI_Get_count( &status, MPI_DOUBLE, &msgcnt[1] );
		  recv_req[1] = MPI_REQUEST_NULL;
		} else {
		  msgcnt[1] = msgcntsU[look_id][1];
		}
#if ( PROFlevel>=1 )
		TOC(t2, t1);
		stat->utime[COMM] += t2;
#endif
#if ( DEBUGlevel>=2 )
		printf("(%d) Recv L(:,%4d): lsub %4d, lusup %4d from Pc %2d\n",
		       iam, k, msgcnt[0], msgcnt[1], kcol);
		fflush(stdout);
#endif

#if ( PRNTlevel==3 )
		++total_msg;
		if ( !msgcnt[0] ) ++zero_msg;
#endif
	    } else msgcnt[0] = 0;

	    lsub  = Lsub_buf_2[look_id];
	    lusup = Lval_buf_2[look_id];
	} /* if mycol = Pc(k) */
	/* stat->time1 += SuperLU_timer_()-tt1; */

	scp = &grid->cscp; /* The scope of process column. */

	/* tt1 = SuperLU_timer_(); */
	if ( myrow == krow ) {
	  lk = LBi( k, grid );
	  usub = Ufstnz_br_ptr[lk];
	  uval = Unzval_br_ptr[lk];

	  if ( factoredU[k0] == -1 ) {
	    /* Parallel triangular solve across process row *krow* --
	       U(k,j) = L(k,k) \ A(k,j).  */
	    /* double ttt2 = SuperLU_timer_(); */
#ifdef _CRAY
	    pzgstrs2(n, k0, k, Glu_persist, grid, Llu, stat, ftcs1, ftcs2, ftcs3);
#else
	    pzgstrs2(n, k0, k, Glu_persist, grid, Llu, stat);
#endif
	    /* stat->time8 += SuperLU_timer_() - ttt2; */

	    /* Multicasts U(k,:) to process columns. */
	    if ( usub )	{
		msgcnt[2] = usub[2];
		msgcnt[3] = usub[1];
	    } else {
		msgcnt[2] = msgcnt[3] = 0;
	    }

	    if ( ToSendD[lk] == YES ) {
		for (pi = 0; pi < Pr; ++pi) {
		    if ( pi != myrow ) {
#if ( PROFlevel>=1 )
			TIC(t1);
#endif
			MPI_Send( usub, msgcnt[2], mpi_int_t, pi,
				  SLU_MPI_TAG(2,k0) /* (4*k0+2)%tag_ub */, 
                                  scp->comm );
			MPI_Send( uval, msgcnt[3], SuperLU_MPI_DOUBLE_COMPLEX, pi,
				  SLU_MPI_TAG(3,k0) /* (4*k0+3)%tag_ub */, 
                                  scp->comm );
#if ( PROFlevel>=1 )
			TOC(t2, t1);
			stat->utime[COMM] += t2;
			msg_cnt += 2;
			msg_vol += msgcnt[2]*iword + msgcnt[3]*zword;
#endif
#if ( DEBUGlevel>=2 )
			printf("(%d) Send U(%4d,:) to Pr %2d\n", iam, k, pi);
#endif
		    } /* if pi ... */
		} /* for pi ... */
	    } /* if ToSendD ... */
	  } else {

	    /* =========================================== *
	     * waiting for U(k,:) for outer-product update *
	     * =========================================== */
	    if ( ToSendD[lk] == YES ) {
                for (pi = 0; pi < Pr; ++pi) {
                    if ( pi != myrow ) {
		      MPI_Wait( &send_reqs_u[look_id][pi],    &status );
		      MPI_Wait( &send_reqs_u[look_id][pi+Pr], &status );
		    }
		}
	    }
	    msgcnt[2] = msgcntsU[look_id][2];
	    msgcnt[3] = msgcntsU[look_id][3];
	  }
	  /* stat->time2 += SuperLU_timer_()-tt1; */
	} else { /* myrow != krow */

	   /* ========================================= *
	    * wait for U(k,:) for outer-product updates *
	    * ========================================= */

	    if ( ToRecv[k] == 2 ) { /* Recv block row U(k,:). */
#if ( PROFlevel>=1 )
		TIC(t1);
#endif

#ifdef IRECV_U
		MPI_Wait( &recv_reqs_u[look_id][0], &status );
		MPI_Get_count( &status, mpi_int_t, &msgcnt[2] );
		MPI_Wait( &recv_reqs_u[look_id][1], &status );
		MPI_Get_count( &status, SuperLU_MPI_DOUBLE_COMPLEX, &msgcnt[3] );
#else
		MPI_Recv( Usub_buf, Llu->bufmax[2], mpi_int_t, krow,
		 	  SLU_MPI_TAG(2,k0) /* (4*k0+2)%tag_ub */, 
                          scp->comm, &status );
		MPI_Get_count( &status, mpi_int_t, &msgcnt[2] );
		MPI_Recv( Uval_buf, Llu->bufmax[3], SuperLU_MPI_DOUBLE_COMPLEX, krow, 
		  	  SLU_MPI_TAG(3,k0) /* (4*k0+3)%tag_ub */, 
                          scp->comm, &status );
		MPI_Get_count( &status, SuperLU_MPI_DOUBLE_COMPLEX, &msgcnt[3] );
#endif

#if ( PROFlevel>=1 )
		TOC(t2, t1);
		stat->utime[COMM] += t2;
#endif
		usub = Usub_buf;
		uval = Uval_buf;
#if ( DEBUGlevel>=2 )
		printf("(%d) Recv U(%4d,:) from Pr %2d\n", iam, k, krow);
#endif
#if ( PRNTlevel==3 )
		++total_msg;
		if ( !msgcnt[2] ) ++zero_msg;
#endif
	    } else msgcnt[2] = 0;
	    /* stat->time6 += SuperLU_timer_()-tt1; */
	} /* if myrow == Pr(k) */

	/* 
	 * Parallel rank-k update; pair up blocks L(i,k) and U(k,j).
	 *  for (j = k+1; k < N; ++k) {
	 *     for (i = k+1; i < N; ++i) 
	 *         if ( myrow == PROW( i, grid ) && mycol == PCOL( j, grid )
	 *              && L(i,k) != 0 && U(k,j) != 0 )
	 *             A(i,j) = A(i,j) - L(i,k) * U(k,j);
	 */
	msg0 = msgcnt[0];
	msg2 = msgcnt[2];
	/* tt1 = SuperLU_timer_(); */
	if ( msg0 && msg2 ) { /* L(:,k) and U(k,:) are not empty. */
	  nsupr = lsub[1]; /* LDA of lusup. */
	  if ( myrow == krow ) { /* Skip diagonal block L(k,k). */
		lptr0 = BC_HEADER + LB_DESCRIPTOR + lsub[BC_HEADER+1];
		luptr0 = knsupc;
		nlb = lsub[0] - 1;
	  } else {
		lptr0 = BC_HEADER;
		luptr0 = 0;
		nlb = lsub[0];
	  }
	  iukp = BR_HEADER; /* Skip header; Pointer to index[] of U(k,:) */
	  rukp = 0;         /* Pointer to nzval[] of U(k,:) */
	  nub = usub[0];    /* Number of blocks in the block row U(k,:) */
	  klst = FstBlockC( k+1 );

	    
	  /* --------------------------------------------------------------
	     Update the look-ahead block columns A(:,k+1:k+num_look_ahead).
	     -------------------------------------------------------------- */
	  iukp0 = iukp;
	  rukp0 = rukp;
	  /* reorder the remaining columns in bottome-up */
          for (jj = 0; jj < nub; jj++) {
#ifdef ISORT
            iperm_u[jj] = iperm_c_supno[usub[iukp]];  /* Global block number of block U(k,j). */
            perm_u [jj] = jj;
#else
            perm_u[2*jj]   = iperm_c_supno[usub[iukp]];  /* Global block number of block U(k,j). */
            perm_u[2*jj+1] = jj;
#endif
            jb = usub[iukp];       /* Global block number of block U(k,j). */
            nsupc = SuperSize( jb );
            iukp += UB_DESCRIPTOR; /* Start fstnz of block U(k,j). */
            iukp += nsupc;
          }
	  iukp = iukp0;
#ifdef ISORT
	  isort( nub, iperm_u, perm_u );
#else
	  qsort(perm_u, (size_t)nub, 2*sizeof(int_t), &superlu_sort_perm );
#endif
	  j = jj0 = 0;
#ifdef ISORT
	  while( j < nub && iperm_u[j] <= k0+num_look_aheads ) {
#else
	  while( j < nub && perm_u[2*j] <= k0+num_look_aheads ) {
#endif
            iukp = iukp0;
            rukp = rukp0;
	    /* move the pointer in U-factor to perm_u[j] column */
#ifdef ISORT
            for (jj = 0; jj < perm_u[j]; jj++) {
#else
            for (jj = 0; jj < perm_u[2*j+1]; jj++) {
#endif
              jb = usub[iukp];       /* Global block number of block U(k,j). */
              nsupc = SuperSize( jb );
              iukp += UB_DESCRIPTOR; /* Start fstnz of block U(k,j). */

              rukp += usub[iukp - 1]; /* Move to block U(k,j+1) */
              iukp += nsupc;
            }

	    j ++;
	    jb = usub[iukp];   /* Global block number of block U(k,j). */
	    { /* unlabeled */
		jj0 ++;

		lptr = lptr0;
		luptr = luptr0;
		ljb = LBj( jb, grid ); /* Local block number of U(k,j). */
		nsupc = SuperSize( jb );
		iukp += UB_DESCRIPTOR; /* Start fstnz of block U(k,j). */

		/* Prepare to call DGEMM. */
		jj = iukp;
		while ( usub[jj] == klst ) ++jj;
		ldu = klst - usub[jj++];
		ncols = 1;
		full = 1;
		for (; jj < iukp+nsupc; ++jj) {
		    segsize = klst - usub[jj];
		    if ( segsize ) {
		        ++ncols;
			if ( segsize != ldu ) full = 0;
		        if ( segsize > ldu ) ldu = segsize;
		    }
		}
#if ( DEBUGlevel>=3 )
		++num_update;
#endif
		if ( full ) {
		    tempu = &uval[rukp];
		} else { /* Copy block U(k,j) into tempU2d. */
#if ( DEBUGlevel>=3 )
		  printf("(%d) full=%d,k=%d,jb=%d,ldu=%d,ncols=%d,nsupc=%d\n",
			 iam, full, k, jb, ldu, ncols, nsupc);
		  ++num_copy;
#endif
		    tempu = tempU2d;
		    for (jj = iukp; jj < iukp+nsupc; ++jj) {
		        segsize = klst - usub[jj];
			if ( segsize ) {
			    lead_zero = ldu - segsize;
                            for (i = 0; i < lead_zero; ++i) tempu[i] = zero;
			    tempu += lead_zero;
			    for (i = 0; i < segsize; ++i) {
			        tempu[i] = uval[rukp+i];
			    }
			    rukp += segsize;
			    tempu += segsize;
			}
		    }
		    tempu = tempU2d;
		    rukp -= usub[iukp - 1]; /* Return to start of U(k,j). */
		} /* if full ... */

		for (lb = 0; lb < nlb; ++lb) { 
		    ib = lsub[lptr]; /* Row block L(i,k). */
		    nbrow = lsub[lptr+1];  /* Number of full rows. */
		    lptr += LB_DESCRIPTOR; /* Skip descriptor. */
		    tempv = tempv2d;
#ifdef _CRAY
		    CGEMM(ftcs, ftcs, &nbrow, &ncols, &ldu, &alpha, 
			  &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr, 
			  tempu, &ldu, &beta, tempv, &ldt);
#elif defined (USE_VENDOR_BLAS)
		    zgemm_("N", "N", &nbrow, &ncols, &ldu, &alpha, 
			   &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr, 
			   tempu, &ldu, &beta, tempv, &ldt, 1, 1);
#else
		    zgemm_("N", "N", &nbrow, &ncols, &ldu, &alpha, 
			   &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr, 
			   tempu, &ldu, &beta, tempv, &ldt);
#endif
                    stat->ops[FACT] += 8 * nbrow * ldu * ncols;

		    /* Now gather the result into the destination block. */
		    if ( ib < jb ) { /* A(i,j) is in U. */
			ilst = FstBlockC( ib+1 );
			lib = LBi( ib, grid );
			index = Ufstnz_br_ptr[lib];

			/* reinitialize the pointer to each row of U */
			iuip[lib] = BR_HEADER;
			ruip[lib] = 0;

			ijb = index[iuip[lib]];
			while ( ijb < jb ) { /* Search for dest block. */
			    ruip[lib] += index[iuip[lib]+1];
			    iuip[lib] += UB_DESCRIPTOR + SuperSize( ijb );
			    ijb = index[iuip[lib]];
			}
			iuip[lib] += UB_DESCRIPTOR; /* Skip descriptor. */

			tempv = tempv2d;
			for (jj = 0; jj < nsupc; ++jj) {
			    segsize = klst - usub[iukp + jj];
			    fnz = index[iuip[lib]++];
			    if ( segsize ) { /* Nonzero segment in U(k.j). */
				ucol = &Unzval_br_ptr[lib][ruip[lib]];
				for (i = 0, it = 0; i < nbrow; ++i) {
				    rel = lsub[lptr + i] - fnz;
                                    z_sub(&ucol[rel], &ucol[rel], &tempv[it]);
                                    ++it;
				}
				tempv += ldt;
			    }
			    ruip[lib] += ilst - fnz;
			}
		    } else { /* A(i,j) is in L. */
			index = Lrowind_bc_ptr[ljb];
			ldv = index[1];   /* LDA of the dest lusup. */
			lptrj = BC_HEADER;
			luptrj = 0;
			ijb = index[lptrj];
			while ( ijb != ib ) { /* Search for dest block -- 
						 blocks are not ordered! */
			    luptrj += index[lptrj+1];
			    lptrj += LB_DESCRIPTOR + index[lptrj+1];
			    ijb = index[lptrj];
			}
			/*
			 * Build indirect table. This is needed because the
			 * indices are not sorted.
			 */
			fnz = FstBlockC( ib );
			lptrj += LB_DESCRIPTOR;
			for (i = 0; i < index[lptrj-1]; ++i) {
			    rel = index[lptrj + i] - fnz;
			    indirect[rel] = i;
			}
			nzval = Lnzval_bc_ptr[ljb] + luptrj;
			tempv = tempv2d;
			for (jj = 0; jj < nsupc; ++jj) {
			    segsize = klst - usub[iukp + jj];
			    if ( segsize ) {
/*#pragma _CRI cache_bypass nzval,tempv*/
				for (it = 0, i = 0; i < nbrow; ++i) {
				    rel = lsub[lptr + i] - fnz;
                                    z_sub(&nzval[indirect[rel]],
                                          &nzval[indirect[rel]],
                                          &tempv[it]);
                                    ++it;
				}
				tempv += ldt;
			    }
			    nzval += ldv;
			}
		    } /* if ib < jb ... */
		    lptr += nbrow;
		    luptr += nbrow;
		} /* for lb ... */
		rukp += usub[iukp - 1]; /* Move to block U(k,j+1) */
		iukp += nsupc;
	    }  /* unlabeled */

            /* ==================================== *
	     * == factorize and send if possible == *
             * ==================================== */
            kk = jb;
	    kcol = PCOL( kk, grid );
#ifdef ISORT
            kk0 = iperm_u[j-1];
#else
            kk0 = perm_u[2*(j-1)];
#endif
	    look_id = kk0 % (1+num_look_aheads);

	    if( look_ahead[kk] == k0 && kcol == mycol ) { /* current column is the last dependency */
	      look_id = kk0 % (1+num_look_aheads);

	      /* Factor diagonal and subdiagonal blocks and test for exact
	         singularity.  */
	      factored[kk] = 0;
	      /* double ttt1 = SuperLU_timer_(); */
	      pzgstrf2(options, nsupers, kk0, kk, thresh, Glu_persist, grid, Llu,
                       U_diag_blk_send_req, stat, info);
	      /* stat->time7 += SuperLU_timer_() - ttt1; */

	      /* Process column *kcol+1* multicasts numeric values of L(:,k+1) 
	         to process rows. */
	      send_req = send_reqs[look_id];
	      msgcnt   = msgcnts[look_id];

	      lk = LBj( kk, grid ); /* Local block number. */
	      lsub1 = Lrowind_bc_ptr[lk];
	      lusup1 = Lnzval_bc_ptr[lk];
 	      if ( lsub1 ) {
		  msgcnt[0] = lsub1[1] + BC_HEADER + lsub1[0]*LB_DESCRIPTOR;
		  msgcnt[1] = lsub1[1] * SuperSize( kk );
	      } else {
		  msgcnt[0] = 0;
		  msgcnt[1] = 0;
	      }

	      scp = &grid->rscp; /* The scope of process row. */
	      for (pj = 0; pj < Pc; ++pj) {
		  if ( ToSendR[lk][pj] != EMPTY ) {
#if ( PROFlevel>=1 )
                    TIC(t1);
#endif
		    MPI_Isend( lsub1, msgcnt[0], mpi_int_t, pj,
			       SLU_MPI_TAG(0,kk0) /* (4*kk0)%tag_ub */, 
                               scp->comm, &send_req[pj] );
		    MPI_Isend( lusup1, msgcnt[1], SuperLU_MPI_DOUBLE_COMPLEX, pj,
			       SLU_MPI_TAG(1,kk0) /* (4*kk0+1)%tag_ub */, 
                               scp->comm, &send_req[pj+Pc] );
#if ( PROFlevel>=1 )
                    TOC(t2, t1);
                    stat->utime[COMM] += t2;
                    msg_cnt += 2;
                    msg_vol += msgcnt[0]*iword + msgcnt[1]*zword;
#endif
#if ( DEBUGlevel>=2 )
                    printf("(%d) Send L(:,%4d): lsub %4d, lusup %4d to Pc %2d\n",
                           iam, k+1, msgcnt[0], msgcnt[1], pj);
#endif
		  }
	      } /* for pj ... */
	    }
	  }
	} /* if L(:,k) and U(k,:) not empty */
	/* stat->time3 += SuperLU_timer_()-tt1; */


        /* ================== */
	/* == post receive == */
        /* ================== */
	kk1 = SUPERLU_MIN( k0+num_look_aheads, nsupers-1 );
	for( kk0 = k0+1; kk0<=kk1; kk0++ ) {
	    kk = perm_c_supno[kk0];
	    kcol = PCOL( kk, grid );

            if( look_ahead[kk] == k0 ) {
	      if( mycol != kcol ) {
                if( ToRecv[kk] >= 1 ) {
	          scp = &grid->rscp; /* The scope of process row. */

	          look_id = kk0 % (1+num_look_aheads);
	          recv_req = recv_reqs[look_id];
	          MPI_Irecv( Lsub_buf_2[look_id], Llu->bufmax[0], mpi_int_t, kcol,
			     SLU_MPI_TAG(0,kk0) /* (4*kk0)%tag_ub */, 
                             scp->comm, &recv_req[0]);
	          MPI_Irecv( Lval_buf_2[look_id], Llu->bufmax[1], SuperLU_MPI_DOUBLE_COMPLEX, kcol, 
			     SLU_MPI_TAG(1,kk0) /* (4*kk0+1)%tag_ub */, 
                             scp->comm, &recv_req[1]);
                }
	      } else {
	        lk = LBj( kk, grid ); /* Local block number. */
	        lsub1 = Lrowind_bc_ptr[lk];
	        lusup1 = Lnzval_bc_ptr[lk];
		if( factored[kk] == -1 ) 
		{
	          /* Factor diagonal and subdiagonal blocks and test for exact
	             singularity.  */
	          factored[kk] = 0;
		  /* double ttt1 = SuperLU_timer_(); */
	          pzgstrf2(options, nsupers, kk0, kk, thresh, Glu_persist, grid, Llu,
                           U_diag_blk_send_req, stat, info);
		  /* stat->time7 += SuperLU_timer_() - ttt1; */

	          /* Process column *kcol+1* multicasts numeric values of L(:,k+1) 
	             to process rows. */

	          look_id = kk0 % (1+num_look_aheads);
	          send_req = send_reqs[look_id];
	          msgcnt   = msgcnts[look_id];

 	          if ( lsub1 ) {
		    msgcnt[0] = lsub1[1] + BC_HEADER + lsub1[0]*LB_DESCRIPTOR;
		    msgcnt[1] = lsub1[1] * SuperSize( kk );
	          } else {
		    msgcnt[0] = 0;
		    msgcnt[1] = 0;
	          }

	          scp = &grid->rscp; /* The scope of process row. */
	          for (pj = 0; pj < Pc; ++pj) {
		    if ( ToSendR[lk][pj] != EMPTY ) {
		      MPI_Isend( lsub1, msgcnt[0], mpi_int_t, pj,
			         SLU_MPI_TAG(0,kk0) /* (4*kk0)%tag_ub */, 
                                 scp->comm, &send_req[pj] );
		      MPI_Isend( lusup1, msgcnt[1], SuperLU_MPI_DOUBLE_COMPLEX, pj,
			         SLU_MPI_TAG(1,kk0) /* (4*kk0+1)%tag_ub */, 
                                 scp->comm, &send_req[pj+Pc] );
		    }
		  }
	        } /* for pj ... */
	      }
	    }
        }

	/* tt1 = SuperLU_timer_(); */
	if ( msg0 && msg2 ) { /* L(:,k) and U(k,:) are not empty. */
	    /* -----------------------------------------------
	       Update remaining blocks using block row U(k,:)
	       ----------------------------------------------- */
	    for (j = jj0; j < nub; ++j) { 
                /* == processing each of the remaining columns == */
                iukp = iukp0;
                rukp = rukp0;
#ifdef ISORT
                for (jj = 0; jj < perm_u[j]; jj++) {
#else
                for (jj = 0; jj < perm_u[2*j+1]; jj++) {
#endif
                  /* reinitilize the pointers to the begining of the */
                  /* kth column/row of L/U factors                   */
                  jb = usub[iukp];       /* Global block number of block U(k,j). */
                  nsupc = SuperSize( jb );
                  iukp += UB_DESCRIPTOR; /* Start fstnz of block U(k,j). */

                  rukp += usub[iukp - 1]; /* Move to block U(k,j+1) */
                  iukp += nsupc;
                }

		/* reinitilize the pointers to the begining of the */
		/* kth column/row of L/U factors                   */
		jb = usub[iukp];       /* Global block number of block U(k,j). */
		ljb = LBj( jb, grid ); /* Local block number of U(k,j). */
		nsupc = SuperSize( jb );
		iukp += UB_DESCRIPTOR; /* Start fstnz of block U(k,j). */

		/* Prepare to call DGEMM. */
		jj = iukp;
		while ( usub[jj] == klst ) ++jj;
		ldu = klst - usub[jj++];
		ncols = 1;
		full = 1;
		for (; jj < iukp+nsupc; ++jj) {
		    segsize = klst - usub[jj];
		    if ( segsize ) {
		        ++ncols;
			if ( segsize != ldu ) full = 0;
		        if ( segsize > ldu ) ldu = segsize;
		    }
		}
#if ( DEBUGlevel>=3 )
		printf("(%d) full=%d,k=%d,jb=%d,ldu=%d,ncols=%d,nsupc=%d\n",
		       iam, full, k, jb, ldu, ncols, nsupc);
		++num_update;
#endif
		if ( full ) {
		    tempu = &uval[rukp];
		} else { /* Copy block U(k,j) into tempU2d. */
#if ( DEBUGlevel>=3 )
		    ++num_copy;
#endif
		    tempu = tempU2d;
		    for (jj = iukp; jj < iukp+nsupc; ++jj) {
		        segsize = klst - usub[jj];
			if ( segsize ) {
			    lead_zero = ldu - segsize;
                            for (i = 0; i < lead_zero; ++i) tempu[i] = zero;
			    tempu += lead_zero;
			    for (i = 0; i < segsize; ++i)
			        tempu[i] = uval[rukp+i];
			    rukp += segsize;
			    tempu += segsize;
			}
		    }
		    tempu = tempU2d;
		    rukp -= usub[iukp - 1]; /* Return to start of U(k,j). */
		} /* if full ... */

		/* do update with the kth column of L and (k,j)th block of U */
		lptr = lptr0;
		luptr = luptr0;
		for (lb = 0; lb < nlb; lb++ ) { 
		    ib = lsub[lptr];       /* Row block L(i,k). */
		    nbrow = lsub[lptr+1];  /* Number of full rows. */
		    lptr += LB_DESCRIPTOR; /* Skip descriptor. */
		    tempv = tempv2d;
#ifdef _CRAY
		    CGEMM(ftcs, ftcs, &nbrow, &ncols, &ldu, &alpha, 
			  &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr, 
			  tempu, &ldu, &beta, tempv, &ldt);
#elif defined (USE_VENDOR_BLAS)
		    zgemm_("N", "N", &nbrow, &ncols, &ldu, &alpha, 
			   &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr, 
			   tempu, &ldu, &beta, tempv, &ldt, 1, 1);
#else
		    zgemm_("N", "N", &nbrow, &ncols, &ldu, &alpha, 
			   &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr, 
			   tempu, &ldu, &beta, tempv, &ldt);
#endif
                    stat->ops[FACT] += 8 * nbrow * ldu * ncols;

		    /* Now gather the result into the destination block. */
		    if ( ib < jb ) { /* A(i,j) is in U. */
			ilst = FstBlockC( ib+1 );
			lib = LBi( ib, grid );
			index = Ufstnz_br_ptr[lib];

			/* reinitialize the pointer to each row of U */
			iuip[lib] = BR_HEADER; 
			ruip[lib] = 0;

			ijb = index[iuip[lib]];
			while ( ijb < jb ) { /* Search for dest block. */
			    ruip[lib] += index[iuip[lib]+1];
			    iuip[lib] += UB_DESCRIPTOR + SuperSize( ijb );
			    ijb = index[iuip[lib]];
			}
			/* Skip descriptor.  Now point to fstnz index of 
			   block U(i,j). */
			iuip[lib] += UB_DESCRIPTOR;

			tempv = tempv2d;
			for (jj = 0; jj < nsupc; ++jj) {
			    segsize = klst - usub[iukp + jj];
			    fnz = index[iuip[lib]++];
			    if ( segsize ) { /* Nonzero segment in U(k.j). */
				ucol = &Unzval_br_ptr[lib][ruip[lib]];
				for (i = 0 ; i < nbrow; ++i) {
				    rel = lsub[lptr + i] - fnz;
                                    z_sub(&ucol[rel], &ucol[rel], &tempv[i]);
				}
				tempv += ldt;
			    }
			    ruip[lib] += ilst - fnz;
			}
		    } else { /* A(i,j) is in L. */
			index = Lrowind_bc_ptr[ljb];
			ldv = index[1];   /* LDA of the dest lusup. */
			lptrj = BC_HEADER;
			luptrj = 0;
			ijb = index[lptrj];
			while ( ijb != ib ) { /* Search for dest block -- 
						 blocks are not ordered! */
			    luptrj += index[lptrj+1];
			    lptrj += LB_DESCRIPTOR + index[lptrj+1];

			    ijb = index[lptrj];
			}
			/*
			 * Build indirect table. This is needed because the
			 * indices are not sorted for the L blocks.
			 */
			fnz = FstBlockC( ib );
			lptrj += LB_DESCRIPTOR;
			for (i = 0; i < index[lptrj-1]; ++i) {
			    rel = index[lptrj + i] - fnz;
			    indirect[rel] = i;
			}
			nzval = Lnzval_bc_ptr[ljb] + luptrj;
			tempv = tempv2d;
			for (jj = 0; jj < nsupc; ++jj) {
			    segsize = klst - usub[iukp + jj];
			    if ( segsize ) {
/*#pragma _CRI cache_bypass nzval,tempv*/
				for (i = 0; i < nbrow; ++i) {
				    rel = lsub[lptr + i] - fnz;
                                    z_sub(&nzval[indirect[rel]],
                                          &nzval[indirect[rel]],
                                          &tempv[i]);
				}
				tempv += ldt;
			    }
			    nzval += ldv;
			}
		    } /* if ib < jb ... */
		    lptr += nbrow;
		    luptr += nbrow;
		} /* for lb ... */
	    } /* for j ... */
	} /* if  k L(:,k) and U(k,:) are not empty */

	/* stat->time5 += SuperLU_timer_()-tt1; */
    } /* for k0 = 0, ... */

    /* ------------------------------------------
       END MAIN LOOP: for k = ...
       ------------------------------------------ */
#if ( DEBUGlevel>=2 )
    for (i = 0; i < Pr * Pc; ++i) {
        if ( iam == i ) {
            /* dPrintLblocks(iam, nsupers, grid, Glu_persist, Llu); */
            /* dPrintUblocks(iam, nsupers, grid, Glu_persist, Llu); */
            printf("(%d)\n", iam);
            PrintInt10("Recv", nsupers, Llu->ToRecv);
        }
        MPI_Barrier( grid->comm );
    }
#endif


#if ( VAMPIR>=1 )
    VT_end(100);
    VT_traceoff();
#endif

    if ( Pr*Pc > 1 ) {
	SUPERLU_FREE(Lsub_buf_2[0]); /* also free Lsub_buf_2[1] */
	SUPERLU_FREE(Lval_buf_2[0]); /* also free Lval_buf_2[1] */
	if ( Llu->bufmax[2] != 0 ) SUPERLU_FREE(Usub_buf_2[0]);
	if ( Llu->bufmax[3] != 0 ) SUPERLU_FREE(Uval_buf_2[0]);
 	if ( U_diag_blk_send_req[myrow] ) {
 	    /* wait for last Isend requests to complete, deallocate objects */ 
 	    for (krow = 0; krow < Pr; ++krow)
 		if ( krow != myrow )
                     MPI_Wait(U_diag_blk_send_req + krow, &status);
        }
 	SUPERLU_FREE(U_diag_blk_send_req);
    }
    SUPERLU_FREE(Lsub_buf_2);  
    SUPERLU_FREE(Lval_buf_2);  
    SUPERLU_FREE(Usub_buf_2);
    SUPERLU_FREE(Uval_buf_2);
    SUPERLU_FREE(iperm_c_supno);
    SUPERLU_FREE(perm_c_supno);
    SUPERLU_FREE(perm_u);
#ifdef ISORT
    SUPERLU_FREE(iperm_u);
#endif
    SUPERLU_FREE(look_ahead);
    SUPERLU_FREE(factoredU);
    SUPERLU_FREE(factored);

    for( i=0; i<=num_look_aheads; i++ ) {
       SUPERLU_FREE(send_reqs_u[i]);
       SUPERLU_FREE(recv_reqs_u[i]);
       SUPERLU_FREE(send_reqs[i]);
       SUPERLU_FREE(recv_reqs[i]);
    }
    SUPERLU_FREE(send_reqs_u);
    SUPERLU_FREE(recv_reqs_u);
    SUPERLU_FREE(send_reqs);
    SUPERLU_FREE(recv_reqs);
    for( i=0; i<=num_look_aheads; i++ ) {
       SUPERLU_FREE(msgcnts[i]);
       SUPERLU_FREE(msgcntsU[i]);
    }
    SUPERLU_FREE(msgcnts);
    SUPERLU_FREE(msgcntsU);

    SUPERLU_FREE(Llu->ujrow);
    SUPERLU_FREE(tempv2d);
    SUPERLU_FREE(indirect);
    SUPERLU_FREE(iuip);
    SUPERLU_FREE(ruip);

    /* Prepare error message. */
    if ( *info == 0 ) *info = n + 1;
#if ( PROFlevel>=1 )
    TIC(t1);
#endif
    MPI_Allreduce( info, &iinfo, 1, MPI_INT, MPI_MIN, grid->comm );
#if ( PROFlevel>=1 )
    TOC(t2, t1);
    stat->utime[COMM] += t2;
    {
	float msg_vol_max, msg_vol_sum, msg_cnt_max, msg_cnt_sum;
	
	MPI_Reduce( &msg_cnt, &msg_cnt_sum,
		   1, MPI_FLOAT, MPI_SUM, 0, grid->comm );
	MPI_Reduce( &msg_cnt, &msg_cnt_max,
		   1, MPI_FLOAT, MPI_MAX, 0, grid->comm );
	MPI_Reduce( &msg_vol, &msg_vol_sum,
		   1, MPI_FLOAT, MPI_SUM, 0, grid->comm );
	MPI_Reduce( &msg_vol, &msg_vol_max,
		   1, MPI_FLOAT, MPI_MAX, 0, grid->comm );
	if ( !iam ) {
	    printf("\tPZGSTRF comm stat:"
		   "\tAvg\tMax\t\tAvg\tMax\n"
		   "\t\t\tCount:\t%.0f\t%.0f\tVol(MB)\t%.2f\t%.2f\n",
		   msg_cnt_sum/Pr/Pc, msg_cnt_max,
		   msg_vol_sum/Pr/Pc*1e-6, msg_vol_max*1e-6);
	}
    }
#endif
    if ( iinfo == n + 1 ) *info = 0;
    else *info = iinfo;


#if ( PRNTlevel==3 )
    MPI_Allreduce( &zero_msg, &iinfo, 1, MPI_INT, MPI_SUM, grid->comm );
    if ( !iam ) printf(".. # msg of zero size\t%d\n", iinfo);
    MPI_Allreduce( &total_msg, &iinfo, 1, MPI_INT, MPI_SUM, grid->comm );
    if ( !iam ) printf(".. # total msg\t%d\n", iinfo);
#endif

#if ( DEBUGlevel>=2 )
    for (i = 0; i < Pr * Pc; ++i) {
	if ( iam == i ) {
	    zPrintLblocks(iam, nsupers, grid, Glu_persist, Llu);
	    zPrintUblocks(iam, nsupers, grid, Glu_persist, Llu);
	    printf("(%d)\n", iam);
	    PrintInt10("Recv", nsupers, Llu->ToRecv);
	}
	MPI_Barrier( grid->comm );
    }
#endif

#if ( DEBUGlevel>=3 )
    printf("(%d) num_copy=%d, num_update=%d\n", iam, num_copy, num_update);
#endif
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit pzgstrf()");
#endif

    return 0;
} /* PZGSTRF */


/************************************************************************/
/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *   Panel factorization -- block column k
 *
 *   Factor diagonal and subdiagonal blocks and test for exact singularity.
 *   Only the column processes that own block column *k* participate
 *   in the work.
 * 
 * Arguments
 * =========
 * options (input) superlu_options_t* (global)
 *         The structure defines the input parameters to control
 *         how the LU decomposition will be performed.
 *
 * nsupers (input) int_t (global)
 *         Number of supernodes.
 *
 * k0     (input) int (global)
 *        Counter of the next supernode to be factorized.
 *
 * k      (input) int (global)
 *        The column number of the block column to be factorized.
 *
 * thresh (input) double (global)
 *        The threshold value = s_eps * anorm.
 *
 * Glu_persist (input) Glu_persist_t*
 *        Global data structures (xsup, supno) replicated on all processes.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 *
 * Llu    (input/output) LocalLU_t*
 *        Local data structures to store distributed L and U matrices.
 *
 * U_diag_blk_send_req (input/output) MPI_Request*
 *        List of send requests to send down the diagonal block of U.
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics about the factorization.
 *        See SuperLUStat_t structure defined in util.h.
 *
 * info   (output) int*
 *        = 0: successful exit
 *        < 0: if info = -i, the i-th argument had an illegal value
 *        > 0: if info = i, U(i,i) is exactly zero. The factorization has
 *             been completed, but the factor U is exactly singular,
 *             and division by zero will occur if it is used to solve a
 *             system of equations.
 * </pre>
 */
static void pzgstrf2
(
 superlu_options_t *options, int_t nsupers,
 int_t k0, int_t k, double thresh, Glu_persist_t *Glu_persist, gridinfo_t *grid,
 LocalLU_t *Llu, MPI_Request *U_diag_blk_send_req, 
 SuperLUStat_t *stat, int* info
 )
{
    int    cols_left, iam, l, pkk, pr, tag;
    int    incx = 1, incy = 1;
    int    nsupr; /* number of rows in the block (LDA) */
    int    luptr;
    int_t  i, myrow, krow, j, jfst, jlst, u_diag_cnt;
    int_t  nsupc; /* number of columns in the block */
    int_t  *xsup = Glu_persist->xsup;
    doublecomplex *lusup, temp;
    doublecomplex *ujrow, *ublk_ptr; /* pointer to the U block */
    doublecomplex one = {1.0, 0.0}, alpha = {-1.0, 0.0};
    int_t  Pr;
    MPI_Status status;
    MPI_Comm comm = (grid->cscp).comm;

    /* Quick return. */
    *info = 0;

    /* Initialization. */
    iam   = grid->iam;
    Pr    = grid->nprow;
    myrow = MYROW( iam, grid );
    krow  = PROW( k, grid );
    pkk   = PNUM( PROW(k, grid), PCOL(k, grid), grid );
    j     = LBj( k, grid ); /* Local block number */
    jfst  = FstBlockC( k );
    jlst  = FstBlockC( k+1 );
    lusup = Llu->Lnzval_bc_ptr[j];
    nsupc = SuperSize( k );
    if ( Llu->Lrowind_bc_ptr[j] ) nsupr = Llu->Lrowind_bc_ptr[j][1];
    ublk_ptr = ujrow = Llu->ujrow;

    luptr = 0; /* Point to the diagonal entries. */
    cols_left = nsupc; /* supernode size */
    u_diag_cnt = 0;

    if ( U_diag_blk_send_req && U_diag_blk_send_req[myrow] ) {
        /* There are pending sends - wait for all Isend to complete */
        for (pr = 0; pr < Pr; ++pr)
            if (pr != myrow) {
                MPI_Wait(U_diag_blk_send_req + pr, &status);
	}
    }

    if ( iam == pkk ) { /* diagonal process */

        for (j = 0; j < jlst - jfst; ++j) { /* for each column in panel */

	    /* Diagonal pivot */
	    i = luptr;
            if ( options->ReplaceTinyPivot == YES ) {
                if ( slud_z_abs1(&lusup[i]) < thresh ) { /* Diagonal */
#if ( PRNTlevel>=2 )
		    printf("(%d) .. col %d, tiny pivot %e  ",
			   iam, jfst+j, lusup[i]);
#endif
		    /* Keep the new diagonal entry with the same sign. */
                    if ( lusup[i].r < 0 ) lusup[i].r = -thresh;
                    else lusup[i].r = thresh;
                    lusup[i].i = 0.0;
#if ( PRNTlevel>=2 )
		    printf("replaced by %e\n", lusup[i]);
#endif
		    ++(stat->TinyPivots);
		}
	    } 

	    for (l = 0; l < cols_left; ++l, i += nsupr, ++u_diag_cnt)
                ublk_ptr[u_diag_cnt] = lusup[i]; /* copy one row of U */

            /* Test for singularity. */
            if ( ujrow[0].r == 0.0 && ujrow[0].i == 0.0 ) {
		*info = j+jfst+1;
	    } else { /* Scale the j-th column. */
                slud_z_div(&temp, &one, &ujrow[0]);
                for (i = luptr+1; i < luptr-j+nsupr; ++i)
                    zz_mult(&lusup[i], &lusup[i], &temp);
                stat->ops[FACT] += 6*(nsupr-j-1) + 10;
	    }

	    /* Rank-1 update of the trailing submatrix. */
	    if ( --cols_left ) {
		l = nsupr - j - 1;
#ifdef _CRAY
                CGERU(&l, &cols_left, &alpha, &lusup[luptr+1], &incx,
                      &ujrow[1], &incy, &lusup[luptr+nsupr+1], &nsupr);
#else
                zgeru_(&l, &cols_left, &alpha, &lusup[luptr+1], &incx,
                       &ujrow[1], &incy, &lusup[luptr+nsupr+1], &nsupr);
#endif
                stat->ops[FACT] += 8 * l * cols_left;
	    }
	    ujrow = ublk_ptr + u_diag_cnt;  /* move to next row of U */
	    luptr += nsupr + 1;	            /* move to next column */

	} /* for column j ... */

	if ( U_diag_blk_send_req && iam == pkk ) { /* Send the U block */
	    /** ALWAYS SEND TO ALL OTHERS - TO FIX **/
	    for (pr = 0; pr < Pr; ++pr)
		if (pr != krow) {
		    /* tag = ((k0<<2)+2) % tag_ub;        */
		    /* tag = (4*(nsupers+k0)+2) % tag_ub; */
		    MPI_Isend( ublk_ptr, u_diag_cnt, SuperLU_MPI_DOUBLE_COMPLEX, pr,
			       SLU_MPI_TAG(4,k0) /* tag */, 
                               comm, U_diag_blk_send_req + pr );
		}
	    U_diag_blk_send_req[krow] = 1; /* flag outstanding Isend */
	}

    } else  { /* non-diagonal process */

	/* ================================================ *
	 * Receive the diagonal block of U                  *
         * for panel factorization of L(:,k)                *
	 * note: we block for panel factorization of L(:,k) *
	 * but panel factorization of U(:,k) don't          *
	 * ================================================ */

	/* tag = ((k0<<2)+2) % tag_ub;        */
	/* tag = (4*(nsupers+k0)+2) % tag_ub; */
        MPI_Recv( ublk_ptr, (nsupc*(nsupc+1))>>1, SuperLU_MPI_DOUBLE_COMPLEX, krow, 
                  SLU_MPI_TAG(4,k0) /* tag */, 
                  comm, &status );

	for (j = 0; j < jlst - jfst; ++j) { /* for each column in panel */
	    u_diag_cnt += cols_left;

	    if ( !lusup ) { /* empty block column */
		--cols_left;
                if ( ujrow[0].r == 0.0 && ujrow[0].i == 0.0 ) *info = j+jfst+1;
		continue;
	    }

	    /* Test for singularity. */
            if ( ujrow[0].r == 0.0 && ujrow[0].i == 0.0 ) {
		*info = j+jfst+1;
	    } else {
		/* Scale the j-th column. */
                slud_z_div(&temp, &one, &ujrow[0]);
                for (i = luptr; i < luptr+nsupr; ++i)
                    zz_mult(&lusup[i], &lusup[i], &temp);
                stat->ops[FACT] += 6*nsupr + 10;
	    }

	    /* Rank-1 update of the trailing submatrix. */
	    if ( --cols_left ) {
#ifdef _CRAY
                CGERU(&nsupr, &cols_left, &alpha, &lusup[luptr], &incx,
                      &ujrow[1], &incy, &lusup[luptr+nsupr], &nsupr);
#else
                zgeru_(&nsupr, &cols_left, &alpha, &lusup[luptr], &incx,
                       &ujrow[1], &incy, &lusup[luptr+nsupr], &nsupr);
#endif
                stat->ops[FACT] += 8 * nsupr * cols_left;
	    }

	    ujrow = ublk_ptr + u_diag_cnt; /* move to next row of U */
	    luptr += nsupr;                /* move to next column */

	} /* for column j ... */

    } /* end if pkk ... */

} /* PZGSTRF2 */


/************************************************************************/
static void pzgstrs2
/************************************************************************/
#ifdef _CRAY
(
 int_t m, int_t k0, int_t k, Glu_persist_t *Glu_persist, gridinfo_t *grid,
 LocalLU_t *Llu, SuperLUStat_t *stat, _fcd ftcs1, _fcd ftcs2, _fcd ftcs3
 )
#else
(
 int_t m, int_t k0, int_t k, Glu_persist_t *Glu_persist, gridinfo_t *grid,
 LocalLU_t *Llu, SuperLUStat_t *stat
 )
#endif
/* 
 * Purpose
 * =======
 *   Perform parallel triangular solves
 *           U(k,:) := A(k,:) \ L(k,k). 
 *   Only the process row that owns block row *k* participates
 *   in the work.
 * 
 * Arguments
 * =========
 *
 * m      (input) int (global)
 *        Number of rows in the matrix.
 *
 * k      (input) int (global)
 *        The row number of the block row to be factorized.
 *
 * Glu_persist (input) Glu_persist_t*
 *        Global data structures (xsup, supno) replicated on all processes.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 *
 * Llu    (input/output) LocalLU_t*
 *        Local data structures to store distributed L and U matrices.
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics about the factorization; 
 *        See SuperLUStat_t structure defined in util.h.
 *
 */
{
    int    iam, pkk;
    int    incx = 1;
    int    nsupr; /* number of rows in the block L(:,k) (LDA) */
    int    segsize;
    int_t  nsupc; /* number of columns in the block */
    int_t  luptr, iukp, rukp;
    int_t  b, gb, j, klst, knsupc, lk, nb;
    int_t  *xsup = Glu_persist->xsup;
    int_t  *usub;
    doublecomplex *lusup, *uval;

    /* Quick return. */
    lk = LBi( k, grid ); /* Local block number */
    if ( !Llu->Unzval_br_ptr[lk] ) return;

    /* Initialization. */
    iam  = grid->iam;
    pkk  = PNUM( PROW(k, grid), PCOL(k, grid), grid );
    klst = FstBlockC( k+1 );
    knsupc = SuperSize( k );
    usub = Llu->Ufstnz_br_ptr[lk]; /* index[] of block row U(k,:) */
    uval = Llu->Unzval_br_ptr[lk];
    nb = usub[0];
    iukp = BR_HEADER;
    rukp = 0;
    if ( iam == pkk ) {
	lk = LBj( k, grid );
	nsupr = Llu->Lrowind_bc_ptr[lk][1]; /* LDA of lusup[] */
	lusup = Llu->Lnzval_bc_ptr[lk];
    } else {
	nsupr = Llu->Lsub_buf_2[k0%(1+stat->num_look_aheads)][1]; /* LDA of lusup[] */
	lusup = Llu->Lval_buf_2[k0%(1+stat->num_look_aheads)];
    }

    /* Loop through all the row blocks. */
    for (b = 0; b < nb; ++b) {
	gb = usub[iukp];
	nsupc = SuperSize( gb );
	iukp += UB_DESCRIPTOR;

	/* Loop through all the segments in the block. */
	for (j = 0; j < nsupc; ++j) {
	    segsize = klst - usub[iukp++]; 
	    if ( segsize ) { /* Nonzero segment. */
		luptr = (knsupc - segsize) * (nsupr + 1);
#ifdef _CRAY
		CTRSV(ftcs1, ftcs2, ftcs3, &segsize, &lusup[luptr], &nsupr, 
		      &uval[rukp], &incx);
#elif defined (USE_VENDOR_BLAS)
		ztrsv_("L", "N", "U", &segsize, &lusup[luptr], &nsupr, 
		       &uval[rukp], &incx, 1, 1, 1);
#else
		ztrsv_("L", "N", "U", &segsize, &lusup[luptr], &nsupr, 
		       &uval[rukp], &incx);
#endif
		stat->ops[FACT] += segsize * (segsize + 1);
		rukp += segsize;
	    }
	}
    } /* for b ... */

} /* PDGSTRS2 */

/*
static int
probe_recv(int iam, int source, int tag, MPI_Datatype datatype, MPI_Comm comm,
	   int buf_size)
{
    MPI_Status status;
    int count; 

    MPI_Probe( source, tag, comm, &status );
    MPI_Get_count( &status, datatype, &count );
    if ( count > buf_size ) {
        printf("(%d) Recv'ed count %d > buffer size %d\n",
	       iam, count, buf_size);
	exit(-1);
    }
    return 0;
} */


#undef SLU_MPI_TAG
