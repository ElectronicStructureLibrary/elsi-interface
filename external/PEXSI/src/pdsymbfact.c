/// @file pdsymbfact.c
/// @brief Symbolic factorization routine using real arithmetic.
///
///
/// @date 2013-07-20 modified by Lin Lin.
/// @date 2016-02-23 Compatible with SuperLU_DIST_v4.3
#include <math.h>
#include "superlu_ddefs.h"
#define  PRNTlevel 0

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef Add_
#define lsame_ lsame
#endif

int lsame_(char *a, char *b);

/// @brief pdsymbfact performs symbolic factorization that can be
/// reused.
///
/// pdsymbfact is a simplified version of the driver subroutine pdgssvx
/// from SuperLU_DIST. For its use see SuperLUMatrix for more
/// information.
	void
#ifdef SUPERLU_DIST_MAJOR_VERSION
#if SUPERLU_DIST_MAJOR_VERSION < 5
pdsymbfact(superlu_options_t *options, SuperMatrix *A,
					 ScalePermstruct_t *ScalePermstruct, gridinfo_t *grid,
					 LUstruct_t *LUstruct, SuperLUStat_t *stat,
					 int *numProcSymbFact, int *info, double *totalMemory,
					 double *maxMemory)
#else
pdsymbfact(superlu_dist_options_t *options, SuperMatrix *A,
					 ScalePermstruct_t *ScalePermstruct, gridinfo_t *grid,
					 LUstruct_t *LUstruct, SuperLUStat_t *stat,
					 int *numProcSymbFact, int *info, double *totalMemory,
					 double *maxMemory)
#endif
#else
pdsymbfact(superlu_dist_options_t *options, SuperMatrix *A,
					 ScalePermstruct_t *ScalePermstruct, gridinfo_t *grid,
					 LUstruct_t *LUstruct, SuperLUStat_t *stat,
					 int *numProcSymbFact, int *info, double *totalMemory,
					 double *maxMemory)
#endif
{
	NRformat_loc *Astore;
	SuperMatrix GA;      /* Global A in NC format */
	NCformat *GAstore;
	double   *a_GA;
	SuperMatrix GAC;      /* Global A in NCP format (add n end pointers) */
	NCPformat *GACstore;
	Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
	Glu_freeable_t *Glu_freeable;
	/* The nonzero structures of L and U factors, which are
		 replicated on all processrs.
		 (lsub, xlsub) contains the compressed subscript of
		 supernodes in L.
		 (usub, xusub) contains the compressed subscript of
		 nonzero segments in U.
		 If options->Fact != SamePattern_SameRowPerm, they are
		 computed by SYMBFACT routine, and then used by PDDISTRIBUTE
		 routine. They will be freed after PDDISTRIBUTE routine.
		 If options->Fact == SamePattern_SameRowPerm, these
		 structures are not used.                                  */
	fact_t   Fact;
	double   *a;
	int_t    *colptr, *rowind;
	int_t    *perm_r; /* row permutations from partial pivoting */
	int_t    *perm_c; /* column permutation vector */
	int_t    *etree;  /* elimination tree */
	int_t    *rowptr, *colind;  /* Local A in NR*/
	int_t    colequ, Equil, factored, job, notran, rowequ, need_value;
	int_t    i, iinfo, j, irow, m, n, nnz, permc_spec;
	int_t    nnz_loc, m_loc, fst_row, icol;
	int      iam;
	int      ldx;  /* LDA for matrix X (local). */
	char     equed[1], norm[1];
	double   *C, *R, *C1, *R1, amax, anorm, colcnd, rowcnd;
	double   *X, *b_col, *b_work, *x_col;
	double   t;
	float    GA_mem_use;    /* memory usage by global A */
	float    dist_mem_use; /* memory usage during distribution */
#ifdef SUPERLU_DIST_MAJOR_VERSION
#if SUPERLU_DIST_MAJOR_VERSION < 5
	mem_usage_t num_mem_usage, symb_mem_usage;
#else
	superlu_dist_mem_usage_t num_mem_usage, symb_mem_usage;
#endif
#else
	superlu_dist_mem_usage_t num_mem_usage, symb_mem_usage;
#endif
#if ( PRNTlevel>= 2 )
	double   dmin, dsum, dprod;
#endif

	/* Structures needed for parallel symbolic factorization */
	int_t *sizes, *fstVtxSep, parSymbFact;
	int   noDomains, nprocs_num;
	MPI_Comm symb_comm; /* communicator for symbolic factorization */
	int   col, key; /* parameters for creating a new communicator */
	Pslu_freeable_t Pslu_freeable;
	float  flinfo;

	/* Initialization. */
	m       = A->nrow;
	n       = A->ncol;
	Astore  = (NRformat_loc *) A->Store;
	nnz_loc = Astore->nnz_loc;
	m_loc   = Astore->m_loc;
	fst_row = Astore->fst_row;
	a       = (double *) Astore->nzval;
	rowptr  = Astore->rowptr;
	colind  = Astore->colind;
	sizes   = NULL;
	fstVtxSep = NULL;
	symb_comm = MPI_COMM_NULL;
	symb_mem_usage.total = 0.;

  /* Suggested from Valgrind debugging by Patrick Seewald */
  stat->peak_buffer    = 0.0;


	/* Test the input parameters. */
	*info = 0;
	Fact = options->Fact;
	if ( Fact < 0 || Fact > FACTORED )
		*info = -1;
	else if ( options->RowPerm < 0 || options->RowPerm > MY_PERMR )
		*info = -1;
	else if ( options->ColPerm < 0 || options->ColPerm > MY_PERMC )
		*info = -1;
	else if ( options->IterRefine < 0 || options->IterRefine > SLU_EXTRA )
		*info = -1;
	else if ( options->IterRefine == SLU_EXTRA ) {
		*info = -1;
		fprintf(stderr, "Extra precise iterative refinement yet to support.");
	} else if ( A->nrow != A->ncol || A->nrow < 0 || A->Stype != SLU_NR_loc
							|| A->Dtype != SLU_D || A->Mtype != SLU_GE )
		*info = -2;
	if ( *info ) {
		i = -(*info);
    // pxerbla("pdsymbfact", grid, -*info);
    // SuperLU_DIST v4.3
    pxerr_dist("pdgssvx", grid, -*info);
		return;
	}

	factored = (Fact == FACTORED);
	Equil = (!factored && options->Equil == YES);
	notran = (options->Trans == NOTRANS);
	parSymbFact = options->ParSymbFact;

	iam = grid->iam;
	job = 5;
	if ( factored || (Fact == SamePattern_SameRowPerm && Equil) ) {
		rowequ = (ScalePermstruct->DiagScale == ROW) ||
			(ScalePermstruct->DiagScale == BOTH);
		colequ = (ScalePermstruct->DiagScale == COL) ||
			(ScalePermstruct->DiagScale == BOTH);
	} else rowequ = colequ = FALSE;

#if ( PRNTlevel >= 1 )
	if( !iam ) fprintf(stderr,"Entering pdsymbfact.\n");
#endif

	/* The following arrays are replicated on all processes. */
	perm_r = ScalePermstruct->perm_r;
	perm_c = ScalePermstruct->perm_c;
	etree = LUstruct->etree;
	R = ScalePermstruct->R;
	C = ScalePermstruct->C;
	/********/

#if ( DEBUGlevel>=1 )
	CHECK_MALLOC(iam, "Enter pdsymbfact()");
#endif

	/* Not factored & ask for equilibration */
	if ( Equil && Fact != SamePattern_SameRowPerm ) {
		/* Allocate storage if not done so before. */
		switch ( ScalePermstruct->DiagScale ) {
			case NOEQUIL:
				if ( !(R = (double *) doubleMalloc_dist(m)) )
					ABORT("Malloc fails for R[].");
				if ( !(C = (double *) doubleMalloc_dist(n)) )
					ABORT("Malloc fails for C[].");
				ScalePermstruct->R = R;
				ScalePermstruct->C = C;
				break;
			case ROW:
				if ( !(C = (double *) doubleMalloc_dist(n)) )
					ABORT("Malloc fails for C[].");
				ScalePermstruct->C = C;
				break;
			case COL:
				if ( !(R = (double *) doubleMalloc_dist(m)) )
					ABORT("Malloc fails for R[].");
				ScalePermstruct->R = R;
				break;
		}
	}

	/* ------------------------------------------------------------
		 Diagonal scaling to equilibrate the matrix.
		 ------------------------------------------------------------*/
	if ( Equil ) {
#if ( DEBUGlevel>=1 )
		CHECK_MALLOC(iam, "Enter equil");
#endif
		t = SuperLU_timer_();

		if ( Fact == SamePattern_SameRowPerm ) {
			/* Reuse R and C. */
			switch ( ScalePermstruct->DiagScale ) {
				case NOEQUIL:
					break;
				case ROW:
					irow = fst_row;
					for (j = 0; j < m_loc; ++j) {
						for (i = rowptr[j]; i < rowptr[j+1]; ++i) {
							a[i] *= R[irow];       /* Scale rows. */
						}
						++irow;
					}
					break;
				case COL:
					for (j = 0; j < m_loc; ++j)
						for (i = rowptr[j]; i < rowptr[j+1]; ++i){
							icol = colind[i];
							a[i] *= C[icol];          /* Scale columns. */
						}
					break;
				case BOTH:
					irow = fst_row;
					for (j = 0; j < m_loc; ++j) {
						for (i = rowptr[j]; i < rowptr[j+1]; ++i) {
							icol = colind[i];
							a[i] *= R[irow] * C[icol]; /* Scale rows and cols. */
						}
						++irow;
					}
					break;
			}
		} else { /* Compute R & C from scratch */
			/* Compute the row and column scalings. */
	    pdgsequ(A, R, C, &rowcnd, &colcnd, &amax, &iinfo, grid);

      if ( iinfo > 0 ) {
        if ( iinfo <= m ) {
#if ( PRNTlevel>=1 )
          fprintf(stderr, "The " IFMT "-th row of A is exactly zero\n", iinfo);
#endif
        } else {
#if ( PRNTlevel>=1 )
          fprintf(stderr, "The " IFMT "-th column of A is exactly zero\n", iinfo-n);
#endif
        }
      } else if ( iinfo < 0 ) return;

			/* Equilibrate matrix A if it is badly-scaled. */
			pdlaqgs(A, R, C, rowcnd, colcnd, amax, equed);

			if ( lsame_(equed, "R") ) {
				ScalePermstruct->DiagScale = rowequ = ROW;
        rowequ = ROW;
			} else if ( lsame_(equed, "C") ) {
				ScalePermstruct->DiagScale = colequ = COL;
        colequ = COL;
			} else if ( lsame_(equed, "B") ) {
				ScalePermstruct->DiagScale = BOTH;
				rowequ = ROW;
				colequ = COL;
			} else ScalePermstruct->DiagScale = NOEQUIL;

#if ( PRNTlevel>=1 )
			if ( !iam ) {
				printf(".. equilibrated? *equed = %c\n", *equed);
				/*fflush(stdout);*/
			}
#endif
		} /* if Fact ... */

		stat->utime[EQUIL] = SuperLU_timer_() - t;
#if ( DEBUGlevel>=1 )
		CHECK_MALLOC(iam, "Exit equil");
#endif
	} /* if Equil ... */

	if ( !factored ) { /* Skip this if already factored. */
		/*
		 * Gather A from the distributed compressed row format to
		 * global A in compressed column format.
		 * Numerical values are gathered only when a row permutation
		 * for large diagonal is sought after.
		 */
		if ( Fact != SamePattern_SameRowPerm &&
				 (parSymbFact == NO || options->RowPerm != NO) ) {

#ifdef SUPERLU_DIST_MAJOR_VERSION
#if SUPERLU_DIST_MAJOR_VERSION < 6
			need_value = (options->RowPerm == LargeDiag);
#else
			need_value = (options->RowPerm == LargeDiag_MC64);
#endif
#else
			need_value = (options->RowPerm == LargeDiag);
#endif

			pdCompRow_loc_to_CompCol_global(need_value, A, grid, &GA);

			GAstore = (NCformat *) GA.Store;
			colptr = GAstore->colptr;
			rowind = GAstore->rowind;
			nnz = GAstore->nnz;
			GA_mem_use = (nnz + n + 1) * sizeof(int_t);

			if ( need_value ) {
				a_GA = (double *) GAstore->nzval;
				GA_mem_use += nnz * sizeof(double);
			} else assert(GAstore->nzval == NULL);
		}

		/* ------------------------------------------------------------
			 Find the row permutation for A.
			 ------------------------------------------------------------*/
		if ( options->RowPerm != NO ) {
			t = SuperLU_timer_();
			if ( Fact != SamePattern_SameRowPerm ) {
				if ( options->RowPerm == MY_PERMR ) { /* Use user's perm_r. */
					/* Permute the global matrix GA for symbfact() */
					for (i = 0; i < colptr[n]; ++i) {
						irow = rowind[i];
						rowind[i] = perm_r[irow];
					}
				} else { /* options->RowPerm == LargeDiag */
					/* Get a new perm_r[] */
					if ( job == 5 ) {
						/* Allocate storage for scaling factors. */
						if ( !(R1 = doubleMalloc_dist(m)) )
							ABORT("SUPERLU_MALLOC fails for R1[]");
						if ( !(C1 = doubleMalloc_dist(n)) )
							ABORT("SUPERLU_MALLOC fails for C1[]");
					}

          // SuperLU_DIST v4.3
          if ( !iam ) {
            /* Process 0 finds a row permutation */
            // SuperLU_DIST v4.3
            iinfo = dldperm_dist(job, m, nnz, colptr, rowind, a_GA,
                perm_r, R1, C1);

            MPI_Bcast( &iinfo, 1, mpi_int_t, 0, grid->comm );
            if ( iinfo == 0 ) {
              MPI_Bcast( perm_r, m, mpi_int_t, 0, grid->comm );
              if ( job == 5 && Equil ) {
                MPI_Bcast( R1, m, MPI_DOUBLE, 0, grid->comm );
                MPI_Bcast( C1, n, MPI_DOUBLE, 0, grid->comm );
              }
            }
          } else {
            MPI_Bcast( &iinfo, 1, mpi_int_t, 0, grid->comm );
            if ( iinfo == 0 ) {
              MPI_Bcast( perm_r, m, mpi_int_t, 0, grid->comm );
              if ( job == 5 && Equil ) {
                MPI_Bcast( R1, m, MPI_DOUBLE, 0, grid->comm );
                MPI_Bcast( C1, n, MPI_DOUBLE, 0, grid->comm );
              }
            }
          }

          if ( iinfo && job == 5) {
            SUPERLU_FREE(R1);
            SUPERLU_FREE(C1);
          }
#if ( PRNTlevel>=2 )
          dmin = dmach_dist("Overflow");
          dsum = 0.0;
          dprod = 1.0;
#endif
          if ( iinfo == 0 ) {
            if ( job == 5 ) {
              if ( Equil ) {
                for (i = 0; i < n; ++i) {
                  R1[i] = exp(R1[i]);
                  C1[i] = exp(C1[i]);
                }

                /* Scale the distributed matrix */
                irow = fst_row;
                for (j = 0; j < m_loc; ++j) {
                  for (i = rowptr[j]; i < rowptr[j+1]; ++i) {
                    icol = colind[i];
                    a[i] *= R1[irow] * C1[icol];
#if ( PRNTlevel>=2 )
                    if ( perm_r[irow] == icol ) { /* New diagonal */
                      if ( job == 2 || job == 3 )
                        dmin = SUPERLU_MIN(dmin, fabs(a[i]));
                      else if ( job == 4 )
                        dsum += fabs(a[i]);
                      else if ( job == 5 )
                        dprod *= fabs(a[i]);
                    }
#endif
                  }
                  ++irow;
                }

                /* Multiply together the scaling factors. */
                if ( rowequ ) for (i = 0; i < m; ++i) R[i] *= R1[i];
                else for (i = 0; i < m; ++i) R[i] = R1[i];
                if ( colequ ) for (i = 0; i < n; ++i) C[i] *= C1[i];
                else for (i = 0; i < n; ++i) C[i] = C1[i];

                ScalePermstruct->DiagScale = BOTH;
                rowequ = colequ = 1;

              } /* end Equil */

              /* Now permute global A to prepare for symbfact() */
              for (j = 0; j < n; ++j) {
                for (i = colptr[j]; i < colptr[j+1]; ++i) {
                  irow = rowind[i];
                  rowind[i] = perm_r[irow];
                }
              }
              SUPERLU_FREE (R1);
              SUPERLU_FREE (C1);
            } else { /* job = 2,3,4 */
              for (j = 0; j < n; ++j) {
                for (i = colptr[j]; i < colptr[j+1]; ++i) {
                  irow = rowind[i];
                  rowind[i] = perm_r[irow];
                } /* end for i ... */
              } /* end for j ... */
            } /* end else job ... */
          } else { /* if iinfo != 0 */
            for (i = 0; i < m; ++i) perm_r[i] = i;
          }

#if ( PRNTlevel>=2 )
					if ( job == 2 || job == 3 ) {
						if ( !iam ) printf("\tsmallest diagonal %e\n", dmin);
					} else if ( job == 4 ) {
						if ( !iam ) printf("\tsum of diagonal %e\n", dsum);
					} else if ( job == 5 ) {
						if ( !iam ) printf("\t product of diagonal %e\n", dprod);
					}
#endif

				} /* end if options->RowPerm ... */

				t = SuperLU_timer_() - t;
				stat->utime[ROWPERM] = t;
#if ( PRNTlevel>=1 )
        if ( !iam ) printf(".. LDPERM job " IFMT "\t time: %.2f\n", job, t);
#endif
			} /* end if Fact ... */
		} else { /* options->RowPerm == NOROWPERM */
			for (i = 0; i < m; ++i) perm_r[i] = i;
		}

#if ( DEBUGlevel>=2 )
		if ( !iam ) PrintInt10("perm_r",  m, perm_r);
#endif
	} /* end if (!factored) */

	if ( !factored || options->IterRefine ) {
		/* Compute norm(A), which will be used to adjust small diagonal. */
		if ( notran ) *(unsigned char *)norm = '1';
		else *(unsigned char *)norm = 'I';
		anorm = pdlangs(norm, A, grid);
#if ( PRNTlevel>=1 )
		if ( !iam ) printf(".. anorm %e\n", anorm);
#endif
	}

	/* ------------------------------------------------------------
		 Perform the LU factorization.
		 ------------------------------------------------------------*/
	if ( !factored ) {
		t = SuperLU_timer_();
		/*
		 * Get column permutation vector perm_c[], according to permc_spec:
		 *   permc_spec = NATURAL:  natural ordering
		 *   permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
		 *   permc_spec = MMD_ATA:  minimum degree on structure of A'*A
		 *   permc_spec = METIS_AT_PLUS_A: METIS on structure of A'+A
		 *   permc_spec = PARMETIS: parallel METIS on structure of A'+A
		 *   permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
		 */
		permc_spec = options->ColPerm;


		if ( parSymbFact == YES || permc_spec == PARMETIS ) {
			nprocs_num = grid->nprow * grid->npcol;
			// Default value
			noDomains = (int) ( pow(2, ((int) LOG2( nprocs_num ))));
			if( *numProcSymbFact != 0 ){
				// User-provided value
				noDomains = MIN(noDomains,*numProcSymbFact);
			}
#if ( PRNTlevel >= 1 )
			if( !iam ) fprintf(stderr,"Using %d processors for ParMETIS.\n", noDomains);
#endif

			/* create a new communicator for the first noDomains processes in
				 grid->comm */
			key = iam;
			if (iam < noDomains) col = 0;
			else col = MPI_UNDEFINED;
			MPI_Comm_split (grid->comm, col, key, &symb_comm );

			permc_spec = PARMETIS; /* only works with PARMETIS */
		}

		if ( permc_spec != MY_PERMC && Fact == DOFACT ) {
			if ( permc_spec == PARMETIS ) {
				/* Get column permutation vector in perm_c.                    *
				 * This routine takes as input the distributed input matrix A  *
				 * and does not modify it.  It also allocates memory for       *
				 * sizes[] and fstVtxSep[] arrays, that contain information    *
				 * on the separator tree computed by ParMETIS.                 */
#if ( PRNTlevel >= 1 )
				if( !iam ) fprintf(stderr,"Before get_perm_c_parmetis.");
#endif
				flinfo = get_perm_c_parmetis(A, perm_r, perm_c, nprocs_num,
																		 noDomains, &sizes, &fstVtxSep,
																		 grid, &symb_comm);
#if ( PRNTlevel >= 1 )
				if( !iam ) fprintf(stderr,"After get_perm_c_parmetis.");
#endif
        if (flinfo > 0) {
#if ( PRNTlevel>=1 )
          fprintf(stderr, "Insufficient memory for get_perm_c parmetis\n");
#endif
          *info = flinfo;
          return;
        }
			} else {
				get_perm_c_dist(iam, permc_spec, &GA, perm_c);
			}
		}

		stat->utime[COLPERM] = SuperLU_timer_() - t;

		/* Compute the elimination tree of Pc*(A'+A)*Pc' or Pc*A'*A*Pc'
			 (a.k.a. column etree), depending on the choice of ColPerm.
			 Adjust perm_c[] to be consistent with a postorder of etree.
			 Permute columns of A to form A*Pc'. */
		if ( Fact != SamePattern_SameRowPerm ) {
			if ( parSymbFact == NO ) {
				int_t *GACcolbeg, *GACcolend, *GACrowind;

				sp_colorder(options, &GA, perm_c, etree, &GAC);

				/* Form Pc*A*Pc' to preserve the diagonal of the matrix GAC. */
				GACstore = (NCPformat *) GAC.Store;
				GACcolbeg = GACstore->colbeg;
				GACcolend = GACstore->colend;
				GACrowind = GACstore->rowind;
				for (j = 0; j < n; ++j) {
					for (i = GACcolbeg[j]; i < GACcolend[j]; ++i) {
						irow = GACrowind[i];
						GACrowind[i] = perm_c[irow];
					}
				}

				/* Perform a symbolic factorization on Pc*Pr*A*Pc' and set up
					 the nonzero data structures for L & U. */
#if ( PRNTlevel>=1 )
				if ( !iam )
		  printf(".. symbfact(): relax " IFMT ", maxsuper " IFMT ", fill " IFMT "\n",
								 sp_ienv_dist(2), sp_ienv_dist(3), sp_ienv_dist(6));
#endif
				t = SuperLU_timer_();
				if ( !(Glu_freeable = (Glu_freeable_t *)
							 SUPERLU_MALLOC(sizeof(Glu_freeable_t))) )
					ABORT("Malloc fails for Glu_freeable.");

				/* Every process does this. */
				iinfo = symbfact(options, iam, &GAC, perm_c, etree,
												 Glu_persist, Glu_freeable);

				stat->utime[SYMBFAC] = SuperLU_timer_() - t;
				if ( iinfo <= 0 ) { /* Successful return */
					QuerySpace_dist(n, -iinfo, Glu_freeable, &symb_mem_usage);
#if ( PRNTlevel>=1 )
          if ( !iam ) {
            printf("\tNo of supers %ld\n", (long long) Glu_persist->supno[n-1]+1);
            printf("\tSize of G(L) %ld\n", (long long) Glu_freeable->xlsub[n]);
            printf("\tSize of G(U) %ld\n", (long long) Glu_freeable->xusub[n]);
            printf("\tint %d, short %d, float %d, double %d\n",
                (int) sizeof(int_t), (int) sizeof(short),
                (int) sizeof(float), (int) sizeof(double));
            printf("\tSYMBfact (MB):\tL\\U %.2f\ttotal %.2f\texpansions " IFMT "\n",
                symb_mem_usage.for_lu*1e-6,
                symb_mem_usage.total*1e-6,
                symb_mem_usage.expansions);
          }
#endif
        } else { /* symbfact out of memory */
#if ( PRNTlevel>=1 )
          if ( !iam )
            fprintf(stderr,"symbfact() error returns " IFMT "\n",iinfo);
#endif
          *info = iinfo;
          return;
        }
			} /* end serial symbolic factorization */
			else {  /* parallel symbolic factorization */
				t = SuperLU_timer_();
#if ( PRNTlevel >= 1 )
				if( !iam ) fprintf(stderr,"Before symbfact_dist.");
#endif
				flinfo = symbfact_dist(nprocs_num, noDomains, A, perm_c, perm_r,
															 sizes, fstVtxSep, &Pslu_freeable,
															 &(grid->comm), &symb_comm,
															 &symb_mem_usage);
#if ( PRNTlevel >= 1 )
				if( !iam ) fprintf(stderr,"After symbfact_dist.");
#endif
				stat->utime[SYMBFAC] = SuperLU_timer_() - t;
        if (flinfo > 0) {
#if ( PRNTlevel>=1 )
          fprintf(stderr, "Insufficient memory for parallel symbolic factorization.");
#endif
          *info = flinfo;
          return;
        }
			}

			/* Destroy GA */
			if ( parSymbFact == NO || options->RowPerm != NO )
				Destroy_CompCol_Matrix_dist(&GA);
			if ( parSymbFact == NO )
				Destroy_CompCol_Permuted_dist(&GAC);

		} /* end if Fact ... */

		if (sizes) SUPERLU_FREE (sizes);
		if (fstVtxSep) SUPERLU_FREE (fstVtxSep);
		if (symb_comm != MPI_COMM_NULL)
			MPI_Comm_free (&symb_comm);

		if (parSymbFact == NO || Fact == SamePattern_SameRowPerm) {
			/* Apply column permutation to the original distributed A */
			for (j = 0; j < nnz_loc; ++j) colind[j] = perm_c[colind[j]];

			/* Distribute Pc*Pr*diag(R)*A*diag(C)*Pc' into L and U storage.
NOTE: the row permutation Pc*Pr is applied internally in the
distribution routine. */
			t = SuperLU_timer_();
			dist_mem_use = pddistribute(Fact, n, A, ScalePermstruct,
																	Glu_freeable, LUstruct, grid);
			stat->utime[DIST] = SuperLU_timer_() - t;

			/* Deallocate storage used in symbolic factorization. */
			if ( Fact != SamePattern_SameRowPerm ) {
				iinfo = symbfact_SubFree(Glu_freeable);
				SUPERLU_FREE(Glu_freeable);
			}
		} else {
			/* Distribute Pc*Pr*diag(R)*A*diag(C)*Pc' into L and U storage.
NOTE: the row permutation Pc*Pr is applied internally in the
distribution routine. */
			/* Apply column permutation to the original distributed A */
			for (j = 0; j < nnz_loc; ++j) colind[j] = perm_c[colind[j]];

			t = SuperLU_timer_();
			dist_mem_use = ddist_psymbtonum(Fact, n, A, ScalePermstruct,
																			&Pslu_freeable, LUstruct, grid);
			if (dist_mem_use > 0)
				ABORT ("Not enough memory available for dist_psymbtonum\n");

			stat->utime[DIST] = SuperLU_timer_() - t;
		}

		/*if (!iam) printf ("\tDISTRIBUTE time  %8.2f\n", stat->utime[DIST]);*/

		/* DO NOT Perform numerical factorization in parallel. */
		t = SuperLU_timer_();
		// pdgstrf(options, m, n, anorm, LUstruct, grid, stat, info);
		stat->utime[FACT] = SuperLU_timer_() - t;

		// LL: The memory output is modified.
		{
			int_t TinyPivots;
			float for_lu, total, max, avg, temp;

			dQuerySpace_dist(n, LUstruct, grid, stat, &num_mem_usage);

			if (parSymbFact == TRUE) {
				/* The memory used in the redistribution routine
					 includes the memory used for storing the symbolic
					 structure and the memory allocated for numerical
					 factorization */
				temp = SUPERLU_MAX(symb_mem_usage.total, -dist_mem_use);
				if ( options->RowPerm != NO )
					temp = SUPERLU_MAX(temp, GA_mem_use);
			} else {
				temp = SUPERLU_MAX (
														symb_mem_usage.total + GA_mem_use, /* symbfact step */
														symb_mem_usage.for_lu + dist_mem_use +
														num_mem_usage.for_lu  /* distribution step */
													 );
			}

			temp = SUPERLU_MAX(temp, num_mem_usage.total);

			MPI_Allreduce( &temp, &max,
									1, MPI_FLOAT, MPI_MAX, grid->comm );
			MPI_Allreduce( &temp, &avg,
									1, MPI_FLOAT, MPI_SUM, grid->comm );
			MPI_Allreduce( &stat->TinyPivots, &TinyPivots, 1, mpi_int_t,
										 MPI_SUM, grid->comm );
			stat->TinyPivots = TinyPivots;

			MPI_Allreduce( &num_mem_usage.for_lu, &for_lu,
									1, MPI_FLOAT, MPI_SUM, grid->comm );
			MPI_Allreduce( &num_mem_usage.total, &total,
									1, MPI_FLOAT, MPI_SUM, grid->comm );

			*totalMemory = avg * 1e-6;
			*maxMemory   = max * 1e-6;


			if ( options->PrintStat ) {
        if (!iam) {
          printf("\n** Memory Usage **********************************\n");
          printf("** NUMfact space (MB): (sum-of-all-processes)\n"
              "    L\\U :        %8.2f |  Total : %8.2f\n",
              for_lu * 1e-6, total * 1e-6);
          printf("** Total highmark (MB):\n"
              "    Sum-of-all : %8.2f | Avg : %8.2f  | Max : %8.2f\n",
              avg * 1e-6,
              avg / grid->nprow / grid->npcol * 1e-6,
              max * 1e-6);
          printf("**************************************************\n");
        }
			}
		}

	} /* end if (!factored) */


#if ( PRNTlevel>=1 )
	if ( !iam ) printf(".. DiagScale = %d\n", ScalePermstruct->DiagScale);
#endif

	/* Deallocate R and/or C if it was not used. */
	if ( Equil && Fact != SamePattern_SameRowPerm ) {
		switch ( ScalePermstruct->DiagScale ) {
			case NOEQUIL:
				SUPERLU_FREE(R);
				SUPERLU_FREE(C);
				break;
			case ROW:
				SUPERLU_FREE(C);
				break;
			case COL:
				SUPERLU_FREE(R);
				break;
		}
	}

#if ( DEBUGlevel>=1 )
	CHECK_MALLOC(iam, "Exit pdsymbfact()");
#endif

}
