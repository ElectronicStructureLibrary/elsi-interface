      SUBROUTINE PZBSSOLVER1( N, M, IM, JM, DESCM, LAMBDA, X, IX, JX,
     $                        DESCX, WORK, LWORK, IWORK, LIWORK, INFO )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, IM, JM, IX, JX, LWORK, LIWORK, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            DESCM( * ), DESCX( * ), IWORK( * )
      DOUBLE PRECISION   M( * ), LAMBDA( * ), WORK( * )
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  PZBSSOLVER1() computes all eigenvalues and (both right and
*  left) eigenvectors of 2n-by-2n complex matrix
*
*     H = [       A,        B;
*          -conj(B), -conj(A) ],
*
*  where A is n-by-n Hermitian, and B is n-by-n symmetric.
*
*  On entry, the information of H is provided in the lower triangular
*  part of
*
*     M = [ real(A)+real(B), imag(A)-imag(B);
*          -imag(A)-imag(B), real(A)-real(B) ].
*
*  The matrix M is required to be positive definite.
*
*  The structure of H leads to the following properties.
*
*     1) H is diagonalizable with n pairs of real eigenvalues
*        (lambda_i, -lambda_i).
*
*     2) The eigenvectors of H has the block structure
*
*           X = [ X_1, conj(X_2);     Y = [ X_1, -conj(X_2);
*                 X_2, conj(X_1) ],        -X_2,  conj(X_1) ],
*
*        and satisfy that
*
*           X_1**T * X_2 = X_2**T * X_1,
*           X_1**H * X_1 - X_2**H * X_2 = I,
*           Y**H * X = I,
*           H * X = X * diag(lambda, -lambda),
*           Y**H * H = diag(lambda, -lambda) * Y**H.
*
*  On exit, only the positive eigenvalues and the corresponding right
*  eigenvectors are returned.  The eigenvalues are sorted in ascending
*  order.  The eigenvectors are normalized (i.e., X = [ X_1; X_2 ] with
*  X_1**H * X_1 - X_2**H * X_2 = I).
*
*  M is destroyed on exit.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The number of rows and columns of A and B.
*          The size of M and X are N*N and (2N)*N, respectively.
*          N >= 0.
*
*  M       (local input and output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension
*          (LLD_M, LOCc(JM+2*N-1)).
*          On entry, the symmetric positive definite matrix M. Only the
*          lower triangular part of M is used to define the elements of
*          the symmetric matrix.
*          On exit, all entries of M are destroyed.
*
*  IM      (global input) INTEGER
*          The row index in the global array M indicating the first
*          row of sub( M ).
*
*  JM      (global input) INTEGER
*          The column index in the global array M indicating the
*          first column of sub( M ).
*
*  DESCM   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix M.
*
*  LAMBDA  (global output) DOUBLE PRECISION array, dimension (N)
*          On normal exit LAMBDA contains the positive eigenvalues of H
*          in ascending order.
*
*  X       (local output) COMPLEX*16 array,
*          global dimension (2N, N),
*          local dimension ( LLD_X, LOCc(JX+N-1) )
*          On normal exit X contains the normalized right eigenvectors
*          of H corresponding to the positive eigenvalues.
*
*  IX      (global input) INTEGER
*          X's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*          In this version, only IX = JX = 1 is supported.
*
*  JX      (global input) INTEGER
*          X's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*          In this version, only IX = JX = 1 is supported.
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*          DESCX( CTXT_ ) must equal DESCA( CTXT_ )
*
*  WORK    (local workspace/output) DOUBLE PRECISION array,
*          dimension (LWORK)
*          On output, WORK( 1 ) returns the minimal amount of workspace
*          needed to guarantee completion.
*          If the input parameters are incorrect, WORK( 1 ) may also be
*          incorrect.
*
*  LWORK   (local input) INTEGER
*          The length of the workspace array WORK.
*          If LWORK = -1, the LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          size for the WORK/IWORK array. The required workspace is
*          returned as the first element of WORK/IWORK and no error
*          message is issued by PXERBLA.
*
*  IWORK   (local workspace/output) INTEGER array,
*          dimension (LIWORK)
*          On output, IWORK( 1 ) returns the minimal amount of workspace
*          needed to guarantee completion.
*          If the input parameters are incorrect, IWORK( 1 ) may also be
*          incorrect.
*
*  LIWORK   (local input) INTEGER
*          The length of the workspace array IWORK.
*          If LIWORK = -1, the LIWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          size for the WORK/IWORK array. The required workspace is
*          returned as the first element of WORK/IWORK and no error
*          message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  The eigensolver did not converge.
*
*  Alignment requirements
*  ======================
*
*  This subroutine requires M and X to be distributed consistently in
*  the sense that DESCM( : ) = DESCX( : ) except for the fourth entry,
*  despite that X is complex and M is real.
*  In addition, square blocks ( i.e., MB = NB ) are required.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            ICTXT, NPROCS, NPROW, NPCOL, MYROW, MYCOL, NB,
     $                   TWON, I, LWKOPT, LIWKOPT, LLWORK, LLIWORK,
     $                   LOCALMAT, MROWS, MCOLS, LLDM, DIMV, NSPLIT,
     $                   INDD, INDE, INDTAU, INDU, INDV, INDW,
     $                   INDGAP, INDWORK, INDIBLOCK, INDISPLIT,
     $                   INDIFAIL, INDICLUSTR, INDIWORK, ITMP
      DOUBLE PRECISION   ABSTOL, ORFAC
      DOUBLE PRECISION   T_CHOL, T_FORMW, T_TRIDIAG1, T_TRIDIAG2,
     $                   T_DIAG1, T_DIAG2, T_DIAG3, T_VEC1, T_VEC2,
     $                   T_PREP
*     ..
*     .. Local Arrays ..
      INTEGER            DESCU( DLEN_ ), DESCV( DLEN_ ), DESCW( DLEN_ ),
     $                   ITMP2( 14 ), ITMP3( N )
      DOUBLE PRECISION   DTMP( 1 )
      DOUBLE PRECISION   DTMP2( 7 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DSQRT
*     ..
*     .. External Functions ..
      EXTERNAL           NUMROC, PDLAMCH, MPI_WTIME
      INTEGER            NUMROC
      DOUBLE PRECISION   PDLAMCH, MPI_WTIME
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDSYRK, PDTRMM, PDTRSM, PDLACPY, PDLASET,
     $                   PDPOTRF, PDSTEBZ, PDSTEIN, PXERBLA,
     $                   BLACS_GRIDINFO, CHK1MAT,
     $                   PDSORTEIG, PZBSAUX1, PDSSTRD, PDORMTR,
     $                   PDGEADD, PZLASCL
*     ..
*     .. Executable Statements ..
*
      T_PREP = MPI_WTIME()
      INFO = 0
      TWON = 2*N
      ICTXT = DESCM( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
      IF ( NPROW .EQ. -1 ) THEN
         INFO = -( 500+CTXT_ )
      END IF
*
*     Test the input arguments.
*
      IF ( INFO .EQ. 0 .AND. N .LT. 0 ) THEN
         INFO = -1
      END IF
      IF ( INFO .EQ. 0 )
     $   CALL CHK1MAT( TWON, 1, TWON, 1, IM, JM, DESCM, 5, INFO )
      IF ( INFO .EQ. 0 .AND. DESCM( MB_ ) .NE. DESCM( NB_ ) )
     $   INFO = -( 500+MB_ )
      IF ( INFO .EQ. 0 )
     $   CALL CHK1MAT( TWON, 1, N, 1, IX, JX, DESCX, 10, INFO )
      IF ( INFO .EQ. 0 .AND. DESCX( MB_ ) .NE. DESCX( NB_ ) )
     $   INFO = -( 1000+MB_ )
*
*     Compute required workspace.
*
      IF ( INFO .EQ. 0 ) THEN
*
*        Set up local indices for the workspace.
*
         LQUERY = LWORK .EQ. -1 .OR. LIWORK .EQ. -1
         ABSTOL = TWO*PDLAMCH( ICTXT, 'S' )
         ORFAC = PDLAMCH( ICTXT, 'P' )
         NB = DESCM( NB_ )
         LLDM = DESCM( LLD_ )
         MROWS = NUMROC( TWON, NB, MYROW, 0, NPROW )
         MCOLS = NUMROC( TWON, NB, MYCOL, 0, NPCOL )
         LOCALMAT = MAX( NB, LLDM )*MCOLS
         DESCW( 1:DLEN_ ) = DESCM( 1:DLEN_ )
         DESCU( 1:DLEN_ ) = DESCM( 1:DLEN_ )
         DESCV( 1:DLEN_ ) = DESCM( 1:DLEN_ )
         INDE = 1
         INDW = INDE + TWON
         INDU = INDW + MAX( TWON+NPROCS, LOCALMAT )
         INDV = INDU + LOCALMAT
         INDWORK = INDV + MAX( TWON, LOCALMAT )
         INDD = INDW
         INDTAU = INDV
         INDGAP = INDW + TWON
         LLWORK = LWORK - INDWORK + 1
         INDIBLOCK = 1
         INDISPLIT = INDIBLOCK + TWON
         INDIFAIL = INDISPLIT + TWON
         INDICLUSTR = INDIFAIL + N
         INDIWORK = INDICLUSTR + 2*NPROCS
         LLIWORK = LIWORK - INDIWORK + 1
*
*        Estimate the workspace required by external subroutines.
*
         CALL PDSSTRD( 'L', TWON, DTMP, 1, 1, DESCW, DTMP, DTMP, WORK,
     $        -1, ITMP )
         LWKOPT = INT( WORK( 1 ) )
         CALL PDORMTR( 'R', 'L', 'N', TWON, TWON, DTMP, 1, 1, DESCW,
     $        DTMP, DTMP, 1, 1, DESCU, WORK, -1, ITMP )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         CALL PDSTEBZ( ICTXT, 'I', 'B', TWON, ZERO, ZERO, N+1, TWON,
     $        ABSTOL, DTMP, DTMP, DIMV, NSPLIT, LAMBDA, ITMP, ITMP,
     $        DTMP2, -1, ITMP2, -1, ITMP )
         LWKOPT = MAX( LWKOPT, INT( DTMP2( 1 ) ) )
         LIWKOPT = ITMP2( 1 )
         DO I = 1, N
            ITMP3(I) = I
         END DO
         CALL PDSTEIN( TWON, DTMP, DTMP, N, LAMBDA, ITMP3, ITMP,
     $        ORFAC, DTMP, 1, 1, DESCU, WORK, -1, IWORK, -1, ITMP,
     $        ITMP, DTMP, ITMP )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         LIWKOPT = MAX( LIWKOPT, IWORK( 1 ) )
*
         LWKOPT = INDWORK - 1 + MAX( LWKOPT, LOCALMAT )
         LIWKOPT = INDIWORK - 1 + LIWKOPT
         IF ( .NOT. LQUERY ) THEN
            IF ( LWORK .LT. LWKOPT ) THEN
               INFO = -12
            ELSE IF ( LIWORK .LT. LIWKOPT ) THEN
               INFO = -14
            END IF
         END IF
      END IF
*
      IF ( INFO .NE. 0 ) THEN
         CALL PXERBLA( ICTXT, 'PZBSSOLVER1', -INFO )
         RETURN
      END IF
      WORK( 1 ) = DBLE( LWKOPT )
      IWORK( 1 ) = LIWKOPT
      IF ( LQUERY )
     $   RETURN
*
*     Quick return if possible.
*
      IF ( N .EQ. 0 )
     $   RETURN
      T_PREP = MPI_WTIME() - T_PREP
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 't_prep = ', T_PREP, ';'
*
*     Compute the Cholesky factorization M = L * L**T.
*     In case of failure, a general solver is needed.
*
      T_CHOL = MPI_WTIME()
      CALL PDPOTRF( 'L', TWON, M, 1, 1, DESCM, ITMP )
      T_CHOL = MPI_WTIME() - T_CHOL
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 't_chol = ', T_CHOL, ';'
!      CALL PDLAPRNT( TWON, TWON, M, 1, 1, DESCM, 0, 0, 'L', 6,
!     $     WORK( INDWORK ) )
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) "L=tril(L);"
      IF ( ITMP .NE. 0 ) THEN
         INFO = -2
         RETURN
      END IF
*
*     Explicitly formulate W = L**T * J * L,
*     where
*
*        J = [   0, I_n;
*             -I_n,   0 ]
*
*     is the standard symplectic matrix.
*     Only the lower triangular part of W is filled by the correct
*     values.
*
      T_FORMW = MPI_WTIME()
      CALL PDLASET( 'A', TWON, TWON, ZERO, ZERO, WORK( INDW ), 1, 1,
     $     DESCW )
      CALL PDGEADD( 'T', N, N, -ONE, M, N+1, 1, DESCM, ONE,
     $     WORK( INDW ), 1, 1, DESCW )
      CALL PDTRADD( 'U', 'T', N, N, -ONE, M, N+1, N+1, DESCM, ONE,
     $     WORK( INDW ), N+1, 1, DESCW )
      CALL PDTRMM( 'R', 'L', 'N', 'N', TWON, N, ONE, M, 1, 1, DESCM,
     $     WORK( INDW ), 1, 1, DESCW )
      CALL PDLACPY( 'A', N, N, WORK( INDW ), 1, 1, DESCW,
     $     WORK( INDV ), 1, 1, DESCV )
      CALL PDGEADD( 'T', N, N, -ONE, WORK( INDV ), 1, 1, DESCV,
     $     ONE, WORK( INDW ), 1, 1, DESCW )
      T_FORMW = MPI_WTIME() - T_FORMW
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 't_formw = ', T_FORMW, ';'
!      CALL PDLAPRNT( TWON, TWON, WORK( INDW ), 1, 1, DESCW, 0, 0, 'W',
!     $     6, WORK( INDWORK ) )
*
*     Tridiagonalization: U**T * W * U = T.
*
      T_TRIDIAG1 = MPI_WTIME()
      CALL PDSSTRD( 'L', TWON, WORK( INDW ), 1, 1, DESCW, WORK( INDE ),
     $     WORK( INDTAU ), WORK( INDWORK ), LLWORK, ITMP )
!      CALL PDLAPRNT( TWON, TWON, WORK( INDW ), 1, 1, DESCW, 0, 0, 'T',
!     $     6, WORK( INDWORK ) )
      DO I = 1, TWON-1
         CALL PDELGET( 'A', ' ', WORK( INDE + I-1 ), WORK( INDW ), I+1,
     $        I, DESCW )
         WORK( INDE + I-1 ) = -WORK( INDE + I-1 )
!         IF ( MYROW+MYCOL .EQ. 0 )
!     $      WRITE( *, * ) 'E(', I, ', 1) =', WORK( INDE + I-1 ), ';'
      END DO
      T_TRIDIAG1 = MPI_WTIME() - T_TRIDIAG1
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 't_tridiag1 = ', T_TRIDIAG1, ';'
*
*     Formulate sqrt(1/2) * L * U.
*
      T_TRIDIAG2 = MPI_WTIME()
      CALL PDLASET( 'A', TWON, TWON, ZERO, ZERO, WORK( INDU ), 1, 1,
     $     DESCU )
      CALL PDTRADD( 'L', 'N', TWON, TWON, DSQRT( ONE/TWO ), M, 1, 1,
     $     DESCM, ZERO, WORK( INDU ), 1, 1, DESCU )
      CALL PDORMTR( 'R', 'L', 'N', TWON, TWON, WORK( INDW ), 1, 1,
     $     DESCW, WORK( INDTAU ), WORK( INDU ), 1, 1, DESCU,
     $     WORK( INDWORK ), LLWORK, ITMP )
      T_TRIDIAG2 = MPI_WTIME() - T_TRIDIAG2
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 't_tridiag2 = ', T_TRIDIAG2, ';'
!      CALL PDLAPRNT( TWON, TWON, WORK( INDU ), 1, 1, DESCU, 0, 0, 'U',
!     $     6, WORK( INDWORK ) )
*
*     Diagonalize the tridiagonal matrix
*
*        D**H * (-i*T) * D = V * diag(lambda) * V**T,
*
*     where D = diag{1,i,i^2,i^3,...}.
*     Only the positive half of the eigenpairs are computed.
*
      DO I = 0, TWON-1
         WORK( INDD + I ) = ZERO
      END DO
*
      T_DIAG1 = MPI_WTIME()
      CALL PDSTEBZ( ICTXT, 'I', 'B', TWON, ZERO, ZERO, N+1, TWON,
     $     ABSTOL, WORK( INDD ), WORK( INDE ), DIMV, NSPLIT, LAMBDA,
     $     IWORK( INDIBLOCK ), IWORK( INDISPLIT ), WORK( INDWORK ),
     $     LLWORK, IWORK( INDIWORK ), LLIWORK, ITMP )
      T_DIAG1 = MPI_WTIME() - T_DIAG1
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 't_diag1 = ', T_DIAG1, ';'
      IF ( ITMP .NE. 0 ) THEN
         INFO = ITMP
         WRITE( *, * ) '% PDSTEBZ fails with INFO =', INFO
         RETURN
      END IF
!      IF ( MYROW+MYCOL .EQ. 0 ) THEN
!         DO I = 1, TWON
!            WRITE( *, * ) 'lambda0(', I, ', 1) =', LAMBDA( I ), ';'
!         END DO
!      END IF
*
      T_DIAG2 = MPI_WTIME()
      CALL PDSTEIN( TWON, WORK( INDD ), WORK( INDE ), DIMV, LAMBDA,
     $     IWORK( INDIBLOCK ), IWORK( INDISPLIT ), ORFAC,
     $     WORK( INDV ), 1, 1, DESCV, WORK( INDWORK ), LLWORK,
     $     IWORK( INDIWORK ), LLIWORK, IWORK( INDIFAIL ),
     $     IWORK( INDICLUSTR ), WORK( INDGAP ), ITMP )
      IF ( ITMP .NE. 0 ) THEN
         INFO = ITMP
         WRITE( *, * ) '% PDSTEIN fails with INFO =', INFO
         RETURN
      END IF
      T_DIAG2 = MPI_WTIME() - T_DIAG2
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 't_diag2 = ', T_DIAG2, ';'
      T_DIAG3 = MPI_WTIME()
      CALL PDSORTEIG( TWON, N, LAMBDA, WORK( INDV ), 1, 1, DESCV, 0 )
*
*     Orthogonalize the eigenvectors.
*
      CALL PDSYRK( 'L', 'T', N, TWON, ONE, WORK( INDV ), 1, 1, DESCV,
     $     ZERO, WORK( INDW ), 1, 1, DESCW )
      CALL PDPOTRF( 'L', N, WORK( INDW ), 1, 1, DESCW, ITMP )
      CALL PDTRSM( 'R', 'L', 'T', 'N', TWON, N, ONE, WORK( INDW ),
     $     1, 1, DESCW, WORK( INDV ), 1, 1, DESCV )
      T_DIAG3 = MPI_WTIME() - T_DIAG3
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 't_diag3 = ', T_DIAG3, ';'
!      CALL PDLAPRNT( TWON, N, WORK( INDV ), 1, 1, DESCV, 0, 0, 'V0',
!     $     6, WORK( INDWORK ) )
*
*     Restore the eigenvectors of H:
*
*        X = Q * L**{-T} * U * D * V * Lambda**{1/2},
*        Y = Q * L * U * D * V * Lambda**{-1/2},
*
*     where
*
*        Q = sqrt(1/2) * [ I_n, -i*I_n;
*                          I_n,  i*I_n ].
*
      T_VEC1 = MPI_WTIME()
*
*     Scale V by Lambda**{-1/2}.
*
      DO I = 1, N
         CALL PDSCAL( TWON, ONE/DSQRT( LAMBDA( I ) ),
     $        WORK( INDV ), 1, I, DESCV, 1 )
      END DO
      T_VEC1 = MPI_WTIME() - T_VEC1
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 't_vec1 = ', T_VEC1, ';'
!      CALL PDLAPRNT( TWON, TWON, WORK( INDU ), 1, 1, DESCU, 0, 0, 'Z',
!     $     6, WORK( INDWORK ) )
*
*     Formulate
*
*        Y = Q * L * U * D * V * Lambda**{-1/2},
*        X = Q * L**{-T} * U * D * V * Lambda**{1/2},
*
*     Once Y is calculated, X is constructed from Y through
*
*        Y = [ Y1; Y2 ], X = [ Y1; -Y2 ].
*
*     Use WORK and M as workspace for real and imaginary parts,
*     respectively.
*
      T_VEC2 = MPI_WTIME()
      CALL PZBSAUX1( N, WORK( INDV ), 1, 1, DESCV, X, IX, JX, DESCX,
     $     WORK( INDU ), 1, 1, DESCU, WORK( INDWORK ), 1, 1, DESCW,
     $     M, IM, JM, DESCM )
      CALL PZLASCL( 'G', -ONE, ONE, N, N, X, IX+N, JX, DESCX, ITMP )
      T_VEC2 = MPI_WTIME() - T_VEC2
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 't_vec2 = ', T_VEC2, ';'
!      CALL PZLAPRNT( TWON, N, X, IX, JX, DESCX, 0, 0, 'X', 6,
!     $     WORK( INDWORK ) )
*
      WORK( 1 ) = DBLE( LWKOPT )
      IWORK( 1 ) = LIWKOPT
*
      RETURN
*
*     End of PZBSSOLVER1().
*
      END
