      SUBROUTINE PZBSEIG( SOLVER, N, A, IA, JA, DESCA, B, IB, JB, DESCB,
     $                    LAMBDA, X, IX, JX, DESCX, WORK, LWORK, RWORK,
     $                    LRWORK, IWORK, LIWORK, INFO )
*
      IMPLICIT NONE
      INCLUDE 'solver.f'
*
*     .. Scalar Arguments ..
      INTEGER            SOLVER, N, IA, JA, IB, JB, IX, JX, LWORK,
     $                   LRWORK, LIWORK, INFO
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( * ), B( * ), X( * ), WORK( * )
      DOUBLE PRECISION   LAMBDA( * ), RWORK( * )
      INTEGER            DESCA( * ), DESCB( * ), DESCX( * ), IWORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZBSEIG() computes all eigenvalues and eigenvectors (both right and
*  left) of a 2n-by-2n complex matrix
*
*     H = [       A,        B;
*          -conj(B), -conj(A) ],
*
*  where A is n-by-n Hermitian, and B is n-by-n symmetric.
*  In addition, the matrix
*
*     [      A,       B;
*      conj(B), conj(A) ]
*
*  is required to be positive definite.
*
*  Only the lower triangular parts of A and B are referenced.
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
*  This routine requires square block decomposition ( MB_A = NB_A ).
*
*  Arguments
*  =========
*
*  SOLVER  (global input) INTEGER
*          Currently supported solvers:
*             BSE_FULLBSE
*             BSE_TDA
*             BSE_TDA + BSE_LAPACK_HEEVR
*             BSE_TDA + BSE_LAPACK_HEEV
*             BSE_TDA + BSE_LAPACK_HEEVD
*             BSE_TDA + BSE_LAPACK_HEEVX
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrices sub( A ) and sub( B ).
*          N >= 0.
*
*  A       (local input) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          This array contains the local pieces of the N-by-N Hermitian
*          distributed matrix sub( A ). The leading N-by-N lower
*          triangular part of sub( A ) contains the lower triangular
*          part of the distributed matrix, and its strictly upper
*          triangular part is not referenced.
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  B       (local input) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_B, LOCc(JB+N-1)).
*          This array contains the local pieces of the N-by-N symmetric
*          distributed matrix sub( B ). The leading N-by-N lower
*          triangular part of sub( B ) contains the lower triangular
*          part of the distributed matrix, and its strictly upper
*          triangular part is not referenced.
*          B is not referenced if TDA is used.
*
*  IB      (global input) INTEGER
*          The row index in the global array B indicating the first
*          row of sub( B ).
*
*  JB      (global input) INTEGER
*          The column index in the global array B indicating the
*          first column of sub( B ).
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*
*  LAMBDA  (global output) DOUBLE PRECISION array, dimension (N)
*          On normal exit LAMBDA contains the positive eigenvalues of H
*          in ascending order.
*
*  X       (local output) COMPLEX*16 array,
*          global dimension (2N, N) for full BSE, and (N, N) for TDA,
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
*  This subroutine requires square blocks ( i.e., MB_X = NB_X ).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, NB,
     $                   LWKOPT, LRWKOPT, LIWKOPT, LLRWORK, LLIWORK,
     $                   ACOLS, LLDA, MCOLS, LLDM, INDM, INDRWORK,
     $                   INDGAP, INDIFAIL, INDICLUSTER, INDIWORK, ITMP1,
     $                   ITMP2
      DOUBLE PRECISION   T_EMBED, T_SLV
*     ..
*     .. Local Arrays ..
      INTEGER            DESCM( DLEN_ ), DESCW( DLEN_ )
      DOUBLE PRECISION   DDUM( 3 )
      COMPLEX*16         ZDUM( 3 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, IAND, INT
*     ..
*     .. External Functions ..
      EXTERNAL           NUMROC, MPI_WTIME
      INTEGER            NUMROC
      DOUBLE PRECISION   MPI_WTIME
*     ..
*     .. External Subroutines ..
      EXTERNAL           PZEMBED1, PZBSSOLVER1, PZHEEV, PZHEEVD,
     $                   PZHEEVX, PZHEEVR, PXERBLA, BLACS_GRIDINFO,
     $                   CHK1MAT, PCHK2MAT, DESCSET
!      EXTERNAL           PDLAPRNT, PZLAPRNT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      IF ( NPROW .EQ. -1 ) THEN
         INFO = -( 600+CTXT_ )
      END IF
*
*     Test the input arguments.
*
      IF ( INFO .EQ. 0 .AND. N .LT. 0 )
     $   INFO = -1
      IF ( INFO .EQ. 0 )
     $   CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
      IF ( INFO .EQ. 0 )
     $   CALL CHK1MAT( N, 2, N, 2, IB, JB, DESCB, 10, INFO )
c      IF ( INFO .EQ. 0 )
c     $   CALL CHK1MAT( 2*N, 2, N, 2, IX, JX, DESCX, 15, INFO )
      IF ( INFO .EQ. 0 .AND. DESCX( MB_ ) .NE. DESCX( NB_ ) )
     $   INFO = -( 1500+MB_ )
*
*     Compute required workspace.
*
      IF ( INFO .EQ. 0 ) THEN
         LQUERY = LWORK .EQ. -1 .OR. LIWORK .EQ. -1
         NB = DESCX( NB_ )
*
         IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_FULLBSE ) THEN
            LLDA = DESCA( LLD_ )
            ACOLS = NUMROC( N, NB, MYCOL, 0, NPCOL )
            LLDM = DESCX( LLD_ )
            MCOLS = NUMROC( 2*N, NB, MYCOL, 0, NPCOL )
            CALL DESCSET( DESCM, 2*N, 2*N, NB, NB, DESCX( RSRC_ ),
     $           DESCX( CSRC_ ), ICTXT, LLDM )
            INDM = 1
            INDRWORK = INDM + LLDM*MCOLS
            LLRWORK = LRWORK - INDRWORK + 1
            CALL PZBSSOLVER1( N, DDUM, 1, 1, DESCM, DDUM, X, IX, JX,
     $           DESCX, RWORK, -1, IWORK, -1, INFO )
            LWKOPT = 1
            LRWKOPT = INT( RWORK( 1 ) )
            LIWKOPT = IWORK( 1 )
            LRWKOPT = MAX( LRWKOPT, LLDA*ACOLS ) + INDRWORK - 1
         ELSE IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_TDA ) THEN
            CALL DESCSET( DESCW, N, N, NB, NB, DESCX( RSRC_ ),
     $           DESCX( CSRC_ ), ICTXT, DESCX( LLD_ ) )
            IF ( IAND( SOLVER, BSE_VARIANT ) .EQ. BSE_LAPACK_HEEV )
     $           THEN
               CALL PZHEEV( 'V', 'L', N, A, IA, JA, DESCA, LAMBDA,
     $              X, IX, JX, DESCW, WORK, -1, RWORK, -1, INFO )
               LWKOPT = INT( WORK( 1 ) )
               LRWKOPT = INT( RWORK( 1 ) )
               LIWKOPT = 1
            ELSE IF ( IAND( SOLVER, BSE_VARIANT ) .EQ.
     $           BSE_LAPACK_HEEVD ) THEN
               CALL PZHEEVD( 'V', 'L', N, A, IA, JA, DESCA, LAMBDA,
     $              X, IX, JX, DESCW, WORK, -1, RWORK, -1, IWORK, -1,
     $              INFO )
               LWKOPT = INT( WORK( 1 ) )
               LRWKOPT = INT( RWORK( 1 ) )
               LIWKOPT = IWORK( 1 )
            ELSE IF ( IAND( SOLVER, BSE_VARIANT ) .EQ.
     $           BSE_LAPACK_HEEVX ) THEN
               INDGAP = 1
               INDRWORK = INDGAP + NPROW*NPCOL
               INDIFAIL = 1
               INDICLUSTER = INDIFAIL + N
               INDIWORK = INDICLUSTER + 2*NPROW*NPCOL
               CALL PZHEEVX( 'V', 'A', 'L', N, A, IA, JA, DESCA,
     $              ZERO, ZERO, 0, 0, ZERO, ITMP1, ITMP2, LAMBDA, -ONE,
     $              X, IX, JX, DESCW, ZDUM, -1, DDUM, -1, IWORK, -1,
     $              IWORK, IWORK, RWORK, INFO )
               LWKOPT = INT( ZDUM( 1 ) )
               LRWKOPT = INT( DDUM( 1 ) ) + INDRWORK - 1
               LIWKOPT = IWORK( 1 ) + INDIWORK - 1
            ELSE IF ( IAND( SOLVER, BSE_VARIANT ) .EQ.
     $           BSE_LAPACK_HEEVR ) THEN
               CALL PZHEEVR( 'V', 'A', 'L', N, A, IA, JA, DESCA,
     $              ZERO, ZERO, 0, 0, ITMP1, ITMP2, LAMBDA,
     $              X, IX, JX, DESCW, ZDUM, -1, DDUM, -1, IWORK, -1,
     $              INFO )
               LWKOPT = INT( ZDUM( 1 ) )
               LRWKOPT = INT( DDUM( 1 ) )
               LIWKOPT = IWORK( 1 )
            END IF
         END IF
*
         IF ( .NOT. LQUERY ) THEN
            IF ( LWORK .LT. LWKOPT )
     $         INFO = -17
            IF ( INFO .EQ. 0 .AND. LIWORK .LT. LIWKOPT )
     $         INFO = -19
         END IF
      END IF
*
      IF ( INFO .NE. 0 ) THEN
         CALL PXERBLA( ICTXT, 'PZBSEIG', -INFO )
         RETURN
      END IF
*
      WORK( 1 ) = DCMPLX( DBLE( LWKOPT ), ZERO )
      RWORK( 1 ) = DBLE( LRWKOPT )
      IWORK( 1 ) = LIWKOPT
      IF ( LQUERY )
     $   RETURN
*
*     Quick return if possible.
*
      IF ( N .EQ. 0 )
     $   RETURN
*
*     Call the eigensolver.
*
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 'SOLVER = ', SOLVER, ';'
      IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_FULLBSE ) THEN
*
*        Full BSE.
*
         T_EMBED = MPI_WTIME()
         CALL PZEMBED1( N, A, IA, JA, DESCA, B, IB, JB, DESCB,
     $        RWORK( INDM ), 1, 1, DESCM, RWORK( INDRWORK ) )
         T_EMBED = MPI_WTIME() - T_EMBED
!         IF ( MYROW+MYCOL .EQ. 0 )
!     $      WRITE( *, * ) 't_embed = ', T_EMBED, ';'
!         CALL PDLAPRNT( 2*N, 2*N, RWORK( INDM ), 1, 1, DESCM, 0, 0, 'M',
!     $        6, RWORK( INDRWORK ) )
!         IF ( MYROW+MYCOL .EQ. 0 )
!     $      WRITE( *, * ) "MM=tril(M)+tril(M)'-diag(diag(M));"
         T_SLV = MPI_WTIME()
         CALL PZBSSOLVER1( N, RWORK( INDM ), 1, 1, DESCM, LAMBDA,
     $        X, IX, JX, DESCX, RWORK( INDRWORK ), LLRWORK, IWORK,
     $        LIWORK, INFO )
         T_SLV = MPI_WTIME() - T_SLV
!         IF ( MYROW+MYCOL .EQ. 0 )
!     $      WRITE( *, * ) 't_slv = ', T_SLV, ';'
*
      ELSE IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_TDA ) THEN
*
*        Tamm--Dancoff approximation (TDA).
*
!
!        PZHEEV does not work yet.
!
         IF ( IAND( SOLVER, BSE_VARIANT ) .EQ. BSE_LAPACK_HEEV ) THEN
            CALL PZHEEV( 'V', 'L', N, A, IA, JA, DESCA, LAMBDA,
     $           X, IX, JX, DESCW, WORK, LWORK, RWORK, LRWORK, INFO )
         ELSE IF ( IAND( SOLVER, BSE_VARIANT ) .EQ.
     $        BSE_LAPACK_HEEVD ) THEN
            CALL PZHEEVD( 'V', 'L', N, A, IA, JA, DESCA, LAMBDA,
     $           X, IX, JX, DESCW, WORK, LWORK, RWORK, LRWORK,
     $           IWORK, LIWORK, INFO )
         ELSE IF ( IAND( SOLVER, BSE_VARIANT ) .EQ.
     $        BSE_LAPACK_HEEVX ) THEN
            LLRWORK = LRWORK - INDRWORK + 1
            LLIWORK = LIWORK - INDIWORK + 1
            CALL PZHEEVX( 'V', 'A', 'L', N, A, IA, JA, DESCA,
     $           ZERO, ZERO, 0, 0, ZERO, ITMP1, ITMP2, LAMBDA, -ONE,
     $           X, IX, JX, DESCW, WORK, LWORK, RWORK( INDRWORK ),
     $           LLRWORK, IWORK( INDIWORK ), LLIWORK,
     $           IWORK( INDIFAIL ), IWORK( INDICLUSTER ),
     $           RWORK( INDGAP ), INFO )
         ELSE IF ( IAND( SOLVER, BSE_VARIANT ) .EQ.
     $        BSE_LAPACK_HEEVR ) THEN
            CALL PZHEEVR( 'V', 'A', 'L', N, A, IA, JA, DESCA,
     $           ZERO, ZERO, 0, 0, ITMP1, ITMP2, LAMBDA,
     $           X, IX, JX, DESCW, WORK, LWORK, RWORK, LRWORK,
     $           IWORK, LIWORK, INFO )
         END IF
*
      END IF
*
      WORK( 1 ) = DCMPLX( DBLE( LWKOPT ), ZERO )
      RWORK( 1 ) = DBLE( LRWKOPT )
      IWORK( 1 ) = LIWKOPT
*
      RETURN
*
*     End of PZBSEIG().
*
      END
