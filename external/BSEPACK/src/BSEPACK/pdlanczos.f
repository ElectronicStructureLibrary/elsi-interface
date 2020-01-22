      SUBROUTINE PDLANCZOS( N, STEPS, A, IA, JA, DESCA, X, IX, JX,
     $                      DESCX, ALPHA, BETA, SCL, INFO )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, STEPS, IA, JA, IX, JX, INFO
      DOUBLE PRECISION   SCL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( * ), X( * ), ALPHA( * ), BETA( * )
      INTEGER            DESCA( * ), DESCX( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLANCZOS() computes the Lanczos decomposition
*
*     A * X = X * T + [rank 1],
*
*  where A is symmetric and X**T * X = I.
*
*  No argument check is performed, i.e.,
*  all arguments are assumed to be valid.
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
*          N is the number of rows and columns of A.
*          N >= 0.
*
*  STEPS   (global input) INTEGER
*          The number of Lanczos steps.
*
*  A       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          This array contains the local pieces of the N-by-N symmetric
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
*  X       (local input/output) DOUBLE PRECISION pointer into the local
*          memory to an array of dimension (LLD_X, LOCc(JX+1)).
*          On entry, X( 1:N, 1 ) is used to store the starting vector x,
*          which does not need to be normalized.
*          On exit, X( 1:N, j ) is used to store X( :, j ).
*
*  IX      (global input) INTEGER
*          The row index in the global array X indicating the first
*          row of sub( X ).
*          In this version, only IX = JX = 1 is supported.
*
*  JX      (global input) INTEGER
*          The column index in the global array X indicating the
*          first column of sub( X ).
*          In this version, only IX = JX = 1 is supported.
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*  ALPHA   (global output) DOUBLE PRECISION array, dimension (STEPS)
*          The diagonal of the tridiagonal matrix T.
*
*  BETA    (global output) DOUBLE PRECISION array, dimension (STEPS)
*          The subdiagonal of the tridiagonal matrix T.
*
*  SCL     (global output) DOUBLE PRECISION
*          SCL = x**T * x,
*          where x is the unnormalized starting vector on entry.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          > 0:  Lanczos process breaks down after INFO steps.
*
*  Alignment requirements
*  ======================
*
*  The distributed matrices A and X must satisfy the following
*  properties:
*
*     DESCA( MB_ ) = DESCA( NB_ )
*     DESCA( MB_ ) = DESCX( MB_ )
*     DESCA( NB_ ) = DESCX( NB_ )
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
      INTEGER            ICTXT, J
      DOUBLE PRECISION   DTMP, TMAX, TOL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DSQRT, MAX
*     ..
*     .. External Functions ..
      EXTERNAL           PDLAMCH
      DOUBLE PRECISION   PDLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDAXPY, PDSCAL, PDSYMM, PDDOT_ALL
*
      INFO = 0
      ICTXT = DESCX( CTXT_ )
      TOL = 10*N*PDLAMCH( ICTXT, 'E' )
      TMAX = ZERO
*
*     scl = X(:, 1)**T * X(:, 1)
*     X(:, 1) = X(:, 1)/sqrt(scl)
*
      CALL PDDOT_ALL( N, SCL, X, IX, JX, DESCX, 1, X, IX, JX, DESCX, 1 )
      CALL PDSCAL( N, DSQRT( ONE/SCL ), X, IX, JX, DESCX, 1 )
*
*     Lanczos process.
*
      DO J = 1, STEPS
*
*        Matrix vector multiplication.
*
*        X(:, j+1) = A*X(:, j)
*
         CALL PDSYMM( 'L', 'L', N, 1, ONE, A, IA, JA, DESCA,
     $        X, IX, JX+J-1, DESCX, ZERO, X, IX, JX+J, DESCX )
*
*        Gram--Schmidt orthogonalization.
*
*        alpha(j) = X(:, j+1)**H * X(:, j)
*        X(:, j+1) = X(:, j+1) - alpha(j)*X(:, j) - beta(j-1)*X(:, j-1)
*
         CALL PDDOT_ALL( N, ALPHA( J ), X, IX, JX+J-1, DESCX, 1,
     $        X, IX, JX+J, DESCX, 1 )
         TMAX = MAX( TMAX, ALPHA( J ) )
         CALL PDAXPY( N, -ALPHA( J ), X, IX, JX+J-1, DESCX, 1,
     $        X, IX, JX+J, DESCX, 1 )
         IF ( J .GT. 1 ) THEN
            CALL PDAXPY( N, -BETA( J-1 ), X, IX, JX+J-2, DESCX, 1,
     $           X, IX, JX+J, DESCX, 1 )
         END IF
*
*        Normalization.
*
*        beta(j) = sqrt(real(X(:, j+1)**H * X(:, j+1)))
*        X(:, j+1) = X(:, j+1)/beta(j)
*
         CALL PDDOT_ALL( N, DTMP, X, IX, JX+J, DESCX, 1,
     $        X, IX, JX+J, DESCX, 1 )
         BETA( J ) = DSQRT( MAX( ZERO, DTMP ) )
         IF ( BETA( J ) .GT. TOL*TMAX ) THEN
            CALL PDSCAL( N, DSQRT( ONE/DTMP ), X, IX, JX+J, DESCX, 1 )
            TMAX = MAX( TMAX, BETA( J ) )
         ELSE
*
*           Breakdown of Lanczos process.
*
            BETA( J ) = ZERO
            INFO = J
            EXIT
         END IF
      END DO
*
      RETURN
*
*     End of PDLANCZOS().
*
      END
