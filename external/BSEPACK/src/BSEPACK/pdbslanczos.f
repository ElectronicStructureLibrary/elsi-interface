      SUBROUTINE PDBSLANCZOS( N, STEPS, M, IM, JM, DESCM, K, IK, JK,
     $                        DESCK, X, IX, JX, DESCX, ALPHA, BETA, SCL,
     $                        INFO )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, STEPS, IM, JM, IK, JK, IX, JX, INFO
      DOUBLE PRECISION   SCL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   M( * ), K( * ), X( * ), ALPHA( * ), BETA( * )
      INTEGER            DESCM( * ), DESCK( * ), DESCX( * )
*     ..
*
*  Purpose
*  =======
*
*  PDBSLANCZOS() computes the structure preserving Lanczos decomposition
*
*     K * M * U = U * T + [rank 1],
*
*  where both M and K are positive definite, and U satisfies the
*  orthogonal property
*
*     U**T * V = I,
*
*  where
*
*     V = M * U.
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
*  Such a global array has an associated description vector DESCM.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCM( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCM( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCM( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCM( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCM( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCM( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCM( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCM( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCM( LLD_ )  The leading dimension of the local
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
*          The number of rows and columns of M and K.
*          N >= 0.
*
*  STEPS   (global input) INTEGER
*          The number of Lanczos steps.
*
*  M       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_M, LOCc(JM+N-1)).
*          This array contains the local pieces of the N-by-N Hermitian
*          distributed matrix sub( M ). The leading N-by-N lower
*          triangular part of sub( M ) contains the lower triangular
*          part of the distributed matrix, and its strictly upper
*          triangular part is not referenced.
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
*  K       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_K, LOCc(JK+N-1)).
*          This array contains the local pieces of the N-by-N symmetric
*          distributed matrix sub( K ). The leading N-by-N lower
*          triangular part of sub( K ) contains the lower triangular
*          part of the distributed matrix, and its strictly upper
*          triangular part is not referenced.
*
*  IK      (global input) INTEGER
*          The row index in the global array K indicating the first
*          row of sub( K ).
*
*  JK      (global input) INTEGER
*          The column index in the global array B indicating the
*          first column of sub( K ).
*
*  DESCK   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix K.
*
*  X       (local input/output) DOUBLE PRECISION pointer into the local
*          memory to an array of dimension (LLD_X, LOCc(JX+1)).
*          On entry, X( 1:N, 1 ) is used to store the starting vector u,
*          which does not need to be normalized.
*          On exit, X( 1:N, j ) is used to store U( :, j ), while
*          X( N+1:2*N, j ) is used to store V( :, j ).
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
*  ALPHA   (global output) DOUBLE PRECISION array, dimension (K)
*          The diagonal of the tridiagonal matrix T.
*
*  BETA    (global output) DOUBLE PRECISION array, dimension (K)
*          The subdiagonal of the tridiagonal matrix T.
*
*  SCL     (global output) DOUBLE PRECISION
*          SCL = u**T * M * u,
*          where u is the unnormalized starting vector on entry.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          > 0:  Lanczos process breaks down after INFO steps.
*
*  Alignment requirements
*  ======================
*
*  The distributed matrices M, K, and X must satisfy the following
*  properties:
*
*     DESCM( MB_ ) = DESCM( NB_ )
*     DESCK( MB_ ) = DESCK( NB_ )
*     DESCM( MB_ ) = DESCK( MB_ ) = DESCX( MB_ )
*     DESCM( NB_ ) = DESCK( NB_ ) = DESCX( NB_ )
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
*     V(:, 1) = M*U(:, 1)
*     scl = U(:, 1)**T * V(:, 1)
*     U(:, 1) = U(:, 1)/sqrt(scl)
*     V(:, 1) = V(:, 1)/sqrt(scl)
*
      CALL PDSYMM( 'L', 'L', N, 1, ONE, M, IM, JM, DESCM,
     $     X, IX, JX, DESCX, ZERO, X, IX+N, JX, DESCX )
      CALL PDDOT_ALL( N, SCL, X, IX, JX, DESCX, 1,
     $     X, IX+N, JX, DESCX, 1 )
      CALL PDSCAL( 2*N, DSQRT( ONE/SCL ), X, IX, JX, DESCX, 1 )
*
*     Lanczos process.
*
      DO J = 1, STEPS
*
*        First matrix vector multiplication.
*
*        U(:, j+1) = K*V(:, j)
*
         CALL PDSYMM( 'L', 'L', N, 1, ONE, K, IK, JK, DESCK,
     $        X, IX+N, JX+J-1, DESCX, ZERO, X, IX, JX+J, DESCX )
*
*        Gram--Schmidt.
*
         IF ( J .GT. 1 ) THEN
*
*           U(:, j+1) = U(:, j+1) - beta(j-1)*U(:, j-1)
*
            CALL PDAXPY( N, -BETA( J-1 ), X, IX, JX+J-2, DESCX, 1,
     $           X, IX, JX+J, DESCX, 1 )
         END IF
*
*        alpha(j) = V(:, j)**T * U(:, j+1)
*        U(:, j+1) = U(:, j+1) - alpha(j)*U(:, j)
*
         CALL PDDOT_ALL( N, ALPHA( J ), X, IX+N, JX+J-1, DESCX, 1,
     $        X, IX, JX+J, DESCX, 1 )
         TMAX = MAX( TMAX, ALPHA( J ) )
         CALL PDAXPY( N, -ALPHA( J ), X, IX, JX+J-1, DESCX, 1,
     $        X, IX, JX+J, DESCX, 1 )
*
*        Second matrix vector multiplication.
*
*        V(:, j+1) = M*U(:, j+1)
*
         CALL PDSYMM( 'L', 'L', N, 1, ONE, M, IM, JM, DESCM,
     $        X, IX, JX+J, DESCX, ZERO, X, IX+N, JX+J, DESCX )
*
*        Normalization.
*
*        beta(j) = sqrt(U(:, j+1)**T * V(:, j+1))
*        U(:, j+1) = U(:, j+1)/beta(j)
*        V(:, j+1) = V(:, j+1)/beta(j)
*
         CALL PDDOT_ALL( N, DTMP, X, IX, JX+J, DESCX, 1,
     $        X, IX+N, JX+J, DESCX, 1 )
         BETA( J ) = DSQRT( MAX( ZERO, DTMP ) )
         IF ( BETA( J ) .GT. TOL*TMAX ) THEN
            CALL PDSCAL( 2*N, DSQRT( ONE/DTMP ), X, IX, JX+J, DESCX, 1 )
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
*     End of PDBSLANCZOS().
*
      END
