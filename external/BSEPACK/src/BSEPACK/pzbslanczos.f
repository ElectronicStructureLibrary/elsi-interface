      SUBROUTINE PZBSLANCZOS( N, STEPS, A, IA, JA, DESCA, B, IB, JB,
     $                        DESCB, X, IX, JX, DESCX, ALPHA, BETA, SCL,
     $                        INFO )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, STEPS, IA, JA, IB, JB, IX, JX, INFO
      DOUBLE PRECISION   SCL
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( * ), B( * ), X( * )
      DOUBLE PRECISION   ALPHA( * ), BETA( * )
      INTEGER            DESCA( * ), DESCB( * ), DESCX( * )
*     ..
*
*  Purpose
*  =======
*
*  PZBSLANCZOS() computes the structure preserving Lanczos decomposition
*
*     H**2 * [ U; conj(U) ] = [ U; conj(U) ] * T + [rank 1],
*
*  where
*
*     H = [       A,        B;
*          -conj(B), -conj(A) ],
*
*  is the Behte--Salpeter Hamiltonian with
*
*     Omega = [      A,       B;
*              conj(B), conj(A) ]
*
*  positive definite, and [ U; conj(U) ] satisfies the orthogonal
*  property
*
*     [ U; conj(U) ]**H * [ V; conj(V) ] = 2 * I,
*
*  where
*
*     [ V; conj(V) ] = Omega * [ U; conj(U) ].
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
*          2*N is the number of rows and columns of H.
*          N >= 0.
*
*  STEPS   (global input) INTEGER
*          The number of Lanczos steps.
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
*  X       (local input/output) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_X, LOCc(JX+1)).
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
*  ALPHA   (global output) DOUBLE PRECISION array, dimension (STEPS)
*          The diagonal of the tridiagonal matrix T.
*
*  BETA    (global output) DOUBLE PRECISION array, dimension (STEPS)
*          The subdiagonal of the tridiagonal matrix T.
*
*  SCL     (global output) DOUBLE PRECISION
*          SCL = [ u; conj(u) ]**H * Omega * [ u; conj(u) ],
*          where u is the unnormalized starting vector on entry.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          > 0:  Lanczos process breaks down after INFO steps.
*
*  Alignment requirements
*  ======================
*
*  The distributed matrices A, B, and X must satisfy the following
*  properties:
*
*     DESCA( MB_ ) = DESCA( NB_ )
*     DESCB( MB_ ) = DESCB( NB_ )
*     DESCA( MB_ ) = DESCB( MB_ ) = DESCX( MB_ )
*     DESCA( NB_ ) = DESCB( NB_ ) = DESCX( NB_ )
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
      COMPLEX*16         CPLX_ZERO, CPLX_ONE
      PARAMETER          ( CPLX_ZERO = ( 0.0D+0, 0.0D+0 ),
     $                     CPLX_ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            ICTXT, J
      DOUBLE PRECISION   DTMP, TMAX, TOL
      COMPLEX*16         ZTMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DSQRT, MAX
*     ..
*     .. External Functions ..
      EXTERNAL           PDLAMCH
      DOUBLE PRECISION   PDLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           PZAXPY, PZDSCAL, PZHEMM, PZSYMM, PZDOTC_ALL,
     $                   PZLACONJ
*
      INFO = 0
      ICTXT = DESCX( CTXT_ )
      TOL = 10*N*PDLAMCH( ICTXT, 'E' )
      TMAX = ZERO
*
*     V(:, 1) = A*U(:, 1) + B*conj(U(:, 1))
*     scl = real(U(:, 1)**H * V(:, 1))
*     U(:, 1) = U(:, 1)/sqrt(scl)
*     V(:, 1) = V(:, 1)/sqrt(scl)
*
      CALL PZHEMM( 'L', 'L', N, 1, CPLX_ONE, A, IA, JA, DESCA,
     $     X, IX, JX, DESCX, CPLX_ZERO, X, IX+N, JX, DESCX )
      CALL PZLACONJ( N, 1, X, IX, JX, DESCX )
      CALL PZSYMM( 'L', 'L', N, 1, CPLX_ONE, B, IB, JB, DESCB,
     $     X, IX, JX, DESCX, CPLX_ONE, X, IX+N, JX, DESCX )
      CALL PZLACONJ( N, 1, X, IX, JX, DESCX )
      CALL PZDOTC_ALL( N, ZTMP, X, IX, JX, DESCX, 1,
     $     X, IX+N, JX, DESCX, 1 )
      SCL = DBLE( ZTMP )
      CALL PZDSCAL( 2*N, DSQRT( ONE/SCL ), X, IX, JX, DESCX, 1 )
*
*     Lanczos process.
*
      DO J = 1, STEPS
*
*        First matrix vector multiplication.
*
*        U(:, j+1) = A*V(:, j) - B*conj(V(:, j))
*
         CALL PZHEMM( 'L', 'L', N, 1, CPLX_ONE, A, IA, JA, DESCA,
     $        X, IX+N, JX+J-1, DESCX, CPLX_ZERO, X, IX, JX+J, DESCX )
         CALL PZLACONJ( N, 1, X, IX+N, JX+J-1, DESCX )
         CALL PZSYMM( 'L', 'L', N, 1, -CPLX_ONE, B, IB, JB, DESCB,
     $        X, IX+N, JX+J-1, DESCX, CPLX_ONE, X, IX, JX+J, DESCX )
         CALL PZLACONJ( N, 1, X, IX+N, JX+J-1, DESCX )
*
*        Gram--Schmidt.
*
         IF ( J .GT. 1 ) THEN
*
*           U(:, j+1) = U(:, j+1) - beta(j-1)*U(:, j-1)
*
            CALL PZAXPY( N, DCMPLX( -BETA( J-1 ), ZERO ),
     $           X, IX, JX+J-2, DESCX, 1, X, IX, JX+J, DESCX, 1 )
         END IF
*
*        alpha(j) = real(V(:, j)**H * U(:, j+1))
*        U(:, j+1) = U(:, j+1) - alpha(j)*U(:, j)
*
         CALL PZDOTC_ALL( N, ZTMP, X, IX+N, JX+J-1, DESCX, 1,
     $        X, IX, JX+J, DESCX, 1 )
         ALPHA( J ) = DBLE( ZTMP )
         TMAX = MAX( TMAX, ALPHA( J ) )
         CALL PZAXPY( N, DCMPLX( -ALPHA( J ), ZERO ),
     $        X, IX, JX+J-1, DESCX, 1, X, IX, JX+J, DESCX, 1 )
*
*        Second matrix vector multiplication.
*
*        V(:, j+1) = A*U(:, j+1) + B*conj(U(:, j+1))
*
         CALL PZHEMM( 'L', 'L', N, 1, CPLX_ONE, A, IA, JA, DESCA,
     $        X, IX, JX+J, DESCX, CPLX_ZERO, X, IX+N, JX+J, DESCX )
         CALL PZLACONJ( N, 1, X, IX, JX+J, DESCX )
         CALL PZSYMM( 'L', 'L', N, 1, CPLX_ONE, B, IB, JB, DESCB,
     $        X, IX, JX+J, DESCX, CPLX_ONE, X, IX+N, JX+J, DESCX )
         CALL PZLACONJ( N, 1, X, IX, JX+J, DESCX )
*
*        Normalization.
*
*        beta(j) = sqrt(real(U(:, j+1)**H * V(:, j+1)))
*        U(:, j+1) = U(:, j+1)/beta(j)
*        V(:, j+1) = V(:, j+1)/beta(j)
*
         CALL PZDOTC_ALL( N, ZTMP, X, IX, JX+J, DESCX, 1,
     $        X, IX+N, JX+J, DESCX, 1 )
         DTMP = MAX( ZERO, DBLE( ZTMP ) )
         BETA( J ) = DSQRT( DTMP )
         IF ( BETA( J ) .GT. TOL*TMAX ) THEN
            CALL PZDSCAL( 2*N, DSQRT( ONE/DTMP ), X, IX, JX+J, DESCX,
     $           1 )
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
*     End of PZBSLANCZOS().
*
      END
