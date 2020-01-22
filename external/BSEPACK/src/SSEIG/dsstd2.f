      SUBROUTINE DSSTD2( UPLO, N, A, LDA, E, TAU, INFO )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * )
*     ..
*
*  Purpose
*  =======
*
*  DSSTD2 reduces a real skew-symmetric matrix A to skew-symmetric
*  tridiagonal form T by an orthogonal similarity transformation:
*     Q**T * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          skew-symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the skew-symmetric matrix A.  If UPLO = 'U', the
*          leading n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, if UPLO = 'U', the first superdiagonal of A is
*          overwritten by the corresponding elements of the tridiagonal
*          matrix T, and the elements above the first superdiagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of elementary reflectors; if UPLO = 'L', the first
*          subdiagonal of A is overwritten by the corresponding elements
*          of the tridiagonal matrix T, and the elements below the first
*          subdiagonal, with the array TAU, represent the orthogonal
*          matrix Q as a product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The superdiagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U',
*          E(i) = -A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v**T,
*
*  where tau is a real scalar, and v is a real vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v**T,
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  0   e   v2  v3  v4 )              (  0                  )
*    (      0   e   v3  v4 )              (  -e  0              )
*    (          0   e   v4 )              (  v1  -e  0          )
*    (              0   e  )              (  v1  v2  -e  0      )
*    (                  0  )              (  v1  v2  v3  -e  0  )
*
*  where e denotes super-diagonal elements of T, and vi denotes an
*  element of the vector defining H(i).
*
*  =====================================================================
*
*  Level 2 LAPACK-like routine.
*
*  DSSTD2 is modified from DSYTD2 in LAPACK version 3.5.0.
*
*  Written by Meiyue Shao, Lawrence Berkeley National Laboratory.
*  Last change: October 2014
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      DOUBLE PRECISION   TAUI
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFG, DSSMV, DSSR2, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF ( .NOT. UPPER .AND. .NOT. LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF ( N .LT. 0 ) THEN
         INFO = -2
      ELSE IF ( LDA .LT. MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'DSSTD2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF ( N .LE. 1 ) RETURN
*
      DO I = 1, N-1
         TAU( I ) = ZERO
      END DO
      IF ( UPPER ) THEN
*
*        Reduce the upper triangle of A.
*
         DO I = N-1, 1, -1
*
*           Generate elementary reflector H(i) = I - tau * v * v**T
*           to annihilate A(1:i-1,i+1).
*
            CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
            E( I ) = A( I, I+1 )
*
            IF ( TAUI .NE. ZERO ) THEN
*
*              Apply H(i) from both sides to A(1:i,1:i).
*
               A( I, I+1 ) = ONE
*
*              Compute x := tau * A * v, storing x in TAU(1:i).
*
               CALL DSSMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
     $              TAU, 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A + v * x**T - x * v**T.
*
               CALL DSSR2( UPLO, I, ONE, A( 1, I+1 ), 1, TAU, 1, A,
     $              LDA )
*
               A( I, I+1 ) = E( I )
            END IF
            TAU( I ) = TAUI
         END DO
      ELSE
*
*        Reduce the lower triangle of A.
*
         DO I = 1, N-1
*
*           Generate elementary reflector H(i) = I - tau * v * v**T
*           to annihilate A(i+2:n,i).
*
            CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $           TAUI )
            E( I ) = -A( I+1, I )
*
            IF ( TAUI .NE. ZERO ) THEN
*
*              Apply H(i) from both sides to A(i+1:n,i+1:n).
*
               A( I+1, I ) = ONE
*
*              Compute x := tau * A * v, storing y in TAU(i:n-1).
*
               CALL DSSMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
     $              A( I+1, I ), 1, ZERO, TAU( I ), 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A + v * x**T - x * v**T.
*
               CALL DSSR2( UPLO, N-I, ONE, A( I+1, I ), 1, TAU( I ), 1,
     $              A( I+1, I+1 ), LDA )
*
               A( I+1, I ) = -E( I )
            END IF
            TAU( I ) = TAUI
         END DO
      END IF
*
      RETURN
*
*     End of DSSTD2.
*
      END
