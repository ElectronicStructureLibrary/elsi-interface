      SUBROUTINE DLASTD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASTD reduces NB rows and columns of a real skew-symmetric matrix A
*  to skew-symmetric tridiagonal form by an orthogonal similarity
*  transformation Q**T * A * Q, and returns the matrices V and W which
*  are needed to apply the transformation to the unreduced part of A.
*
*  If UPLO = 'U', DLASTD reduces the last NB rows and columns of a
*  matrix, of which the upper triangle is supplied;
*  if UPLO = 'L', DLASTD reduces the first NB rows and columns of a
*  matrix, of which the lower triangle is supplied.
*
*  This is an auxiliary routine called by DSSTRD.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          skew-symmetric matrix A is stored:
*          = 'U': Upper triangular
*          = 'L': Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.
*
*  NB      (input) INTEGER
*          The number of rows and columns to be reduced.
*
*  A       (input/output) DOUBLE PRECISION   array, dimension (LDA,N)
*          On entry, the skew-symmetric matrix A.
*          If UPLO = 'U', the leading n-by-n upper triangular part of A
*          contains the upper triangular part of the matrix A, and the
*          strictly lower triangular part of A is not referenced.
*          If UPLO = 'L', the leading n-by-n lower triangular part of A
*          contains the lower triangular part of the matrix A, and the
*          strictly upper triangular part of A is not referenced.
*          On exit:
*          if UPLO = 'U', the last NB columns have been reduced to
*            tridiagonal form; the elements above the diagonal with the
*            array TAU, represent the orthogonal matrix Q as a product
*            of elementary reflectors;
*          if UPLO = 'L', the first NB columns have been reduced to
*            tridiagonal form; the elements below the diagonal with the
*            array TAU, represent the  orthogonal matrix Q as a product
*            of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= (1,N).
*
*  E       (output) DOUBLE PRECISION   array, dimension (N-1)
*          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
*          elements of the last NB columns of the reduced matrix;
*          if UPLO = 'L', E(1:nb) contains the superdiagonal elements of
*          the first NB columns of the reduced matrix.
*
*  TAU     (output) DOUBLE PRECISION   array, dimension (N-1)
*          The scalar factors of the elementary reflectors, stored in
*          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
*          See Further Details.
*
*  W       (output) DOUBLE PRECISION   array, dimension (LDW,NB)
*          The n-by-nb matrix W required to update the unreduced part
*          of A.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array W. LDW >= max(1,N).
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n) H(n-1) . . . H(n-nb+1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v**T
*
*  where tau is a real scalar, and v is a real vector with
*  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
*  and tau in TAU(i-1).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(nb).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v**T
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
*  and tau in TAU(i).
*
*  The elements of the vectors v together form the n-by-nb matrix V
*  which is needed, with W, to apply the transformation to the unreduced
*  part of the matrix, using a skew-symmetric rank-2k update of the
*  form:
*     A := A + V*W**T - W*V**T.
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5 and nb = 2:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  0   a   a   v4  v5 )              (  0                  )
*    (      0   a   v4  v5 )              (  1   0              )
*    (          0   1   v5 )              (  v1  1   0          )
*    (              0   1  )              (  v1  v2  a   0      )
*    (                  0  )              (  v1  v2  a   a   0  )
*
*  where a denotes an element of the original matrix that is unchanged,
*  and vi denotes an element of the vector defining H(i).
*
*  =====================================================================
*
*  LAPACK-like auxiliary routine.
*
*  DLASTD is modified from DLATRD in LAPACK version 3.5.0.
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
      INTEGER            I, IW
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DLARFG, DSCAL, DSSMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      IF ( N .LE. 0 )
     $   RETURN
*
      IF ( LSAME( UPLO, 'U' ) ) THEN
*
*        Reduce last NB columns of upper triangle.
*
         DO I = N, N-NB+1, -1
            IW = I-N+NB
            IF ( I .GT. 1 ) THEN
*
*              Update A(1:i-1,i).
*
               IF ( I .LT. N ) THEN
                  CALL DGEMV( 'N', I-1, N-I, ONE, A( 1, I+1 ), LDA,
     $                 W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
                  CALL DGEMV( 'N', I-1, N-I, -ONE, W( 1, IW+1 ), LDW,
     $                 A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
               END IF
*
*              Generate elementary reflector H(i) to annihilate
*              A(1:i-2,i).
*
               CALL DLARFG( I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = A( I-1, I )
               A( I-1, I ) = ONE
*
*              Compute W(1:i-1,i).
*
               CALL DSSMV( 'U', I-1, ONE, A, LDA, A( 1, I ), 1, ZERO,
     $              W( 1, IW ), 1 )
               IF ( I .LT. N ) THEN
                  CALL DGEMV( 'T', I-1, N-I, ONE, W( 1, IW+1 ), LDW,
     $                 A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'N', I-1, N-I, ONE, A( 1, I+1 ), LDA,
     $                 W( I+1, IW ), 1, ONE, W( 1, IW ), 1 )
                  CALL DGEMV( 'T', I-1, N-I, ONE, A( 1, I+1 ), LDA,
     $                 A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'N', I-1, N-I, -ONE, W( 1, IW+1 ), LDW,
     $                 W( I+1, IW ), 1, ONE, W( 1, IW ), 1 )
               END IF
               CALL DSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
            END IF
*
         END DO
*
      ELSE
*
*        Reduce first NB columns of lower triangle.
*
         DO I = 1, NB
            IF ( I .LT. N ) THEN
*
*              Update A(i+1:n,i).
*
               IF ( I .GT. 1 ) THEN
                  CALL DGEMV( 'N', N-I, I-1, ONE, A( I+1, 1 ), LDA,
     $                 W( I, 1 ), LDW, ONE, A( I+1, I ), 1 )
                  CALL DGEMV( 'N', N-I, I-1, -ONE, W( I+1, 1 ), LDW,
     $                 A( I, 1 ), LDA, ONE, A( I+1, I ), 1 )
               END IF
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:n,i).
*
               CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $              TAU( I ) )
               E( I ) = -A( I+1, I )
               A( I+1, I ) = ONE
*
*              Compute W(i+1:n,i).
*
               CALL DSSMV( 'L', N-I, ONE, A( I+1, I+1 ), LDA,
     $              A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL DGEMV( 'T', N-I, I-1, ONE, W( I+1, 1 ), LDW,
     $              A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'N', N-I, I-1, ONE, A( I+1, 1 ), LDA,
     $              W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DGEMV( 'T', N-I, I-1, ONE, A( I+1, 1 ), LDA,
     $              A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'N', N-I, I-1, -ONE, W( I+1, 1 ), LDW,
     $              W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
            END IF
*
         END DO
*
      END IF
*
      RETURN
*
*     End of DLASTD.
*
      END
