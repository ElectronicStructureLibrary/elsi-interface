      SUBROUTINE DSSMM( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA,
     $                  C, LDC )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DSSMM performs one of the matrix-matrix operations
*
*     C := alpha*A*B + beta*C,
*
*  or
*
*     C := alpha*B*A + beta*C,
*
*  where alpha and beta are scalars, A is a skew-symmetric matrix and
*  B and C are m by n matrices.
*
*  Arguments
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether the skew-symmetric matrix A
*           appears on the left or right in the operation as follows:
*
*              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
*
*              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the skew-symmetric matrix A is to be
*           referenced as follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of the
*                                  skew-symmetric matrix is to be
*                                  referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of the
*                                  skew-symmetric matrix is to be
*                                  referenced.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix C.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix C.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           m when SIDE = 'L' or 'l' and is n otherwise.
*           Before entry with SIDE = 'L' or 'l', the m by m part of
*           the array A must contain the skew-symmetric matrix, such
*           that when UPLO = 'U' or 'u', the leading m by m upper
*           triangular part of the array A must contain the upper
*           triangular part of the skew-symmetric matrix and the
*           strictly lower triangular part of A is not referenced, and
*           when UPLO = 'L' or 'l', the leading m by m lower triangular
*           part of the array A must contain the lower triangular part
*           of the skew-symmetric matrix and the strictly upper
*           triangular part of A is not referenced.
*           Before entry with SIDE = 'R' or 'r', the n by n part of
*           the array A must contain the skew-symmetric matrix, such
*           that when UPLO = 'U' or 'u', the leading n by n upper
*           triangular part of the array A must contain the upper
*           triangular part of the skew-symmetric matrix and the
*           strictly lower triangular part of A is not referenced, and
*           when UPLO = 'L' or 'l', the leading n by n lower triangular
*           part of the array A must contain the lower triangular part
*           of the skew-symmetric matrix and the strictly upper
*           triangular part of A is not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When SIDE = 'L' or 'l' then
*           LDA must be at least max( 1, m ), otherwise LDA must be at
*           least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry, the leading m by n part of the array B must
*           contain the matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           inthe calling (sub) program. LDB must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading m by n part of the array C must
*           contain the matrix C, except when beta is zero, in which
*           case C need not be set on entry.
*           On exit, the array C is overwritten by the m by n updated
*           matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in the calling (sub) program. LDC must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  Further Details
*  ===============
*
*  Level 3 BLAS-like routine.
*
*  Written by Meiyue Shao, Lawrence Berkeley National Laboratory.
*  Last change: October 2014
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      INTEGER            NB
      PARAMETER          ( NB = 64 )
*
*     .. Local Scalars ..
      INTEGER            II, JJ, KK, IC, JC, KC, NROWA, INFO
      LOGICAL            LEFT, UPPER
      DOUBLE PRECISION   TEMP
*     .. Local Arrays ..
      DOUBLE PRECISION   WORK( M, N )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSSMM_SM, DGEMM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Set NROWA as the number of rows of A.
*
      LEFT = LSAME( SIDE, 'L' )
      IF ( LEFT ) THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      UPPER = LSAME( UPLO, 'U' )
*
*     Test the input parameters.
*
      INFO = 0
      IF ( .NOT. LEFT .AND. .NOT. LSAME( SIDE, 'R' ) ) THEN
         INFO = 1
      ELSE IF ( .NOT. UPPER .AND. .NOT. LSAME( UPLO, 'L' ) ) THEN
         INFO = 2
      ELSE IF ( M .LT. 0 ) THEN
         INFO = 3
      ELSE IF ( N .LT. 0 ) THEN
         INFO = 4
      ELSE IF ( LDA .LT. MAX( 1, NROWA ) ) THEN
         INFO = 7
      ELSE IF ( LDB .LT. MAX( 1, M ) ) THEN
         INFO = 9
      ELSE IF ( LDC .LT. MAX( 1, M ) ) THEN
         INFO = 12
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'DSSMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF ( ( M .EQ. 0 ) .OR. ( N .EQ. 0 ) .OR.
     $     ( ( ALPHA .EQ. ZERO ) .AND. ( BETA .EQ. ONE ) ) ) RETURN
*
*     And when alpha equals zero.
*
      IF ( ALPHA .EQ. ZERO ) THEN
         IF ( BETA .EQ. ZERO) THEN
            DO JJ = 1, N
               DO II = 1, M
                  C( II, JJ ) = ZERO
               END DO
            END DO
         ELSE
            DO JJ = 1, N
               DO II = 1, M
                  C( II, JJ ) = BETA * C( II, JJ )
               END DO
            END DO
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF ( LEFT ) THEN
*
*        Form  C := alpha*A*B + beta*C.
*
         DO KK = 1, N, NB
            KC = MIN( NB, N-KK+1 )
            TEMP = BETA
            DO JJ = 1, M, NB
               JC = MIN( NB, M-JJ+1 )
               DO II = 1, M, NB
                  IC = MIN( NB, M-II+1 )
*
*                 Call DSSMM_SM for diagonal blocks,
*                 and DGEMM for off-diagonal blocks.
*
                  IF ( UPPER ) THEN
                     IF ( II .GT. JJ ) THEN
                        CALL DGEMM( 'T', 'N', IC, KC, JC, -ALPHA,
     $                       A( JJ, II ), LDA, B( JJ, KK ), LDB, TEMP,
     $                       C( II, KK ), LDC )
                     ELSE IF ( II .LT. JJ ) THEN
                        CALL DGEMM( 'N', 'N', IC, KC, JC, ALPHA,
     $                       A( II, JJ ), LDA, B( JJ, KK ), LDB, TEMP,
     $                       C( II, KK ), LDC )
                     ELSE
                        CALL DSSMM_SM( LEFT, UPPER, JC, KC, ALPHA,
     $                       A( JJ, JJ ), LDA, B( JJ, KK ), LDB, TEMP,
     $                       C( JJ, KK ), LDC, WORK )
                     END IF
                  ELSE
                     IF ( II .LT. JJ ) THEN
                        CALL DGEMM( 'T', 'N', IC, KC, JC, -ALPHA,
     $                       A( JJ, II ), LDA, B( JJ, KK ), LDB, TEMP,
     $                       C( II, KK ), LDC )
                     ELSE IF ( II .GT. JJ ) THEN
                        CALL DGEMM( 'N', 'N', IC, KC, JC, ALPHA,
     $                       A( II, JJ ), LDA, B( JJ, KK ), LDB, TEMP,
     $                       C( II, KK ), LDC )
                     ELSE
                        CALL DSSMM_SM( LEFT, UPPER, JC, KC, ALPHA,
     $                       A( JJ, JJ ), LDA, B( JJ, KK ), LDB, TEMP,
     $                       C( JJ, KK ), LDC, WORK )
                     END IF
                  END IF
               END DO
               TEMP = ONE
            END DO
         END DO
      ELSE
*
*        Form  C := alpha*B*A + beta*C.
*
         DO KK = 1, M, NB
            KC = MIN( NB, M-KK+1 )
            DO JJ = 1, N, NB
               JC = MIN( NB, N-JJ+1 )
               TEMP = BETA
               DO II = 1, N, NB
                  IC = MIN( NB, N-II+1 )
*
*                 Call DSSMM_SM for diagonal blocks,
*                 and DGEMM for off-diagonal blocks.
*
                  IF ( UPPER ) THEN
                     IF ( II .GT. JJ ) THEN
                        CALL DGEMM( 'N', 'T', KC, JC, IC, -ALPHA,
     $                       B( KK, II ), LDB, A( JJ, II ), LDA, TEMP,
     $                       C( KK, JJ ), LDC )
                     ELSE IF ( II .LT. JJ ) THEN
                        CALL DGEMM( 'N', 'N', KC, JC, IC, ALPHA,
     $                       B( KK, II ), LDB, A( II, JJ ), LDA, TEMP,
     $                       C( KK, JJ ), LDC )
                     ELSE
                        CALL DSSMM_SM( LEFT, UPPER, KC, JC, ALPHA,
     $                       A( JJ, JJ ), LDA, B( KK, JJ ), LDB, TEMP,
     $                       C( KK, JJ ), LDC, WORK )
                     END IF
                  ELSE
                     IF ( II .LT. JJ ) THEN
                        CALL DGEMM( 'N', 'T', KC, JC, IC, -ALPHA,
     $                       B( KK, II ), LDB, A( JJ, II ), LDA, TEMP,
     $                       C( KK, JJ ), LDC )
                     ELSE IF ( II .GT. JJ ) THEN
                        CALL DGEMM( 'N', 'N', KC, JC, IC, ALPHA,
     $                       B( KK, II ), LDB, A( II, JJ ), LDA, TEMP,
     $                       C( KK, JJ ), LDC )
                     ELSE
                        CALL DSSMM_SM( LEFT, UPPER, KC, JC, ALPHA,
     $                       A( JJ, JJ ), LDA, B( KK, JJ ), LDB, TEMP,
     $                       C( KK, JJ ), LDC, WORK )
                     END IF
                  END IF
                  TEMP = ONE
               END DO
            END DO
         END DO
      END IF
*
      RETURN
*
*     End of DSSMM.
*
      END
*
*
*
      SUBROUTINE DSSMM_SM( LEFT, UPPER, M, N, ALPHA, A, LDA, B, LDB,
     $                     BETA, C, LDC, WORK )
*
*     The dimension of WORK is at least (m,n).
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      LOGICAL            LEFT, UPPER
      INTEGER            M, N, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   WORK( M, * )
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
*     .. Local Scalars ..
      INTEGER            I, J, NROWA
      CHARACTER          SIDE, UPLO
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACPY, DTRMM
*     ..
*     .. Executable Statements ..
*
*     Set NROWA as the number of rows of A.
*
      IF ( LEFT ) THEN
         SIDE = 'L'
         NROWA = M
      ELSE
         SIDE = 'R'
         NROWA = N
      END IF
      IF ( UPPER ) THEN
         UPLO = 'U'
      ELSE
         UPLO = 'L'
      END IF
*
*     Avoid unexpected behavior when C contains NaN entries.
*
      IF ( BETA .EQ. ZERO ) THEN
         DO J = 1, N
            DO I = 1, M
               C( I, J ) = ZERO
            END DO
         END DO
      END IF
*
*     Quick return when alpha equals zero.
*
      IF ( ALPHA .EQ. ZERO ) THEN
         IF ( BETA .NE. ZERO ) THEN
            DO J = 1, N
               DO I = 1, M
                  C( I, J ) = BETA * C( I, J )
               END DO
            END DO
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      CALL DLACPY( 'A', M, N, B, LDB, WORK, M )
      CALL DTRMM( SIDE, UPLO, 'N', 'N', M, N, ALPHA, A, LDA, WORK, M )
      DO J = 1, N
         DO I = 1, M
            C( I, J ) = BETA * C( I, J ) + WORK( I, J )
         END DO
      END DO
      CALL DLACPY( 'A', M, N, B, LDB, WORK, M )
      CALL DTRMM( SIDE, UPLO, 'T', 'N', M, N, ALPHA, A, LDA, WORK, M )
      DO J = 1, N
         DO I = 1, M
            C( I, J ) = C( I, J ) - WORK( I, J )
         END DO
      END DO
*
      RETURN
*
*     End of DSSMM_SM.
*
      END
