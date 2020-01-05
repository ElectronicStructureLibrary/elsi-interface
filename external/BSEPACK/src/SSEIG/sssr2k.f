      SUBROUTINE SSSR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA,
     $                   C, LDC )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      REAL               ALPHA, BETA
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  SSSR2K performs one of the skew-symmetric rank 2k operations
*
*     C := alpha*A*B**T - alpha*B*A**T + beta*C,
*
*  or
*
*     C := alpha*A**T*B - alpha*B**T*A + beta*C,
*
*  where alpha and beta are scalars, C is an n by n skew-symmetric
*  matrix and A and B are n by k matrices in the first case and k by n
*  matrices in the second case.
*
*  Arguments
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array C is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   C := alpha*A*B**T - alpha*B*A**T +
*                                        beta*C.
*
*              TRANS = 'T' or 't'   C := alpha*A**T*B - alpha*B**T*A +
*                                        beta*C.
*
*              TRANS = 'C' or 'c'   C := alpha*A**T*B - alpha*B**T*A +
*                                        beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with TRANS = 'N' or 'n', K specifies the number
*           of columns of the matrices A and B, and on entry with
*           TRANS = 'T' or 't', K specifies the number of rows of the
*           matrices A and B. K must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL array of DIMENSION ( LDA, ka ), where ka is
*           k when TRANS = 'N' or 'n', and is n otherwise.
*           Before entry with TRANS = 'N' or 'n', the leading n by k
*           part of the array A must contain the matrix A, otherwise
*           the leading k by n part of the array A must contain the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When TRANS = 'N' or 'n'
*           then LDA must be at least max( 1, n ), otherwise LDA must
*           be at least max( 1, k ).
*           Unchanged on exit.
*
*  B      - REAL array of DIMENSION ( LDB, kb ), where kb is
*           k when TRANS = 'N' or 'n', and is n otherwise.
*           Before entry with TRANS = 'N' or 'n', the leading n by k
*           part of the array B must contain the matrix B, otherwise
*           the leading k by n part of the array B must contain the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When TRANS = 'N' or 'n'
*           then LDB must be at least max( 1, n ), otherwise LDB must
*           be at least max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - REAL.
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - REAL array of DIMENSION ( LDC, n ).
*           Before entry with UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array C must contain the upper
*           triangular part of the skew-symmetric matrix and the
*           strictly lower triangular part of C is not referenced.
*           On exit, the upper triangular part of the array C is
*           overwritten by the upper triangular part of the updated
*           matrix.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array C must contain the lower
*           triangular part of the skew-symmetric matrix and the
*           strictly upper triangular part of C is not referenced.
*           On exit, the lower triangular part of the array C is
*           overwritten by the lower triangular part of the updated
*           matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in the calling (sub) program. LDC must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  Further Details
*  ===============
*
*  Level 3 BLAS-like routine.
*
*  SSSR2K_SM is modified from SSYR2K in reference BLAS version 3.5.0.
*
*  Written by Meiyue Shao, Lawrence Berkeley National Laboratory.
*  Last change: October 2014
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      INTEGER            NB
      PARAMETER          ( NB = 64 )
*     ..
*     .. Local Scalars ..
      INTEGER            II, JJ, KK, IC, JC, KC, NROWA, INFO
      LOGICAL            NOTRANS, UPPER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SSSR2K_SM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*
*     Test the input parameters.
*
      NOTRANS = LSAME( TRANS, 'N' )
      UPPER = LSAME( UPLO, 'U' )
      IF ( NOTRANS ) THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
*
      INFO = 0
      IF ( .NOT. UPPER .AND. .NOT. LSAME( UPLO, 'L' ) ) THEN
         INFO = 1
      ELSE IF ( .NOT. LSAME( TRANS, 'N' ) .AND.
     $     .NOT. LSAME( TRANS, 'T' ) .AND.
     $     .NOT. LSAME( TRANS, 'C' ) ) THEN
         INFO = 2
      ELSE IF ( N .LT. 0 ) THEN
         INFO = 3
      ELSE IF ( K .LT. 0 ) THEN
         INFO = 4
      ELSE IF ( LDA .LT. MAX( 1, NROWA ) ) THEN
         INFO = 7
      ELSE IF ( LDB .LT. MAX( 1, NROWA ) ) THEN
         INFO = 9
      ELSE IF ( LDC .LT. MAX( 1, N ) ) THEN
         INFO = 12
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'SSSR2K', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF ( ( N .EQ. 0 ) .OR. ( ( ( ALPHA .EQ. ZERO ) .OR.
     &     ( K .EQ. 0 ) ) .AND. ( BETA .EQ. ONE ) ) ) RETURN
*
*     Form  C := beta*C.
*
      IF ( UPPER ) THEN
         IF ( BETA .EQ. ZERO ) THEN
            DO JJ = 1, N
               DO II = 1, JJ
                  C( II, JJ ) = ZERO
               END DO
            END DO
         ELSE
            DO JJ = 1, N
               DO II = 1, JJ
                  C( II, JJ ) = BETA*C( II, JJ )
               END DO
            END DO
         END IF
      ELSE
         IF ( BETA .EQ. ZERO ) THEN
            DO JJ = 1, N
               DO II = JJ, N
                  C( II, JJ ) = ZERO
               END DO
            END DO
         ELSE
            DO JJ = 1, N
               DO II = JJ, N
                  C( II, JJ ) = BETA*C( II, JJ )
               END DO
            END DO
         END IF
      END IF
*
*     Quick return when alpha equals zero.
*
      IF ( ALPHA .EQ. ZERO ) RETURN
*
*     Start the operations.
*
      IF ( NOTRANS ) THEN
*
*        Form  C := alpha*A*B**T - alpha*B*A**T + C.
*
         IF ( UPPER ) THEN
            DO KK = 1, K, NB
               KC = MIN( NB, K-KK+1 )
               DO JJ = 1, N, NB
                  JC = MIN( NB, N-JJ+1 )
*
*                 Call SSSR2K_SM for diagonal blocks,
*                 and SGEMM for off-diagonal blocks.
*
                  DO II = 1, JJ-NB, NB
                     CALL SGEMM( 'N', 'T', NB, JC, KC, ALPHA,
     $                    A( II, KK ), LDA, B( JJ, KK ), LDB, ONE,
     $                    C( II, JJ ), LDC )
                     CALL SGEMM( 'N', 'T', NB, JC, KC, -ALPHA,
     $                    B( II, KK ), LDB, A( JJ, KK ), LDA, ONE,
     $                    C( II, JJ ), LDC )
                  END DO
                  CALL SSSR2K_SM( UPPER, NOTRANS, JC, KC, ALPHA,
     $                 A( JJ, KK ), LDA, B( JJ, KK ), LDB, ONE,
     $                 C( JJ, JJ ), LDC )
               END DO
            END DO
         ELSE
            DO KK = 1, K, NB
               KC = MIN( NB, K-KK+1 )
               DO JJ = 1, N, NB
                  JC = MIN( NB, N-JJ+1 )
*
*                 Call SSSR2K_SM for diagonal blocks,
*                 and SGEMM for off-diagonal blocks.
*
                  CALL SSSR2K_SM( UPPER, NOTRANS, JC, KC, ALPHA,
     $                 A( JJ, KK ), LDA, B( JJ, KK ), LDB, ONE,
     $                 C( JJ, JJ ), LDC )
                  DO II = JJ+NB, N, NB
                     IC = MIN( NB, N-II+1 )
                     CALL SGEMM( 'N', 'T', IC, NB, KC, ALPHA,
     $                    A( II, KK ), LDA, B( JJ, KK ), LDB, ONE,
     $                    C( II, JJ ), LDC )
                     CALL SGEMM( 'N', 'T', IC, NB, KC, -ALPHA,
     $                    B( II, KK ), LDB, A( JJ, KK ), LDA, ONE,
     $                    C( II, JJ ), LDC )
                  END DO
               END DO
            END DO
         END IF
      ELSE
*
*        Form  C := alpha*A**T*B - alpha*B**T*A + C.
*
         IF ( UPPER ) THEN
            DO KK = 1, K, NB
               KC = MIN( NB, K-KK+1 )
               DO JJ = 1, N, NB
                  JC = MIN( NB, N-JJ+1 )
*
*                 Call SSSR2K_SM for diagonal blocks,
*                 and SGEMM for off-diagonal blocks.
*
                  DO II = 1, JJ-NB, NB
                     CALL SGEMM( 'T', 'N', NB, JC, KC, ALPHA,
     $                    A( KK, II ), LDA, B( KK, JJ ), LDB, ONE,
     $                    C( II, JJ ), LDC )
                     CALL SGEMM( 'T', 'N', NB, JC, KC, -ALPHA,
     $                    B( KK, II ), LDB, A( KK, JJ ), LDA, ONE,
     $                    C( II, JJ ), LDC )
                  END DO
                  CALL SSSR2K_SM( UPPER, NOTRANS, JC, KC, ALPHA,
     $                 A( KK, JJ ), LDA, B( KK, JJ ), LDB, ONE,
     $                 C( JJ, JJ ), LDC )
               END DO
            END DO
         ELSE
            DO KK = 1, K, NB
               KC = MIN( NB, K-KK+1 )
               DO JJ = 1, N, NB
                  JC = MIN( NB, N-JJ+1 )
*
*                 Call SSSR2K_SM for diagonal blocks,
*                 and SGEMM for off-diagonal blocks.
*
                  CALL SSSR2K_SM( UPPER, NOTRANS, JC, KC, ALPHA,
     $                 A( KK, JJ ), LDA, B( KK, JJ ), LDB, ONE,
     $                 C( JJ, JJ ), LDC )
                  DO II = JJ+NB, N, NB
                     IC = MIN( NB, N-II+1 )
                     CALL SGEMM( 'T', 'N', IC, NB, KC, ALPHA,
     $                    A( KK, II ), LDA, B( KK, JJ ), LDB, ONE,
     $                    C( II, JJ ), LDC )
                     CALL SGEMM( 'T', 'N', IC, NB, KC, -ALPHA,
     $                    B( KK, II ), LDB, A( KK, JJ ), LDA, ONE,
     $                    C( II, JJ ), LDC )
                  END DO
               END DO
            END DO
         END IF
      END IF
*
      RETURN
*
*     End of SSSR2K.
*
      END
*
*
*
      SUBROUTINE SSSR2K_SM( UPPER, NOTRANS, N, K, ALPHA, A, LDA, B, LDB,
     $                      BETA, C, LDC )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      LOGICAL            UPPER, NOTRANS
      INTEGER            N, K, LDA, LDB, LDC
      REAL               ALPHA, BETA
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*     .. Local Scalars ..
      REAL               TEMP1, TEMP2
      INTEGER            I, J, L, NROWA
*     ..
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Executable Statements ..
*
      IF ( NOTRANS ) THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
*
      IF ( NOTRANS ) THEN
*
*        Form  C := alpha*A*B**T - alpha*B*A**T + beta*C.
*
         IF ( UPPER ) THEN
            DO J = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  DO I = 1, J
                     C( I, J ) = ZERO
                  END DO
               ELSE IF ( BETA .NE. ONE ) THEN
                  DO I = 1 , J
                     C( I, J ) = BETA*C( I, J )
                  END DO
               END IF
               DO L = 1, K
                  IF ( ( A (J, L ) .NE. ZERO ) .OR.
     $                 ( B(J, L ) .NE. ZERO ) ) THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO I = 1, J
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1
     $                       - B( I, L )*TEMP2
                     END DO
                  END IF
               END DO
            END DO
         ELSE
            DO J = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  DO I = J, N
                     C( I, J ) = ZERO
                  END DO
               ELSE IF ( BETA .NE. ONE ) THEN
                  DO I = J, N
                     C( I, J ) = BETA*C( I, J )
                  END DO
               END IF
               DO L = 1, K
                  IF ( ( A( J, L ) .NE. ZERO ) .OR.
     $                 ( B( J, L ) .NE. ZERO ) ) THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO I = J, N
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1
     +                       - B( I, L )*TEMP2
                     END DO
                  END IF
               END DO
            END DO
         END IF
      ELSE
*
*        Form  C := alpha*A**T*B - alpha*B**T*A + beta*C.
*
         IF ( UPPER ) THEN
            DO J = 1, N
               DO I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
                  END DO
                  IF ( BETA .EQ. ZERO ) THEN
                     C( I, J ) = ALPHA*TEMP1 - ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA*C( I, J ) + ALPHA*TEMP1
     $                    - ALPHA*TEMP2
                  END IF
               END DO
            END DO
         ELSE
            DO J = 1, N
               DO I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
                  END DO
                  IF ( BETA .EQ. ZERO ) THEN
                     C( I, J ) = ALPHA*TEMP1 - ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA*C( I, J ) + ALPHA*TEMP1
     $                    - ALPHA*TEMP2
                  END IF
               END DO
            END DO
         END IF
      END IF
*
      RETURN
*
*     End of SSSR2K_SM.
*
      END
