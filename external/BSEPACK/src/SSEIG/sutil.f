      SUBROUTINE SGENMAT( M, N, A, LDA, ISEED, AMAX )
*
      IMPLICIT NONE
*
      INTEGER            M, N, LDA
      REAL               AMAX
      INTEGER            ISEED( 4 )
      REAL               A( LDA, * )
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           SLARND
      REAL               SLARND
*
      INTEGER            I, J
*
      AMAX = ZERO
      DO J = 1, N
         DO I = 1, M
            A( I, J ) = SLARND( 3, ISEED )
            AMAX = MAX( AMAX, ABS( A( I, J ) ) )
         END DO
      END DO
*
      END
*
*
*
      SUBROUTINE SGENSSMAT( N, A, LDA, ISEED, AMAX )
*
      IMPLICIT NONE
*
      INTEGER            N, LDA
      REAL               AMAX
      INTEGER            ISEED( 4 )
      REAL               A( LDA, * )
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           SLARND
      REAL               SLARND
*
      INTEGER            I, J
*
      DO J = 1, N
         DO I = J+1, N
            A( I, J ) = SLARND( 3, ISEED )
         END DO
      END DO
      AMAX = ZERO
      DO J = 1, N
         DO I = 1, J-1
            A( I, J ) = -A( J, I )
            AMAX = MAX( AMAX, ABS( A( I, J ) ) )
         END DO
      END DO
      DO I = 1, N
         A( I, I ) = ZERO
      END DO
*
      END
*
*
*
      SUBROUTINE SGENVEC( N, X, INCX, ISEED, XMAX )
*
      IMPLICIT NONE
*
      INTEGER            N, INCX
      REAL               XMAX
      INTEGER            ISEED( 4 )
      REAL               X( * )
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           SLARND
      REAL               SLARND
*
      INTEGER            I, IX
*
      IX = 1
      IF ( INCX .LT. 0 ) IX = 1 - (N-1)*INCX
      XMAX = ZERO
      DO I = 1, N
         X( IX ) = SLARND( 3, ISEED )
         XMAX = MAX( XMAX, ABS( X( IX ) ) )
         IX = IX + INCX
      END DO
*
      END
*
*
*
      REAL FUNCTION SERRMAT( UPLO, M, N, A, LDA, B, LDB )
*
*     SERRMAT returns || A - B ||_F.
*     When UPLO = 'U'/'L', A and B are assumed to be upper/lower
*     trangular.
*
      IMPLICIT NONE
*
      CHARACTER          UPLO
      INTEGER            M, N, LDA, LDB
      REAL               A( LDA, * ), B( LDB, * )
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*
      INTEGER            I, J
      REAL               TMP
      INTRINSIC          SQRT
      EXTERNAL           LSAME
      LOGICAL            LSAME
*
      TMP = ZERO
      IF ( LSAME( UPLO, 'L' ) ) THEN
         DO J = 1, N
            DO I = J, M
               TMP = TMP + ( A( I, J ) - B( I, J ) )**2
            END DO
         END DO
      ELSE IF ( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            DO I = 1, J
               TMP = TMP + ( A( I, J ) - B( I, J ) )**2
            END DO
         END DO
      ELSE
         DO J = 1, N
            DO I = 1, M
               TMP = TMP + ( A( I, J ) - B( I, J ) )**2
            END DO
         END DO
      END IF
      SERRMAT = SQRT( TMP )
      RETURN
*
      END
*
*
*
      REAL FUNCTION SERRVEC( N, X, INCX, ALPHA, Y, INCY )
*
*     SERRVEC returns || X - ALPHA*Y ||_2.
*
      IMPLICIT NONE
*
      INTEGER            N, INCX, INCY
      REAL               ALPHA
      REAL               X( * ), Y( * )
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*
      INTEGER            I, IX, IY
      REAL               TMP
      INTRINSIC          SQRT
*
      TMP = ZERO
      IX = 1
      IY = 1
      IF ( INCX .LT. 0 ) IX = 1 - (N-1)*INCX
      IF ( INCY .LT. 0 ) IY = 1 - (N-1)*INCY
      DO I = 1, N
         TMP = TMP + ( X( IX ) - ALPHA*Y( IY ) )**2
         IX = IX + INCX
         IY = IY + INCY
      END DO
      SERRVEC = SQRT( TMP )
      RETURN
*
      END
