      SUBROUTINE DGENMAT( M, N, A, LDA, ISEED, AMAX )
*
      IMPLICIT NONE
*
      INTEGER            M, N, LDA
      DOUBLE PRECISION   AMAX
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * )
*
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           DLARND
      DOUBLE PRECISION   DLARND
*
      INTEGER            I, J
*
      AMAX = ZERO
      DO J = 1, N
         DO I = 1, M
            A( I, J ) = DLARND( 3, ISEED )
            AMAX = MAX( AMAX, ABS( A( I, J ) ) )
         END DO
      END DO
*
      END
*
*
*
      SUBROUTINE DGENSSMAT( N, A, LDA, ISEED, AMAX )
*
      IMPLICIT NONE
*
      INTEGER            N, LDA
      DOUBLE PRECISION   AMAX
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * )
*
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           DLARND
      DOUBLE PRECISION   DLARND
*
      INTEGER            I, J
*
      DO J = 1, N
         DO I = J+1, N
            A( I, J ) = DLARND( 3, ISEED )
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
      SUBROUTINE DGENVEC( N, X, INCX, ISEED, XMAX )
*
      IMPLICIT NONE
*
      INTEGER            N, INCX
      DOUBLE PRECISION   XMAX
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( * )
*
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           DLARND
      DOUBLE PRECISION   DLARND
*
      INTEGER            I, IX
*
      IX = 1
      IF ( INCX .LT. 0 ) IX = 1 - (N-1)*INCX
      XMAX = ZERO
      DO I = 1, N
         X( IX ) = DLARND( 3, ISEED )
         XMAX = MAX( XMAX, ABS( X( IX ) ) )
         IX = IX + INCX
      END DO
*
      END
*
*
*
      DOUBLE PRECISION FUNCTION DERRMAT( UPLO, M, N, A, LDA, B, LDB )
*
*     DERRMAT returns || A - B ||_F.
*     When UPLO = 'U'/'L', A and B are assumed to be upper/lower
*     trangular.
*
      IMPLICIT NONE
*
      CHARACTER          UPLO
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
      INTEGER            I, J
      DOUBLE PRECISION   TMP
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
      DERRMAT = SQRT( TMP )
      RETURN
*
      END
*
*
*
      DOUBLE PRECISION FUNCTION DERRVEC( N, X, INCX, ALPHA, Y, INCY )
*
*     DERRVEC returns || X - ALPHA*Y ||_2.
*
      IMPLICIT NONE
*
      INTEGER            N, INCX, INCY
      DOUBLE PRECISION   ALPHA
      DOUBLE PRECISION   X( * ), Y( * )
*
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
      INTEGER            I, IX, IY
      DOUBLE PRECISION   TMP
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
      DERRVEC = SQRT( TMP )
      RETURN
*
      END
