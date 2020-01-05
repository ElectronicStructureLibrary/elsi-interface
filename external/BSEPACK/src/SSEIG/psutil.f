      SUBROUTINE PSGENMAT( M, N, A, IA, JA, DESCA, ISEED, AMAX )
*
      IMPLICIT NONE
*
      INTEGER            M, N, IA, JA
      REAL               AMAX
      INTEGER            DESCA( 9 ), ISEED( 4 )
      REAL               A( * )
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           SLARND, PSELSET
      REAL               SLARND
*
      INTEGER            I, J
      REAL               TMP
*
      AMAX = ZERO
      DO J = 0, N-1
         DO I = 0, M-1
            TMP = SLARND( 3, ISEED )
            CALL PSELSET( A, IA+I, JA+J, DESCA, TMP )
            AMAX = MAX( AMAX, ABS( TMP ) )
         END DO
      END DO
*
      END
*
*
*
      SUBROUTINE PSGENSSMAT( N, A, IA, JA, DESCA, ISEED, AMAX )
*
      IMPLICIT NONE
*
      INTEGER            N, IA, JA
      REAL               AMAX
      INTEGER            DESCA( 9 ), ISEED( 4 )
      REAL               A( * )
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           SLARND, PSELSET
      REAL               SLARND
*
      INTEGER            I, J
      REAL               TMP
*
      AMAX = ZERO
      DO J = 0, N-1
         DO I = J+1, N-1
            TMP = SLARND( 3, ISEED )
            CALL PSELSET( A, IA+I, JA+J, DESCA, TMP )
            CALL PSELSET( A, IA+J, JA+I, DESCA, -TMP )
            AMAX = MAX( AMAX, ABS( TMP ) )
         END DO
      END DO
      DO I = 0, N-1
         CALL PSELSET( A, IA+I, JA+I, DESCA, ZERO )
      END DO
*
      END
*
*
*
      SUBROUTINE PSGENVEC( N, X, IX, JX, DESCX, INCX, ISEED, XMAX )
*
      IMPLICIT NONE
*
      INTEGER            N, IX, JX, INCX
      REAL               XMAX
      INTEGER            DESCX( 9 ), ISEED( 4 )
      REAL               X( * )
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           SLARND, PSELSET
      REAL               SLARND
*
      INTEGER            I
      REAL               TMP
*
      XMAX = ZERO
      DO I = 0, N-1
         TMP = SLARND( 3, ISEED )
         IF ( INCX .EQ. 1 ) THEN
            CALL PSELSET( X, IX+I, JX, DESCX, TMP )
         ELSE
            CALL PSELSET( X, IX, JX+I, DESCX, TMP )
         END IF
         XMAX = MAX( XMAX, ABS( TMP ) )
      END DO
*
      END
*
*
*
      REAL FUNCTION PSERRMAT( UPLO, M, N, A, IA, JA, DESCA,
     $     B, IB, JB, DESCB )
*
*     SERRMAT returns || A - B ||_F.
*     When UPLO = 'U'/'L', A and B are assumed to be upper/lower
*     trangular.
*
      IMPLICIT NONE
*
      CHARACTER          UPLO
      INTEGER            M, N, IA, JA, IB, JB
      INTEGER            DESCA( 9 ), DESCB( 9 )
      REAL               A( * ), B( * )
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*
      INTEGER            I, J
      REAL               TMP, TMP1, TMP2
      INTRINSIC          SQRT
      EXTERNAL           LSAME, PSELGET
      LOGICAL            LSAME
*
      TMP = ZERO
      IF ( LSAME( UPLO, 'L' ) ) THEN
         DO J = 0, N-1
            DO I = J, M-1
               CALL PSELGET( 'A', ' ', TMP1, A, IA+I, JA+J, DESCA )
               CALL PSELGET( 'A', ' ', TMP2, B, IB+I, JB+J, DESCB )
               TMP = TMP + ( TMP1 - TMP2 )**2
            END DO
         END DO
      ELSE IF ( LSAME( UPLO, 'U' ) ) THEN
         DO J = 0, N-1
            DO I = 0, J-1
               CALL PSELGET( 'A', ' ', TMP1, A, IA+I, JA+J, DESCA )
               CALL PSELGET( 'A', ' ', TMP2, B, IB+I, JB+J, DESCB )
               TMP = TMP + ( TMP1 - TMP2 )**2
            END DO
         END DO
      ELSE
         DO J = 0, N-1
            DO I = 0, M-1
               CALL PSELGET( 'A', ' ', TMP1, A, IA+I, JA+J, DESCA )
               CALL PSELGET( 'A', ' ', TMP2, B, IB+I, JB+J, DESCB )
               TMP = TMP + ( TMP1 - TMP2 )**2
            END DO
         END DO
      END IF
      PSERRMAT = SQRT( TMP )
      RETURN
*
      END
*
*
*
      REAL FUNCTION PSERRVEC( N, X, IX, JX, DESCX, INCX,
     $     ALPHA, Y, IY, JY, DESCY, INCY )
*
*     SERRVEC returns || X - ALPHA*Y ||_2.
*
      IMPLICIT NONE
*
      INTEGER            N, IX, JX, INCX, IY, JY, INCY
      REAL               ALPHA
      INTEGER            DESCX( 9 ), DESCY( 9 )
      REAL               X( * ), Y( * )
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*
      INTEGER            I
      REAL               TMP, TMP1, TMP2
      INTRINSIC          SQRT
*
      TMP = ZERO
      DO I = 0, N-1
         IF ( INCX .EQ. 1 ) THEN
            CALL PSELGET( 'A', ' ', TMP1, X, IX+I, JX, DESCX )
         ELSE
            CALL PSELGET( 'A', ' ', TMP1, X, IX, JX+I, DESCX )
         END IF
         IF ( INCY .EQ. 1 ) THEN
            CALL PSELGET( 'A', ' ', TMP2, Y, IY+I, JY, DESCY )
         ELSE
            CALL PSELGET( 'A', ' ', TMP2, Y, IY, JY+I, DESCY )
         END IF
         TMP = TMP + ( TMP1 - ALPHA*TMP2 )**2
      END DO
      PSERRVEC = SQRT( TMP )
      RETURN
*
      END
