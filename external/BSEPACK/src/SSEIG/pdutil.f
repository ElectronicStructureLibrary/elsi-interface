      SUBROUTINE PDGENMAT( M, N, A, IA, JA, DESCA, ISEED, AMAX )
*
      IMPLICIT NONE
*
      INTEGER            M, N, IA, JA
      DOUBLE PRECISION   AMAX
      INTEGER            DESCA( 9 ), ISEED( 4 )
      DOUBLE PRECISION   A( * )
*
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           DLARND, PDELSET
      DOUBLE PRECISION   DLARND
*
      INTEGER            I, J
      DOUBLE PRECISION   TMP
*
      AMAX = ZERO
      DO J = 0, N-1
         DO I = 0, M-1
            TMP = DLARND( 3, ISEED )
            CALL PDELSET( A, IA+I, JA+J, DESCA, TMP )
            AMAX = MAX( AMAX, ABS( TMP ) )
         END DO
      END DO
*
      END
*
*
*
      SUBROUTINE PDGENSSMAT( N, A, IA, JA, DESCA, ISEED, AMAX )
*
      IMPLICIT NONE
*
      INTEGER            N, IA, JA
      DOUBLE PRECISION   AMAX
      INTEGER            DESCA( 9 ), ISEED( 4 )
      DOUBLE PRECISION   A( * )
*
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           DLARND, PDELSET
      DOUBLE PRECISION   DLARND
*
      INTEGER            I, J
      DOUBLE PRECISION   TMP
*
      AMAX = ZERO
      DO J = 0, N-1
         DO I = J+1, N-1
            TMP = DLARND( 3, ISEED )
            CALL PDELSET( A, IA+I, JA+J, DESCA, TMP )
            CALL PDELSET( A, IA+J, JA+I, DESCA, -TMP )
            AMAX = MAX( AMAX, ABS( TMP ) )
         END DO
      END DO
      DO I = 0, N-1
         CALL PDELSET( A, IA+I, JA+I, DESCA, ZERO )
      END DO
*
      END
*
*
*
      SUBROUTINE PDGENVEC( N, X, IX, JX, DESCX, INCX, ISEED, XMAX )
*
      IMPLICIT NONE
*
      INTEGER            N, IX, JX, INCX
      DOUBLE PRECISION   XMAX
      INTEGER            DESCX( 9 ), ISEED( 4 )
      DOUBLE PRECISION   X( * )
*
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
      INTRINSIC          ABS, MAX
      EXTERNAL           DLARND, PDELSET
      DOUBLE PRECISION   DLARND
*
      INTEGER            I
      DOUBLE PRECISION   TMP
*
      XMAX = ZERO
      DO I = 0, N-1
         TMP = DLARND( 3, ISEED )
         IF ( INCX .EQ. 1 ) THEN
            CALL PDELSET( X, IX+I, JX, DESCX, TMP )
         ELSE
            CALL PDELSET( X, IX, JX+I, DESCX, TMP )
         END IF
         XMAX = MAX( XMAX, ABS( TMP ) )
      END DO
*
      END
*
*
*
      DOUBLE PRECISION FUNCTION PDERRMAT( UPLO, M, N, A, IA, JA, DESCA,
     $     B, IB, JB, DESCB )
*
*     DERRMAT returns || A - B ||_F.
*     When UPLO = 'U'/'L', A and B are assumed to be upper/lower
*     trangular.
*
      IMPLICIT NONE
*
      CHARACTER          UPLO
      INTEGER            M, N, IA, JA, IB, JB
      INTEGER            DESCA( 9 ), DESCB( 9 )
      DOUBLE PRECISION   A( * ), B( * )
*
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
      INTEGER            I, J
      DOUBLE PRECISION   TMP, TMP1, TMP2
      INTRINSIC          SQRT
      EXTERNAL           LSAME, PDELGET
      LOGICAL            LSAME
*
      TMP = ZERO
      IF ( LSAME( UPLO, 'L' ) ) THEN
         DO J = 0, N-1
            DO I = J, M-1
               CALL PDELGET( 'A', ' ', TMP1, A, IA+I, JA+J, DESCA )
               CALL PDELGET( 'A', ' ', TMP2, B, IB+I, JB+J, DESCB )
               TMP = TMP + ( TMP1 - TMP2 )**2
            END DO
         END DO
      ELSE IF ( LSAME( UPLO, 'U' ) ) THEN
         DO J = 0, N-1
            DO I = 0, J-1
               CALL PDELGET( 'A', ' ', TMP1, A, IA+I, JA+J, DESCA )
               CALL PDELGET( 'A', ' ', TMP2, B, IB+I, JB+J, DESCB )
               TMP = TMP + ( TMP1 - TMP2 )**2
            END DO
         END DO
      ELSE
         DO J = 0, N-1
            DO I = 0, M-1
               CALL PDELGET( 'A', ' ', TMP1, A, IA+I, JA+J, DESCA )
               CALL PDELGET( 'A', ' ', TMP2, B, IB+I, JB+J, DESCB )
               TMP = TMP + ( TMP1 - TMP2 )**2
            END DO
         END DO
      END IF
      PDERRMAT = SQRT( TMP )
      RETURN
*
      END
*
*
*
      DOUBLE PRECISION FUNCTION PDERRVEC( N, X, IX, JX, DESCX, INCX,
     $     ALPHA, Y, IY, JY, DESCY, INCY )
*
*     DERRVEC returns || X - ALPHA*Y ||_2.
*
      IMPLICIT NONE
*
      INTEGER            N, IX, JX, INCX, IY, JY, INCY
      DOUBLE PRECISION   ALPHA
      INTEGER            DESCX( 9 ), DESCY( 9 )
      DOUBLE PRECISION   X( * ), Y( * )
*
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
      INTEGER            I
      DOUBLE PRECISION   TMP, TMP1, TMP2
      INTRINSIC          SQRT
*
      TMP = ZERO
      DO I = 0, N-1
         IF ( INCX .EQ. 1 ) THEN
            CALL PDELGET( 'A', ' ', TMP1, X, IX+I, JX, DESCX )
         ELSE
            CALL PDELGET( 'A', ' ', TMP1, X, IX, JX+I, DESCX )
         END IF
         IF ( INCY .EQ. 1 ) THEN
            CALL PDELGET( 'A', ' ', TMP2, Y, IY+I, JY, DESCY )
         ELSE
            CALL PDELGET( 'A', ' ', TMP2, Y, IY, JY+I, DESCY )
         END IF
         TMP = TMP + ( TMP1 - ALPHA*TMP2 )**2
      END DO
      PDERRVEC = SQRT( TMP )
      RETURN
*
      END
