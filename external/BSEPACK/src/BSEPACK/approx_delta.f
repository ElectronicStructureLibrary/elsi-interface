      DOUBLE PRECISION FUNCTION APPROX_DELTA( APPROX, X, SIGMA )
*
      IMPLICIT NONE
      INCLUDE 'solver.f'
*
*     .. Scalar Arguments ..
      INTEGER            APPROX
      DOUBLE PRECISION   X, SIGMA
*     ..
*
*  Purpose
*  =======
*
*  APPROX_DELTA() computes an approximate value of the delta function.
*   returns the value of the Gaussian function
*
*  Arguments
*  =========
*
*  APPROX  (input) INTEGER
*          Indicate how the delta function is approximated.
*          See solver.f for details.
*
*  X       (input) DOUBLE PRECISION
*          Value of the variable x.
*
*  SIGMA   (input) DOUBLE PRECISION
*          Broadening factor in the approximation.
*          SIGMA > 0.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, NEG_HALF, INV_SQRT_TWOPI, INV_PI
      PARAMETER          ( ZERO = 0.0D+0, NEG_HALF = -0.5D+0,
     $                     INV_SQRT_TWOPI = 0.39894228040143268D+0,
     $                     INV_PI = 0.31830988618379067D+0 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          EXP, IAND
*     ..
*     .. Executable Statements ..
*
      IF ( IAND( APPROX, BSE_DELTA ) .EQ. BSE_GAUSSIAN ) THEN
*
*        Gaussian:
*           exp( -x**2/( 2*sigma**2 ) )/( sqrt( 2*pi )*sigma )
*
         APPROX_DELTA =
     $        INV_SQRT_TWOPI/SIGMA*EXP( NEG_HALF*( X/SIGMA )**2 )
      ELSE IF ( IAND( APPROX, BSE_DELTA ) .EQ. BSE_LORENTZIAN ) THEN
*
*        Lorentzian:
*           1/pi * sigma / ( x**2 + sigma**2 )
*
         APPROX_DELTA =
     $        INV_PI*SIGMA/( X**2 + SIGMA**2 )
      ELSE
*
*        Return zero for unspecified approximation.
*
         APPROX_DELTA = ZERO
      END IF
*
      RETURN
*
*     End of APPROX_DELTA().
*
      END
