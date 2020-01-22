      SUBROUTINE DBSABSP_AUX( APPROX, K, NPTS, SIGMA, OMEGA, EPS,
     $                        NORMD2, X, INCX, MU, DIV )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            APPROX, K, NPTS, INCX, DIV
      DOUBLE PRECISION   SIGMA, NORMD2
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   OMEGA( * ), EPS( * ), X( INCX, * ), MU( * )
*     ..
*
*  Purpose
*  =======
*
*  DBSABSP_AUX() estimates the optical absorption spectra
*
*     epsilon( omega ) =
*        ||d||**2 * sum_i |x_i(1)|^2 * delta( omega - mu_i ) / mu_i**div
*
*  using the normalized eigenpairs (mu_j,x_j) of the k-by-k symmetric
*  tridiagonal matrix.
*  Here the norm is not necessarily the 2-norm.
*  The value ||d||**2 is an input argument.
*
*  No argument check is performed, i.e.,
*  all arguments are assumed to be valid.
*
*  Arguments
*  =========
*
*  APPROX  (input) INTEGER
*          Indicate how the delta function is approximated.
*          See solver.f for details.
*
*  K       (input) INTEGER
*          The number of rows and columns of the tridiagonal matrix.
*          K >= 0.
*
*  NPTS    (input) INTEGER
*          NPTS is the number of sampling points in OMEGA.
*          NPTS >= 0.
*
*  SIGMA   (input) DOUBLE PRECISION
*          Standard deviation of the Gaussian function.
*          SIGMA > 0.
*
*  OMEGA   (input) DOUBLE PRECISION array, dimension (NPTS)
*          Sampling points of omega.
*          All entries of OMEGA are nonnegative.
*
*  EPS     (output) DOUBLE PRECISION array, dimension (NPTS)
*          Sampling points of epsilon, i.e.,
*          EPS( I ) = epsilon( OMEGA( I ) ).
*
*  NORMD2  (input) DOUBLE PRECISION
*          The squared norm of the first block, d_1, of the dipole
*          vector d, i.e., ||d_1||**2.
*
*  X       (input) DOUBLE PRECISION array, dimension 1+(N-1)*INCX.
*          The first row from normalized eigenvectors of the
*          tridiagonal matrix obtained by the Lanczos process.
*
*  INCX    (input) INGEGER
*          The increment between elements of X. INCX >= 1.
*
*  MU      (input) DOUBLE PRECISION array, dimension (K)
*          The eigenvalues of the tridiagonal matrix.
*
*  DIV     (input) INTEGER
*          If DIV = 1, the delta function is divided by MU;
*          if DIV = 0, use the delta function without further dividing.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   DTMP
*     ..
*     .. External Functions ..
      EXTERNAL           APPROX_DELTA
      DOUBLE PRECISION   APPROX_DELTA
*     ..
*     .. Executable Statements ..
*
      DO J = 1, NPTS
         EPS( J ) = ZERO
      END DO
      DO I = 1, K
         IF ( MU( I ) .GT. ZERO ) THEN
            IF ( DIV .EQ. 0 ) THEN
               DTMP = NORMD2*X( 1, I )*X( 1, I )
            ELSE
               DTMP = NORMD2*X( 1, I )*X( 1, I )/MU( I )
            END IF
            DO J = 1, NPTS
               EPS( J ) = EPS( J ) + DTMP
     $              *(
     $              APPROX_DELTA( APPROX, OMEGA( J ) - MU( I ), SIGMA )
     $              -APPROX_DELTA( APPROX, OMEGA( J ) + MU( I ), SIGMA )
     $              )
            END DO
         END IF
      END DO
*
      RETURN
*
*     End of DBSABSP_AUX().
*
      END
