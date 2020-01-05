      SUBROUTINE PDSSMV( UPLO, N, ALPHA, A, IA, JA, DESCA, X, IX, JX,
     $                   DESCX, INCX, BETA, Y, IY, JY, DESCY, INCY,
     $                   WORK )
*
      IMPLICIT NONE
*
*  -- PBLAS-like routine --
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            N, IA, JA, IX, JX, INCX, IY, JY, INCY
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCX( * ), DESCY( * )
      DOUBLE PRECISION   A( * ), X( * ), Y( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDSSMV performs the matrix-vector operation
*
*     sub( Y ) := alpha*sub( A )*sub( X ) + beta*sub( Y ),
*
*  where
*
*     sub( A ) denotes A(IA:IA+M-1,JA:JA+N-1),
*
*     sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                      X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X, and,
*
*     sub( Y ) denotes Y(IY,JY:JY+N-1) if INCY = M_Y,
*                      Y(IY:IY+N-1,JY) if INCY = 1 and INCY <> M_Y.
*
*  Alpha and beta are scalars, sub( X ) and sub( Y ) are n element
*  subvectors and sub( A ) is an n by n skew-symmetric submatrix.
*
*  Unlike PBLAS routines, this subroutine requires some workspace.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  UPLO    (global input) CHARACTER*1
*          On entry, UPLO specifies whether the local pieces of the
*          array A containing the upper or lower triangular part of the
*          skew-symmetric submatrix sub( A ) are to be referenced as
*          follows:
*
*             UPLO = 'U' or 'u'   Only the local pieces corresponding to
*                                 the upper triangular part of the skew-
*                                 symmetric submatrix sub( A ) are to be
*                                 referenced,
*
*             UPLO = 'L' or 'l'   Only the local pieces corresponding to
*                                 the lower triangular part of the skew-
*                                 symmetric submatrix sub( A ) are to be
*                                 referenced.
*
*  N       (global input) INTEGER
*          On entry, N specifies the order of the submatrix sub( A ).
*          N must be at least zero.
*
*  ALPHA   (global input) DOUBLE PRECISION
*          On entry, ALPHA specifies the scalar alpha.  When ALPHA is
*          supplied as zero then the local entries of the arrays A and
*          X corresponding to the entries of the submatrix sub( A ) and
*          the subvector sub( X ) need not be set on input.
*
*  A       (local input) DOUBLE PRECISION array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before entry, this array contains
*          the local entries of the matrix A.
*          Before entry with UPLO = 'U' or 'u', this array contains the
*          local entries of the upper triangular part of the
*          skew-symmetric submatrix sub( A ), and the local entries of
*          the strictly lower triangular of sub( A ) are not referenced.
*          Before entry with UPLO = 'L' or 'l', this array contains the
*          local entries of the lower triangular part of the
*          skew-symmetric submatrix sub( A ), and the local entries of
*          the strictly upper triangular of sub( A ) are not referenced.
*
*  IA      (global input) INTEGER
*          On entry, IA specifies A's global row index, which points to
*          the beginning of the submatrix sub( A ).
*
*  JA      (global input) INTEGER
*          On entry, JA specifies A's global column index, which points
*          to the beginning of the submatrix sub( A ).
*
*  DESCA   (global and local input) INTEGER array
*          On entry, DESCA is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix A.
*
*  X       (local input) DOUBLE PRECISION array
*          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
*          is at least MAX( 1, Lr( 1, IX ) ) when INCX = M_X and
*          MAX( 1, Lr( 1, IX+N-1 ) ) otherwise, and, Kx is at least
*          Lc( 1, JX+N-1 ) when INCX = M_X and Lc( 1, JX ) otherwise.
*          Before entry, this array contains the local entries of the
*          matrix X.
*
*  IX      (global input) INTEGER
*          On entry, IX specifies X's global row index, which points to
*          the beginning of the submatrix sub( X ).
*
*  JX      (global input) INTEGER
*          On entry, JX specifies X's global column index, which points
*          to the beginning of the submatrix sub( X ).
*
*  DESCX   (global and local input) INTEGER array
*          On entry, DESCX is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix X.
*
*  INCX    (global input) INTEGER
*          On entry, INCX  specifies the global increment for the
*          elements of X.  Only two values of INCX are supported in this
*          version, namely 1 and M_X.  INCX must not be zero.
*
*  BETA    (global input) DOUBLE PRECISION
*          On entry, BETA specifies the scalar beta.  When BETA is
*          supplied as zero then the local entries of the array Y
*          corresponding to the entries of the subvector sub( Y ) need
*          not be set on input.
*
*  Y       (local input/local output) DOUBLE PRECISION array
*          On entry, Y is an array of dimension (LLD_Y, Ky), where LLD_Y
*          is at least MAX( 1, Lr( 1, IY ) ) when INCY = M_Y and
*          MAX( 1, Lr( 1, IY+N-1 ) ) otherwise, and, Ky is at least
*          Lc( 1, JY+N-1 ) when INCY = M_Y and Lc( 1, JY ) otherwise.
*          Before entry, this array contains the local entries of the
*          matrix Y.  On exit, sub( Y ) is overwritten by the updated
*          subvector.
*
*  IY      (global input) INTEGER
*          On entry, IY specifies Y's global row index, which points to
*          the beginning of the submatrix sub( Y ).
*
*  JY      (global input) INTEGER
*          On entry, JY specifies Y's global column index, which points
*          to the beginning of the submatrix sub( Y ).
*
*  DESCY   (global and local input) INTEGER array
*          On entry, DESCY is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix Y.
*
*  INCY    (global input) INTEGER
*          On entry, INCY specifies the global increment for the
*          elements of Y.  Only two values of INCY are supported in this
*          version, namely 1 and M_Y.  INCY must not be zero.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array
*          WORK needs to sufficiently long to hold a full copy of the
*          local part of sub( X ).
*
*  =====================================================================
*
*  Level 2 PBLAS-like routine.
*
*  Written by Meiyue Shao, Lawrence Berkeley National Laboratory.
*  Last change: October 2014
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICTXT, MB, NB, LLD
*     ..
*     .. Local Arrays ..
      INTEGER            DESCWK( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DESCSET, PDAXPY, PDCOPY, PDSCAL, PDTRMV
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      IF ( N .LE. 0 ) RETURN
*
      ICTXT = DESCX( CTXT_ )
      MB = DESCX( MB_ )
      NB = DESCX( NB_ )
      LLD = DESCX( LLD_ )
      CALL DESCSET( DESCWK, N, 1, MB, NB, 0, 0, ICTXT, LLD )
*
*     Use two TRMV calls to achieve SSMV.
*
      CALL PDSCAL( N, BETA, Y, IY, JY, DESCY, INCY )
*
      IF ( ALPHA .NE. ZERO ) THEN
         CALL PDCOPY( N, X, IX, JX, DESCX, INCX, WORK, 1, 1, DESCWK, 1 )
         CALL PDTRMV( UPLO, 'N', 'N', N, A, IA, JA, DESCA, WORK, 1, 1,
     $        DESCWK, 1 )
         CALL PDAXPY( N, ALPHA, WORK, 1, 1, DESCWK, 1, Y, IY, JY, DESCY,
     $        INCY )
         CALL PDCOPY( N, X, IX, JX, DESCX, INCX, WORK, 1, 1, DESCWK, 1 )
         CALL PDTRMV( UPLO, 'T', 'N', N, A, IA, JA, DESCA, WORK, 1, 1,
     $        DESCWK, 1 )
         CALL PDAXPY( N, -ALPHA, WORK, 1, 1, DESCWK, 1, Y, IY, JY,
     $        DESCY, INCY )
      END IF
*
      RETURN
*
*     End of PDSSMV.
*
      END
