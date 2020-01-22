      SUBROUTINE PSSSR2K( UPLO, TRANS, N, K, ALPHA, A, IA, JA, DESCA, B,
     $                    IB, JB, DESCB, BETA, C, IC, JC, DESCC, WORK )
*
      IMPLICIT NONE
*
*  -- PBLAS-like routine --
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, TRANS
      INTEGER            N, K, IA, JA, IB, JB, IC, JC
      REAL               ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCC( * )
      REAL               A( * ), B( * ), C( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSSSR2K performs one of the skew-symmetric rank 2k operations
*
*     sub( C ) := alpha*sub( A )*sub( B )**T +
*                 alpha*sub( B )*sub( A )**T +
*                 beta*sub( C ),
*
*  or
*
*     sub( C ) := alpha*sub( A )**T*sub( B ) +
*                 alpha*sub( B )**T*sub( A ) +
*                 beta*sub( C ),
*
*  where
*
*     sub( C ) denotes C(IC:IC+N-1,JC:JC+N-1),
*
*     sub( A ) denotes A(IA:IA+N-1,JA:JA+K-1) if TRANS = 'N',
*                      A(IA:IA+K-1,JA:JA+N-1) otherwise, and,
*
*     sub( B ) denotes B(IB:IB+N-1,JB:JB+K-1) if TRANS = 'N',
*                      B(IB:IB+K-1,JB:JB+N-1) otherwise.
*
*  Alpha and beta are scalars, sub( C ) is an n by n skew-symmetric
*  submatrix and sub( A ) and sub( B ) are n by k submatrices in the
*  first case and k by n submatrices in the second case.
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
*
*  TRANS   (global input) CHARACTER*1
*          On entry, TRANS specifies the operation to be performed as
*          follows:
*
*             TRANS = 'N' or 'n'
*               sub( C ) := alpha*sub( A )*sub( B )**T +
*                           alpha*sub( B )*sub( A )**T +
*                           beta*sub( C ),
*
*             TRANS = 'T' or 't'
*               sub( C ) := alpha*sub( B )**T*sub( A ) +
*                           alpha*sub( A )**T*sub( B ) +
*                           beta*sub( C ),
*
*             TRANS = 'C' or 'c'
*               sub( C ) := alpha*sub( B )**T*sub( A ) +
*                           alpha*sub( A )**T*sub( B ) +
*                           beta*sub( C ).
*
*  N       (global input) INTEGER
*          On entry, N specifies the order of the submatrix sub( C ).
*          N must be at least zero.
*
*  K       (global input) INTEGER
*          On entry with TRANS = 'N' or 'n', K specifies the number of
*          columns of the submatrix sub( A ) and sub( B ), and on entry
*          with TRANS = 'T' or 't' or 'C' or 'c', K specifies the number
*          of rows of the submatrices sub( A ) and sub( B ).
*          K must be at least zero.
*
*  ALPHA   (global input) REAL
*          On entry, ALPHA specifies the scalar alpha.  When ALPHA is
*          supplied as zero then the local entries of the arrays A and
*          B corresponding to the entries of the submatrices sub( A )
*          and sub( B ) respectively need not be set on input.
*
*  A       (local input) REAL array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+K-1 ) when TRANS = 'N' or 'n', and is at
*          least Lc( 1, JA+N-1 ) otherwise.  Before entry, this array
*          contains the local entries of the matrix A.
*          Before entry with TRANS = 'N' or 'n', this array contains the
*          local entries corresponding to the entries of the n by k
*          submatrix sub( A ), otherwise the local entries corresponding
*          to the entries of the k by n submatrix sub( A ).
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
*  B       (local input) REAL array
*          On entry, B is an array of dimension (LLD_B, Kb), where Kb is
*          at least Lc( 1, JB+K-1 ) when TRANS = 'N' or 'n', and is at
*          least Lc( 1, JB+N-1 ) otherwise.  Before entry, this array
*          contains the local entries of the matrix B.
*          Before entry with TRANS = 'N' or 'n', this array contains the
*          local entries corresponding to the entries of the n by k
*          submatrix sub( B ), otherwise the local entries corresponding
*          to the entries of the k by n submatrix sub( B ).
*
*  IB      (global input) INTEGER
*          On entry, IB specifies B's global row index, which points to
*          the beginning of the submatrix sub( B ).
*
*  JB      (global input) INTEGER
*          On entry, JB specifies B's global column index, which points
*          to the beginning of the submatrix sub( B ).
*
*  DESCB   (global and local input) INTEGER array
*          On entry, DESCB is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix B.
*
*  BETA    (global input) REAL
*          On entry, BETA specifies the scalar beta.  When BETA is
*          supplied as zero then the local entries of the array C
*          corresponding to the entries of the submatrix sub( C ) need
*          not be set on input.
*
*  C       (local input) REAL array
*          On entry, C is an array of dimension (LLD_C, Kc), where Kc is
*          at least Lc( 1, JC+N-1 ).  Before entry, this array contains
*          the local entries of the matrix C.
*          Before entry with UPLO = 'U' or 'u', this array contains the
*          local entries of the upper triangular part of the
*          skew-symmetric submatrix sub( C ), and the local entries of
*          the strictly lower triangular of sub( C ) are not referenced.
*          On exit, the upper triangular part of sub( C ) is overwritten
*          by the upper triangular part of the updated submatrix.
*          Before entry with UPLO = 'L' or 'l', this array contains the
*          local entries of the lower triangular part of the
*          skew-symmetric submatrix sub( C ), and the local entries of
*          the strictly upper triangular of sub( C ) are not referenced.
*          On exit, the lower triangular part of sub( C ) is overwritten
*          by the lower triangular part of the updated submatrix.
*
*  IC      (global input) INTEGER
*          On entry, IC specifies C's global row index, which points to
*          the beginning of the submatrix sub( C ).
*
*  JC      (global input) INTEGER
*          On entry, JC specifies C's global column index, which points
*          to the beginning of the submatrix sub( C ).
*
*  DESCC   (global and local input) INTEGER array
*          On entry, DESCC is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix C.
*
*  WORK    (local workspace/local output) REAL array
*          WORK needs to sufficiently long to hold a full copy of the
*          local part of sub( C ).
*
*  =====================================================================
*
*  Level 3 PBLAS-like routine.
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
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICTXT, MB, NB, LLD
*     ..
*     .. Local Arrays ..
      INTEGER            DESCWK( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DESCSET, PSGEMM, PSTRADD
*     ..
*     .. External Functions ..
      EXTERNAL           LSAME
      LOGICAL            LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      IF ( N .LE. 0 ) RETURN
*
      ICTXT = DESCC( CTXT_ )
      MB = DESCC( MB_ )
      NB = DESCC( NB_ )
      LLD = DESCC( LLD_ )
      CALL DESCSET( DESCWK, N, N, MB, NB, 0, 0, ICTXT, LLD )
*
*     Use one GEMM and two TRADD calls to achieve SSR2K.
*
      IF ( LSAME( TRANS, 'N' ) ) THEN
         CALL PSGEMM( 'N', 'T', N, N, K, ALPHA, A, IA, JA, DESCA,
     $        B, IB, JB, DESCB, ZERO, WORK, 1, 1, DESCWK )
      ELSE
         CALL PSGEMM( 'T', 'N', N, N, K, ALPHA, A, IA, JA, DESCA,
     $        B, IB, JB, DESCB, ZERO, WORK, 1, 1, DESCWK )
      END IF
      CALL PSTRADD( UPLO, 'N', N, N, ONE, WORK, 1, 1, DESCWK, BETA,
     $     C, IC, JC, DESCC )
      IF ( ALPHA .NE. ZERO ) THEN
         CALL PSTRADD( UPLO, 'T', N, N, -ONE, WORK, 1, 1, DESCWK, ONE,
     $        C, IC, JC, DESCC )
      END IF
      RETURN
*
*     End of PSSSR2K.
*
      END
