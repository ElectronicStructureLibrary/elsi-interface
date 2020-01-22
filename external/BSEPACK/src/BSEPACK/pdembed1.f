      SUBROUTINE PDEMBED1( N, A, IA, JA, DESCA, B, IB, JB, DESCB, M, IM,
     $                     JM, DESCM, K, IK, JK, DESCK )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, IA, JA, IB, JB, IM, JM, IK, JK
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCM( * ), DESCK( * )
      DOUBLE PRECISION   A( * ), B( * ), M( * ), K( * )
*     ..
*
*  Purpose
*  =======
*
*  PDEMBED1() is an auxiliary routine called by PDBSEIG().
*  It generates two n-by-n real symmetric matrices M = A+B and K = A-B,
*  where both A and B are n-by-n symmetric.
*  Only the lower triangular parts of these matrices are referenced.
*
*  No argument check is performed, i.e.,
*  all arguments are assumed to be valid.
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
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrices sub( A ) and sub( B ).
*          N >= 0.
*
*  A       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          This array contains the local pieces of the N-by-N symmetric
*          distributed matrix sub( A ). The leading N-by-N lower
*          triangular part of sub( A ) contains the lower triangular
*          part of the distributed matrix, and its strictly upper
*          triangular part is not referenced.
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  B       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_B, LOCc(JB+N-1)).
*          This array contains the local pieces of the N-by-N symmetric
*          distributed matrix sub( B ). The leading N-by-N lower
*          triangular part of sub( B ) contains the lower triangular
*          part of the distributed matrix, and its strictly upper
*          triangular part is not referenced.
*
*  IB      (global input) INTEGER
*          The row index in the global array B indicating the first
*          row of sub( B ).
*
*  JB      (global input) INTEGER
*          The column index in the global array B indicating the
*          first column of sub( B ).
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*
*  M       (local output) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_M, LOCc(JM+N-1)).
*          On exit, this array contains the local pieces of the N-by-N
*          symmetric distributed matrix sub( M ). The leading N-by-N
*          lower triangular part of sub( M ) contains the lower
*          triangular part of the distributed matrix, and its strictly
*          upper triangular part is not referenced.
*
*  IM      (global input) INTEGER
*          The row index in the global array M indicating the first
*          row of sub( M ).
*
*  JM      (global input) INTEGER
*          The column index in the global array M indicating the
*          first column of sub( M ).
*
*  K       (local output) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_K, LOCc(JK+N-1)).
*          On exit, this array contains the local pieces of the N-by-N
*          symmetric distributed matrix sub( K ). The leading N-by-N
*          lower triangular part of sub( K ) contains the lower
*          triangular part of the distributed matrix, and its strictly
*          upper triangular part is not referenced.
*
*  IK      (global input) INTEGER
*          The row index in the global array K indicating the first
*          row of sub( K ).
*
*  JK      (global input) INTEGER
*          The column index in the global array K indicating the
*          first column of sub( K ).
*
*  DESCK   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix K.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDTRADD, PDLACPY
*     ..
*     .. Executable Statements ..
*
      CALL PDLACPY( 'L', N, N, A, IA, JA, DESCA, M, IM, JM, DESCM )
      CALL PDLACPY( 'L', N, N, A, IA, JA, DESCA, K, IK, JK, DESCK )
      CALL PDTRADD( 'L', 'N', N, N, ONE, B, IB, JB, DESCB, ONE,
     $     M, IM, JM, DESCM )
      CALL PDTRADD( 'L', 'N', N, N, -ONE, B, IB, JB, DESCB, ONE,
     $     K, IK, JK, DESCK )
*
      RETURN
*
*     End of PDEMBED1().
*
      END
