      SUBROUTINE PZEMBED1( N, A, IA, JA, DESCA, B, IB, JB, DESCB, M, IM,
     $                     JM, DESCM, WORK )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, IA, JA, IB, JB, IM, JM
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCM( * )
      DOUBLE PRECISION   M( * ), WORK( * )
      COMPLEX*16         A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PZEMBED1() is an auxiliary routine called by PZBSEIG().
*  It generates a 2n-by-2n real symmetric matrix
*
*     M = [ real(A)+real(B), imag(A)-imag(B);
*          -imag(A)-imag(B), real(A)-real(B) ],
*
*  where A is n-by-n Hermitian, and B is n-by-n symmetric.
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
*  A       (local input) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          This array contains the local pieces of the N-by-N Hermitian
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
*  B       (local input) COMPLEX*16 pointer into the local memory
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
*          to an array of dimension (LLD_M, LOCc(JM+2*N-1)).
*          On exit, this array contains the local pieces of the 2N-by-2N
*          symmetric distributed matrix sub( M ). The leading 2N-by-2N
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
*  DESCM   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix M.
*
*  WORK    (local output) DOUBLE PRECISION array.
*          The workspace needs to be at least as big as A and B.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, MBA, MBB,
     $                   NBA, NBB, AROWS, ACOLS, LLDA, BROWS, BCOLS,
     $                   LLDB, RSRC, CSRC, I, J, I0, J0
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DIMAG
*     ..
*     .. External Functions ..
      EXTERNAL           NUMROC
      INTEGER            NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, PDTRADD,
     $                   PDELSET, PZELGET
*     ..
*     .. Executable Statements ..
*
      ICTXT = DESCA( CTXT_ )
      MBA = DESCA( MB_ )
      NBA = DESCA( NB_ )
      MBB = DESCB( MB_ )
      NBB = DESCB( NB_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF ( ICTXT .EQ. DESCB( CTXT_ ) .AND. MBA .EQ. NBA .AND.
     $     MBB .EQ. NBB .AND. IA .EQ. 1 .AND. JA .EQ. 1 .AND.
     $     IB .EQ. 1 .AND. JB .EQ. 1 ) THEN
*
*        A quick version for almost perfect alignment.
*
         LLDB = DESCB( LLD_ )
         BROWS = NUMROC( DESCB( M_ ), MBB, MYROW, 0, NPROW )
         BCOLS = NUMROC( DESCB( N_ ), NBB, MYCOL, 0, NPCOL )
         CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL, I0,
     $        J0, RSRC, CSRC )
         DO J = J0, BCOLS
            DO I = I0, BROWS
               WORK( I + (J-1)*LLDB ) = DBLE( B( I + (J-1)*LLDB ) )
            END DO
         END DO
         CALL PDTRADD( 'L', 'N', N, N, ONE, WORK, IB, JB, DESCB, ZERO,
     $        M, IM, JM, DESCM )
         CALL PDTRADD( 'L', 'N', N, N, -ONE, WORK, IB, JB, DESCB, ZERO,
     $        M, IM + N, JM + N, DESCM )
*
         DO J = J0, BCOLS
            DO I = I0, BROWS
               WORK( I + (J-1)*LLDB ) = DIMAG( B( I + (J-1)*LLDB ) )
            END DO
         END DO
*
*        The diagonal of imag(B) is accumulated only once.
*
         CALL PDTRADD( 'L', 'N', N, N, -ONE, WORK, IB, JB, DESCB, ZERO,
     $        M, IM + N, JM, DESCM )
         CALL PDTRADD( 'U', 'T', N, N, -ONE, WORK, IB, JB, DESCB, ZERO,
     $        M, IM + N, JM, DESCM )
*
         LLDA = DESCA( LLD_ )
         AROWS = NUMROC( DESCA( M_ ), MBA, MYROW, 0, NPROW )
         ACOLS = NUMROC( DESCA( N_ ), NBA, MYCOL, 0, NPCOL )
         CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, I0,
     $        J0, RSRC, CSRC )
*
         DO J = J0, ACOLS
            DO I = I0, AROWS
               WORK( I + (J-1)*LLDA ) = DBLE( A( I + (J-1)*LLDA ) )
            END DO
         END DO
         CALL PDTRADD( 'L', 'N', N, N, ONE, WORK, IA, JA, DESCA, ONE,
     $        M, IM, JM, DESCM )
         CALL PDTRADD( 'L', 'N', N, N, ONE, WORK, IA, JA, DESCA, ONE,
     $        M, IM + N, JM + N, DESCM )
*
         DO J = J0, ACOLS
            DO I = I0, AROWS
               WORK( I + (J-1)*LLDA ) = DIMAG( A( I + (J-1)*LLDA ) )
            END DO
         END DO
*
*        The diagonal of imag(A) is accumulated twice.
*        This is safe because these entries are assumed to be zero.
*
         CALL PDTRADD( 'L', 'N', N, N, -ONE, WORK, IA, JA, DESCA, ONE,
     $        M, IM + N, JM, DESCM )
         CALL PDTRADD( 'U', 'T', N, N, ONE, WORK, IA, JA, DESCA, ONE,
     $        M, IM + N, JM, DESCM )
*
      ELSE
*
*        A slow version that handles the general case.
*
         DO J = 0, N-1
            DO I = J, N-1
               CALL PZELGET( 'A', ' ', ALPHA, A, IA + I, JA + J, DESCA )
               CALL PZELGET( 'A', ' ', BETA, B, IB + I, JB + J, DESCB )
               CALL PDELSET( M, IM + I, JM + J, DESCM,
     $              DBLE( ALPHA ) + DBLE( BETA ) )
               CALL PDELSET( M, IM + N+I, JM + N+J, DESCM,
     $              DBLE( ALPHA ) - DBLE( BETA ) )
               CALL PDELSET( M, IM + N+J, JM + I, DESCM,
     $              DIMAG( ALPHA ) - DIMAG( BETA ) )
               CALL PDELSET( M, IM + N+I, JM + J, DESCM,
     $              -DIMAG( ALPHA ) - DIMAG( BETA ) )
            END DO
         END DO
      END IF
*
      RETURN
*
*     End of PZEMBED1().
*
      END
