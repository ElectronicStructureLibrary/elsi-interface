      SUBROUTINE PZEMBED2( N, A, IA, JA, DESCA, B, IB, JB, DESCB, H, IH,
     $                     JH, DESCH, WORK )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, IA, JA, IB, JB, IH, JH
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCH( * )
      COMPLEX*16         A( * ), B( * ), H( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZEMBED2() is an auxiliary routine called by PZBSEIG().
*  It generates a 2n-by-2n complex matrix
*
*     H = [       A,        B;
*          -conj(B), -conj(A) ],
*
*  where A is n-by-n Hermitian, and B is n-by-n symmetric.
*  Only the lower triangular parts of A and B are referenced,
*  but all entries of H are explicitly formed.
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
*  H       (local output) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_M, LOCc(JH+2*N-1)).
*          On exit, this array contains the local pieces of the 2N-by-2N
*          distributed matrix sub( H ).
*
*  IH      (global input) INTEGER
*          The row index in the global array H indicating the first
*          row of sub( H ).
*
*  JH      (global input) INTEGER
*          The column index in the global array H indicating the
*          first column of sub( H ).
*
*  DESCH   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix H.
*
*  WORK    (local output) COMPLEX*16 array.
*          The workspace needs to be at least as long as A or B.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                     ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      COMPLEX*16         ZTMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
*     ..
*     .. External Subroutines ..
      EXTERNAL           PZGEADD, PZTRADD, PZLACPY, PZLASET,
     $                   PZELGET, PZELSET
*     ..
*     .. Executable Statements ..
*
      IF ( IA .EQ. 1 .AND. JA .EQ. 1 .AND. IB .EQ. 1 .AND. JB .EQ. 1 )
     $     THEN
*
*        A quick version for IA=JA=IB=JB=1.
*
         CALL PZLASET( 'L', N, N, ZERO, ZERO, WORK, IA, JA, DESCA )
         CALL PZLACPY( 'L', N-1, N-1, A, IA + 1, JA, DESCA, WORK,
     $        IA + 1, JA, DESCA )
         CALL PZTRADD( 'U', 'C', N, N, ONE, WORK, IA, JA, DESCA, ZERO,
     $        H, IH, JH, DESCH )
         CALL PZTRADD( 'L', 'N', N, N, ONE, A, IA, JA, DESCA, ZERO,
     $        H, IH, JH, DESCH )
         CALL PZGEADD( 'T', N, N, -ONE, H, IH, JH, DESCH, ZERO,
     $        H, IH + N, JH + N, DESCH )
*
         CALL PZLASET( 'L', N, N, ZERO, ZERO, WORK, IB, JB, DESCB )
         CALL PZLACPY( 'L', N-1, N-1, B, IB + 1, JB, DESCA, WORK,
     $        IB + 1, JB, DESCB )
         CALL PZTRADD( 'U', 'T', N, N, ONE, WORK, IB, JB, DESCB, ZERO,
     $        H, IH, JH + N, DESCH )
         CALL PZTRADD( 'L', 'N', N, N, ONE, B, IB, JB, DESCB, ZERO,
     $        H, IH, JH + N, DESCH )
         CALL PZGEADD( 'C', N, N, -ONE, H, IH, JH + N, DESCH, ZERO,
     $        H, IH + N, JH, DESCH )
*
      ELSE
*
*        A slow version that handles the general case.
*
         DO J = 0, N-1
            DO I = J, N-1
               CALL PZELGET( 'A', ' ', ZTMP, A, IA + I, JA + J, DESCA )
               CALL PZELSET( H, IH + I, JH + J, DESCH, ZTMP )
               CALL PZELSET( H, IH + J, JH + I, DESCH, DCONJG( ZTMP ) )
               CALL PZELSET( H, IH + N+I, JH + N+J, DESCH,
     $              -DCONJG( ZTMP ) )
               CALL PZELSET( H, IH + N+J, JH + N+I, DESCH, -ZTMP )
               CALL PZELGET( 'A', ' ', ZTMP, B, IB + I, JB + J, DESCB )
               CALL PZELSET( H, IH + I, JH + N+J, DESCH, ZTMP )
               CALL PZELSET( H, IH + J, JH + N+I, DESCH, ZTMP )
               CALL PZELSET( H, IH + N+I, JH + J, DESCH,
     $              -DCONJG( ZTMP ) )
               CALL PZELSET( H, IH + N+J, JH + I, DESCH,
     $              -DCONJG( ZTMP ) )
            END DO
         END DO
      END IF
*
      RETURN
*
*     End of PZEMBED2().
*
      END
