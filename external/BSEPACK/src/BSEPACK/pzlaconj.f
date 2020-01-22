      SUBROUTINE PZLACONJ( M, N, A, IA, JA, DESCA )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            M, N, IA, JA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( * )
      INTEGER            DESCA( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLACONJ() converts A to its complex conjugate conj(A).
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
*  M       (global input) INTEGER
*          The number of rows of the distributed submatrix sub( A ).
*          M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns of the distributed submatrix sub( A ).
*          N >= 0.
*
*  A       (local iutput/output) COMPLEX*16 pointer into the local
*          memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On exit, all entries of A are converted to the complex
*          conjugate.
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
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, II, JJ,
     $                   MP, NQ, LLDA, IAROW, IACOL, IIA, JJA, IOFFA,
     $                   IROFFA, ICOFFA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MOD
*     ..
*     .. External Functions ..
      EXTERNAL           NUMROC
      INTEGER            NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L
*
      IF ( M .EQ. 0 .OR. N .EQ. 0 ) RETURN
*
      ICTXT = DESCA( CTXT_ )
      LLDA = DESCA( LLD_ )
      IROFFA = MOD( IA-1, DESCA( MB_ ) )
      ICOFFA = MOD( JA-1, DESCA( NB_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $     IAROW, IACOL)
      MP = NUMROC( M+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IF ( MYROW .EQ. IAROW ) MP = MP - IROFFA
      IF ( MYCOL .EQ. IACOL ) NQ = NQ - ICOFFA
*
      IOFFA = ( JJA - 1 ) * LLDA
      DO JJ = JJA, JJA+NQ-1
         DO II = IIA, IIA+MP-1
            A( IOFFA+II ) = DCONJG( A( IOFFA+II ) )
         END DO
         IOFFA = IOFFA + LLDA
      END DO
*
      RETURN
*
*     End of PZLACONJ().
*
      END
