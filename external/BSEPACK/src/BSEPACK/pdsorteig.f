      SUBROUTINE PDSORTEIG( N, M, LAMBDA, V, IV, JV, DESCV, MODE )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, M, IV, JV, MODE
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   LAMBDA( * ), V( * )
      INTEGER            DESCV( * )
*     ..
*
*  Purpose
*  =======
*
*  PDSORTEIG() is an auxiliary routine.
*  It sorts the eigenvalues stored in lambda and the corresponding
*  eigenvectors stored in V.
*
*  On exit,
*  if MODE = 0, then
*     lambda( 1 ) <= lambda( 2 ) <= ... <= lambda( m );
*  if MODE = 1, then
*     lambda( 1 ) >= lambda( 2 ) >= ... >= lambda( m ).
*
*  V( :, j ) always contains the eigenvector corresponds to lambda( j ).
*
*  No argument check is performed in this subroutine, i.e.,
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
*          The number of rows of the distributed submatrix sub( V ).
*          N >= 0.
*
*  M       (global input) INTEGER
*          The number of columns of the distributed submatrix sub( V ),
*          as well as the length oF lambda.
*          M >= 0.
*
*  LAMBDA  (global input and output) DOUBLE PRECISION array,
*          dimension (M).
*          The entries of lambda must be consistent among all
*          processors.
*
*  V       (local input and output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_V, LOCc(JV+N-1)).
*          This array contains the local pieces of the N-by-M
*          distributed eigenvector matrix sub( V ) to be sorted
*          according to the value of lambda.
*
*  IV      (global input) INTEGER
*          The row index in the global array V indicating the first
*          row of sub( V ).
*
*  JV      (global input) INTEGER
*          The column index in the global array V indicating the
*          first column of sub( V ).
*
*  DESCV   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix V.
*
*  MODE    (global input) INTEGER
*          = 0:  LAMBDA is sorted in ascending order;
*          = 1:  LAMBDA is sorted in descending order.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J, JJ
      DOUBLE PRECISION   TMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDSWAP
*
      IF ( MODE .EQ. 0 ) THEN
         DO J = 1, M
            I = J
            TMP = LAMBDA( J )
            DO JJ = J + 1, M
               IF ( LAMBDA( JJ ) .LT. TMP ) I = JJ
            END DO
            IF ( I .NE. J ) THEN
               LAMBDA( J ) = LAMBDA( I )
               LAMBDA( I ) = TMP
               CALL PDSWAP( N, V, IV, JV-1+J, DESCV, 1, V, IV, JV-1+I,
     $              DESCV, 1)
            END IF
         END DO
      ELSE
         DO J = 1, M
            I = J
            TMP = LAMBDA( J )
            DO JJ = J + 1, M
               IF ( LAMBDA( JJ ) .GT. TMP ) I = JJ
            END DO
            IF ( I .NE. J ) THEN
               LAMBDA( J ) = LAMBDA( I )
               LAMBDA( I ) = TMP
               CALL PDSWAP( N, V, IV, JV-1+J, DESCV, 1, V, IV, JV-1+I,
     $              DESCV, 1)
            END IF
         END DO
      END IF
*
      RETURN
*
*     End of PDSORTEIG().
*
      END
