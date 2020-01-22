      SUBROUTINE PZDOTC_ALL( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY,
     $                       DESCY, INCY )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, IX, JX, INCX, IY, JY, INCY
      COMPLEX*16         DOT
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * ), DESCY( * )
      COMPLEX*16         X( * ), Y( * )
*
*  Purpose
*  =======
*
*  PZDOTC_ALL forms the dot product of two subvectors,
*
*     DOT := sub( X )**H * sub( Y ),
*
*  where
*
*     sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                      X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X, and,
*
*     sub( Y ) denotes Y(IY,JY:JY+N-1) if INCY = M_Y,
*                      Y(IY:IY+N-1,JY) if INCY = 1 and INCY <> M_Y.
*
*  The result is available on all processes.
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
*          On entry,  N  specifies the  length of the  subvectors to  be
*          multiplied. N must be at least zero.
*
*  DOT     (local output) DOUBLE PRECISION array
*          On exit, DOT  specifies the dot product of the two subvectors
*          sub( X ) and sub( Y ).
*
*  X       (local input) COMPLEX* array
*          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
*          is   at  least  MAX( 1, Lr( 1, IX ) )  when  INCX = M_X   and
*          MAX( 1, Lr( 1, IX+N-1 ) )  otherwise,  and,  Kx  is  at least
*          Lc( 1, JX+N-1 )  when  INCX = M_X  and Lc( 1, JX ) otherwise.
*          Before  entry,  this  array contains the local entries of the
*          matrix X.
*
*  IX      (global input) INTEGER
*          On entry, IX  specifies X's global row index, which points to
*          the beginning of the submatrix sub( X ).
*
*  JX      (global input) INTEGER
*          On entry, JX  specifies X's global column index, which points
*          to the beginning of the submatrix sub( X ).
*
*  DESCX   (global and local input) INTEGER array
*          On entry, DESCX  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix X.
*
*  INCX    (global input) INTEGER
*          On entry,  INCX   specifies  the  global  increment  for  the
*          elements of  X.  Only two values of  INCX   are  supported in
*          this version, namely 1 and M_X. INCX  must not be zero.
*
*  Y       (local input) COMPLEX*16 array
*          On entry, Y is an array of dimension (LLD_Y, Ky), where LLD_Y
*          is   at  least  MAX( 1, Lr( 1, IY ) )  when  INCY = M_Y   and
*          MAX( 1, Lr( 1, IY+N-1 ) )  otherwise,  and,  Ky  is  at least
*          Lc( 1, JY+N-1 )  when  INCY = M_Y  and Lc( 1, JY ) otherwise.
*          Before  entry,  this array  contains the local entries of the
*          matrix Y.
*
*  IY      (global input) INTEGER
*          On entry, IY  specifies Y's global row index, which points to
*          the beginning of the submatrix sub( Y ).
*
*  JY      (global input) INTEGER
*          On entry, JY  specifies Y's global column index, which points
*          to the beginning of the submatrix sub( Y ).
*
*  DESCY   (global and local input) INTEGER array
*          On entry, DESCY  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix Y.
*
*  INCY    (global input) INTEGER
*          On entry,  INCY   specifies  the  global  increment  for  the
*          elements of  Y.  Only two values of  INCY   are  supported in
*          this version, namely 1 and M_Y. INCY  must not be zero.
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
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, R, C, RSRC,
     $                   CSRC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, ZGEBR2D, ZGEBS2D, INFOG2L,
     $                   PZDOTC
*     ..
*     .. Executable Statements ..
*
      ICTXT = DESCX( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Call PZDOTC().
*
      CALL PZDOTC( N, DOT, X, IX, JX, DESCX, INCX,
     $     Y, IY, JY, DESCY, INCY )
*
*     Broadcast the result to all processes.
*
      IF ( NPROW*NPCOL .GT. 1 ) THEN
*
*        Find the first entry of X.
*
         CALL INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL, R, C,
     $        RSRC, CSRC )
         IF ( MYROW .EQ. RSRC .AND. MYCOL .EQ. CSRC ) THEN
            CALL ZGEBS2D( ICTXT, 'A', ' ', 1, 1, DOT, 1, RSRC, CSRC )
         ELSE
            CALL ZGEBR2D( ICTXT, 'A', ' ', 1, 1, DOT, 1, RSRC, CSRC )
         END IF
      END IF
*
      RETURN
*
*     End of PZDOTC_ALL().
*
      END
