      SUBROUTINE PZBSAUX1( N, V, IV, JV, DESCV, X, IX, JX, DESCX, Z, IZ,
     $                     JZ, DESCZ, ZR, IZR, JZR, DESCZR, ZI, IZI,
     $                     JZI, DESCZI )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, IV, JV, IX, JX, IZ, JZ, IZR, JZR, IZI, JZI
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   V( * ), Z( * ), ZR( * ), ZI( * )
      COMPLEX*16         X( * )
      INTEGER            DESCV( * ), DESCX( * ), DESCZ( * ),
     $                   DESCZR( * ), DESCZI( * )
*     ..
*
*  Purpose
*  =======
*
*  PZBSAUX1() is an auxiliary routine called by PZBSSOLVER1().
*  It applies three 2n-by-2n or 2n-by-n unitary transformations
*
*     X = Q * Z * D * V,
*
*  where D = diag{1,i,i^2,i^3,...},
*
*        Q = 1/sqrt(2) * [ I_n, -i*I_n;
*                          I_n,  i*I_n ],
*
*  and V is a given 2n-by-n real orthogonal matrix.
*  On entry, Z has already been prescaled by the factor 1/sqrt(2) in Q.
*
*  ZR and ZI are provided as workspace to store the real and imaginary
*  parts of Q * Z * D, respectively, i.e.,
*
*     Q * Z * D = ZR + i*ZI.
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
*          2*N is the number of rows and columns to be operated on, i.e.
*          the order of the distributed submatrices sub( V ), sub( X ),
*          sub( XR ), and sub( XI ).
*          N >= 0.
*
*  V       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_V, LOCc(JV+N-1)).
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
*  X       (local output) COMPLEX*16 pointer into the local memory to an
*          array of dimension (LLD_X, LOCc(JX+N-1)).
*
*  IX      (global input) INTEGER
*          The row index in the global array X indicating the first
*          row of sub( X ).
*          In this version, only IX = JX = 1 is supported.
*
*  JX      (global input) INTEGER
*          The column index in the global array X indicating the
*          first column of sub( X ).
*          In this version, only IX = JX = 1 is supported.
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*  Z       (local input and output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension
*          (LLD_Z, LOCc(JZ+2*N-1)).
*
*  IZ      (global input) INTEGER
*          The row index in the global array Z indicating the first
*          row of sub( Z ).
*          In this version, only IZ = JZ = 1 is supported.
*
*  JZ      (global input) INTEGER
*          The column index in the global array Z indicating the
*          first column of sub( Z ).
*          In this version, only IZ = JZ = 1 is supported.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*  ZR      (local output) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_ZR, LOCc(JZR+2*N-1)).
*
*  IZR     (global input) INTEGER
*          The row index in the global array ZR indicating the first
*          row of sub( ZR ).
*
*  JZR     (global input) INTEGER
*          The column index in the global array ZR indicating the
*          first column of sub( ZR ).
*
*  DESCZR  (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix ZR.
*
*  ZI      (local output) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_ZI, LOCc(JZI+2*N-1)).
*
*  IZI     (global input) INTEGER
*          The row index in the global array ZI indicating the first
*          row of sub( ZI ).
*
*  JZI     (global input) INTEGER
*          The column index in the global array ZI indicating the
*          first column of sub( ZI ).
*
*  DESCZI  (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix ZI.
*
*  Alignment requirements
*  ======================
*
*  This subroutine requires X and Z to be distributed identically in the
*  sense that DESCX( : ) = DESCZ( : ), despite that X is complex and Z
*  is real.
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
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, NB,
     $                   I, J, XROWS, XCOLS, LLDX, IX0, JX0, RSRC, CSRC
      DOUBLE PRECISION   T, T_VEC_O1, T_VEC_O2, T_VEC_O3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MOD
*     ..
*     .. External Functions ..
      EXTERNAL           NUMROC, MPI_WTIME
      INTEGER            NUMROC
      DOUBLE PRECISION   MPI_WTIME
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, PDCOPY, PDGEMM, PDSCAL
*     ..
*     .. Executable Statements ..
*
*     Apply Q from the left and D from the right.
*
      ICTXT = DESCX( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NB = DESCX( NB_ )
      LLDX = DESCX( LLD_ )
      XROWS = NUMROC( 2*N, NB, MYROW, 0, NPROW )
      XCOLS = NUMROC( N, NB, MYCOL, 0, NPCOL )
      CALL INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL, IX0, JX0,
     $     RSRC, CSRC )
*
      T_VEC_O2 = MPI_WTIME()
      DO J = 0, 2*N-1
         IF ( MOD( J, 4 ) .EQ. 0 ) THEN
            CALL PDCOPY( N, Z, IZ, JZ + J, DESCZ, 1, ZR, IZR, JZR + J,
     $           DESCZR, 1 )
            CALL PDCOPY( N, Z, IZ, JZ + J, DESCZ, 1, ZR, IZR + N,
     $           JZR + J, DESCZR, 1 )
            CALL PDCOPY( N, Z, IZ + N, JZ + J, DESCZ, 1, ZI, IZI,
     $           JZI + J, DESCZI, 1 )
            CALL PDSCAL( N, -ONE, ZI, IZI, JZI + J, DESCZI, 1 )
            CALL PDCOPY( N, Z, IZ + N, JZ + J, DESCZ, 1, ZI, IZI + N,
     $           JZI + J, DESCZI, 1 )
         ELSE IF ( MOD( J, 4 ) .EQ. 1 ) THEN
            CALL PDCOPY( N, Z, IZ, JZ + J, DESCZ, 1, ZI, IZI, JZI + J,
     $           DESCZI, 1 )
            CALL PDCOPY( N, Z, IZ, JZ + J, DESCZ, 1, ZI, IZI + N,
     $           JZI + J, DESCZI, 1 )
            CALL PDCOPY( N, Z, IZ + N, JZ + J, DESCZ, 1, ZR, IZR,
     $           JZR + J, DESCZR, 1 )
            CALL PDCOPY( N, Z, IZ + N, JZ + J, DESCZ, 1, ZR, IZR + N,
     $           JZR + J, DESCZR, 1 )
            CALL PDSCAL( N, -ONE, ZR, IZR + N, JZR + J, DESCZR, 1 )
         ELSE IF ( MOD( J, 4 ) .EQ. 2 ) THEN
            CALL PDCOPY( N, Z, IZ, JZ + J, DESCZ, 1, ZR, IZR, JZR + J,
     $           DESCZR, 1 )
            CALL PDSCAL( N, -ONE, ZR, IZR, JZR + J, DESCZR, 1 )
            CALL PDCOPY( N, Z, IZ, JZ + J, DESCZ, 1, ZR, IZR + N,
     $           JZR + J, DESCZR, 1 )
            CALL PDSCAL( N, -ONE, ZR, IZR + N, JZR + J, DESCZR, 1 )
            CALL PDCOPY( N, Z, IZ + N, JZ + J, DESCZ, 1, ZI, IZI,
     $           JZI + J, DESCZI, 1 )
            CALL PDCOPY( N, Z, IZ + N, JZ + J, DESCZ, 1, ZI, IZI + N,
     $           JZI + J, DESCZI, 1 )
            CALL PDSCAL( N, -ONE, ZI, IZI + N, JZI + J, DESCZI, 1 )
         ELSE !IF ( MOD( J, 4 ) .EQ. 3 ) THEN
            CALL PDCOPY( N, Z, IZ, JZ + J, DESCZ, 1, ZI, IZI, JZI + J,
     $           DESCZI, 1 )
            CALL PDSCAL( N, -ONE, ZI, IZI, JZI + J, DESCZI, 1 )
            CALL PDCOPY( N, Z, IZ, JZ + J, DESCZ, 1, ZI, IZI + N,
     $           JZI + J, DESCZI, 1 )
            CALL PDSCAL( N, -ONE, ZI, IZI + N, JZI + J, DESCZI, 1 )
            CALL PDCOPY( N, Z, IZ + N, JZ + J, DESCZ, 1, ZR, IZR,
     $           JZR + J, DESCZR, 1 )
            CALL PDSCAL( N, -ONE, ZR, IZR, JZR + J, DESCZR, 1 )
            CALL PDCOPY( N, Z, IZ + N, JZ + J, DESCZ, 1, ZR, IZR + N,
     $           JZR + J, DESCZR, 1 )
         END IF
      END DO
      T_VEC_O2 = MPI_WTIME() - T_VEC_O2
*
*     Apply V from the right.
*
      T = MPI_WTIME()
      CALL PDGEMM( 'N', 'N', 2*N, N, 2*N, ONE, ZR, IZR, JZR, DESCZR,
     $     V, IV, JV, DESCV, ZERO, Z, IZ, JZ, DESCZ )
      T_VEC_O3 = MPI_WTIME() - T
      T = MPI_WTIME()
      DO J = IX0, XCOLS
         DO I = JX0, XROWS
            X( I + (J-1)*LLDX ) = DCMPLX( Z( I + (J-1)*LLDX ) )
         END DO
      END DO
      T_VEC_O1 = MPI_WTIME() - T
*
      T = MPI_WTIME()
      CALL PDGEMM( 'N', 'N', 2*N, N, 2*N, ONE, ZI, IZI, JZI, DESCZI,
     $     V, IV, JV, DESCV, ZERO, Z, IZ, JZ, DESCZ )
      T_VEC_O3 = T_VEC_O3 + MPI_WTIME() - T
      T = MPI_WTIME()
      DO J = IX0, XCOLS
         DO I = JX0, XROWS
            X( I + (J-1)*LLDX ) = DCMPLX( DBLE( X( I + (J-1)*LLDX ) ),
     $           Z( I + (J-1)*LLDX ) )
         END DO
      END DO
      T_VEC_O1 = T_VEC_O1 + MPI_WTIME() - T
!      IF ( MYROW+MYCOL .EQ. 0 ) THEN
!         WRITE( *, * ) 't_vec_o1 =', T_VEC_O1, ';'
!         WRITE( *, * ) 't_vec_o2 =', T_VEC_O2, ';'
!         WRITE( *, * ) 't_vec_o3 =', T_VEC_O3, ';'
!      END IF
*
      RETURN
*
*     End of ZBSAUX1().
*
      END
