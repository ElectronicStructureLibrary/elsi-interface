      SUBROUTINE PDBSABSP( SOLVER, N, NPTS, SIGMA, OMEGA, EPS, A, IA,
     $                     JA, DESCA, B, IB, JB, DESCB, LAMBDA, X, IX,
     $                     JX, DESCX, D, ID, JD, DESCD, ALPHA, BETA,
     $                     RESUME, ITMAX, WORK, LWORK, IWORK, LIWORK,
     $                     INFO )
*
      IMPLICIT NONE
      INCLUDE 'solver.f'
*
*     .. Scalar Arguments ..
      INTEGER            SOLVER, N, NPTS, IA, JA, IB, JB, IX, JX, ID,
     $                   JD, RESUME, ITMAX, LWORK, LIWORK, INFO
      DOUBLE PRECISION   SIGMA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   OMEGA( * ), EPS( * ), LAMBDA( * ), A( * ),
     $                   B( * ), X( * ), D( * ), ALPHA( * ), BETA( * ),
     $                   WORK( * )
      INTEGER            DESCA( * ), DESCB( * ), DESCX( * ), DESCD( * ),
     $                   IWORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDBSABSP() computes the optical absorption spectrum
*
*     epsilon( omega ) = d_r**T * delta( omega * I - H ) * d_l,
*
*  where
*
*     H = [ A,  B;
*          -B, -A ],
*
*  is the Behte--Salpeter Hamiltonian,
*
*     d_r = [ d_1**T, -d_1**T ]**T
*
*  and
*
*     d_l = [ d_1**T, d_1**T ]**T
*
*  are right and left dipole vectors, respectively.
*
*  The matrix
*
*     [ A, B;
*       B, A ]
*
*  is assumed to be positive definite.
*
*  The delta function is replaced by, e.g., the Gaussian function
*
*     g( x ) = exp( -x**2/( 2*sigma**2 ) )/( sqrt( 2*pi )*sigma ).
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
*  SOLVER  (global input) INTEGER
*          SOLVER specifies the solver to be applied;
*          see solver.f for details.
*          Currently supported solvers:
*          SOLVER = S1 + S2 + S3 + S4, where
*             S1 \in { BSE_FULL, BSE_TDA },
*             S2 \in { BSE_DIRECT, BSE_LANCZOS },
*             S3 \in { BSE_GAUSSIAN, BSE_LORENTZIAN },
*             if S2 = BSE_DIRECT, see PDBSEIG for supported variants,
*             if S2 = BSE_LANCZOS, S4 \in { 0, BSE_QUADAVGGAUSS }.
*
*  N       (global input) INTEGER
*          The number of rows and columns of A and B.
*          N >= 0.
*
*  NPTS    (global input) INTEGER
*          NPTS is the number of sampling points in OMEGA.
*          NPTS >= 0.
*
*  SIGMA   (global input) DOUBLE PRECISION
*          Broadening factor in the approximation of the delta function.
*          SIGMA > 0.
*
*  OMEGA   (global input) DOUBLE PRECISION array, dimension (NPTS)
*          Sampling points of omega.
*          All entries of OMEGA are nonnegative.
*
*  EPS     (global output) DOUBLE PRECISION array, dimension (NPTS)
*          Sampling points of epsilon, i.e.,
*          EPS( I ) = epsilon( OMEGA( I ) ).
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
*          B is not referenced if TDA is used.
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
*  LAMBDA  (global output) DOUBLE PRECISION array
*          LAMBDA is used to store the eigenvalues or Ritz values.
*          If IAND( SOLVER, BSE_ALGORITHM ) = BSE_DIRECT,
*             LAMBDA contains N positive eigenvalues of H;
*          if IAND( SOLVER, BSE_ALGORITHM ) = BSE_LANCZOS,
*             LAMBDA contains MIN( N, ITMAX ) Ritz values.
*
*  X       (local input/output) DOUBLE PRECISION pointer into the local
*          memory to an array of dimension (LLD_X, LOCc(JX+1)).
*          X is used to store the eigenvectors or Lanczos vectors.
*          X has 2N rows for full BSE, and has N rows for TDA.
*          If IAND( SOLVER, BSE_ALGORITHM ) = BSE_DIRECT,
*             X has 2*N columns;
*          if IAND( SOLVER, BSE_ALGORITHM ) = BSE_LANCZOS,
*             X has ITER + 1 columns, where ITER is the number of
*             Lanczos steps.
*          If m = RESUME > 0, the first m columns of X stores the
*          previously calculated eigenvectors or Lanczos vectors.
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
*  D       (local iutput) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_D, LOCc(JD+1)).
*          Only the first block (i.e., d_1) is stored in D.
*
*  ID      (global input) INTEGER
*          The row index in the global array D indicating the first
*          row of sub( D ).
*
*  JD      (global input) INTEGER
*          The column index in the global array D indicating the
*          first column of sub( D ).
*
*  DESCD   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix D.
*
*  ALPHA   (global input/output) DOUBLE PRECISION array,
*          dimension (ITMAX)
*          If IAND( SOLVER, BSE_ALGORITHM ) = BSE_DIRECT,
*             ALPHA is not referenced.
*          if IAND( SOLVER, BSE_ALGORITHM ) = BSE_LANCZOS,
*             On exit, ALPHA contains the diagonal entries of the
*             tridiagonal matrix.
*             If m = RESUME > 0, the first m entries of ALPHA stores the
*             previously calculated results on entry.
*
*  BETA    (global input/output) DOUBLE PRECISION array,
*          dimension (ITMAX)
*          If IAND( SOLVER, BSE_ALGORITHM ) = BSE_DIRECT,
*             BETA is not referenced;
*          if IAND( SOLVER, BSE_ALGORITHM ) = BSE_LANCZOS,
*             On exit, BETA contains the subdiagonal entries of the
*             tridiagonal matrix.
*             If m = RESUME > 0, the first m entries of BETA stores the
*             previously calculated results on entry.
*
*  RESUME  (global input) INTEGER
*          If m = resume > 0, the solver reuses the previously
*          calculated spectral decomposition (for full diagonalization)
*          or m-step Lanczos decomposition (for Lanczos).
*          In this version only m = 0 is supported.
*
*  ITMAX   (global input) INTEGER
*          If IAND( SOLVER, BSE_ALGORITHM ) = BSE_DIRECT,
*             BETA is not referenced;
*          if IAND( SOLVER, BSE_ALGORITHM ) = BSE_LANCZOS,
*            The maximum number of Lanczos steps.
*
*  WORK    (local workspace/output) DOUBLE PRECISION array,
*          dimension (LWORK)
*          On output, WORK( 1 ) returns the minimal amount of workspace
*          needed to guarantee completion.
*          If the input parameters are incorrect, WORK( 1 ) may also be
*          incorrect.
*
*  LWORK   (local input) INTEGER
*          The length of the workspace array WORK.
*          If LWORK = -1, LWORK is global input and a workspace query is
*          assumed; the routine only calculates the minimum size for the
*          WORK/IWORK array. The required workspace is returned as the
*          first element of WORK/IWORK and no error message is issued by
*          PXERBLA.
*
*  IWORK   (local workspace/output) INTEGER array,
*          dimension (LIWORK)
*          On output, IWORK( 1 ) returns the minimal amount of workspace
*          needed to guarantee completion.
*          If the input parameters are incorrect, IWORK( 1 ) may also be
*          incorrect.
*
*  LIWORK  (local input) INTEGER
*          The length of the workspace array IWORK.
*          If LIWORK = -1, LIWORK is global input and a workspace query
*          is assumed; the routine only calculates the minimum size for
*          the WORK/IWORK array. The required workspace is returned as
*          the first element of WORK/IWORK and no error message is
*          issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  Lanczos process breaks down after INFO steps.
*
*  Alignment requirements
*  ======================
*
*  The distributed matrices A, B, X, and D must satisfy the following
*  properties:
*
*     DESCA( MB_ ) = DESCA( NB_ )
*     DESCB( MB_ ) = DESCB( NB_ )
*     DESCA( MB_ ) = DESCB( MB_ ) = DESCX( MB_ ) = DESCD( MB_ )
*     DESCA( NB_ ) = DESCB( NB_ ) = DESCX( NB_ )
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
      LOGICAL            LQUERY
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, NB, LLDW,
     $                   WROWS, LWKOPT, LIWKOPT, DIV, I, J, KK,
     $                   STEPS, IERR, MROWS, MCOLS, LLDM, INDM, INDK,
     $                   INDALPHA, INDBETA, INDZ, INDWORK
      DOUBLE PRECISION   NORMD2
*     ..
*     .. Local Arrays ..
      INTEGER            DESCM( DLEN_ ), DESCK( DLEN_ ), DESCW( DLEN_ )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DSQRT, IAND, MAX, MIN
*     ..
*     .. External Functions ..
      EXTERNAL           NUMROC, APPROX_DELTA, PDLAMCH
      INTEGER            NUMROC
      DOUBLE PRECISION   APPROX_DELTA, PDLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PXERBLA, DESCSET,
     $                   PDEMBED1, PDBSEIG, PDGEMV, PDELGET,
     $                   PDCOPY, PDBSLANCZOS, PDLANCZOS,
     $                   DSTEQR, DBSABSP_AUX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      ICTXT = DESCD( CTXT_ )
      NB = DESCD( MB_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      IF ( NPROW .EQ. -1 ) THEN
         INFO = -( 1000+CTXT_ )
      END IF
*
*     Test the input arguments.
*
      LQUERY = LWORK .EQ. -1 .OR. LWORK .EQ. -1
      IF ( INFO .EQ. 0 .AND. SOLVER .LT. 0 )
     $   INFO = -1
      IF ( INFO .EQ. 0 .AND. N .LT. 0 )
     $   INFO = -2
      IF ( INFO .EQ. 0 .AND. NPTS .LT. 0 )
     $   INFO = -3
      IF ( INFO .EQ. 0 .AND. .NOT. SIGMA .GT. 0 )
     $   INFO = -4
      IF ( INFO .EQ. 0 )
     $   CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 10, INFO )
      IF ( INFO .EQ. 0 .AND. DESCA( MB_ ) .NE. DESCA( NB_ ) )
     $   INFO = -( 1000+MB_ )
      IF ( INFO .EQ. 0 )
     $   CALL CHK1MAT( N, 2, N, 2, IB, JB, DESCB, 14, INFO )
      IF ( INFO .EQ. 0 .AND. DESCB( MB_ ) .NE. DESCB( NB_ ) )
     $   INFO = -( 1400+MB_ )
c      IF ( INFO .EQ. 0 )
c     $   CALL CHK1MAT( 2*N, 2, 1, 0, IX, JX, DESCX, 19, INFO )
      IF ( INFO .EQ. 0 )
     $   CALL CHK1MAT( N, 2, 1, 0, ID, JD, DESCD, 23, INFO )
c      IF ( INFO .EQ. 0 .AND. RESUME .LT. 0 )
c     $   INFO = -26
      IF ( INFO .EQ. 0 .AND.
     $     IAND( SOLVER, BSE_ALGORITHM ) .NE. BSE_DIRECT
     $     .AND. ITMAX .LT. 0 )
     $   INFO = -27
*
*     Compute required workspace.
*
      IF ( INFO .EQ. 0 ) THEN
         IF ( IAND( SOLVER, BSE_ALGORITHM ) .EQ. BSE_DIRECT ) THEN
            CALL PDBSEIG( SOLVER, N, A, IA, JA, DESCA, B, IB, JB, DESCB,
     $           LAMBDA, X, IX, JX, DESCX, WORK, -1, IWORK, -1, INFO )
            LWKOPT = INT( WORK( 1 ) )
            LIWKOPT = IWORK( 1 )
            WROWS = NUMROC( N, NB, MYROW, 0, NPROW )
            LLDW = MAX( WROWS, 1 )
            CALL DESCSET( DESCW, N, N, NB, NB, 0, 0, ICTXT, LLDW )
            INDWORK = N + 1
            LWKOPT = MAX( LWKOPT, N + LLDW )
         ELSE IF ( IAND( SOLVER, BSE_ALGORITHM ) .EQ. BSE_LANCZOS ) THEN
            LWKOPT = 1
            LIWKOPT = 1
            STEPS = MIN( N, ITMAX )
            IF ( IAND( SOLVER, BSE_QUADAVGGAUSS ) .NE. 0 ) THEN
               KK = 2*STEPS - 1
            ELSE
               KK = STEPS
            END IF
            MROWS = NUMROC( N, NB, MYROW, 0, NPROW )
            MCOLS = NUMROC( N, NB, MYCOL, 0, NPCOL )
            LLDM = MAX( 1, MROWS )
            CALL DESCSET( DESCM, N, N, NB, NB, 0, 0, ICTXT, LLDM )
            CALL DESCSET( DESCK, N, N, NB, NB, 0, 0, ICTXT, LLDM )
            INDALPHA = 1
            INDBETA = INDALPHA + KK
            INDM = INDBETA + KK
            INDK = INDM + LLDM*MCOLS
            LWKOPT = INDK + LLDM*MCOLS - 1
            INDZ = INDM
            INDWORK = INDZ + KK*KK
            LWKOPT = MAX( LWKOPT, MAX( 1, 2*KK - 2 ) + INDWORK - 1 )
         END IF
         IF ( .NOT. LQUERY .AND. LWORK .LT. LWKOPT )
     $      INFO = -29
         IF ( .NOT. LQUERY .AND. LIWORK .LT. LIWKOPT )
     $      INFO = -31
      END IF
*
      IF ( INFO .NE. 0 ) THEN
         CALL PXERBLA( ICTXT, 'PDBSABSP', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      WORK( 1 ) = DBLE( LWKOPT )
      IWORK( 1 ) = LIWKOPT
      IF ( LQUERY ) RETURN
*
      IF ( N .EQ. 0 ) RETURN
      IF ( NPTS .EQ. 0 .AND. RESUME .EQ. 0 ) RETURN
*
*     Compute the optical absorption spectrum.
*
      IF ( IAND( SOLVER, BSE_ALGORITHM ) .EQ. BSE_DIRECT ) THEN
*
*        Fully diagonalize the BSE Hamiltonian.
*
         CALL PDBSEIG( SOLVER, N, A, IA, JA, DESCA, B, IB, JB, DESCB,
     $        LAMBDA, X, IX, JX, DESCX, WORK, LWORK, IWORK, LIWORK,
     $        INFO )
*
*        Compute the absorption spectrum.
*
         IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_FULLBSE ) THEN
            CALL PDGEMV( 'T', N, N, ONE, X, IX+N, JX, DESCX, D, ID,
     $           JD, DESCD, 1, ZERO, WORK( INDWORK ), 1, 1, DESCW, 1 )
            CALL PDGEMV( 'T', N, N, ONE, X, IX, JX, DESCX, D, ID, JD,
     $           DESCD, 1, -ONE, WORK( INDWORK ), 1, 1, DESCW, 1 )
         ELSE IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_TDA ) THEN
            CALL PDGEMV( 'T', N, N, ONE, X, IX, JX, DESCX, D, ID, JD,
     $           DESCD, 1, ZERO, WORK( INDWORK ), 1, 1, DESCW, 1 )
         END IF
         DO I = 1, N
            CALL PDELGET( 'A', ' ', WORK( I ), WORK( INDWORK ), I, 1,
     $           DESCW )
         END DO
         CALL DBSABSP_AUX( IAND( SOLVER, BSE_DELTA ), N, NPTS,
     $        SIGMA, OMEGA, EPS, ONE, WORK, 1, LAMBDA, 0 )
*
      ELSE IF ( IAND( SOLVER, BSE_ALGORITHM ) .EQ. BSE_LANCZOS ) THEN
*
*        Lanczos algorithm.
*
         CALL PDCOPY( N, D, ID, JD, DESCD, 1, X, IX, JX, DESCX, 1 )
         IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_FULLBSE ) THEN
*
*           Lanczos for full BSE.
*
            CALL PDEMBED1( N, A, IA, JA, DESCA, B, IB, JB, DESCB,
     $           WORK( INDM ), 1, 1, DESCM, WORK( INDK ), 1, 1, DESCK )
            CALL PDBSLANCZOS( N, STEPS, WORK( INDM ), 1, 1, DESCM,
     $           WORK( INDK ), 1, 1, DESCK, X, IX, JX, DESCX,
     $           ALPHA, BETA, NORMD2, INFO )
*
         ELSE IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_TDA ) THEN
*
*           Lanczos for TDA.
*
            CALL PDLANCZOS( N, STEPS, A, IA, JA, DESCA, X, IX, JX,
     $           DESCX, ALPHA, BETA, NORMD2, INFO )
*
         END IF
         IF ( INFO .GT. 0 ) STEPS = INFO
*
*        Expand the tridiagonal matrix for averaged Gauss quadrature.
*
         DO I = 1, STEPS
            WORK( INDALPHA + I - 1 ) = ALPHA( I )
            WORK( INDBETA + I - 1 ) = BETA( I )
         END DO
         IF ( IAND( SOLVER, BSE_QUADAVGGAUSS ) .NE. 0 ) THEN
            KK = 2*STEPS - 1
            DO J = 0, STEPS - 2
               WORK( INDALPHA + 2*STEPS - J - 2 )
     $              = WORK( INDALPHA + J )
            END DO
            DO J = 0, STEPS - 3
               WORK( INDBETA + 2*STEPS - J - 3 )
     $              = WORK( INDBETA + J )
            END DO
         ELSE
            KK = STEPS
         END IF
*
*        Diagonalize the projected tridiagonal matrix.
*
         CALL DSTEQR( 'I', KK, WORK( INDALPHA ), WORK( INDBETA ),
     $        WORK( INDZ ), KK, WORK( INDWORK ), IERR )
         IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_FULLBSE ) THEN
            DO J = 0, KK - 1
               WORK( INDALPHA + J ) = DSQRT( WORK( INDALPHA + J ) )
            END DO
         END IF
         IF ( IAND( SOLVER, BSE_QUADAVGGAUSS ) .NE. 0 ) THEN
            DO J = 0, STEPS - 1
               LAMBDA( J + 1 ) = WORK( INDALPHA + 2*J )
            END DO
         ELSE
            DO J = 0, STEPS - 1
               LAMBDA( J + 1 ) = WORK( INDALPHA + J )
            END DO
         END IF
*
*        Compute the absorption spectrum.
*
         IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_FULLBSE ) THEN
            DIV = 1
         ELSE IF ( IAND( SOLVER, BSE_OFFDIAG ) .EQ. BSE_TDA ) THEN
            DIV = 0
         END IF
         CALL DBSABSP_AUX( IAND( SOLVER, BSE_DELTA ), KK, NPTS,
     $        SIGMA, OMEGA, EPS, NORMD2, WORK( INDZ ), KK,
     $        WORK( INDALPHA ), DIV )
*
      END IF
*
      WORK( 1 ) = DBLE( LWKOPT )
      IWORK( 1 ) = LIWKOPT
*
      RETURN
*
*     End of PDBSABSP().
*
      END
