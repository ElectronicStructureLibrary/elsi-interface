!
!     The binary representation of SOLVER is of the form
!        SOLVER = AAAABBBBCCCCDDDD
!     AAAA = 0000: Full Bethe--Salpeter solver
!     AAAA = 0001: Tamm--Dancoff approximation (TDA)
!     BBBB = 0000: Full diagonalization
!     BBBB = 0001: Lanczos algorithm
!     CCCC = 0000: delta function is replaced by Gaussian
!     CCCC = 0001: delta function is replaced by Lorentzian
!     DDDD:        Algorithmic variant
!
      INTEGER            BSE_OFFDIAG, BSE_OFFDIAG_SHIFT
      PARAMETER          ( BSE_OFFDIAG_SHIFT = 4096 )
      PARAMETER          ( BSE_OFFDIAG = BSE_OFFDIAG_SHIFT*15 )
      INTEGER            BSE_FULLBSE, BSE_TDA
      PARAMETER          ( BSE_FULLBSE = BSE_OFFDIAG_SHIFT*0 )
      PARAMETER          ( BSE_TDA = BSE_OFFDIAG_SHIFT*1 )
!
      INTEGER            BSE_ALGORITHM, BSE_ALGORITHM_SHIFT
      PARAMETER          ( BSE_ALGORITHM_SHIFT = 256 )
      PARAMETER          ( BSE_ALGORITHM = BSE_ALGORITHM_SHIFT*15 )
      INTEGER            BSE_DIRECT, BSE_LANCZOS
      PARAMETER          ( BSE_DIRECT = BSE_ALGORITHM_SHIFT*0 )
      PARAMETER          ( BSE_LANCZOS = BSE_ALGORITHM_SHIFT*1 )
!
      INTEGER            BSE_DELTA, BSE_DELTA_SHIFT
      PARAMETER          ( BSE_DELTA_SHIFT = 16 )
      PARAMETER          ( BSE_DELTA = BSE_DELTA_SHIFT*15 )
      INTEGER            BSE_GAUSSIAN, BSE_LORENTZIAN
      PARAMETER          ( BSE_GAUSSIAN = BSE_DELTA_SHIFT*0 )
      PARAMETER          ( BSE_LORENTZIAN = BSE_DELTA_SHIFT*1 )
!
      INTEGER            BSE_VARIANT, BSE_VARIANT_SHIFT
      PARAMETER          ( BSE_VARIANT_SHIFT = 1 )
      PARAMETER          ( BSE_VARIANT = BSE_VARIANT_SHIFT*15 )
!     BSE_PRODUCT only works with BSE_FULLBSE+BSE_DIRECT
      INTEGER            BSE_PRODUCT
      PARAMETER          ( BSE_PRODUCT = BSE_VARIANT_SHIFT*0 )
!     BSE_SVD only works with BSE_FULLBSE+BSE_DIRECT
      INTEGER            BSE_SVD
      PARAMETER          ( BSE_SVD = BSE_VARIANT_SHIFT*9 )
!     BSE_LAPACK_HEEVR only works with BSE_TDA+BSE_DIRECT
      INTEGER            BSE_LAPACK_HEEVR
      PARAMETER          ( BSE_LAPACK_HEEVR = BSE_VARIANT_SHIFT*0 )
!     BSE_LAPACK_HEEV only works with BSE_TDA+BSE_DIRECT
      INTEGER            BSE_LAPACK_HEEV
      PARAMETER          ( BSE_LAPACK_HEEV = BSE_VARIANT_SHIFT*1 )
!     BSE_LAPACK_HEEVD only works with BSE_TDA+BSE_DIRECT
      INTEGER            BSE_LAPACK_HEEVD
      PARAMETER          ( BSE_LAPACK_HEEVD = BSE_VARIANT_SHIFT*2 )
!     BSE_LAPACK_HEEVX only works with BSE_TDA+BSE_DIRECT
      INTEGER            BSE_LAPACK_HEEVX
      PARAMETER          ( BSE_LAPACK_HEEVX = BSE_VARIANT_SHIFT*3 )
!     BSE_QUADAVGGAUSS only works with BSE_LANCZOS
      INTEGER            BSE_QUADAVGGAUSS
      PARAMETER          ( BSE_QUADAVGGAUSS = BSE_VARIANT_SHIFT*1 )
!
      INTEGER            BSE_UNDEFINED
      PARAMETER          ( BSE_UNDEFINED = -1 )
