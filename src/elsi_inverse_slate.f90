! The wrapper to call slate (potentially w/ GPU acceleration) for matrix inversion
! complex 64 format
subroutine elsi_inverse_slate_c64 &
      (n, n_bb_row, n_bb_col, A, lda, nb, p, q, mpi_comm_ )
      use iso_c_binding
      use slate, only:slate_HermitianMatrix_create_fromScaLAPACK_c64, &
                      slate_chol_factor_c64, &
                      slate_chol_inverse_using_factor_c64, &
                      slate_HermitianMatrix_destroy_c64
      implicit none
      integer(kind=c_int64_t)   :: n
      integer(kind=c_int64_t)   :: n_bb_row
      integer(kind=c_int64_t)   :: n_bb_col
      complex*16, intent(INOUT) :: A(n_bb_row, n_bb_col)
      integer(kind=c_int64_t)   :: lda
      integer(kind=c_int64_t)   :: nb
      integer(kind=c_int)       :: p, q
      integer(kind=c_int)       :: mpi_comm_
      type(c_ptr)               :: slate_A, opts

      slate_A =  slate_HermitianMatrix_create_fromScaLAPACK_c64( &
                         'U', n, A, lda, &
                         nb, p, q, mpi_comm_ )
      call slate_chol_factor_c64( slate_A, 0, opts )
      call slate_chol_inverse_using_factor_c64( slate_A, 0, opts )
      call slate_HermitianMatrix_destroy_c64( slate_A )
end subroutine elsi_inverse_slate_c64

! The wrapper to call slate (potentially w/ GPU acceleration) for matrix inversion
! complex 32 format
subroutine elsi_inverse_slate_c32 &
      ( n, n_bb_row, n_bb_col, A, lda, nb, p, q, mpi_comm_ )
      use iso_c_binding
      use slate, only:slate_HermitianMatrix_create_fromScaLAPACK_c32, &
                      slate_chol_factor_c32, &
                      slate_chol_inverse_using_factor_c32, &
                      slate_HermitianMatrix_destroy_c32
      implicit none
      integer(kind=c_int64_t)   :: n
      integer(kind=c_int64_t)   :: n_bb_row
      integer(kind=c_int64_t)   :: n_bb_col
      complex*8, intent(INOUT) :: A(n_bb_row, n_bb_col)
      integer(kind=c_int64_t)   :: lda
      integer(kind=c_int64_t)   :: nb
      integer(kind=c_int)       :: p, q
      integer(kind=c_int)       :: mpi_comm_
      type(c_ptr)               :: slate_A, opts

      slate_A =  slate_HermitianMatrix_create_fromScaLAPACK_c32( &
                         'U', n, A, lda, &
                         nb, p, q, mpi_comm_ )
      call slate_chol_factor_c32( slate_A, 0, opts )
      call slate_chol_inverse_using_factor_c32( slate_A, 0, opts )
      call slate_HermitianMatrix_destroy_c32( slate_A )
end subroutine elsi_inverse_slate_c32
