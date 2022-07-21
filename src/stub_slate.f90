module slate

contains

   function slate_HermitianMatrix_create_fromScaLAPACK_c64(uplo, n, A, lda, nb, p, q, mpi_comm_) result(ret_val)
        use iso_c_binding
        implicit none
        character(kind=c_char), value :: uplo
        integer(kind=c_int64_t) :: n
        complex(kind=c_double_complex), target :: A(*)
        integer(kind=c_int64_t) :: lda
        integer(kind=c_int64_t) :: nb
        integer(kind=c_int) :: p
        integer(kind=c_int) :: q
        integer(kind=c_int) :: mpi_comm_
        type(c_ptr) :: ret_val


        write(*,"(A)") "**Error! A SLATE stub routine was called"
        stop

   end function

   subroutine slate_chol_factor_c64(A, num_opts, opts)
        use iso_c_binding
        implicit none
        type(c_ptr) :: A
        integer(kind=c_int) :: num_opts
        type(c_ptr), target :: opts


        write(*,"(A)") "**Error! A SLATE stub routine was called"
        stop

   end subroutine

   subroutine slate_chol_inverse_using_factor_c64(A, num_opts, opts)
       use iso_c_binding
       implicit none
       type(c_ptr) :: A
       integer(kind=c_int) :: num_opts
       type(c_ptr), target :: opts

       write(*,"(A)") "**Error! A SLATE stub routine was called"
       stop

   end subroutine


   subroutine slate_HermitianMatrix_destroy_c64(A)
       use iso_c_binding
       implicit none
       type(c_ptr) :: A

       write(*,"(A)") "**Error! A SLATE stub routine was called"
       stop

   end subroutine


   function slate_HermitianMatrix_create_fromScaLAPACK_c32(uplo, n, A, lda, nb, p, q, mpi_comm_) result(ret_val)
      use iso_c_binding
      implicit none
      character(kind=c_char), value :: uplo
      integer(kind=c_int64_t) :: n
      complex(kind=c_float_complex), target :: A(*)
      integer(kind=c_int64_t) :: lda
      integer(kind=c_int64_t) :: nb
      integer(kind=c_int) :: p
      integer(kind=c_int) :: q
      integer(kind=c_int) :: mpi_comm_
      type(c_ptr) :: ret_val

      write(*,"(A)") "**Error! A SLATE stub routine was called"
      stop

   end function

   subroutine slate_chol_factor_c32(A, num_opts, opts)
       use iso_c_binding
       implicit none
       type(c_ptr) :: A
       integer(kind=c_int) :: num_opts
       type(c_ptr), target :: opts

       write(*,"(A)") "**Error! A SLATE stub routine was called"
       stop

   end subroutine

   subroutine slate_chol_inverse_using_factor_c32(A, num_opts, opts)
       use iso_c_binding
       implicit none
       type(c_ptr) :: A
       integer(kind=c_int) :: num_opts
       type(c_ptr), target :: opts

       write(*,"(A)") "**Error! A SLATE stub routine was called"
       stop

   end subroutine


   subroutine slate_HermitianMatrix_destroy_c32(A)
       use iso_c_binding
       implicit none
       type(c_ptr) :: A

       write(*,"(A)") "**Error! A SLATE stub routine was called"
       stop

   end subroutine

end module slate
