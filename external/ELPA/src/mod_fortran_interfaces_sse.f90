module kernel_interfaces

  implicit none

  interface
    subroutine double_hh_trafo_double(q, hh, nb, nq, ldq, ldh) bind(C,name="double_hh_trafo_double")
      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: nb, nq, ldq, ldh
      type(c_ptr), value  :: q
      real(kind=c_double) :: hh(nb,6)
    end subroutine
  end interface

  interface
    subroutine single_hh_trafo_complex_double(q, hh, nb, nq, ldq) bind(C,name="single_hh_trafo_complex_double")
      use, intrinsic :: iso_c_binding
      integer(kind=c_int)    :: nb, nq, ldq
      complex(kind=c_double) :: q(*)
      complex(kind=c_double) :: hh(nb,2)
    end subroutine
  end interface

end module
