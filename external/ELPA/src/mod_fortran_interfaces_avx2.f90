module kernel_interfaces

  implicit none

  interface
    subroutine single_hh_trafo_complex_avx_avx2_1hv_double(q, hh, pnb, pnq, pldq) &
      bind(C, name="single_hh_trafo_complex_avx_avx2_1hv_double")

      use, intrinsic :: iso_c_binding

      integer(kind=c_int) :: pnb, pnq, pldq
      complex(kind=c_double) :: q(*)
      complex(kind=c_double) :: hh(pnb,2)

    end subroutine
  end interface

  interface
    subroutine double_hh_trafo_real_avx_avx2_2hv_double(q, hh, pnb, pnq, pldq, pldh) &
      bind(C, name="double_hh_trafo_real_avx_avx2_2hv_double")

      use, intrinsic :: iso_c_binding

      integer(kind=c_int) :: pnb, pnq, pldq, pldh
      type(c_ptr), value :: q
      real(kind=c_double) :: hh(pnb,6)

    end subroutine
  end interface

end module
