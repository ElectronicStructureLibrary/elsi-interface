#if defined(HAVE_AVX) || defined(HAVE_AVX2)
 interface
   subroutine single_hh_trafo_complex_avx_avx2_1hv_double(q, hh, pnb, pnq, pldq) &
                             bind(C, name="single_hh_trafo_complex_avx_avx2_1hv_double")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq
     complex(kind=c_double)     :: q(*)
     complex(kind=c_double)     :: hh(pnb,2)
   end subroutine
 end interface
#endif
#if defined(HAVE_AVX) || defined(HAVE_AVX2)
 interface
   subroutine single_hh_trafo_complex_avx_avx2_1hv_single(q, hh, pnb, pnq, pldq) &
                             bind(C, name="single_hh_trafo_complex_avx_avx2_1hv_single")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq
     complex(kind=c_float)   :: q(*)
     complex(kind=c_float)   :: hh(pnb,2)
   end subroutine
 end interface
#endif
#if defined(HAVE_AVX) || defined(HAVE_AVX2)
 interface
   subroutine double_hh_trafo_complex_avx_avx2_2hv_double(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="double_hh_trafo_complex_avx_avx2_2hv_double")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     complex(kind=c_double)     :: q(*)
     complex(kind=c_double)     :: hh(pnb,2)
   end subroutine
 end interface
#endif
#if defined(HAVE_AVX) || defined(HAVE_AVX2)
 interface
   subroutine double_hh_trafo_complex_avx_avx2_2hv_single(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="double_hh_trafo_complex_avx_avx2_2hv_single")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     complex(kind=c_float)   :: q(*)
     complex(kind=c_float)   :: hh(pnb,2)
   end subroutine
 end interface
#endif
#ifdef HAVE_SSE_INTRINSICS
 interface
   subroutine single_hh_trafo_complex_sse_1hv_double(q, hh, pnb, pnq, pldq) &
                             bind(C, name="single_hh_trafo_complex_sse_1hv_double")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq
     complex(kind=c_double)     :: q(*)
     complex(kind=c_double)     :: hh(pnb,2)
   end subroutine
 end interface
#endif
#ifdef HAVE_SSE_INTRINSICS
 interface
   subroutine single_hh_trafo_complex_sse_1hv_single(q, hh, pnb, pnq, pldq) &
                             bind(C, name="single_hh_trafo_complex_sse_1hv_single")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq
     complex(kind=c_float)   :: q(*)
     complex(kind=c_float)   :: hh(pnb,2)
   end subroutine
 end interface
#endif
#ifdef HAVE_SSE_INTRINSICS
 interface
   subroutine double_hh_trafo_complex_sse_2hv_double(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="double_hh_trafo_complex_sse_2hv_double")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     complex(kind=c_double)     :: q(*)
     complex(kind=c_double)     :: hh(pnb,2)
   end subroutine
 end interface
#endif
#ifdef HAVE_SSE_INTRINSICS
 interface
   subroutine double_hh_trafo_complex_sse_2hv_single(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="double_hh_trafo_complex_sse_2hv_single")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     complex(kind=c_float)   :: q(*)
     complex(kind=c_float)   :: hh(pnb,2)
   end subroutine
 end interface
#endif
#if defined(HAVE_AVX) || defined(HAVE_AVX2)
 interface
   subroutine double_hh_trafo_real_avx_avx2_2hv_double(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="double_hh_trafo_real_avx_avx2_2hv_double")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     type(c_ptr), value      :: q
     real(kind=c_double)     :: hh(pnb,6)
   end subroutine
 end interface
#endif
#if defined(HAVE_AVX) || defined(HAVE_AVX2)
 interface
   subroutine double_hh_trafo_real_avx_avx2_2hv_single(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="double_hh_trafo_real_avx_avx2_2hv_single")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     type(c_ptr), value      :: q
     real(kind=c_float)      :: hh(pnb,6)
   end subroutine
 end interface
#endif
#if defined(HAVE_AVX) || defined(HAVE_AVX2)
 interface
   subroutine quad_hh_trafo_real_avx_avx2_4hv_double(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="quad_hh_trafo_real_avx_avx2_4hv_double")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     type(c_ptr), value      :: q
     real(kind=c_double)     :: hh(pnb,6)
   end subroutine
 end interface
#endif
#if defined(HAVE_AVX) || defined(HAVE_AVX2)
 interface
   subroutine quad_hh_trafo_real_avx_avx2_4hv_single(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="quad_hh_trafo_real_avx_avx2_4hv_single")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     type(c_ptr), value      :: q
     real(kind=c_float)      :: hh(pnb,6)
   end subroutine
 end interface
#endif
#if defined(HAVE_AVX) || defined(HAVE_AVX2)
 interface
   subroutine hexa_hh_trafo_real_avx_avx2_6hv_double(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="hexa_hh_trafo_real_avx_avx2_6hv_double")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     type(c_ptr), value      :: q
     real(kind=c_double)     :: hh(pnb,6)
   end subroutine
 end interface
#endif
#if defined(HAVE_AVX) || defined(HAVE_AVX2)
 interface
   subroutine hexa_hh_trafo_real_avx_avx2_6hv_single(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="hexa_hh_trafo_real_avx_avx2_6hv_single")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     type(c_ptr), value      :: q
     real(kind=c_float)      :: hh(pnb,6)
   end subroutine
 end interface
#endif
#ifdef HAVE_SSE_INTRINSICS
 interface
   subroutine double_hh_trafo_real_sse_2hv_double(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="double_hh_trafo_real_sse_2hv_double")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int) :: pnb, pnq, pldq, pldh
     type(c_ptr), value  :: q
     real(kind=c_double) :: hh(pnb,6)
   end subroutine
 end interface
#endif
#ifdef HAVE_SSE_INTRINSICS
 interface
   subroutine double_hh_trafo_real_sse_2hv_single(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="double_hh_trafo_real_sse_2hv_single")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     type(c_ptr), value      :: q
     real(kind=c_float)      :: hh(pnb,6)
   end subroutine
 end interface
#endif
#ifdef HAVE_SSE_INTRINSICS
 interface
   subroutine quad_hh_trafo_real_sse_4hv_double(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="quad_hh_trafo_real_sse_4hv_double")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     type(c_ptr), value      :: q
     real(kind=c_double)     :: hh(pnb,6)
   end subroutine
 end interface
#endif
#ifdef HAVE_SSE_INTRINSICS
 interface
   subroutine quad_hh_trafo_real_sse_4hv_single(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="quad_hh_trafo_real_sse_4hv_single")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     type(c_ptr), value      :: q
     real(kind=c_float)      :: hh(pnb,6)
   end subroutine
 end interface
#endif
#ifdef HAVE_SSE_INTRINSICS
 interface
   subroutine hexa_hh_trafo_real_sse_6hv_double(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="hexa_hh_trafo_real_sse_6hv_double")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     type(c_ptr), value      :: q
     real(kind=c_double)     :: hh(pnb,6)
   end subroutine
 end interface
#endif
#ifdef HAVE_SSE_INTRINSICS
 interface
   subroutine hexa_hh_trafo_real_sse_6hv_single(q, hh, pnb, pnq, pldq, pldh) &
                             bind(C, name="hexa_hh_trafo_real_sse_6hv_single")
     use, intrinsic :: iso_c_binding
     integer(kind=c_int)     :: pnb, pnq, pldq, pldh
     type(c_ptr), value      :: q
     real(kind=c_float)      :: hh(pnb,6)
   end subroutine
 end interface
#endif
#ifdef WITH_REAL_SSE_ASSEMBLY_KERNEL
  interface
    subroutine double_hh_trafo_double(q, hh, nb, nq, ldq, ldh) bind(C,name="double_hh_trafo_double")
      use, intrinsic :: iso_c_binding
      integer(kind=c_int)  :: nb, nq, ldq, ldh
      type(c_ptr), value   :: q
      real(kind=c_double)  :: hh(nb,6)
    end subroutine
  end interface
#endif
#ifdef WITH_COMPLEX_SSE_ASSEMBLY_KERNEL
  interface
    subroutine single_hh_trafo_complex_double(q, hh, nb, nq, ldq) bind(C,name="single_hh_trafo_complex_double")
      use, intrinsic :: iso_c_binding
      integer(kind=c_int)     :: nb, nq, ldq
      complex(kind=c_double)  :: q(*)
      complex(kind=c_double)  :: hh(nb,2)
    end subroutine
  end interface
#endif
#ifdef WITH_REAL_SSE_ASSEMBLY_KERNEL
#ifdef WANT_SINGLE_PRECISION_REAL
  interface
    subroutine double_hh_trafo_single(q, hh, nb, nq, ldq, ldh) bind(C,name="double_hh_trafo_single")
      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: nb, nq, ldq, ldh
      type(c_ptr), value  :: q
      real(kind=c_float)  :: hh(nb,6)
    end subroutine
  end interface
#endif
#endif
#ifdef WITH_COMPLEX_SSE_ASSEMBLY_KERNEL
#ifdef WANT_SINGLE_PRECISION_COMPLEX
  interface
    subroutine single_hh_trafo_complex_single(q, hh, nb, nq, ldq) bind(C,name="single_hh_trafo_complex_single")
      use, intrinsic :: iso_c_binding
      integer(kind=c_int)   :: nb, nq, ldq
      complex(kind=c_float) :: q(*)
      complex(kind=c_float) :: hh(nb,2)
    end subroutine
  end interface
#endif
#endif
