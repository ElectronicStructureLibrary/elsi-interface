!
!    Copyright 2017, L. Hüdepohl and A. Marek, MPCDF
!
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!

!#ifdef HAVE_64BIT_INTEGER_MATH_SUPPORT
!#define FORTRAN_INT_TYPE c_int64_t
!#else
!#define FORTRAN_INT_TYPE c_int
!#endif

module elpa_generated_fortran_interfaces
   use, intrinsic :: iso_c_binding
   implicit none

   interface
      subroutine single_hh_trafo_complex_avx512_1hv_double(q, hh, pnb, pnq, pldq) &
         bind(C, name="single_hh_trafo_complex_avx512_1hv_double")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)     :: pnb, pnq, pldq
         ! complex(kind=c_double_complex)     :: q(*)
         type(c_ptr), value                 :: q
         complex(kind=c_double_complex)     :: hh(pnb,2)
      end subroutine
   end interface
   interface
      subroutine single_hh_trafo_complex_avx512_1hv_single(q, hh, pnb, pnq, pldq) &
         bind(C, name="single_hh_trafo_complex_avx512_1hv_single")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)     :: pnb, pnq, pldq
         ! complex(kind=c_float_complex)     :: q(*)
         type(c_ptr), value                  :: q
         complex(kind=c_float_complex)     :: hh(pnb,2)
      end subroutine
   end interface
   interface
      subroutine double_hh_trafo_complex_avx512_2hv_double(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="double_hh_trafo_complex_avx512_2hv_double")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)     :: pnb, pnq, pldq, pldh
         ! complex(kind=c_double_complex)     :: q(*)
         type(c_ptr), value                   :: q
         complex(kind=c_double_complex)     :: hh(pnb,2)
      end subroutine
   end interface
   interface
      subroutine double_hh_trafo_complex_avx512_2hv_single(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="double_hh_trafo_complex_avx512_2hv_single")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)     :: pnb, pnq, pldq, pldh
         ! complex(kind=c_float_complex)     :: q(*)
         type(c_ptr), value                  :: q
         complex(kind=c_float_complex)     :: hh(pnb,2)
      end subroutine
   end interface
   interface
      subroutine double_hh_trafo_real_VSX_2hv_double(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="double_hh_trafo_real_VSX_2hv_double")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int) :: pnb, pnq, pldq, pldh
         type(c_ptr), value  :: q
         real(kind=c_double) :: hh(pnb,6)
      end subroutine
   end interface
   interface
      subroutine double_hh_trafo_real_VSX_2hv_single(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="double_hh_trafo_real_VSX_2hv_single")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)        :: pnb, pnq, pldq, pldh
         type(c_ptr), value        :: q
         real(kind=c_float)        :: hh(pnb,6)
      end subroutine
   end interface
   interface
      subroutine double_hh_trafo_real_avx512_2hv_double(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="double_hh_trafo_real_avx512_2hv_double")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)     :: pnb, pnq, pldq, pldh
         type(c_ptr), value      :: q
         real(kind=c_double)     :: hh(pnb,6)
      end subroutine
   end interface
   interface
      subroutine double_hh_trafo_real_avx512_2hv_single(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="double_hh_trafo_real_avx512_2hv_single")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)     :: pnb, pnq, pldq, pldh
         type(c_ptr), value      :: q
         real(kind=c_float)      :: hh(pnb,6)
      end subroutine
   end interface
   interface
      subroutine quad_hh_trafo_real_VSX_4hv_double(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="quad_hh_trafo_real_VSX_4hv_double")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)        :: pnb, pnq, pldq, pldh
         type(c_ptr), value        :: q
         real(kind=c_double)        :: hh(pnb,6)
      end subroutine
   end interface
   interface
      subroutine quad_hh_trafo_real_VSX_4hv_single(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="quad_hh_trafo_real_VSX_4hv_single")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)        :: pnb, pnq, pldq, pldh
         type(c_ptr), value        :: q
         real(kind=c_float)        :: hh(pnb,6)
      end subroutine
   end interface
   interface
      subroutine quad_hh_trafo_real_avx512_4hv_double(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="quad_hh_trafo_real_avx512_4hv_double")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)     :: pnb, pnq, pldq, pldh
         type(c_ptr), value      :: q
         real(kind=c_double)     :: hh(pnb,6)
      end subroutine
   end interface
   interface
      subroutine quad_hh_trafo_real_avx512_4hv_single(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="quad_hh_trafo_real_avx512_4hv_single")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)     :: pnb, pnq, pldq, pldh
         type(c_ptr), value      :: q
         real(kind=c_float)      :: hh(pnb,6)
      end subroutine
   end interface
   interface
      subroutine hexa_hh_trafo_real_VSX_6hv_double(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="hexa_hh_trafo_real_VSX_6hv_double")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)        :: pnb, pnq, pldq, pldh
         type(c_ptr), value        :: q
         real(kind=c_double)        :: hh(pnb,6)
      end subroutine
   end interface
   interface
      subroutine hexa_hh_trafo_real_VSX_6hv_single(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="hexa_hh_trafo_real_VSX_6hv_single")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)        :: pnb, pnq, pldq, pldh
         type(c_ptr), value        :: q
         real(kind=c_float)        :: hh(pnb,6)
      end subroutine
   end interface
   interface
      subroutine hexa_hh_trafo_real_avx512_6hv_double(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="hexa_hh_trafo_real_avx512_6hv_double")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)     :: pnb, pnq, pldq, pldh
         type(c_ptr), value      :: q
         real(kind=c_double)     :: hh(pnb,6)
      end subroutine
   end interface
   interface
      subroutine hexa_hh_trafo_real_avx512_6hv_single(q, hh, pnb, pnq, pldq, pldh) &
         bind(C, name="hexa_hh_trafo_real_avx512_6hv_single")
         use, intrinsic :: iso_c_binding
         integer(kind=c_int)     :: pnb, pnq, pldq, pldh
         type(c_ptr), value      :: q
         real(kind=c_float)      :: hh(pnb,6)
      end subroutine
   end interface
   interface
      function elpa_index_instance_c() result(index) bind(C, name="elpa_index_instance")
         import c_ptr
         type(c_ptr) :: index
      end function
   end interface
   interface
      subroutine elpa_index_free_c(index) bind(C, name="elpa_index_free")
         import c_ptr
         type(c_ptr), value :: index
      end subroutine
   end interface
   interface
      function elpa_index_get_int_value_c(index, name, success) result(value) &
         bind(C, name="elpa_index_get_int_value")
         import c_ptr, c_int, c_char
         type(c_ptr), value                         :: index
         character(kind=c_char), intent(in)         :: name(*)
         integer(kind=c_int), intent(out)           :: success
         integer(kind=c_int)                        :: value
      end function
   end interface
   interface
      function elpa_index_set_int_value_c(index, name, value) result(success) &
         bind(C, name="elpa_index_set_int_value")
         import c_ptr, c_int, c_char
         type(c_ptr), value                    :: index
         character(kind=c_char), intent(in)    :: name(*)
         integer(kind=c_int),intent(in), value :: value
         integer(kind=c_int)                   :: success
      end function
   end interface
   interface
      function elpa_index_int_value_is_set_c(index, name) result(success) bind(C, name="elpa_index_int_value_is_set")
         import c_ptr, c_int, c_char
         type(c_ptr), value                    :: index
         character(kind=c_char), intent(in)    :: name(*)
         integer(kind=c_int)                   :: success
      end function
   end interface
   interface
      function elpa_index_get_int_loc_c(index, name) result(loc) bind(C, name="elpa_index_get_int_loc")
         import c_ptr, c_char
         type(c_ptr), value                 :: index
         character(kind=c_char), intent(in) :: name(*)
         type(c_ptr)                        :: loc
      end function
   end interface
   interface
      function elpa_index_get_double_value_c(index, name, success) result(value) bind(C, name="elpa_index_get_double_value")
         import c_ptr, c_int, c_double, c_char
         type(c_ptr), value                              :: index
         character(kind=c_char), intent(in)              :: name(*)
         integer(kind=c_int), intent(out)                :: success
         real(kind=c_double)                             :: value
      end function
   end interface
   interface
      function elpa_index_set_double_value_c(index, name, value) result(success) &
         bind(C, name="elpa_index_set_double_value")
         import c_ptr, c_int, c_double, c_char
         type(c_ptr), value                    :: index
         character(kind=c_char), intent(in)    :: name(*)
         real(kind=c_double),intent(in), value :: value
         integer(kind=c_int)                   :: success
      end function
   end interface
   interface
      function elpa_index_double_value_is_set_c(index, name) result(success) &
         bind(C, name="elpa_index_double_value_is_set")
         import c_ptr, c_int, c_char
         type(c_ptr), value                    :: index
         character(kind=c_char), intent(in)    :: name(*)
         integer(kind=c_int)                   :: success
      end function
   end interface
   interface
      function elpa_index_get_double_loc_c(index, name) result(loc) bind(C, name="elpa_index_get_double_loc")
         import c_ptr, c_char
         type(c_ptr), value                 :: index
         character(kind=c_char), intent(in) :: name(*)
         type(c_ptr)                        :: loc
      end function
   end interface
   interface
      function elpa_index_value_is_set_c(index, name) result(success) bind(C, name="elpa_index_value_is_set")
         import c_ptr, c_int, c_char
         type(c_ptr), value                    :: index
         character(kind=c_char), intent(in)    :: name(*)
         integer(kind=c_int)                   :: success
      end function
   end interface
   interface
      pure function elpa_index_int_value_to_strlen_c(index, name) &
         result(length) bind(C, name="elpa_index_int_value_to_strlen")
         import c_int, c_ptr, c_char
         type(c_ptr), intent(in), value :: index
         character(kind=c_char), intent(in) :: name(*)
         integer(kind=c_int) :: length
      end function
   end interface

   interface
      function elpa_int_string_to_value_c(name, string, value) result(error) bind(C, name="elpa_int_string_to_value")
         import c_int, c_ptr, c_char
         character(kind=c_char), intent(in) :: name(*)
         character(kind=c_char), intent(in) :: string(*)
         integer(kind=c_int), intent(out) :: value
         integer(kind=c_int) :: error
      end function
   end interface

   interface
      function elpa_option_cardinality_c(name) result(n) bind(C, name="elpa_option_cardinality")
         import c_int, c_char
         character(kind=c_char), intent(in) :: name(*)
         integer(kind=c_int) :: n
      end function
   end interface

   interface
      function elpa_option_enumerate_c(name, i) result(value) bind(C, name="elpa_option_enumerate")
         import c_int, c_char
         character(kind=c_char), intent(in) :: name(*)
         integer(kind=c_int), intent(in), value :: i
         integer(kind=c_int) :: value
      end function
   end interface

   interface
      function elpa_index_int_is_valid_c(index, name, new_value) result(success) &
         bind(C, name="elpa_index_int_is_valid")
         import c_int, c_ptr, c_char
         type(c_ptr), intent(in), value :: index
         character(kind=c_char), intent(in) :: name(*)
         integer(kind=c_int), intent(in), value :: new_value
         integer(kind=c_int) :: success
      end function
   end interface

   interface
      function elpa_index_autotune_cardinality_c(index, autotune_level, autotune_domain) result(n) &
         bind(C, name="elpa_index_autotune_cardinality")
         import c_int, c_ptr, c_char
         type(c_ptr), intent(in), value :: index
         integer(kind=c_int), intent(in), value :: autotune_level, autotune_domain
         integer(kind=c_int) :: n
      end function
   end interface

   interface
      function elpa_index_set_autotune_parameters_c(index, autotune_level, autotune_domain, n) result(success) &
         bind(C, name="elpa_index_set_autotune_parameters")
         import c_int, c_ptr, c_char
         type(c_ptr), intent(in), value :: index
         integer(kind=c_int), intent(in), value :: autotune_level, autotune_domain, n
         integer(kind=c_int) :: success
      end function
   end interface

   interface
      function elpa_index_print_autotune_parameters_c(index, autotune_level, autotune_domain) result(success) &
         bind(C, name="elpa_index_print_autotune_parameters")
         import c_int, c_ptr, c_char
         type(c_ptr), intent(in), value :: index
         integer(kind=c_int), intent(in), value :: autotune_level, autotune_domain
         integer(kind=c_int) :: success
      end function
   end interface

   interface
      function elpa_index_print_settings_c(index, file_name) result(success) &
         bind(C, name="elpa_index_print_settings")
         import c_int, c_ptr, c_char
         type(c_ptr), intent(in), value :: index
         character(kind=c_char), intent(in)     :: file_name(*)
         integer(kind=c_int) :: success
      end function
   end interface

   interface
      function elpa_index_load_settings_c(index, file_name) result(success) &
         bind(C, name="elpa_index_load_settings")
         import c_int, c_ptr, c_char
         type(c_ptr), intent(in), value :: index
         character(kind=c_char), intent(in)     :: file_name(*)
         integer(kind=c_int) :: success
      end function
   end interface

   interface
      function elpa_index_print_autotune_state_c(index, autotune_level, autotune_domain, min_loc, &
         min_val, current, cardinality, file_name) result(success) &
         bind(C, name="elpa_index_print_autotune_state")
         import c_int, c_ptr, c_char, c_double
         type(c_ptr), intent(in), value :: index
         integer(kind=c_int), intent(in), value :: autotune_level, autotune_domain, min_loc, current, cardinality
         real(kind=c_double), intent(in), value :: min_val
         character(kind=c_char), intent(in)     :: file_name(*)
         integer(kind=c_int) :: success
      end function
   end interface

   interface
      function elpa_index_load_autotune_state_c(index, autotune_level, autotune_domain, min_loc, &
         min_val, current, cardinality, file_name) result(success) &
         bind(C, name="elpa_index_load_autotune_state")
         import c_int, c_ptr, c_char, c_double
         type(c_ptr), intent(in), value :: index
         integer(kind=c_int), intent(in) :: autotune_level, autotune_domain, min_loc, current, cardinality
         real(kind=c_double), intent(in) :: min_val
         character(kind=c_char), intent(in)     :: file_name(*)
         integer(kind=c_int) :: success
      end function
   end interface

   interface
      subroutine cannons_reduction_d(A, U, local_rowsCast, local_colsCast, a_desc, Res, toStore, row_comm, col_comm) &
         bind(C, name="cannons_reduction_c_d")
         use precision
         real(c_double)                :: A(local_rowsCast, local_colsCast), U(local_rowsCast, local_colsCast)
         real(c_double)                :: Res(local_rowsCast, local_colsCast)
         integer(kind=BLAS_KIND)       :: a_desc(9)
         integer(kind=c_int),value     :: local_rowsCast, local_colsCast
         integer(kind=MPI_KIND),value  :: row_comm, col_comm, ToStore
      end subroutine
   end interface
   interface
      subroutine cannons_triang_rectangular_d(U, B, local_rowsCast, local_colsCast, u_desc, b_desc, Res, row_comm, col_comm) &
         bind(C, name="cannons_triang_rectangular_c_d")
         use precision
         real(c_double)                :: U(local_rowsCast, local_colsCast), B(local_rowsCast, local_colsCast)
         real(c_double)                :: Res(local_rowsCast, local_colsCast)
         integer(kind=BLAS_KIND)       :: u_desc(9), b_desc(9)
         integer(kind=c_int),value     :: local_rowsCast, local_colsCast
         integer(kind=MPI_KIND),value  :: row_comm, col_comm
      end subroutine
   end interface
   interface
      subroutine cannons_reduction_f(A, U, local_rowsCast, local_colsCast, a_desc, Res, toStore, row_comm, col_comm) &
         bind(C, name="cannons_reduction_c_f")
         use precision
         real(c_float)                 :: A(local_rowsCast, local_colsCast), U(local_rowsCast, local_colsCast)
         real(c_float)                 :: Res(local_rowsCast, local_colsCast)
         !type(c_ptr), value           :: A, U, Res
         integer(kind=BLAS_KIND)       :: a_desc(9)
         integer(kind=c_int),value     :: local_rowsCast, local_colsCast
         integer(kind=MPI_KIND),value  :: row_comm, col_comm, ToStore
      end subroutine
   end interface
   interface
      subroutine cannons_triang_rectangular_f(U, B, local_rowsCast, local_colsCast, u_desc, b_desc, Res, row_comm, col_comm) &
         bind(C, name="cannons_triang_rectangular_c_f")
         use precision
         real(c_float)                        :: U(local_rowsCast, local_colsCast), B(local_rowsCast, local_colsCast)
         real(c_float)                        :: Res(local_rowsCast, local_colsCast)
         integer(kind=BLAS_KIND)              :: u_desc(9), b_desc(9)
         integer(kind=c_int),value            :: local_rowsCast, local_colsCast
         integer(kind=MPI_KIND),value         :: row_comm, col_comm
      end subroutine
   end interface
   interface
      subroutine cannons_reduction_dc(A, U, local_rowsCast, local_colsCast, a_desc, Res, toStore, row_comm, col_comm) &
         bind(C, name="cannons_reduction_c_dc")
         use precision
         complex(c_double)             :: A(local_rowsCast, local_colsCast), U(local_rowsCast, local_colsCast)
         complex(c_double)             :: Res(local_rowsCast, local_colsCast)
         integer(kind=BLAS_KIND)       :: a_desc(9)
         integer(kind=c_int),value     :: local_rowsCast, local_colsCast
         integer(kind=MPI_KIND),value  :: row_comm, col_comm, ToStore
      end subroutine
   end interface
   interface
      subroutine cannons_triang_rectangular_dc(U, B, local_rowsCast, local_colsCast, u_desc, b_desc, Res, row_comm, col_comm) &
         bind(C, name="cannons_triang_rectangular_c_dc")
         use precision
         complex(c_double)             :: U(local_rowsCast, local_colsCast), B(local_rowsCast, local_colsCast)
         complex(c_double)             :: Res(local_rowsCast, local_colsCast)
         integer(kind=BLAS_KIND)       :: u_desc(9), b_desc(9)
         integer(kind=c_int),value     :: local_rowsCast, local_colsCast
         integer(kind=MPI_KIND),value  :: row_comm, col_comm
      end subroutine
   end interface
   interface
      subroutine cannons_reduction_fc(A, U, local_rowsCast, local_colsCast, a_desc, Res, toStore, row_comm, col_comm) &
         bind(C, name="cannons_reduction_c_fc")
         use precision
         complex(c_float)              :: A(local_rowsCast, local_colsCast), U(local_rowsCast, local_colsCast)
         complex(c_float)              :: Res(local_rowsCast, local_colsCast)
         !type(c_ptr), value           :: A, U, Res
         integer(kind=BLAS_KIND)       :: a_desc(9)
         integer(kind=c_int),value     :: local_rowsCast, local_colsCast
         integer(kind=MPI_KIND),value  :: row_comm, col_comm, ToStore
      end subroutine
   end interface
   interface
      subroutine cannons_triang_rectangular_fc(U, B, local_rowsCast, local_colsCast, u_desc, b_desc, Res, row_comm, col_comm) &
         bind(C, name="cannons_triang_rectangular_c_fc")
         use precision
         complex(c_float)              :: U(local_rowsCast, local_colsCast), B(local_rowsCast, local_colsCast)
         complex(c_float)              :: Res(local_rowsCast, local_colsCast)
         integer(kind=BLAS_KIND)       :: u_desc(9), b_desc(9)
         integer(kind=c_int),value     :: local_rowsCast, local_colsCast
         integer(kind=MPI_KIND),value  :: row_comm, col_comm
      end subroutine
   end interface
   interface
      subroutine test_c_bindings(A, local_rows, local_cols, np_rows, np_cols, my_prow, my_pcol, a_desc, &
         Res, row_comm, col_comm) &
         bind(C, name="d_test_c_bindings_c")
         use, intrinsic :: iso_c_binding
         real(c_double)             :: A(local_rows, local_cols), Res(local_rows, local_cols)
         !type(c_ptr), value        :: A, Res
         integer(kind=c_int)        :: a_desc(9)
         integer(kind=c_int),value  :: local_rows, local_cols
         integer(kind=c_int),value  :: np_rows, np_cols, my_prow, my_pcol, row_comm, col_comm
      end subroutine
   end interface

end module
