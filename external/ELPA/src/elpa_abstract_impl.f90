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
!> \brief Fortran module to provide an abstract definition of the implementation. Do not use directly. Use the module "elpa"
module elpa_abstract_impl
   use elpa_api
   use elpa_generated_fortran_interfaces

   use ftimings

   implicit none

   ! The reason to have this additional layer is to allow for members (here the
   ! 'timer' object) that can be used internally but are not exposed to the
   ! public API. This cannot be done via 'private' members, as the scope of
   ! 'private' is per-file.
   !
   ! Thus, other sub-types or suplementary routines cannot use these members
   ! (unless they would all be implemented in one giant file)
   !
   type, abstract, extends(elpa_t) :: elpa_abstract_impl_t
      type(timer_t) :: timer
      type(timer_t) :: autotune_timer
      type(c_ptr)         :: index = C_NULL_PTR
      logical             :: eigenvalues_only
   contains
      procedure, public :: elpa_set_integer                      !< private methods to implement the setting of an integer/double key/value pair
      procedure, public :: elpa_set_double

      procedure, public :: elpa_get_integer                      !< private methods to implement the querry of an integer/double key/value pair
      procedure, public :: elpa_get_double

   end type

contains

   !> \brief internal subroutine to set an integer key/value pair
   !> Parameters
   !> \param   self       the allocated ELPA object
   !> \param   name       string, the key
   !> \param   value      integer, the value to be set
   !> \result  error      integer, the error code
   subroutine elpa_set_integer(self, name, value, error)
      use, intrinsic :: iso_c_binding
      use elpa_utilities, only : error_unit
      class(elpa_abstract_impl_t)     :: self
      character(*), intent(in)        :: name
      integer(kind=c_int), intent(in) :: value
      integer                         :: error
      integer                         :: actual_error

      actual_error = elpa_index_set_int_value_c(self%index, name // c_null_char, value)

      error = actual_error
   end subroutine

   !> \brief internal subroutine to get an integer key/value pair
   !> Parameters
   !> \param   self       the allocated ELPA object
   !> \param   name       string, the key
   !> \param   value      integer, the value of the key/vaue pair
   !> \param   error      integer, optional, to store an error code
   subroutine elpa_get_integer(self, name, value, error)
      use, intrinsic :: iso_c_binding
      use elpa_utilities, only : error_unit
      class(elpa_abstract_impl_t)    :: self
      character(*), intent(in)       :: name
      integer(kind=c_int)            :: value
      integer, intent(out)           :: error
      integer                        :: actual_error

      value = elpa_index_get_int_value_c(self%index, name // c_null_char, actual_error)

      error = actual_error
   end subroutine

   !> \brief internal subroutine to set a double key/value pair
   !> Parameters
   !> \param   self       the allocated ELPA object
   !> \param   name       string, the key
   !> \param   value      double, the value to be set
   !> \result  error      integer, the error code
   subroutine elpa_set_double(self, name, value, error)
      use, intrinsic :: iso_c_binding
      use elpa_utilities, only : error_unit
      class(elpa_abstract_impl_t)     :: self
      character(*), intent(in)        :: name
      real(kind=c_double), intent(in) :: value
      integer                         :: actual_error

      integer                         :: error
      actual_error = elpa_index_set_double_value_c(self%index, name // c_null_char, value)

      error = actual_error
   end subroutine

   !> \brief internal subroutine to get an double key/value pair
   !> Parameters
   !> \param   self       the allocated ELPA object
   !> \param   name       string, the key
   !> \param   value      double, the value of the key/vaue pair
   !> \param   error      integer, optional, to store an error code
   subroutine elpa_get_double(self, name, value, error)
      use, intrinsic :: iso_c_binding
      use elpa_utilities, only : error_unit
      class(elpa_abstract_impl_t)    :: self
      character(*), intent(in)       :: name
      real(kind=c_double)            :: value
      integer, intent(out)           :: error
      integer                        :: actual_error

      value = elpa_index_get_double_value_c(self%index, name // c_null_char, actual_error)
      error = actual_error
   end subroutine

end module
