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
!> \brief Fortran module which provides the definition of the ELPA API. Do not use directly! Use the module "elpa"

module elpa_api
   use elpa_constants
   use, intrinsic :: iso_c_binding
   implicit none

   interface
      pure function elpa_strerr_c(elpa_error) result(string) bind(C, name="elpa_strerr")
         import c_int, c_ptr
         integer(kind=c_int), intent(in), value :: elpa_error
         type(c_ptr) :: string
      end function
   end interface
   interface
      function elpa_int_value_to_string_c(name, value, string) &
         result(error) bind(C, name="elpa_int_value_to_string")
         import c_int, c_ptr, c_char
         character(kind=c_char), intent(in) :: name(*)
         integer(kind=c_int), intent(in), value :: value
         type(c_ptr), intent(out) :: string
         integer(kind=c_int) :: error
      end function
   end interface

   interface
      pure function elpa_int_value_to_strlen_c(name, value) &
         result(length) bind(C, name="elpa_int_value_to_strlen")
         import c_int, c_ptr, c_char
         character(kind=c_char), intent(in) :: name(*)
         integer(kind=c_int), intent(in), value :: value
         integer(kind=c_int) :: length
      end function
   end interface

   integer, private, parameter :: earliest_api_version = 20170403 !< Definition of the earliest API version supported
   !< with the current release
   integer, private, parameter :: current_api_version  = 20200417  !< Definition of the current API version

   integer, private, parameter :: earliest_autotune_version = 20171201 !< Definition of the earliest API version
   !< which supports autotuning
   integer, private            :: api_version_set
   logical, private            :: initDone = .false.

   public :: elpa_t, &
      c_int, &
      c_double, c_double_complex, &
      c_float, c_float_complex

   !> \brief Abstract definition of the elpa_t type
   type, abstract :: elpa_t
      private

      !< these have to be public for proper bounds checking, sadly
      integer(kind=c_int), public, pointer :: na => NULL()
      integer(kind=c_int), public, pointer :: nev => NULL()
      integer(kind=c_int), public, pointer :: local_nrows => NULL()
      integer(kind=c_int), public, pointer :: local_ncols => NULL()
      integer(kind=c_int), public, pointer :: nblk => NULL()

   contains
      ! general
      procedure(elpa_setup_i),   deferred, public :: setup          !< method to setup an ELPA object
      procedure(elpa_destroy_i), deferred, public :: destroy        !< method to destroy an ELPA object

      ! key/value store
      generic, public :: set => &                                   !< export a method to set integer/double/float key/values
         elpa_set_integer, &
         elpa_set_float, &
         elpa_set_double

      generic, public :: get => &                                   !< export a method to get integer/double/float key/values
         elpa_get_integer, &
         elpa_get_float, &
         elpa_get_double

      procedure(elpa_is_set_i),  deferred, public :: is_set         !< method to check whether key/value is set
      procedure(elpa_can_set_i), deferred, public :: can_set        !< method to check whether key/value can be set

      ! call before setup if created from the legacy api
      ! remove this function completely after the legacy api is dropped
      procedure(elpa_creating_from_legacy_api_i), deferred, public :: creating_from_legacy_api

      ! Timer
      procedure(elpa_get_time_i), deferred, public :: get_time        !< method to get the times from the timer object
      procedure(elpa_print_times_i), deferred, public :: print_times  !< method to print the timings tree
      procedure(elpa_timer_start_i), deferred, public :: timer_start  !< method to start a time measurement
      procedure(elpa_timer_stop_i), deferred, public :: timer_stop    !< method to stop a time measurement

      ! Actual math routines
      generic, public :: eigenvectors => &                          !< method eigenvectors for solving the full eigenvalue problem
         elpa_eigenvectors_d, &                                    !< the eigenvalues and (parts of) the eigenvectors are computed
         elpa_eigenvectors_f, &                                    !< for symmetric real valued / hermitian complex valued matrices
         elpa_eigenvectors_dc, &
         elpa_eigenvectors_fc

      generic, public :: eigenvalues => &                           !< method eigenvalues for solving the eigenvalue problem
         elpa_eigenvalues_d, &                                     !< only the eigenvalues are computed
         elpa_eigenvalues_f, &                                     !< for symmetric real valued / hermitian complex valued matrices
         elpa_eigenvalues_dc, &
         elpa_eigenvalues_fc

      generic, public :: skew_eigenvectors => &                     !< method skew_eigenvectors for solving the full skew-symmetric eigenvalue problem
         elpa_skew_eigenvectors_d, &                               !< the eigenvalues and (parts of) the eigenvectors are computed
         elpa_skew_eigenvectors_f                                  !< for symmetric real valued skew-symmetric matrices

      generic, public :: skew_eigenvalues => &                      !< method skew_eigenvalues for solving the skew-symmetric eigenvalue problem
         elpa_skew_eigenvalues_d, &                                !< only the eigenvalues are computed
         elpa_skew_eigenvalues_f                                   !< for symmetric real valued skew-symmetric matrices

      generic, public :: generalized_eigenvectors => &              !< method eigenvectors for solving the full generalized eigenvalue problem
         elpa_generalized_eigenvectors_d, &                        !< the eigenvalues and (parts of) the eigenvectors are computed
         elpa_generalized_eigenvectors_f, &                        !< for symmetric real valued / hermitian complex valued matrices
         elpa_generalized_eigenvectors_dc, &
         elpa_generalized_eigenvectors_fc

      generic, public :: generalized_eigenvalues => &              !< method eigenvectors for solving the full generalized eigenvalue problem
         elpa_generalized_eigenvalues_d, &                        !< only the eigenvalues
         elpa_generalized_eigenvalues_f, &                        !< for symmetric real valued / hermitian complex valued matrices
         elpa_generalized_eigenvalues_dc, &
         elpa_generalized_eigenvalues_fc

      generic, public :: hermitian_multiply => &                    !< method for a "hermitian" multiplication of matrices a and b
         elpa_hermitian_multiply_d, &                              !< for real valued matrices:   a**T * b
         elpa_hermitian_multiply_dc, &                             !< for complex valued matrices a**H * b
         elpa_hermitian_multiply_f, &
         elpa_hermitian_multiply_fc

      generic, public :: cholesky => &                              !< method for the cholesky factorisation of matrix a
         elpa_cholesky_d, &
         elpa_cholesky_f, &
         elpa_cholesky_dc, &
         elpa_cholesky_fc

      generic, public :: invert_triangular => &                     !< method to invert a upper triangular matrix a
         elpa_invert_trm_d, &
         elpa_invert_trm_f, &
         elpa_invert_trm_dc, &
         elpa_invert_trm_fc

      generic, public :: solve_tridiagonal => &                      !< method to solve the eigenvalue problem for a tridiagonal
         elpa_solve_tridiagonal_d, &                                !< matrix
         elpa_solve_tridiagonal_f

      procedure(print_settings_i), deferred, public :: print_settings !< method to print all parameters
      procedure(store_settings_i), deferred, public :: store_settings !< method to save all parameters
      procedure(load_settings_i), deferred, public :: load_settings !< method to save all parameters
      ! Auto-tune
      procedure(elpa_autotune_setup_i), deferred, public :: autotune_setup       !< method to prepare the ELPA autotuning
      procedure(elpa_autotune_step_i), deferred, public :: autotune_step         !< method to do an autotuning step
      procedure(elpa_autotune_set_best_i), deferred, public :: autotune_set_best !< method to set the best options
      procedure(elpa_autotune_print_best_i), deferred, public :: autotune_print_best !< method to print the best options
      procedure(elpa_autotune_print_state_i), deferred, public :: autotune_print_state !< method to print the state
      procedure(elpa_autotune_save_state_i), deferred, public :: autotune_save_state !< method to save the state
      procedure(elpa_autotune_load_state_i), deferred, public :: autotune_load_state !< method to load the state

      !> \brief These method have to be public, in order to be overrideable in the extension types
      procedure(elpa_set_integer_i), deferred, public :: elpa_set_integer
      procedure(elpa_set_float_i),  deferred, public :: elpa_set_float
      procedure(elpa_set_double_i),  deferred, public :: elpa_set_double

      procedure(elpa_get_integer_i), deferred, public :: elpa_get_integer
      procedure(elpa_get_float_i),  deferred, public :: elpa_get_float
      procedure(elpa_get_double_i),  deferred, public :: elpa_get_double

      procedure(elpa_eigenvectors_d_i),    deferred, public :: elpa_eigenvectors_d
      procedure(elpa_eigenvectors_f_i),    deferred, public :: elpa_eigenvectors_f
      procedure(elpa_eigenvectors_dc_i), deferred, public :: elpa_eigenvectors_dc
      procedure(elpa_eigenvectors_fc_i), deferred, public :: elpa_eigenvectors_fc

      procedure(elpa_eigenvalues_d_i),    deferred, public :: elpa_eigenvalues_d
      procedure(elpa_eigenvalues_f_i),    deferred, public :: elpa_eigenvalues_f
      procedure(elpa_eigenvalues_dc_i), deferred, public :: elpa_eigenvalues_dc
      procedure(elpa_eigenvalues_fc_i), deferred, public :: elpa_eigenvalues_fc

      procedure(elpa_skew_eigenvectors_d_i),    deferred, public :: elpa_skew_eigenvectors_d
      procedure(elpa_skew_eigenvectors_f_i),    deferred, public :: elpa_skew_eigenvectors_f

      procedure(elpa_skew_eigenvalues_d_i),    deferred, public :: elpa_skew_eigenvalues_d
      procedure(elpa_skew_eigenvalues_f_i),    deferred, public :: elpa_skew_eigenvalues_f

      procedure(elpa_generalized_eigenvectors_d_i),    deferred, public :: elpa_generalized_eigenvectors_d
      procedure(elpa_generalized_eigenvectors_f_i),    deferred, public :: elpa_generalized_eigenvectors_f
      procedure(elpa_generalized_eigenvectors_dc_i), deferred, public :: elpa_generalized_eigenvectors_dc
      procedure(elpa_generalized_eigenvectors_fc_i), deferred, public :: elpa_generalized_eigenvectors_fc

      procedure(elpa_generalized_eigenvalues_d_i),    deferred, public :: elpa_generalized_eigenvalues_d
      procedure(elpa_generalized_eigenvalues_f_i),    deferred, public :: elpa_generalized_eigenvalues_f
      procedure(elpa_generalized_eigenvalues_dc_i), deferred, public :: elpa_generalized_eigenvalues_dc
      procedure(elpa_generalized_eigenvalues_fc_i), deferred, public :: elpa_generalized_eigenvalues_fc

      procedure(elpa_hermitian_multiply_d_i),  deferred, public :: elpa_hermitian_multiply_d
      procedure(elpa_hermitian_multiply_f_i),  deferred, public :: elpa_hermitian_multiply_f
      procedure(elpa_hermitian_multiply_dc_i), deferred, public :: elpa_hermitian_multiply_dc
      procedure(elpa_hermitian_multiply_fc_i), deferred, public :: elpa_hermitian_multiply_fc

      procedure(elpa_cholesky_d_i),    deferred, public :: elpa_cholesky_d
      procedure(elpa_cholesky_f_i),    deferred, public :: elpa_cholesky_f
      procedure(elpa_cholesky_dc_i), deferred, public :: elpa_cholesky_dc
      procedure(elpa_cholesky_fc_i), deferred, public :: elpa_cholesky_fc

      procedure(elpa_invert_trm_d_i),    deferred, public :: elpa_invert_trm_d
      procedure(elpa_invert_trm_f_i),    deferred, public :: elpa_invert_trm_f
      procedure(elpa_invert_trm_dc_i), deferred, public :: elpa_invert_trm_dc
      procedure(elpa_invert_trm_fc_i), deferred, public :: elpa_invert_trm_fc

      procedure(elpa_solve_tridiagonal_d_i), deferred, public :: elpa_solve_tridiagonal_d
      procedure(elpa_solve_tridiagonal_f_i), deferred, public :: elpa_solve_tridiagonal_f
   end type elpa_t

   !> \brief Abstract definition of the elpa_autotune type
   type, abstract :: elpa_autotune_t
      private
   contains
      procedure(elpa_autotune_destroy_i), deferred, public :: destroy
      procedure(elpa_autotune_print_i), deferred, public :: print
   end type

   !> \brief definition of helper function to get C strlen
   !> Parameters
   !> \details
   !> \param   ptr         type(c_ptr) : pointer to string
   !> \result  size        integer(kind=c_size_t) : length of string
   interface
      pure function elpa_strlen_c(ptr) result(size) bind(c, name="strlen")
         use, intrinsic :: iso_c_binding
         implicit none
         type(c_ptr), intent(in), value :: ptr
         integer(kind=c_size_t) :: size
      end function
   end interface

   !> \brief abstract definition of the ELPA setup method
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \result  error       integer : error code, which can be queried with elpa_strerr()
   abstract interface
      function elpa_setup_i(self) result(error)
         import elpa_t
         implicit none
         class(elpa_t), intent(inout)   :: self
         integer                        :: error
      end function
   end interface

   !> \brief abstract definition of the print_settings method
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   error       integer, optional
   !> Prints all the elpa parameters
   abstract interface
      subroutine print_settings_i(self, error)
         import elpa_t
         implicit none
         class(elpa_t), intent(inout)   :: self
         integer, intent(out)           :: error
      end subroutine
   end interface

   !> \brief abstract definition of the store_settings method
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   file_name   string, the name of the file where to save the parameters
   !> \param   error       integer, optional
   !> Saves all the elpa parameters
   abstract interface
      subroutine store_settings_i(self, file_name, error)
         import elpa_t
         implicit none
         class(elpa_t), intent(inout)  :: self
         character(*), intent(in)      :: file_name
         integer, intent(out)          :: error
      end subroutine
   end interface

   !> \brief abstract definition of the load_settings method
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   file_name   string, the name of the file from which to load the parameters
   !> \param   error       integer, optional
   !> Loads all the elpa parameters
   abstract interface
      subroutine load_settings_i(self, file_name, error)
         import elpa_t
         implicit none
         class(elpa_t), intent(inout)   :: self
         character(*), intent(in)       :: file_name
         integer, intent(out)           :: error
      end subroutine
   end interface

   !> \brief abstract definition of the autotune setup method
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object, which should be tuned
   !> \param   level       integer: the level of "thoroughness" of the tuning steps
   !> \param   domain      integer: domain (real/complex) which should be tuned
   !> \result  tune_state  class(elpa_autotune_t): the autotuning object
   abstract interface
      function elpa_autotune_setup_i(self, level, domain, error) result(tune_state)
         import elpa_t, elpa_autotune_t
         implicit none
         class(elpa_t), intent(inout), target :: self
         integer, intent(in)                  :: level, domain
         class(elpa_autotune_t), pointer      :: tune_state
         integer                              :: error
      end function
   end interface

   !> \brief abstract definition of the autotune step method
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object, which should be tuned
   !> \param   tune_state  class(elpa_autotune_t): the autotuning object
   !> \param   unfinished  logical: state whether tuning is unfinished or not
   !> \param   error       integer, optional
   abstract interface
      function elpa_autotune_step_i(self, tune_state, error) result(unfinished)
         import elpa_t, elpa_autotune_t
         implicit none
         class(elpa_t), intent(inout)                  :: self
         class(elpa_autotune_t), intent(inout), target :: tune_state
         logical                                       :: unfinished
         integer, intent(out)                          :: error
      end function
   end interface

   !> \brief abstract definition of the autotune set_best method
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object, which should be tuned
   !> \param   tune_state  class(elpa_autotune_t): the autotuning object
   !> \param   error       integer, optional
   !> Sets the best combination of ELPA options
   abstract interface
      subroutine elpa_autotune_set_best_i(self, tune_state, error)
         import elpa_t, elpa_autotune_t
         implicit none
         class(elpa_t), intent(inout)               :: self
         class(elpa_autotune_t), intent(in), target :: tune_state
         integer, intent(out)                       :: error

      end subroutine
   end interface

   !> \brief abstract definition of the autotune print best method
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object, which should be tuned
   !> \param   tune_state  class(elpa_autotune_t): the autotuning object
   !> \param   error       integer, optional
   !> Prints the best combination of ELPA options
   abstract interface
      subroutine elpa_autotune_print_best_i(self, tune_state, error)
         import elpa_t, elpa_autotune_t
         implicit none
         class(elpa_t), intent(inout)               :: self
         class(elpa_autotune_t), intent(in), target :: tune_state
         integer, intent(out)                       :: error

      end subroutine
   end interface

   !> \brief abstract definition of the autotune print state method
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object, which should be tuned
   !> \param   tune_state  class(elpa_autotune_t): the autotuning object
   !> \param   error       integer, optional
   !> Prints the autotuning state
   abstract interface
      subroutine elpa_autotune_print_state_i(self, tune_state, error)
         import elpa_t, elpa_autotune_t
         implicit none
         class(elpa_t), intent(inout)               :: self
         class(elpa_autotune_t), intent(in), target :: tune_state
         integer,  intent(out)                      :: error
      end subroutine
   end interface

   !> \brief abstract definition of the autotune save state method
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object, which should be tuned
   !> \param   tune_state  class(elpa_autotune_t): the autotuning object
   !> \param   file_name   string, the name of the file where to save the state
   !> \param   error       integer, optional
   !> Saves the autotuning state
   abstract interface
      subroutine elpa_autotune_save_state_i(self, tune_state, file_name, error)
         import elpa_t, elpa_autotune_t
         implicit none
         class(elpa_t), intent(inout)               :: self
         class(elpa_autotune_t), intent(in), target :: tune_state
         character(*), intent(in)                   :: file_name
         integer, intent(out)                       :: error
      end subroutine
   end interface

   !> \brief abstract definition of the autotune load state method
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object, which is being tuned
   !> \param   tune_state  class(elpa_autotune_t): the autotuning object
   !> \param   file_name   string, the name of the file from which to load the autotuning state
   !> \param   error       integer, optional
   !> Loads all the elpa parameters
   abstract interface
      subroutine elpa_autotune_load_state_i(self, tune_state, file_name, error)
         import elpa_t, elpa_autotune_t
         implicit none
         class(elpa_t), intent(inout)               :: self
         class(elpa_autotune_t), intent(in), target :: tune_state
         character(*), intent(in)                   :: file_name
         integer, intent(out)                       :: error
      end subroutine
   end interface

   !> \brief abstract definition of set method for integer values
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name        string: the name of the key
   !> \param   value       integer : the value to set for the key
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr()
   abstract interface
      subroutine elpa_set_integer_i(self, name, value, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         character(*), intent(in)        :: name
         integer(kind=c_int), intent(in) :: value
         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of get method for integer values
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name        string: the name of the key
   !> \param   value       integer : the value corresponding to the key
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr()
   abstract interface
      subroutine elpa_get_integer_i(self, name, value, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                  :: self
         character(*), intent(in)       :: name
         integer(kind=c_int)            :: value
         integer, intent(out)           :: error
      end subroutine
   end interface

   !> \brief abstract definition of is_set method for integer values
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name        string: the name of the key
   !> \result  state       integer : 1 is set, 0 if not, else a negativ error code
   !>                                                    which can be queried with elpa_strerr
   abstract interface
      function elpa_is_set_i(self, name) result(state)
         import elpa_t
         implicit none
         class(elpa_t)            :: self
         character(*), intent(in) :: name
         integer                  :: state
      end function
   end interface

   !> \brief abstract definition of can_set method for integer values
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name        string: the name of the key
   !> \param   value       integer: the valye to associate with the key
   !> \result  state       integer : 1 is set, 0 if not, else a negativ error code
   !>                                                    which can be queried with elpa_strerr
   abstract interface
      function elpa_can_set_i(self, name, value) result(state)
         import elpa_t, c_int
         implicit none
         class(elpa_t)                   :: self
         character(*), intent(in)        :: name
         integer(kind=c_int), intent(in) :: value
         integer                         :: state
      end function
   end interface

   !> \brief abstract definition of set method for float values
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name        string: the name of the key
   !? \param   value       float: the value to associate with the key
   !> \param   error       integer. optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_set_float_i(self, name, value, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         character(*), intent(in)        :: name
         real(kind=c_float), intent(in) :: value
         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of get method for float values
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name        string: the name of the key
   !> \param   value       float: the value associated with the key
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_get_float_i(self, name, value, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                  :: self
         character(*), intent(in)       :: name
         real(kind=c_float)            :: value
         integer, intent(out)           :: error
      end subroutine
   end interface

   !> \brief abstract definition of set method for double values
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name        string: the name of the key
   !? \param   value       double: the value to associate with the key
   !> \param   error       integer. optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_set_double_i(self, name, value, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         character(*), intent(in)        :: name
         real(kind=c_double), intent(in) :: value
         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of get method for double values
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name        string: the name of the key
   !> \param   value       double: the value associated with the key
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_get_double_i(self, name, value, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                  :: self
         character(*), intent(in)       :: name
         real(kind=c_double)            :: value
         integer, intent(out)           :: error
      end subroutine
   end interface

   !> \brief abstract definition of associate method for integer pointers
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name        string: the name of the key
   !> \result  value       integer, pointer: the value associated with the key
   abstract interface
      function elpa_associate_int_i(self, name) result(value)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                  :: self
         character(*), intent(in)       :: name
         integer(kind=c_int), pointer   :: value
      end function
   end interface

   ! Timer routines

   !> \brief abstract definition of get_time method to querry the timer
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name1..6    string: the name of the timer entry, supports up to 6 levels
   !> \result  s           double: the time for the entry name1..6
   abstract interface
      function elpa_get_time_i(self, name1, name2, name3, name4, name5, name6) result(s)
         import elpa_t, c_double
         implicit none
         class(elpa_t), intent(in) :: self
         ! this is clunky, but what can you do..
         character(len=*), intent(in), optional :: name1, name2, name3, name4, name5, name6
         real(kind=c_double) :: s
      end function
   end interface

   !> \brief abstract definition of print method for timer
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   abstract interface
      subroutine elpa_print_times_i(self, name1, name2, name3, name4)
         import elpa_t
         implicit none
         class(elpa_t), intent(in) :: self
         character(len=*), intent(in), optional :: name1, name2, name3, name4
      end subroutine
   end interface

   !> \brief abstract definition of the start method for timer
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name        character(len=*) the name of the entry int the timer tree
   abstract interface
      subroutine elpa_timer_start_i(self, name)
         import elpa_t
         implicit none
         class(elpa_t), intent(inout) :: self
         character(len=*), intent(in) :: name
      end subroutine
   end interface

   !> \brief abstract definition of the stop method for timer
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t): the ELPA object
   !> \param   name        character(len=*) the name of the entry int the timer tree
   abstract interface
      subroutine elpa_timer_stop_i(self, name)
         import elpa_t
         implicit none
         class(elpa_t), intent(inout) :: self
         character(len=*), intent(in) :: name
      end subroutine
   end interface

   ! Actual math routines

   !> \brief abstract definition of interface to solve double real eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvectors
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvectors
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double complex matrix q: on output stores the eigenvectors
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single complex matrix q: on output stores the eigenvectors
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_eigenvectors_&
      &d&
      &_i(self, a, ev, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self

         real(kind=c_double) :: a(self%local_nrows, *), q(self%local_nrows,*)
         real(kind=c_double) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve double real skew-symmetric eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvectors
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvectors
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_skew_eigenvectors_&
      &d&
      &_i(self, a, ev, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self

         real(kind=c_double) :: a(self%local_nrows, *), q(self%local_nrows,*)
         real(kind=c_double) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_eigenvalues_&
      &d&
      &_i(self, a, ev, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         real(kind=c_double) :: a(self%local_nrows, *)
         real(kind=c_double) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a skew-symmetric eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_skew_eigenvalues_&
      &d&
      &_i(self, a, ev, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         real(kind=c_double) :: a(self%local_nrows, *)
         real(kind=c_double) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a generalized eigenvalue problem
   !>
   !>  The dimensions of the matrix a and b (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   b           double real matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   b           single real matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   b           double complex matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double complex matrix q: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   b           single complex matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single complex matrix q: on output stores the eigenvalues

   !> \param   is_already_decomposed   logical, input: is it repeated call with the same b (decomposed in the fist call)?
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_generalized_eigenvectors_&
      &d&
      &_i(self, a, b, ev, q, is_already_decomposed, error)
         use, intrinsic :: iso_c_binding
         use elpa_constants
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         real(kind=c_double) :: a(self%local_nrows, *), b(self%local_nrows, *), q(self%local_nrows, *)
         real(kind=c_double) :: ev(self%na)

         logical             :: is_already_decomposed
         integer, optional   :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a generalized eigenvalue problem
   !>
   !>  The dimensions of the matrix a and b (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   b           double real matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   b           single real matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   b           double complex matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   b           single complex matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues

   !> \param   is_already_decomposed   logical, input: is it repeated call with the same b (decomposed in the fist call)?
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_generalized_eigenvalues_&
      &d&
      &_i(self, a, b, ev, is_already_decomposed, error)
         use, intrinsic :: iso_c_binding
         use elpa_constants
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         real(kind=c_double) :: a(self%local_nrows, *), b(self%local_nrows, *)
         real(kind=c_double) :: ev(self%na)

         logical             :: is_already_decomposed
         integer, optional   :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to compute C : = A**T * B
   !>         where   A is a square matrix (self%a,self%na) which is optionally upper or lower triangular
   !>                 B is a (self%na,ncb) matrix
   !>                 C is a (self%na,ncb) matrix where optionally only the upper or lower
   !>                   triangle may be computed
   !>
   !> the MPI commicators are already known to the type. Thus the class method "setup" must be called
   !> BEFORE this method is used
   !> \details
   !>
   !> \param   self                class(elpa_t), the ELPA object
   !> \param  uplo_a               'U' if A is upper triangular
   !>                              'L' if A is lower triangular
   !>                              anything else if A is a full matrix
   !>                              Please note: This pertains to the original A (as set in the calling program)
   !>                                           whereas the transpose of A is used for calculations
   !>                              If uplo_a is 'U' or 'L', the other triangle is not used at all,
   !>                              i.e. it may contain arbitrary numbers
   !> \param uplo_c                'U' if only the upper diagonal part of C is needed
   !>                              'L' if only the upper diagonal part of C is needed
   !>                              anything else if the full matrix C is needed
   !>                              Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
   !>                                            written to a certain extent, i.e. one shouldn't rely on the content there!
   !> \param ncb                   Number of columns  of global matrices B and C
   !> \param a                     matrix a
   !> \param self%local_nrows      number of rows of local (sub) matrix a, set with method set("local_nrows,value")
   !> \param self%local_ncols      number of columns of local (sub) matrix a, set with method set("local_ncols,value")
   !> \param b                     matrix b
   !> \param nrows_b               number of rows of local (sub) matrix b
   !> \param ncols_b               number of columns of local (sub) matrix b
   !> \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
   !> \param c                     matrix c
   !> \param nrows_c               number of rows of local (sub) matrix c
   !> \param ncols_c               number of columns of local (sub) matrix c
   !> \param error                 optional argument, error code which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_hermitian_multiply_&
      &d&
      &_i (self,uplo_a, uplo_c, ncb, a, b, nrows_b, ncols_b, &
         c, nrows_c, ncols_c, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         character*1                     :: uplo_a, uplo_c
         integer(kind=c_int), intent(in) :: nrows_b, ncols_b, nrows_c, ncols_c, ncb
         real(kind=c_double) :: a(self%local_nrows,*), b(nrows_b,*), c(nrows_c,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to do a cholesky decomposition of a matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix: the matrix to be decomposed
   !> \param   a           single real matrix: the matrix to be decomposed
   !> \param   a           double complex matrix: the matrix to be decomposed
   !> \param   a           single complex matrix: the matrix to be decomposed
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_cholesky_&
      &d&
      &_i (self, a, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         real(kind=c_double) :: a(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to invert a triangular matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix: the matrix to be inverted
   !> \param   a           single real matrix: the matrix to be inverted
   !> \param   a           double complex matrix: the matrix to be inverted
   !> \param   a           single complex matrix: the matrix to be inverted

   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_invert_trm_&
      &d&
      &_i (self, a, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         real(kind=c_double) :: a(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve the eigenvalue problem for a valued tridiangular matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   d           double real 1d array: the diagonal elements of a matrix defined in setup, on output the eigenvalues
   !>                      in ascending order
   !> \param   e           double real 1d array: the subdiagonal elements of a matrix defined in setup
   !> \param   q           double real matrix: on output contains the eigenvectors
   !> \param   d           single real 1d array: the diagonal elements of a matrix defined in setup, on output the eigenvalues
   !>                      in ascending order
   !> \param   e           single real 1d array: the subdiagonal elements of a matrix defined in setup
   !> \param   q           single real matrix: on output contains the eigenvectors
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_solve_tridiagonal_&
      &d&
      &_i (self, d, e, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         real(kind=c_double)        :: d(self%na), e(self%na)
         real(kind=c_double)        :: q(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve double real eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvectors
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvectors
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double complex matrix q: on output stores the eigenvectors
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single complex matrix q: on output stores the eigenvectors
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_eigenvectors_&
      &f&
      &_i(self, a, ev, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self

         real(kind=c_float) :: a(self%local_nrows, *), q(self%local_nrows,*)
         real(kind=c_float) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve double real skew-symmetric eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvectors
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvectors
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_skew_eigenvectors_&
      &f&
      &_i(self, a, ev, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self

         real(kind=c_float) :: a(self%local_nrows, *), q(self%local_nrows,*)
         real(kind=c_float) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_eigenvalues_&
      &f&
      &_i(self, a, ev, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         real(kind=c_float) :: a(self%local_nrows, *)
         real(kind=c_float) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a skew-symmetric eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_skew_eigenvalues_&
      &f&
      &_i(self, a, ev, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         real(kind=c_float) :: a(self%local_nrows, *)
         real(kind=c_float) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a generalized eigenvalue problem
   !>
   !>  The dimensions of the matrix a and b (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   b           double real matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   b           single real matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   b           double complex matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double complex matrix q: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   b           single complex matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single complex matrix q: on output stores the eigenvalues

   !> \param   is_already_decomposed   logical, input: is it repeated call with the same b (decomposed in the fist call)?
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_generalized_eigenvectors_&
      &f&
      &_i(self, a, b, ev, q, is_already_decomposed, error)
         use, intrinsic :: iso_c_binding
         use elpa_constants
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         real(kind=c_float) :: a(self%local_nrows, *), b(self%local_nrows, *), q(self%local_nrows, *)
         real(kind=c_float) :: ev(self%na)

         logical             :: is_already_decomposed
         integer, optional   :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a generalized eigenvalue problem
   !>
   !>  The dimensions of the matrix a and b (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   b           double real matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   b           single real matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   b           double complex matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   b           single complex matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues

   !> \param   is_already_decomposed   logical, input: is it repeated call with the same b (decomposed in the fist call)?
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_generalized_eigenvalues_&
      &f&
      &_i(self, a, b, ev, is_already_decomposed, error)
         use, intrinsic :: iso_c_binding
         use elpa_constants
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         real(kind=c_float) :: a(self%local_nrows, *), b(self%local_nrows, *)
         real(kind=c_float) :: ev(self%na)

         logical             :: is_already_decomposed
         integer, optional   :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to compute C : = A**T * B
   !>         where   A is a square matrix (self%a,self%na) which is optionally upper or lower triangular
   !>                 B is a (self%na,ncb) matrix
   !>                 C is a (self%na,ncb) matrix where optionally only the upper or lower
   !>                   triangle may be computed
   !>
   !> the MPI commicators are already known to the type. Thus the class method "setup" must be called
   !> BEFORE this method is used
   !> \details
   !>
   !> \param   self                class(elpa_t), the ELPA object
   !> \param  uplo_a               'U' if A is upper triangular
   !>                              'L' if A is lower triangular
   !>                              anything else if A is a full matrix
   !>                              Please note: This pertains to the original A (as set in the calling program)
   !>                                           whereas the transpose of A is used for calculations
   !>                              If uplo_a is 'U' or 'L', the other triangle is not used at all,
   !>                              i.e. it may contain arbitrary numbers
   !> \param uplo_c                'U' if only the upper diagonal part of C is needed
   !>                              'L' if only the upper diagonal part of C is needed
   !>                              anything else if the full matrix C is needed
   !>                              Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
   !>                                            written to a certain extent, i.e. one shouldn't rely on the content there!
   !> \param ncb                   Number of columns  of global matrices B and C
   !> \param a                     matrix a
   !> \param self%local_nrows      number of rows of local (sub) matrix a, set with method set("local_nrows,value")
   !> \param self%local_ncols      number of columns of local (sub) matrix a, set with method set("local_ncols,value")
   !> \param b                     matrix b
   !> \param nrows_b               number of rows of local (sub) matrix b
   !> \param ncols_b               number of columns of local (sub) matrix b
   !> \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
   !> \param c                     matrix c
   !> \param nrows_c               number of rows of local (sub) matrix c
   !> \param ncols_c               number of columns of local (sub) matrix c
   !> \param error                 optional argument, error code which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_hermitian_multiply_&
      &f&
      &_i (self,uplo_a, uplo_c, ncb, a, b, nrows_b, ncols_b, &
         c, nrows_c, ncols_c, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         character*1                     :: uplo_a, uplo_c
         integer(kind=c_int), intent(in) :: nrows_b, ncols_b, nrows_c, ncols_c, ncb
         real(kind=c_float) :: a(self%local_nrows,*), b(nrows_b,*), c(nrows_c,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to do a cholesky decomposition of a matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix: the matrix to be decomposed
   !> \param   a           single real matrix: the matrix to be decomposed
   !> \param   a           double complex matrix: the matrix to be decomposed
   !> \param   a           single complex matrix: the matrix to be decomposed
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_cholesky_&
      &f&
      &_i (self, a, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         real(kind=c_float) :: a(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to invert a triangular matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix: the matrix to be inverted
   !> \param   a           single real matrix: the matrix to be inverted
   !> \param   a           double complex matrix: the matrix to be inverted
   !> \param   a           single complex matrix: the matrix to be inverted

   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_invert_trm_&
      &f&
      &_i (self, a, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         real(kind=c_float) :: a(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve the eigenvalue problem for a valued tridiangular matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   d           double real 1d array: the diagonal elements of a matrix defined in setup, on output the eigenvalues
   !>                      in ascending order
   !> \param   e           double real 1d array: the subdiagonal elements of a matrix defined in setup
   !> \param   q           double real matrix: on output contains the eigenvectors
   !> \param   d           single real 1d array: the diagonal elements of a matrix defined in setup, on output the eigenvalues
   !>                      in ascending order
   !> \param   e           single real 1d array: the subdiagonal elements of a matrix defined in setup
   !> \param   q           single real matrix: on output contains the eigenvectors
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_solve_tridiagonal_&
      &f&
      &_i (self, d, e, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         real(kind=c_float)        :: d(self%na), e(self%na)
         real(kind=c_float)        :: q(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve double real eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvectors
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvectors
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double complex matrix q: on output stores the eigenvectors
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single complex matrix q: on output stores the eigenvectors
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_eigenvectors_&
      &dc&
      &_i(self, a, ev, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self

         complex(kind=c_double) :: a(self%local_nrows, *), q(self%local_nrows,*)
         real(kind=c_double) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve double real skew-symmetric eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvectors
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvectors
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_skew_eigenvectors_&
      &dc&
      &_i(self, a, ev, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self

         complex(kind=c_double) :: a(self%local_nrows, *), q(self%local_nrows,*)
         real(kind=c_double) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_eigenvalues_&
      &dc&
      &_i(self, a, ev, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         complex(kind=c_double) :: a(self%local_nrows, *)
         real(kind=c_double) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a skew-symmetric eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_skew_eigenvalues_&
      &dc&
      &_i(self, a, ev, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         complex(kind=c_double) :: a(self%local_nrows, *)
         real(kind=c_double) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a generalized eigenvalue problem
   !>
   !>  The dimensions of the matrix a and b (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   b           double real matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   b           single real matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   b           double complex matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double complex matrix q: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   b           single complex matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single complex matrix q: on output stores the eigenvalues

   !> \param   is_already_decomposed   logical, input: is it repeated call with the same b (decomposed in the fist call)?
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_generalized_eigenvectors_&
      &dc&
      &_i(self, a, b, ev, q, is_already_decomposed, error)
         use, intrinsic :: iso_c_binding
         use elpa_constants
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         complex(kind=c_double) :: a(self%local_nrows, *), b(self%local_nrows, *), q(self%local_nrows, *)
         real(kind=c_double) :: ev(self%na)

         logical             :: is_already_decomposed
         integer, optional   :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a generalized eigenvalue problem
   !>
   !>  The dimensions of the matrix a and b (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   b           double real matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   b           single real matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   b           double complex matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   b           single complex matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues

   !> \param   is_already_decomposed   logical, input: is it repeated call with the same b (decomposed in the fist call)?
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_generalized_eigenvalues_&
      &dc&
      &_i(self, a, b, ev, is_already_decomposed, error)
         use, intrinsic :: iso_c_binding
         use elpa_constants
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         complex(kind=c_double) :: a(self%local_nrows, *), b(self%local_nrows, *)
         real(kind=c_double) :: ev(self%na)

         logical             :: is_already_decomposed
         integer, optional   :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to compute C : = A**T * B
   !>         where   A is a square matrix (self%a,self%na) which is optionally upper or lower triangular
   !>                 B is a (self%na,ncb) matrix
   !>                 C is a (self%na,ncb) matrix where optionally only the upper or lower
   !>                   triangle may be computed
   !>
   !> the MPI commicators are already known to the type. Thus the class method "setup" must be called
   !> BEFORE this method is used
   !> \details
   !>
   !> \param   self                class(elpa_t), the ELPA object
   !> \param  uplo_a               'U' if A is upper triangular
   !>                              'L' if A is lower triangular
   !>                              anything else if A is a full matrix
   !>                              Please note: This pertains to the original A (as set in the calling program)
   !>                                           whereas the transpose of A is used for calculations
   !>                              If uplo_a is 'U' or 'L', the other triangle is not used at all,
   !>                              i.e. it may contain arbitrary numbers
   !> \param uplo_c                'U' if only the upper diagonal part of C is needed
   !>                              'L' if only the upper diagonal part of C is needed
   !>                              anything else if the full matrix C is needed
   !>                              Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
   !>                                            written to a certain extent, i.e. one shouldn't rely on the content there!
   !> \param ncb                   Number of columns  of global matrices B and C
   !> \param a                     matrix a
   !> \param self%local_nrows      number of rows of local (sub) matrix a, set with method set("local_nrows,value")
   !> \param self%local_ncols      number of columns of local (sub) matrix a, set with method set("local_ncols,value")
   !> \param b                     matrix b
   !> \param nrows_b               number of rows of local (sub) matrix b
   !> \param ncols_b               number of columns of local (sub) matrix b
   !> \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
   !> \param c                     matrix c
   !> \param nrows_c               number of rows of local (sub) matrix c
   !> \param ncols_c               number of columns of local (sub) matrix c
   !> \param error                 optional argument, error code which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_hermitian_multiply_&
      &dc&
      &_i (self,uplo_a, uplo_c, ncb, a, b, nrows_b, ncols_b, &
         c, nrows_c, ncols_c, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         character*1                     :: uplo_a, uplo_c
         integer(kind=c_int), intent(in) :: nrows_b, ncols_b, nrows_c, ncols_c, ncb
         complex(kind=c_double) :: a(self%local_nrows,*), b(nrows_b,*), c(nrows_c,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to do a cholesky decomposition of a matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix: the matrix to be decomposed
   !> \param   a           single real matrix: the matrix to be decomposed
   !> \param   a           double complex matrix: the matrix to be decomposed
   !> \param   a           single complex matrix: the matrix to be decomposed
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_cholesky_&
      &dc&
      &_i (self, a, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         complex(kind=c_double) :: a(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to invert a triangular matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix: the matrix to be inverted
   !> \param   a           single real matrix: the matrix to be inverted
   !> \param   a           double complex matrix: the matrix to be inverted
   !> \param   a           single complex matrix: the matrix to be inverted

   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_invert_trm_&
      &dc&
      &_i (self, a, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         complex(kind=c_double) :: a(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve the eigenvalue problem for a valued tridiangular matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   d           double real 1d array: the diagonal elements of a matrix defined in setup, on output the eigenvalues
   !>                      in ascending order
   !> \param   e           double real 1d array: the subdiagonal elements of a matrix defined in setup
   !> \param   q           double real matrix: on output contains the eigenvectors
   !> \param   d           single real 1d array: the diagonal elements of a matrix defined in setup, on output the eigenvalues
   !>                      in ascending order
   !> \param   e           single real 1d array: the subdiagonal elements of a matrix defined in setup
   !> \param   q           single real matrix: on output contains the eigenvectors
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_solve_tridiagonal_&
      &dc&
      &_i (self, d, e, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         real(kind=c_double)        :: d(self%na), e(self%na)
         real(kind=c_double)        :: q(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve double real eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvectors
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvectors
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double complex matrix q: on output stores the eigenvectors
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single complex matrix q: on output stores the eigenvectors
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_eigenvectors_&
      &fc&
      &_i(self, a, ev, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self

         complex(kind=c_float) :: a(self%local_nrows, *), q(self%local_nrows,*)
         real(kind=c_float) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve double real skew-symmetric eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvectors
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvectors
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_skew_eigenvectors_&
      &fc&
      &_i(self, a, ev, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self

         complex(kind=c_float) :: a(self%local_nrows, *), q(self%local_nrows,*)
         real(kind=c_float) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_eigenvalues_&
      &fc&
      &_i(self, a, ev, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         complex(kind=c_float) :: a(self%local_nrows, *)
         real(kind=c_float) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a skew-symmetric eigenvalue problem
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_skew_eigenvalues_&
      &fc&
      &_i(self, a, ev, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         complex(kind=c_float) :: a(self%local_nrows, *)
         real(kind=c_float) :: ev(self%na)

         integer             :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a generalized eigenvalue problem
   !>
   !>  The dimensions of the matrix a and b (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   b           double real matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double real matrix q: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   b           single real matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single real matrix q: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   b           double complex matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   q           double complex matrix q: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   b           single complex matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   q           single complex matrix q: on output stores the eigenvalues

   !> \param   is_already_decomposed   logical, input: is it repeated call with the same b (decomposed in the fist call)?
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_generalized_eigenvectors_&
      &fc&
      &_i(self, a, b, ev, q, is_already_decomposed, error)
         use, intrinsic :: iso_c_binding
         use elpa_constants
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         complex(kind=c_float) :: a(self%local_nrows, *), b(self%local_nrows, *), q(self%local_nrows, *)
         real(kind=c_float) :: ev(self%na)

         logical             :: is_already_decomposed
         integer, optional   :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve a generalized eigenvalue problem
   !>
   !>  The dimensions of the matrix a and b (locally ditributed and global), the block-cyclic distribution
   !>  blocksize, the number of eigenvectors
   !>  to be computed and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !>  It is possible to change the behaviour of the method by setting tunable parameters with the
   !>  class method "set"
   !> Parameters
   !> \details
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix a: defines the problem to solve
   !> \param   b           double real matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single real matrix a: defines the problem to solve
   !> \param   b           single real matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues
   !> \param   a           double complex matrix a: defines the problem to solve
   !> \param   b           double complex matrix b: defines the problem to solve
   !> \param   ev          double real: on output stores the eigenvalues
   !> \param   a           single complex matrix a: defines the problem to solve
   !> \param   b           single complex matrix b: defines the problem to solve
   !> \param   ev          single real: on output stores the eigenvalues

   !> \param   is_already_decomposed   logical, input: is it repeated call with the same b (decomposed in the fist call)?
   !> \result  error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_generalized_eigenvalues_&
      &fc&
      &_i(self, a, b, ev, is_already_decomposed, error)
         use, intrinsic :: iso_c_binding
         use elpa_constants
         import elpa_t
         implicit none
         class(elpa_t)       :: self
         complex(kind=c_float) :: a(self%local_nrows, *), b(self%local_nrows, *)
         real(kind=c_float) :: ev(self%na)

         logical             :: is_already_decomposed
         integer, optional   :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to compute C : = A**T * B
   !>         where   A is a square matrix (self%a,self%na) which is optionally upper or lower triangular
   !>                 B is a (self%na,ncb) matrix
   !>                 C is a (self%na,ncb) matrix where optionally only the upper or lower
   !>                   triangle may be computed
   !>
   !> the MPI commicators are already known to the type. Thus the class method "setup" must be called
   !> BEFORE this method is used
   !> \details
   !>
   !> \param   self                class(elpa_t), the ELPA object
   !> \param  uplo_a               'U' if A is upper triangular
   !>                              'L' if A is lower triangular
   !>                              anything else if A is a full matrix
   !>                              Please note: This pertains to the original A (as set in the calling program)
   !>                                           whereas the transpose of A is used for calculations
   !>                              If uplo_a is 'U' or 'L', the other triangle is not used at all,
   !>                              i.e. it may contain arbitrary numbers
   !> \param uplo_c                'U' if only the upper diagonal part of C is needed
   !>                              'L' if only the upper diagonal part of C is needed
   !>                              anything else if the full matrix C is needed
   !>                              Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
   !>                                            written to a certain extent, i.e. one shouldn't rely on the content there!
   !> \param ncb                   Number of columns  of global matrices B and C
   !> \param a                     matrix a
   !> \param self%local_nrows      number of rows of local (sub) matrix a, set with method set("local_nrows,value")
   !> \param self%local_ncols      number of columns of local (sub) matrix a, set with method set("local_ncols,value")
   !> \param b                     matrix b
   !> \param nrows_b               number of rows of local (sub) matrix b
   !> \param ncols_b               number of columns of local (sub) matrix b
   !> \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
   !> \param c                     matrix c
   !> \param nrows_c               number of rows of local (sub) matrix c
   !> \param ncols_c               number of columns of local (sub) matrix c
   !> \param error                 optional argument, error code which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_hermitian_multiply_&
      &fc&
      &_i (self,uplo_a, uplo_c, ncb, a, b, nrows_b, ncols_b, &
         c, nrows_c, ncols_c, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         character*1                     :: uplo_a, uplo_c
         integer(kind=c_int), intent(in) :: nrows_b, ncols_b, nrows_c, ncols_c, ncb
         complex(kind=c_float) :: a(self%local_nrows,*), b(nrows_b,*), c(nrows_c,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to do a cholesky decomposition of a matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix: the matrix to be decomposed
   !> \param   a           single real matrix: the matrix to be decomposed
   !> \param   a           double complex matrix: the matrix to be decomposed
   !> \param   a           single complex matrix: the matrix to be decomposed
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_cholesky_&
      &fc&
      &_i (self, a, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         complex(kind=c_float) :: a(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to invert a triangular matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   a           double real matrix: the matrix to be inverted
   !> \param   a           single real matrix: the matrix to be inverted
   !> \param   a           double complex matrix: the matrix to be inverted
   !> \param   a           single complex matrix: the matrix to be inverted

   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_invert_trm_&
      &fc&
      &_i (self, a, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         complex(kind=c_float) :: a(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to solve the eigenvalue problem for a valued tridiangular matrix
   !>
   !>  The dimensions of the matrix a (locally ditributed and global), the block-cylic-distribution
   !>  block size, and the MPI communicators are already known to the object and MUST be set BEFORE
   !>  with the class method "setup"
   !>
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   d           double real 1d array: the diagonal elements of a matrix defined in setup, on output the eigenvalues
   !>                      in ascending order
   !> \param   e           double real 1d array: the subdiagonal elements of a matrix defined in setup
   !> \param   q           double real matrix: on output contains the eigenvectors
   !> \param   d           single real 1d array: the diagonal elements of a matrix defined in setup, on output the eigenvalues
   !>                      in ascending order
   !> \param   e           single real 1d array: the subdiagonal elements of a matrix defined in setup
   !> \param   q           single real matrix: on output contains the eigenvectors
   !> \param   error       integer, optional : error code, which can be queried with elpa_strerr
   abstract interface
      subroutine elpa_solve_tridiagonal_&
      &fc&
      &_i (self, d, e, q, error)
         use, intrinsic :: iso_c_binding
         import elpa_t
         implicit none
         class(elpa_t)                   :: self
         real(kind=c_float)        :: d(self%na), e(self%na)
         real(kind=c_float)        :: q(self%local_nrows,*)

         integer                         :: error
      end subroutine
   end interface

! end of math routines

   !> \brief abstract definition of interface to destroy an ELPA object
   !> Parameters
   !> \param   self        class(elpa_t), the ELPA object
   !> \param   error       integer, optional, the error code
   abstract interface
      subroutine elpa_destroy_i(self, error)
         import elpa_t
         implicit none
         class(elpa_t)                  :: self
         integer, intent(out)           :: error
      end subroutine
   end interface

   !> \brief abstract definition of interface to print the autotuning state
   !> Parameters
   !> \param   self        class(elpa_autotune_t): the ELPA autotune object
   abstract interface
      subroutine elpa_autotune_print_i(self, error)
         import elpa_autotune_t
         implicit none
         class(elpa_autotune_t), intent(in) :: self
         integer, intent(out)               :: error

      end subroutine
   end interface

   !> \brief abstract definition of interface to destroy the autotuning state
   !> Parameters
   !> \param   self        class(elpa_autotune_t): the ELPA autotune object
   abstract interface
      subroutine elpa_autotune_destroy_i(self, error)
         import elpa_autotune_t
         implicit none
         class(elpa_autotune_t), intent(inout) :: self
         integer, intent(out)                  :: error
      end subroutine
   end interface

   abstract interface
      subroutine elpa_creating_from_legacy_api_i(self)
         import elpa_t
         implicit none
         class(elpa_t), intent(inout)          :: self
      end subroutine
   end interface

contains

   !> \brief function to intialize the ELPA library
   !> Parameters
   !> \param   api_version integer: api_version that ELPA should use
   !> \result  error       integer: error code, which can be queried with elpa_strerr
   !
   !c> int elpa_init(int api_version);
   function elpa_init(api_version) result(error) bind(C, name="elpa_init")
      use elpa_utilities, only : error_unit
      use, intrinsic :: iso_c_binding
      integer(kind=c_int), intent(in), value :: api_version
      integer(kind=c_int)                    :: error

      if (earliest_api_version <= api_version .and. api_version <= current_api_version) then
         initDone = .true.
         api_version_set = api_version
         error = ELPA_OK
      else
         write(error_unit, "(a,i0,a)") "ELPA: Error API version ", api_version," is not supported by this library"
         error = ELPA_ERROR_API_VERSION
      endif

   end function

   !> \brief function to check whether the ELPA library has been correctly initialised
   !> Parameters
   !> \result  state      integer: state is either ELPA_OK or ELPA_ERROR, which can be queried with elpa_strerr
   function elpa_initialized() result(state)
      integer :: state
      if (initDone) then
         state = ELPA_OK
      else
         state = ELPA_ERROR_CRITICAL
      endif
   end function

   function elpa_get_api_version() result(api_version)
      integer :: api_version

      api_version = api_version_set
   end function

   !> \brief subroutine to uninit the ELPA library. Does nothing at the moment. Might do sth. later
   !
   !c_no> #ifndef OPTIONAL_C_ERROR_ARGUMENT
   !c_no> void elpa_uninit(int *error);
   !c_no> #endif
   subroutine elpa_uninit_c(error) bind(C, name="elpa_uninit")
      integer(kind=c_int)        :: error
      call elpa_uninit(error)
   end subroutine

   subroutine elpa_uninit(error)
      integer, intent(out)           :: error

      error = ELPA_OK
   end subroutine
   !> \brief helper function for error strings
   !> Parameters
   !> \param   elpa_error  integer: error code to querry
   !> \result  string      string:  error string
   function elpa_strerr(elpa_error) result(string)
      integer, intent(in) :: elpa_error
      character(kind=C_CHAR, len=elpa_strlen_c(elpa_strerr_c(elpa_error))), pointer :: string
      call c_f_pointer(elpa_strerr_c(elpa_error), string)
   end function

   !> \brief helper function for c strings
   !> Parameters
   !> \param   ptr         type(c_ptr)
   !> \result  string      string
   function elpa_c_string(ptr) result(string)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in) :: ptr
      character(kind=c_char, len=elpa_strlen_c(ptr)), pointer :: string
      call c_f_pointer(ptr, string)
   end function

   !> \brief function to convert an integer in its string representation
   !> Parameters
   !> \param   name        string: the key
   !> \param   value       integer: the value correponding to the key
   !> \param   error       integer, optional: error code, which can be queried with elpa_strerr()
   !> \result  string      string: the string representation
   function elpa_int_value_to_string(name, value, error) result(string)
      use elpa_utilities, only : error_unit
      implicit none
      character(kind=c_char, len=*), intent(in) :: name
      integer(kind=c_int), intent(in) :: value
      integer(kind=c_int), intent(out), optional :: error

      character(kind=c_char, len=elpa_int_value_to_strlen_c(name // C_NULL_CHAR, value)), pointer :: string

      integer(kind=c_int) :: actual_error
      type(c_ptr) :: ptr

      actual_error = elpa_int_value_to_string_c(name // C_NULL_CHAR, value, ptr)
      if (c_associated(ptr)) then
         call c_f_pointer(ptr, string)
      else
         nullify(string)
      endif

      if (present(error)) then
         error = actual_error
      else if (actual_error /= ELPA_OK) then
         write(error_unit,'(a,i0,a)') "ELPA: Error converting value '", value, "' to a string for option '" // &
            name // "' and you did not check for errors: " // elpa_strerr(actual_error)
      endif
   end function

   !> \brief function to convert a string in its integer representation:
   !> Parameters
   !> \param   name        string: the key
   !> \param   string      string: the string whose integer representation should be associated with the key
   !> \param   error       integer, optional: error code, which can be queried with elpa_strerr()
   !> \result  value       integer: the integer representation of the string
   function elpa_int_string_to_value(name, string, error) result(value)
      use elpa_generated_fortran_interfaces
      use elpa_utilities, only : error_unit
      implicit none
      character(kind=c_char, len=*), intent(in)         :: name
      character(kind=c_char, len=*), intent(in), target :: string
      integer(kind=c_int), intent(out), optional        :: error
      integer(kind=c_int)                               :: actual_error

      integer(kind=c_int)                               :: value

      actual_error = elpa_int_string_to_value_c(name // C_NULL_CHAR, string // C_NULL_CHAR, value)

      if (present(error)) then
         error = actual_error
      else if (actual_error /= ELPA_OK) then
         write(error_unit,'(a)') "ELPA: Error converting string '" // string // "' to value for option '" // &
            name // "' and you did not check for errors: " // elpa_strerr(actual_error)
      endif
   end function

   !> \brief function to get the number of possible choices for an option
   !> Parameters
   !> \param   option_name string:   the option
   !> \result  number      integer:  the total number of possible values to be chosen
   function elpa_option_cardinality(option_name) result(number)
      use elpa_generated_fortran_interfaces
      character(kind=c_char, len=*), intent(in) :: option_name
      integer                                   :: number
      number = elpa_option_cardinality_c(option_name // C_NULL_CHAR)
   end function

   !> \brief function to enumerate an option
   !> Parameters
   !> \param   option_name string: the option
   !> \param   i           integer
   !> \result  option      integer
   function elpa_option_enumerate(option_name, i) result(option)
      use elpa_generated_fortran_interfaces
      character(kind=c_char, len=*), intent(in) :: option_name
      integer, intent(in)                       :: i
      integer                                   :: option
      option = elpa_option_enumerate_c(option_name // C_NULL_CHAR, i)
   end function

end module
