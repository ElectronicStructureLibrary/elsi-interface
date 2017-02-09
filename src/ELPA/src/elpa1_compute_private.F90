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
!
! ELPA1 -- Faster replacements for ScaLAPACK symmetric eigenvalue routines
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".

module ELPA1_compute

  use elpa_utilities
  use elpa_mpi

  implicit none

  private

  public :: tridiag_real_double
  public :: tridiag_real
  public :: trans_ev_real_double
  public :: trans_ev_real
  public :: solve_tridi_double

  interface tridiag_real
    module procedure tridiag_real_double
  end interface

  interface trans_ev_real
    module procedure trans_ev_real_double
  end interface

  public :: tridiag_complex_double
  public :: tridiag_complex
  public :: trans_ev_complex_double
  public :: trans_ev_complex

  interface tridiag_complex
    module procedure tridiag_complex_double
  end interface

  interface trans_ev_complex
    module procedure trans_ev_complex_double
  end interface

  public :: local_index           ! Get local index of a block cyclic distributed matrix
  public :: least_common_multiple ! Get least common multiple

  public :: hh_transform_real_double
  public :: hh_transform_real
  public :: elpa_reduce_add_vectors_real_double
  public :: elpa_reduce_add_vectors_real
  public :: elpa_transpose_vectors_real_double
  public :: elpa_transpose_vectors_real

  interface hh_transform_real
    module procedure hh_transform_real_double
  end interface

  interface elpa_reduce_add_vectors_real
    module procedure elpa_reduce_add_vectors_real_double
  end interface

  interface elpa_transpose_vectors_real
    module procedure elpa_transpose_vectors_real_double
  end interface

  public :: hh_transform_complex_double
  public :: hh_transform_complex
  public :: elpa_reduce_add_vectors_complex_double
  public :: elpa_reduce_add_vectors_complex
  public :: elpa_transpose_vectors_complex_double
  public :: elpa_transpose_vectors_complex

  interface hh_transform_complex
    module procedure hh_transform_complex_double
  end interface

  interface elpa_reduce_add_vectors_complex
    module procedure elpa_reduce_add_vectors_complex_double
  end interface

  interface elpa_transpose_vectors_complex
    module procedure elpa_transpose_vectors_complex_double
  end interface

  contains

! real double
#define DOUBLE_PRECISION_REAL 1
#define DATATYPE REAL(kind=rk8)
#define BYTESIZE 8
#define REALCASE 1
#include "elpa_transpose_vectors.X90"
#include "elpa_reduce_add_vectors.X90"
#undef DOUBLE_PRECISION_REAL
#undef DATATYPE
#undef BYTESIZE
#undef REALCASE

! complex double
#define DOUBLE_PRECISION_COMPLEX 1
#define DATATYPE COMPLEX(kind=ck8)
#define BYTESIZE 16
#define COMPLEXCASE 1
#include "elpa_transpose_vectors.X90"
#include "elpa_reduce_add_vectors.X90"
#undef DATATYPE
#undef BYTESIZE
#undef COMPLEXCASE
#undef DOUBLE_PRECISION_COMPLEX

! real double
#define DOUBLE_PRECISION_REAL 1
#define REAL_DATATYPE rk8
#include "elpa1_compute_real_template.X90"
#undef DOUBLE_PRECISION_REAL
#undef REAL_DATATYPE

! complex double
#define DOUBLE_PRECISION_COMPLEX 1
#define REAL_DATATYPE rk8
#define COMPLEX_DATATYPE ck8
#include "elpa1_compute_complex_template.X90"
#undef DOUBLE_PRECISION_COMPLEX
#undef REAL_DATATYPE
#undef COMPLEX_DATATYPE

end module ELPA1_compute
