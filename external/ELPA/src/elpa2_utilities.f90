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

! ELPA2 -- 2-stage solver for ELPA
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".

module ELPA2_utilities

  use ELPA_utilities
  use precision

  implicit none

  private

  public :: REAL_ELPA_KERNEL_GENERIC, REAL_ELPA_KERNEL_GENERIC_SIMPLE, &
            REAL_ELPA_KERNEL_BGP, REAL_ELPA_KERNEL_BGQ, &
            REAL_ELPA_KERNEL_SSE, REAL_ELPA_KERNEL_SSE_BLOCK2, &
            REAL_ELPA_KERNEL_SSE_BLOCK4, REAL_ELPA_KERNEL_SSE_BLOCK6, &
            REAL_ELPA_KERNEL_AVX_BLOCK2, REAL_ELPA_KERNEL_AVX_BLOCK4, &
            REAL_ELPA_KERNEL_AVX_BLOCK6, REAL_ELPA_KERNEL_AVX2_BLOCK2, &
            REAL_ELPA_KERNEL_AVX2_BLOCK4, REAL_ELPA_KERNEL_AVX2_BLOCK6, &
            REAL_ELPA_KERNEL_AVX512_BLOCK2, REAL_ELPA_KERNEL_AVX512_BLOCK4, &
            REAL_ELPA_KERNEL_AVX512_BLOCK6, REAL_ELPA_KERNEL_GPU, &
            DEFAULT_REAL_ELPA_KERNEL

  public :: COMPLEX_ELPA_KERNEL_GENERIC, COMPLEX_ELPA_KERNEL_GENERIC_SIMPLE, &
            COMPLEX_ELPA_KERNEL_BGP, COMPLEX_ELPA_KERNEL_BGQ, &
            COMPLEX_ELPA_KERNEL_SSE, COMPLEX_ELPA_KERNEL_SSE_BLOCK1, &
            COMPLEX_ELPA_KERNEL_SSE_BLOCK2, COMPLEX_ELPA_KERNEL_AVX_BLOCK1, &
            COMPLEX_ELPA_KERNEL_AVX_BLOCK2, COMPLEX_ELPA_KERNEL_AVX2_BLOCK1, &
            COMPLEX_ELPA_KERNEL_AVX2_BLOCK2, &
            COMPLEX_ELPA_KERNEL_AVX512_BLOCK1, &
            COMPLEX_ELPA_KERNEL_AVX512_BLOCK2, COMPLEX_ELPA_KERNEL_GPU, &
            DEFAULT_COMPLEX_ELPA_KERNEL

  integer, parameter :: number_of_real_kernels = 18
  integer, parameter :: REAL_ELPA_KERNEL_GENERIC = 1
  integer, parameter :: REAL_ELPA_KERNEL_GENERIC_SIMPLE = 2
  integer, parameter :: REAL_ELPA_KERNEL_BGP = 3
  integer, parameter :: REAL_ELPA_KERNEL_BGQ = 4
  integer, parameter :: REAL_ELPA_KERNEL_SSE = 5
  integer, parameter :: REAL_ELPA_KERNEL_SSE_BLOCK2 = 6
  integer, parameter :: REAL_ELPA_KERNEL_SSE_BLOCK4 = 7
  integer, parameter :: REAL_ELPA_KERNEL_SSE_BLOCK6 = 8
  integer, parameter :: REAL_ELPA_KERNEL_AVX_BLOCK2 = 9
  integer, parameter :: REAL_ELPA_KERNEL_AVX_BLOCK4 = 10
  integer, parameter :: REAL_ELPA_KERNEL_AVX_BLOCK6 = 11
  integer, parameter :: REAL_ELPA_KERNEL_AVX2_BLOCK2 = 12
  integer, parameter :: REAL_ELPA_KERNEL_AVX2_BLOCK4 = 13
  integer, parameter :: REAL_ELPA_KERNEL_AVX2_BLOCK6 = 14
  integer, parameter :: REAL_ELPA_KERNEL_AVX512_BLOCK2 = 15
  integer, parameter :: REAL_ELPA_KERNEL_AVX512_BLOCK4 = 16
  integer, parameter :: REAL_ELPA_KERNEL_AVX512_BLOCK6 = 17
  integer, parameter :: REAL_ELPA_KERNEL_GPU = 18
  integer, parameter :: DEFAULT_REAL_ELPA_KERNEL = REAL_ELPA_KERNEL_GENERIC

  integer, parameter :: number_of_complex_kernels = 14
  integer, parameter :: COMPLEX_ELPA_KERNEL_GENERIC = 1
  integer, parameter :: COMPLEX_ELPA_KERNEL_GENERIC_SIMPLE = 2
  integer, parameter :: COMPLEX_ELPA_KERNEL_BGP = 3
  integer, parameter :: COMPLEX_ELPA_KERNEL_BGQ = 4
  integer, parameter :: COMPLEX_ELPA_KERNEL_SSE = 5
  integer, parameter :: COMPLEX_ELPA_KERNEL_SSE_BLOCK1 = 6
  integer, parameter :: COMPLEX_ELPA_KERNEL_SSE_BLOCK2 = 7
  integer, parameter :: COMPLEX_ELPA_KERNEL_AVX_BLOCK1 = 8
  integer, parameter :: COMPLEX_ELPA_KERNEL_AVX_BLOCK2 = 9
  integer, parameter :: COMPLEX_ELPA_KERNEL_AVX2_BLOCK1 = 10
  integer, parameter :: COMPLEX_ELPA_KERNEL_AVX2_BLOCK2 = 11
  integer, parameter :: COMPLEX_ELPA_KERNEL_AVX512_BLOCK1 = 12
  integer, parameter :: COMPLEX_ELPA_KERNEL_AVX512_BLOCK2 = 13
  integer, parameter :: COMPLEX_ELPA_KERNEL_GPU = 14
  integer, parameter :: DEFAULT_COMPLEX_ELPA_KERNEL = COMPLEX_ELPA_KERNEL_GENERIC

end module ELPA2_utilities
