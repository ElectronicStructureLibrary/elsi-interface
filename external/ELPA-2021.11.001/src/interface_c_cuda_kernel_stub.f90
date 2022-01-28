










!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
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
!    https://elpa.mpcdf.mpg.de/
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

! This is a module contains all CUDA C Calls
! it was provided by NVIDIA with their ELPA GPU port and
! adapted for an ELPA release by A.Marek, RZG


module cuda_c_kernel

  implicit none

  interface
    subroutine launch_compute_hh_trafo_c_cuda_kernel_real_double(q, hh, hh_tau, nev, nb, ldq, ncols) &
               bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: nev, nb, ldq, ncols
      integer(kind=c_intptr_t), value :: q
      integer(c_intptr_t), value :: hh_tau ,hh
    end subroutine
  end interface

  interface
    subroutine launch_compute_hh_trafo_c_cuda_sm80_kernel_real_double(q, hh, hh_tau, nev, nb, ldq, ncols) &
               bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: nev, nb, ldq, ncols
      integer(kind=c_intptr_t), value :: q
      integer(c_intptr_t), value :: hh_tau ,hh
    end subroutine
  end interface

  interface
    subroutine launch_compute_hh_trafo_c_cuda_kernel_real_single(q, hh, hh_tau, nev, nb, ldq, ncols) &
               bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: nev, nb, ldq, ncols
      integer(kind=c_intptr_t), value :: q
      integer(c_intptr_t), value :: hh_tau ,hh
    end subroutine
  end interface

  interface
    subroutine launch_compute_hh_trafo_c_cuda_sm80_kernel_real_single(q, hh, hh_tau, nev, nb, ldq, ncols) &
               bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: nev, nb, ldq, ncols
      integer(kind=c_intptr_t), value :: q
      integer(c_intptr_t), value :: hh_tau ,hh
    end subroutine
  end interface


  interface
    subroutine launch_compute_hh_trafo_c_cuda_kernel_complex_double(q, hh, hh_tau, nev, nb, ldq, ncols) &
               bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: nev, nb, ldq, ncols
      integer(kind=c_intptr_t), value :: q
      integer(kind=c_intptr_t), value :: hh_tau ,hh
    end subroutine
  end interface

  !interface
  !  subroutine launch_compute_hh_trafo_c_cuda_sm80_kernel_complex_double(q, hh, hh_tau, nev, nb, ldq, ncols) &
  !             bind(c)
  !    use, intrinsic :: iso_c_binding
  !    implicit none
  !    integer(kind=c_int), value :: nev, nb, ldq, ncols
  !    integer(kind=c_intptr_t), value :: q
  !    integer(kind=c_intptr_t), value :: hh_tau ,hh
  !  end subroutine
  !end interface

  interface
    subroutine launch_compute_hh_trafo_c_cuda_kernel_complex_single(q, hh, hh_tau, nev, nb, ldq, ncols) &
               bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: nev, nb, ldq, ncols
      integer(kind=c_intptr_t), value :: q
      integer(kind=c_intptr_t), value :: hh_tau ,hh
    end subroutine
  end interface

  !interface
  !  subroutine launch_compute_hh_trafo_c_cuda_sm80_kernel_complex_single(q, hh, hh_tau, nev, nb, ldq, ncols) &
  !             bind(c)
  !    use, intrinsic :: iso_c_binding
  !    implicit none
  !    integer(kind=c_int), value :: nev, nb, ldq, ncols
  !    integer(kind=c_intptr_t), value :: q
  !    integer(kind=c_intptr_t), value :: hh_tau ,hh
  !  end subroutine
  !end interface


  interface
    subroutine launch_my_unpack_c_cuda_kernel_real_double(row_count, n_offset, max_idx,stripe_width, a_dim2, stripe_count, &
               l_nev,row_group_dev, a_dev) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: row_count
      integer(kind=c_int), value :: n_offset, max_idx,stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t), value :: a_dev, row_group_dev
    end subroutine
  end interface

  interface
    subroutine launch_my_unpack_c_cuda_kernel_real_single(row_count, n_offset, max_idx,stripe_width, a_dim2, stripe_count, &
               l_nev,row_group_dev, a_dev) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: row_count
      integer(kind=c_int), value :: n_offset, max_idx,stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t), value :: a_dev, row_group_dev
    end subroutine
  end interface

  interface
    subroutine launch_my_pack_c_cuda_kernel_real_double(row_count, n_offset, max_idx,stripe_width, a_dim2, &
               stripe_count, l_nev, a_dev, row_group_dev) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t), value :: a_dev
      integer(kind=c_intptr_t), value :: row_group_dev
    end subroutine
  end interface

  interface
    subroutine launch_my_pack_c_cuda_kernel_real_single(row_count, n_offset, max_idx,stripe_width, a_dim2, stripe_count, &
               l_nev, a_dev, row_group_dev) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t), value :: a_dev
      integer(kind=c_intptr_t), value :: row_group_dev
    end subroutine
  end interface

  interface
    subroutine launch_extract_hh_tau_c_cuda_kernel_real_double(hh, hh_tau, nb, n, is_zero) &
               bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t), value :: hh
      integer(kind=c_intptr_t), value :: hh_tau
      integer(kind=c_int), value :: nb, n
      integer(kind=c_int), value :: is_zero
    end subroutine
  end interface

  interface
    subroutine launch_extract_hh_tau_c_cuda_kernel_real_single(hh, hh_tau, nb, n, is_zero) &
               bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t), value :: hh
      integer(kind=c_intptr_t), value :: hh_tau
      integer(kind=c_int), value :: nb, n
      integer(kind=c_int), value :: is_zero
    end subroutine
  end interface

  interface
    subroutine launch_my_unpack_c_cuda_kernel_complex_double(row_count, n_offset, max_idx, stripe_width, a_dim2, &
               stripe_count, l_nev, row_group_dev, a_dev) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: row_count
      integer(kind=c_int), value :: n_offset, max_idx,stripe_width, a_dim2, stripe_count,l_nev
      integer(kind=c_intptr_t), value :: a_dev, row_group_dev
    end subroutine
  end interface

 interface
    subroutine launch_my_unpack_c_cuda_kernel_complex_single(row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, &
               l_nev, row_group_dev, a_dev) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: row_count
      integer(kind=c_int), value :: n_offset, max_idx,stripe_width, a_dim2, stripe_count,l_nev
      integer(kind=c_intptr_t), value :: a_dev, row_group_dev
    end subroutine
  end interface

  interface
    subroutine launch_my_pack_c_cuda_kernel_complex_double(row_count, n_offset, max_idx,stripe_width,a_dim2, &
               stripe_count, l_nev, a_dev, row_group_dev) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: row_count, n_offset, max_idx, stripe_width, a_dim2,stripe_count, l_nev
      integer(kind=c_intptr_t), value :: a_dev
      integer(kind=c_intptr_t), value :: row_group_dev
    end subroutine
  end interface

  interface
    subroutine launch_my_pack_c_cuda_kernel_complex_single(row_count, n_offset, max_idx,stripe_width,a_dim2, &
               stripe_count, l_nev, a_dev, row_group_dev) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: row_count, n_offset, max_idx, stripe_width, a_dim2,stripe_count, l_nev
      integer(kind=c_intptr_t), value :: a_dev
      integer(kind=c_intptr_t), value :: row_group_dev
    end subroutine
  end interface

  interface
    subroutine launch_extract_hh_tau_c_cuda_kernel_complex_double(hh, hh_tau, nb, n, is_zero) &
               bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t), value :: hh
      integer(kind=c_intptr_t), value :: hh_tau
      integer(kind=c_int), value :: nb, n
      integer(kind=c_int), value :: is_zero
    end subroutine
  end interface

  interface
    subroutine launch_extract_hh_tau_c_cuda_kernel_complex_single(hh, hh_tau, nb, n, is_zero) &
               bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t), value :: hh
      integer(kind=c_intptr_t), value :: hh_tau
      integer(kind=c_int), value :: nb, n
      integer(kind=c_int), value :: is_zero
    end subroutine
  end interface

  contains

    subroutine launch_compute_hh_trafo_cuda_kernel_real_double(q, hh, hh_tau, nev, nb, ldq, ncols)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: nev, nb, ldq, ncols
      integer(kind=c_intptr_t) :: q
      integer(c_intptr_t) :: hh_tau ,hh
    end subroutine

    subroutine launch_compute_hh_trafo_cuda_kernel_real_single(q, hh, hh_tau, nev, nb, ldq, ncols)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: nev, nb, ldq, ncols
      integer(kind=c_intptr_t) :: q
      integer(c_intptr_t) :: hh_tau ,hh
    end subroutine

    subroutine launch_compute_hh_trafo_cuda_kernel_complex_double(q, hh, hh_tau, nev, nb, ldq, ncols)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: nev, nb, ldq, ncols
      integer(kind=c_intptr_t) :: q
      integer(kind=c_intptr_t) :: hh_tau ,hh
    end subroutine

    subroutine launch_compute_hh_trafo_cuda_kernel_complex_single(q, hh, hh_tau, nev, nb, ldq, ncols)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: nev, nb, ldq, ncols
      integer(kind=c_intptr_t) :: q
      integer(kind=c_intptr_t) :: hh_tau ,hh
    end subroutine

    subroutine launch_my_unpack_cuda_kernel_real_double(row_count, n_offset, max_idx,stripe_width, a_dim2, stripe_count, &
               l_nev,row_group_dev, a_dev)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: row_count
      integer(kind=c_int) :: n_offset, max_idx,stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t) :: a_dev, row_group_dev
    end subroutine

    subroutine launch_my_unpack_cuda_kernel_real_single(row_count, n_offset, max_idx,stripe_width, a_dim2, stripe_count, &
               l_nev,row_group_dev, a_dev)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: row_count
      integer(kind=c_int) :: n_offset, max_idx,stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t) :: a_dev, row_group_dev
    end subroutine

    subroutine launch_my_pack_cuda_kernel_real_double(row_count, n_offset, max_idx,stripe_width, a_dim2, &
               stripe_count, l_nev, a_dev, row_group_dev)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t) :: a_dev
      integer(kind=c_intptr_t) :: row_group_dev
    end subroutine

    subroutine launch_my_pack_cuda_kernel_real_single(row_count, n_offset, max_idx,stripe_width, &
               a_dim2, stripe_count, l_nev, a_dev, row_group_dev)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t) :: a_dev
      integer(kind=c_intptr_t) :: row_group_dev
    end subroutine

    subroutine launch_extract_hh_tau_cuda_kernel_real_double(hh, hh_tau, nb, n, is_zero)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t) :: hh
      integer(kind=c_intptr_t) :: hh_tau
      integer(kind=c_int) :: nb, n
      integer(kind=c_int) :: is_zero
    end subroutine

    subroutine launch_extract_hh_tau_cuda_kernel_real_single(hh, hh_tau, nb, n, is_zero)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t) :: hh
      integer(kind=c_intptr_t) :: hh_tau
      integer(kind=c_int) :: nb, n
      integer(kind=c_int) :: is_zero
    end subroutine

    subroutine launch_my_unpack_cuda_kernel_complex_double(row_count, n_offset, max_idx, stripe_width, &
               a_dim2, stripe_count, l_nev, row_group_dev, a_dev)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: row_count
      integer(kind=c_int) :: n_offset, max_idx,stripe_width, a_dim2, stripe_count,l_nev
      integer(kind=c_intptr_t) :: a_dev, row_group_dev
    end subroutine

    subroutine launch_my_unpack_cuda_kernel_complex_single(row_count, n_offset, max_idx, stripe_width, &
               a_dim2, stripe_count, l_nev, row_group_dev, a_dev)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: row_count
      integer(kind=c_int) :: n_offset, max_idx,stripe_width, a_dim2, stripe_count,l_nev
      integer(kind=c_intptr_t) :: a_dev, row_group_dev
    end subroutine

    subroutine launch_my_pack_cuda_kernel_complex_double(row_count, n_offset, max_idx,stripe_width,a_dim2, &
               stripe_count, l_nev, a_dev, row_group_dev)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: row_count, n_offset, max_idx, stripe_width, a_dim2,stripe_count, l_nev
      integer(kind=c_intptr_t) :: a_dev
      integer(kind=c_intptr_t) :: row_group_dev
    end subroutine

    subroutine launch_my_pack_cuda_kernel_complex_single(row_count, n_offset, max_idx,stripe_width,a_dim2, &
               stripe_count, l_nev, a_dev, row_group_dev)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: row_count, n_offset, max_idx, stripe_width, a_dim2,stripe_count, l_nev
      integer(kind=c_intptr_t) :: a_dev
      integer(kind=c_intptr_t) :: row_group_dev
    end subroutine

    subroutine launch_extract_hh_tau_cuda_kernel_complex_double(hh, hh_tau, nb, n, is_zero)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t) :: hh
      integer(kind=c_intptr_t) :: hh_tau
      integer(kind=c_int) :: nb, n
      integer(kind=c_int) :: is_zero
    end subroutine

    subroutine launch_extract_hh_tau_cuda_kernel_complex_single(hh, hh_tau, nb, n, is_zero)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t) :: hh
      integer(kind=c_intptr_t) :: hh_tau
      integer(kind=c_int) :: nb, n
      integer(kind=c_int) :: is_zero
    end subroutine

end module cuda_c_kernel
