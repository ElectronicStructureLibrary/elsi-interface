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
!    http://elpa.rzg.mpg.de/
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

!This is a module contains all CUDA C Calls
! it was provided by NVIDIA with their ELPA GPU port and
! adapted for an ELPA release by A.Marek, RZG

module cuda_c_kernel
  implicit none

  interface
    subroutine launch_dot_product_kernel_c_double(hs_dev, hv_new_dev, tau_new, x_dev, h_dev,hv_dev, nr) &
               bind(c,name="launch_dot_product_kernel_double")
      use precision
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value      :: nr
      integer(kind=C_intptr_T), value :: hs_dev ,hv_new_dev,x_dev,h_dev, hv_dev
      complex(kind=ck8),value         :: tau_new

    end subroutine
  end interface



  interface
    subroutine launch_dot_product_kernel_1_c_double(ab_dev, hs_dev, hv_new_dev, x_dev,h_dev,hv_dev,nb, nr, ns) &
               bind(c, name="launch_dot_product_kernel_1_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value        ::  nb, nr, ns
      integer(kind=C_intptr_T), value   :: x_dev,h_dev, hv_dev, ab_dev, hs_dev,hv_new_dev

    end subroutine
  end interface



  interface
    subroutine launch_dot_product_kernel_2_c_double(ab_dev, hs_dev, hv_dev,hd_dev,nb, nr, ne) &
                bind(c,name="launch_dot_product_kernel_2_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value      ::  nb, nr, ne
      integer(kind=C_intptr_T), value :: hd_dev,hv_dev, hs_dev, ab_dev

    end subroutine
  end interface



  interface
    subroutine launch_double_hh_transform_1_c_double(ab_dev, hs_dev,hv_dev,nb,ns) &
               bind(c,name="launch_double_hh_transform_1_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value      ::  nb, ns
      integer(kind=C_intptr_T), value :: hv_dev, ab_dev,hs_dev

    end subroutine
  end interface



  interface
    subroutine launch_double_hh_transform_2_c_double(ab_dev, hd_dev,hv_dev,nc,ns, nb) &
               bind(c,name="launch_double_hh_transform_2_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value      ::  nc, ns, nb
      integer(kind=C_intptr_T), value :: hv_dev, ab_dev,hd_dev
    end subroutine
  end interface



  interface
    subroutine launch_compute_kernel_reduce_c_double(a_dev, lda, n, nbw, h1_dev) &
               bind(c,name="launch_compute_kernel_reduce_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value      :: n,lda,nbw
      integer(kind=C_intptr_T), value :: h1_dev ,a_dev
    end subroutine
  end interface



  interface
    subroutine launch_compute_kernel_reduce_1_c_double(a_dev, lda, n, h1_dev) &
               bind(c,name="launch_compute_kernel_reduce_1_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value      :: n,lda
      integer(kind=C_intptr_T), value :: h1_dev ,a_dev

    end subroutine
  end interface



  interface
    subroutine launch_compute_hh_trafo_c_kernel_real_c_double(q, hh, hh_dot, hh_tau, nev, nb, ldq, off, ncols) &
               bind(c,name="launch_compute_hh_trafo_c_kernel_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value     :: nev, nb, ldq, off, ncols
      integer(kind=c_size_t), value  :: q
      integer(kind=c_size_t), value  :: hh_dot
      integer(C_SIZE_T), value       :: hh_tau ,hh
    end subroutine
  end interface



  interface
    subroutine launch_compute_hh_trafo_c_kernel_complex_c_double(q, hh, hh_tau, nev, nb,ldq,off, ncols) &
               bind(c,name="launch_compute_hh_trafo_c_kernel_complex_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value    :: nev, nb, ldq, off, ncols
      integer(kind=c_size_t), value :: q
      integer(kind=c_size_t), value :: hh_tau ,hh
    end subroutine
  end interface





  interface
    subroutine launch_my_unpack_c_kernel_real_c_double(row_count, n_offset, max_idx,stripe_width, a_dim2, stripe_count, &
                                                l_nev,row_group_dev, a_dev) bind(c,name="launch_my_unpack_c_kernel_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value       :: row_count
      integer(kind=c_int), value       :: n_offset, max_idx,stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t), value  :: a_dev, row_group_dev

    end subroutine
  end interface



  interface
    subroutine launch_my_pack_c_kernel_real_c_double(row_count, n_offset, max_idx,stripe_width, a_dim2, &
                                                     stripe_count, l_nev, a_dev, &
                                                     row_group_dev) bind(c,name="launch_my_pack_c_kernel_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value      :: row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t), value :: a_dev
      integer(kind=c_intptr_t), value :: row_group_dev

    end subroutine
  end interface



  interface
    subroutine launch_compute_hh_dotp_c_kernel_real_c_double(bcast_buffer_dev, hh_dot_dev, nbw, n) &
               bind(c,name="launch_compute_hh_dotp_c_kernel_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_intptr_t), value :: bcast_buffer_dev
      integer(kind=c_intptr_t), value :: hh_dot_dev
      integer(kind=c_int), value      :: nbw, n

    end subroutine
  end interface



  interface
    subroutine launch_extract_hh_tau_c_kernel_real_c_double(hh, hh_tau, nb, n, is_zero) &
               bind(c,NAME="launch_extract_hh_tau_c_kernel_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_size_t), value :: hh
      integer(kind=c_size_t), value :: hh_tau
      integer(kind=c_int), value    :: nb, n
      integer(kind=c_int), value    :: is_zero

    end subroutine
  end interface


  interface
    subroutine launch_my_unpack_c_kernel_complex_c_double(row_count, n_offset, max_idx, stripe_width, a_dim2, &
                                                          stripe_count, l_nev, &
                                                 row_group_dev, a_dev) bind(c,name="launch_my_unpack_c_kernel_complex_double")

      use, intrinsic :: iso_c_binding

      implicit none

      integer(kind=c_int), value       :: row_count
      integer(kind=c_int), value       :: n_offset, max_idx,stripe_width, a_dim2, stripe_count,l_nev
      integer(kind=c_intptr_t), value  :: a_dev, row_group_dev

    end subroutine
  end interface



  interface
    subroutine launch_my_pack_c_kernel_complex_c_double(row_count, n_offset, max_idx,stripe_width,a_dim2, &
                                                        stripe_count, l_nev, a_dev, &
                                               row_group_dev) bind(c,name="launch_my_pack_c_kernel_complex_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int), value      :: row_count, n_offset, max_idx, stripe_width, a_dim2,stripe_count, l_nev
      integer(kind=c_intptr_t), value :: a_dev
      integer(kind=c_intptr_t), value :: row_group_dev

    end subroutine
  end interface



  interface
   subroutine launch_compute_hh_dotp_c_kernel_complex_c_double(bcast_buffer_dev, hh_dot_dev, nbw,n) &
              bind(c,name="launch_compute_hh_dotp_c_kernel_complex_double")

     use, intrinsic :: iso_c_binding

     implicit none
     integer(kind=c_intptr_t), value :: bcast_buffer_dev
     integer(kind=c_intptr_t), value :: hh_dot_dev
     integer(kind=c_int), value      :: nbw, n
   end subroutine
  end interface



  interface
    subroutine launch_extract_hh_tau_c_kernel_complex_c_double(hh, hh_tau, nb, n, is_zero) &
               bind(c,name="launch_extract_hh_tau_c_kernel_complex_double")

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_size_t), value :: hh
      integer(kind=c_size_t), value :: hh_tau
      integer(kind=c_int), value    :: nb, n
      integer(kind=c_int), value    :: is_zero

    end subroutine
  end interface



  contains

    subroutine launch_dot_product_kernel_double(hs_dev, hv_new_dev, tau_new, x_dev, h_dev,hv_dev, nr)

      use, intrinsic :: iso_c_binding
      use precision
      implicit none
      integer(kind=c_int)      :: nr
      integer(kind=C_intptr_T) :: hs_dev ,hv_new_dev,x_dev,h_dev, hv_dev
      complex(kind=ck8)         :: tau_new

    end subroutine



    subroutine launch_dot_product_kernel_1_double(ab_dev, hs_dev, hv_new_dev, x_dev,h_dev,hv_dev,nb, nr, ns)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int)      ::  nb, nr, ns
      integer(kind=C_intptr_T) :: x_dev,h_dev, hv_dev, ab_dev, hs_dev,hv_new_dev

    end subroutine



    subroutine launch_dot_product_kernel_2_double(ab_dev, hs_dev, hv_dev,hd_dev,nb, nr, ne)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int)      ::  nb, nr, ne
      integer(kind=C_intptr_T) :: hd_dev,hv_dev, hs_dev, ab_dev

    end subroutine



    subroutine launch_double_hh_transform_1_double(ab_dev, hs_dev,hv_dev,nb,ns)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int)      ::  nb, ns
      integer(kind=C_intptr_T) :: hv_dev, ab_dev,hs_dev

    end subroutine



    subroutine launch_double_hh_transform_2_double(ab_dev, hd_dev,hv_dev,nc,ns, nb)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int)      ::  nc, ns, nb
      integer(kind=C_intptr_T) :: hv_dev, ab_dev,hd_dev

    end subroutine



    subroutine launch_compute_kernel_reduce_double(a_dev, lda, n, nbw, h1_dev)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int)      :: n,lda,nbw
      integer(kind=C_intptr_T) :: h1_dev ,a_dev

    end subroutine



    subroutine launch_compute_kernel_reduce_1_double(a_dev, lda, n, h1_dev)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int)      :: n,lda
      integer(kind=C_intptr_T) :: h1_dev ,a_dev

    end subroutine



    subroutine launch_compute_hh_trafo_c_kernel_real_double(q, hh, hh_dot, hh_tau, nev, nb, ldq, off, ncols)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int)     :: nev, nb, ldq, off, ncols
      integer(kind=c_size_t)  :: q
      integer(kind=c_size_t)  :: hh_dot
      integer(C_SIZE_T)       :: hh_tau ,hh

    end subroutine



    subroutine launch_compute_hh_trafo_c_kernel_complex_double(q, hh, hh_tau, nev, nb,ldq,off, ncols)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int)    :: nev, nb, ldq, off, ncols
      integer(kind=c_size_t) :: q
      integer(kind=c_size_t) :: hh_tau ,hh

    end subroutine





    subroutine launch_my_unpack_c_kernel_real_double(row_count, n_offset, max_idx,stripe_width, a_dim2, stripe_count, &
                                              l_nev,row_group_dev, a_dev)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int)       :: row_count
      integer(kind=c_int)       :: n_offset, max_idx,stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t)  :: a_dev, row_group_dev


    end subroutine



    subroutine launch_my_pack_c_kernel_real_double(row_count, n_offset, max_idx,stripe_width, a_dim2, stripe_count, l_nev, a_dev, &
                                       row_group_dev)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int)      :: row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev
      integer(kind=c_intptr_t) :: a_dev
      integer(kind=c_intptr_t) :: row_group_dev


    end subroutine



    subroutine launch_compute_hh_dotp_c_kernel_real_double(bcast_buffer_dev, hh_dot_dev, nbw, n)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_intptr_t) :: bcast_buffer_dev
      integer(kind=c_intptr_t) :: hh_dot_dev
      integer(kind=c_int)      :: nbw, n

    end subroutine



    subroutine launch_extract_hh_tau_c_kernel_real_double(hh, hh_tau, nb, n, is_zero)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_size_t) :: hh
      integer(kind=c_size_t) :: hh_tau
      integer(kind=c_int)    :: nb, n
      integer(kind=c_int)    :: is_zero

    end subroutine



    subroutine launch_my_unpack_c_kernel_complex_double(row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev, &
                                                 row_group_dev, a_dev)

      use, intrinsic :: iso_c_binding

      implicit none

      integer(kind=c_int)       :: row_count
      integer(kind=c_int)       :: n_offset, max_idx,stripe_width, a_dim2, stripe_count,l_nev
      integer(kind=c_intptr_t)  :: a_dev, row_group_dev

    end subroutine



    subroutine launch_my_pack_c_kernel_complex_double(row_count, n_offset, max_idx,stripe_width,a_dim2, &
                                                      stripe_count, l_nev, a_dev, &
                                               row_group_dev)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_int)      :: row_count, n_offset, max_idx, stripe_width, a_dim2,stripe_count, l_nev
      integer(kind=c_intptr_t) :: a_dev
      integer(kind=c_intptr_t) :: row_group_dev

    end subroutine



    subroutine launch_compute_hh_dotp_c_kernel_complex_double(bcast_buffer_dev, hh_dot_dev, nbw,n)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_intptr_t) :: bcast_buffer_dev
      integer(kind=c_intptr_t) :: hh_dot_dev
      integer(kind=c_int)      :: nbw, n

    end subroutine



    subroutine launch_extract_hh_tau_c_kernel_complex_double(hh, hh_tau, nb, n, is_zero)

      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=c_size_t) :: hh
      integer(kind=c_size_t) :: hh_tau
      integer(kind=c_int)    :: nb, n
      integer(kind=c_int)    :: is_zero

    end subroutine


end module cuda_c_kernel

