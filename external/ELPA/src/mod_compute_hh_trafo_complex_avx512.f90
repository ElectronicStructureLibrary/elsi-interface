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
! This file was written by A. Marek, MPCDF

module compute_hh_trafo_complex

  use elpa_mpi

  implicit none

  public :: compute_hh_trafo_complex_cpu_double

  contains

  subroutine compute_hh_trafo_complex_cpu_double(a, stripe_width, a_dim2, stripe_count, &
    a_off, nbw, max_blk_size, bcast_buffer, kernel_flops, kernel_time, &
    off, ncols, istripe, last_stripe_width)

    use precision
    use elpa2_utilities
    use, intrinsic :: iso_c_binding
    use kernel_interfaces

    implicit none

    real(kind=c_double), intent(inout) :: kernel_time ! MPI_WTIME always needs double
    integer(kind=lik) :: kernel_flops
    integer(kind=ik), intent(in) :: nbw, max_blk_size
    complex(kind=ck8) :: bcast_buffer(nbw,max_blk_size)
    integer(kind=ik), intent(in) :: a_off

    integer(kind=ik), intent(in) :: stripe_width, a_dim2, stripe_count
    integer(kind=ik), intent(in) :: last_stripe_width
    complex(kind=ck8) :: a(stripe_width,a_dim2,stripe_count)

    integer(kind=ik) :: off, ncols, istripe, j, nl
    real(kind=c_double) :: ttt ! MPI_WTIME always needs double

    nl = merge(stripe_width, last_stripe_width, istripe<stripe_count)

    ttt = mpi_wtime()
    do j = ncols, 1, -1
       call single_hh_trafo_complex_avx512_1hv_double(a(1,j+off+a_off,istripe), &
            bcast_buffer(1,j+off),nbw,nl,stripe_width)
    enddo

    kernel_flops = kernel_flops + 4*4*int(nl,lik)*int(ncols,lik)*int(nbw,lik)
    kernel_time = kernel_time + mpi_wtime()-ttt

  end subroutine compute_hh_trafo_complex_cpu_double

end module
