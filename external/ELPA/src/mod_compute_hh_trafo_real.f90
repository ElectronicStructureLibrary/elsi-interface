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

module compute_hh_trafo_real

  use elpa_mpi

  implicit none

  public :: compute_hh_trafo_real_cpu_double

  contains

  subroutine compute_hh_trafo_real_cpu_double(a, stripe_width, stripe_count, &
    a_off, nbw, max_blk_size, bcast_buffer, kernel_flops, kernel_time, off, ncols, istripe, &
    last_stripe_width)

    use precision
    use, intrinsic :: iso_c_binding
    use elpa2_utilities
    use single_hh_trafo_real

    implicit none

    real(kind=c_double), intent(inout) :: kernel_time ! MPI_WTIME always needs double
    integer(kind=lik) :: kernel_flops
    integer(kind=ik), intent(in) :: nbw, max_blk_size
    real(kind=rk8) :: bcast_buffer(nbw,max_blk_size)
    integer(kind=ik), intent(in) :: a_off
    integer(kind=ik), intent(in) :: stripe_width,stripe_count
    integer(kind=ik), intent(in) :: last_stripe_width
    real(kind=rk8), pointer :: a(:,:,:)

    integer(kind=ik) :: off, ncols, istripe
    integer(kind=ik) :: j, nl
    real(kind=rk8) :: w(nbw,6)
    real(kind=c_double) :: ttt ! MPI_WTIME always needs double

    ttt = mpi_wtime()

    nl = merge(stripe_width, last_stripe_width, istripe<stripe_count)

    do j = ncols, 2, -2
      w(:,1) = bcast_buffer(1:nbw,j+off)
      w(:,2) = bcast_buffer(1:nbw,j+off-1)

      call double_hh_trafo_generic_double(a(1,j+off+a_off-1,istripe),w, &
           nbw, nl, stripe_width, nbw)
    enddo

    if (j==1) then
      call single_hh_trafo_real_cpu_double(a(1:stripe_width,1+off+a_off:1+off+a_off+nbw-1,istripe), &
           bcast_buffer(1:nbw,off+1), nbw, nl, stripe_width)
    endif

    kernel_flops = kernel_flops + 4*int(nl,lik)*int(ncols,lik)*int(nbw,lik)
    kernel_time = kernel_time + mpi_wtime()-ttt

  end subroutine compute_hh_trafo_real_cpu_double

end module
