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

module pack_unpack_complex

  implicit none

  public :: pack_row_complex_cpu_double

  contains

  subroutine pack_row_complex_cpu_double(a, row, n, stripe_width, last_stripe_width, stripe_count)

    use precision

    implicit none

    integer(kind=ik), intent(in) :: stripe_width, last_stripe_width, stripe_count
    complex(kind=ck8), intent(in) :: a(:,:,:)

    complex(kind=ck8) :: row(:)
    integer(kind=ik) :: n, i, noff, nl

    do i=1,stripe_count
      nl = merge(stripe_width, last_stripe_width, i<stripe_count)
      noff = (i-1)*stripe_width
      row(noff+1:noff+nl) = a(1:nl,n,i)
    enddo

  end subroutine pack_row_complex_cpu_double

  subroutine unpack_row_complex_cpu_double(a, row, n, stripe_count, stripe_width, last_stripe_width)

    use precision

    implicit none

    integer(kind=ik), intent(in) :: stripe_count, stripe_width, last_stripe_width, n
    complex(kind=ck8), intent(in) :: row(:)
    complex(kind=ck8) :: a(:,:,:)
    integer(kind=ik) :: i, noff, nl

    do i=1,stripe_count
      nl = merge(stripe_width, last_stripe_width, i<stripe_count)
      noff = (i-1)*stripe_width
      a(1:nl,n,i) = row(noff+1:noff+nl)
    enddo

  end  subroutine unpack_row_complex_cpu_double

end module
