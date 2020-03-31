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

module pack_unpack_cpu
   implicit none

   private

   public pack_row_real_cpu_double, unpack_row_real_cpu_double
   public pack_row_complex_cpu_double, unpack_row_complex_cpu_double

   public pack_row_real_cpu_single, unpack_row_real_cpu_single
   public pack_row_complex_cpu_single, unpack_row_complex_cpu_single

contains

   !real double precision

   subroutine pack_row_&
   &real&
   &_cpu_&
   &double &
      (obj, a, row, n, stripe_width,  &
      last_stripe_width, stripe_count)
      use elpa_abstract_impl
      use precision
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj

      integer(kind=ik), intent(in)               :: n, stripe_count, stripe_width
      integer(kind=ik), intent(in)               :: last_stripe_width
      real(kind=c_double), intent(in)     :: a(:,:,:)

      real(kind=c_double)                 :: row(:)

      integer(kind=ik)                           :: i, noff, nl

      call obj%timer%start("pack_row_&
      &real&
      &_cpu" // &
      &"_double" &
         )

      do i=1,stripe_count
         nl = merge(stripe_width, last_stripe_width, i<stripe_count)
         noff = (i-1)*stripe_width
         row(noff+1:noff+nl) = a(1:nl,n,i)
      enddo

      call obj%timer%stop("pack_row_&
      &real&
      &_cpu" // &
      &"_double" &
         )

   end subroutine

   subroutine unpack_row_&
   &real&
   &_cpu_&
   &double &
      (obj, a, row, n, &
      stripe_count, &
      stripe_width, &
      last_stripe_width)
      use elpa_abstract_impl
      use precision
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik), intent(in)               :: n, stripe_count, stripe_width

      integer(kind=ik), intent(in)               :: last_stripe_width
      real(kind=c_double)                 :: a(:,:,:)

      real(kind=c_double), intent(in)     :: row(:)
      integer(kind=ik)                           :: i, noff, nl

      call obj%timer%start("unpack_row_&
      &real&
      &_cpu" // &
      &"_double" &
         )

      do i=1,stripe_count
         nl = merge(stripe_width, last_stripe_width, i<stripe_count)
         noff = (i-1)*stripe_width
         a(1:nl,n,i) = row(noff+1:noff+nl)

      enddo

      call obj%timer%stop("unpack_row_&
      &real&
      &_cpu" // &
      &"_double" &
         )

   end subroutine

   ! real single precision

   subroutine pack_row_&
   &real&
   &_cpu_&
   &single &
      (obj, a, row, n, stripe_width,  &
      last_stripe_width, stripe_count)
      use elpa_abstract_impl
      use precision
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj

      integer(kind=ik), intent(in)               :: n, stripe_count, stripe_width
      integer(kind=ik), intent(in)               :: last_stripe_width
      real(kind=c_float), intent(in)     :: a(:,:,:)

      real(kind=c_float)                 :: row(:)

      integer(kind=ik)                           :: i, noff, nl

      call obj%timer%start("pack_row_&
      &real&
      &_cpu" // &
      &"_single" &
         )

      do i=1,stripe_count
         nl = merge(stripe_width, last_stripe_width, i<stripe_count)
         noff = (i-1)*stripe_width
         row(noff+1:noff+nl) = a(1:nl,n,i)
      enddo

      call obj%timer%stop("pack_row_&
      &real&
      &_cpu" // &
      &"_single" &
         )

   end subroutine

   subroutine unpack_row_&
   &real&
   &_cpu_&
   &single &
      (obj, a, row, n, &
      stripe_count, &
      stripe_width, &
      last_stripe_width)
      use elpa_abstract_impl
      use precision
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik), intent(in)               :: n, stripe_count, stripe_width

      integer(kind=ik), intent(in)               :: last_stripe_width
      real(kind=c_float)                 :: a(:,:,:)

      real(kind=c_float), intent(in)     :: row(:)
      integer(kind=ik)                           :: i, noff, nl

      call obj%timer%start("unpack_row_&
      &real&
      &_cpu" // &
      &"_single" &
         )

      do i=1,stripe_count
         nl = merge(stripe_width, last_stripe_width, i<stripe_count)
         noff = (i-1)*stripe_width
         a(1:nl,n,i) = row(noff+1:noff+nl)

      enddo

      call obj%timer%stop("unpack_row_&
      &real&
      &_cpu" // &
      &"_single" &
         )

   end subroutine

   !complex double precision

   subroutine pack_row_&
   &complex&
   &_cpu_&
   &double &
      (obj, a, row, n, stripe_width,  &
      last_stripe_width, stripe_count)
      use elpa_abstract_impl
      use precision
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj

      integer(kind=ik), intent(in)               :: n, stripe_count, stripe_width
      integer(kind=ik), intent(in)               :: last_stripe_width
      complex(kind=c_double), intent(in)  :: a(:,:,:)

      complex(kind=c_double)              :: row(:)

      integer(kind=ik)                           :: i, noff, nl

      call obj%timer%start("pack_row_&
      &complex&
      &_cpu" // &
      &"_double" &
         )

      do i=1,stripe_count
         nl = merge(stripe_width, last_stripe_width, i<stripe_count)
         noff = (i-1)*stripe_width
         row(noff+1:noff+nl) = a(1:nl,n,i)
      enddo

      call obj%timer%stop("pack_row_&
      &complex&
      &_cpu" // &
      &"_double" &
         )

   end subroutine

   subroutine unpack_row_&
   &complex&
   &_cpu_&
   &double &
      (obj, a, row, n, &
      stripe_count, &
      stripe_width, &
      last_stripe_width)
      use elpa_abstract_impl
      use precision
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik), intent(in)               :: n, stripe_count, stripe_width

      integer(kind=ik), intent(in)               :: last_stripe_width
      complex(kind=c_double)              :: a(:,:,:)

      complex(kind=c_double), intent(in)  :: row(:)
      integer(kind=ik)                           :: i, noff, nl

      call obj%timer%start("unpack_row_&
      &complex&
      &_cpu" // &
      &"_double" &
         )

      do i=1,stripe_count
         nl = merge(stripe_width, last_stripe_width, i<stripe_count)
         noff = (i-1)*stripe_width
         a(1:nl,n,i) = row(noff+1:noff+nl)

      enddo

      call obj%timer%stop("unpack_row_&
      &complex&
      &_cpu" // &
      &"_double" &
         )

   end subroutine

   ! complex single precision

   subroutine pack_row_&
   &complex&
   &_cpu_&
   &single &
      (obj, a, row, n, stripe_width,  &
      last_stripe_width, stripe_count)
      use elpa_abstract_impl
      use precision
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj

      integer(kind=ik), intent(in)               :: n, stripe_count, stripe_width
      integer(kind=ik), intent(in)               :: last_stripe_width
      complex(kind=c_float), intent(in)  :: a(:,:,:)

      complex(kind=c_float)              :: row(:)

      integer(kind=ik)                           :: i, noff, nl

      call obj%timer%start("pack_row_&
      &complex&
      &_cpu" // &
      &"_single" &
         )

      do i=1,stripe_count
         nl = merge(stripe_width, last_stripe_width, i<stripe_count)
         noff = (i-1)*stripe_width
         row(noff+1:noff+nl) = a(1:nl,n,i)
      enddo

      call obj%timer%stop("pack_row_&
      &complex&
      &_cpu" // &
      &"_single" &
         )

   end subroutine

   subroutine unpack_row_&
   &complex&
   &_cpu_&
   &single &
      (obj, a, row, n, &
      stripe_count, &
      stripe_width, &
      last_stripe_width)
      use elpa_abstract_impl
      use precision
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik), intent(in)               :: n, stripe_count, stripe_width

      integer(kind=ik), intent(in)               :: last_stripe_width
      complex(kind=c_float)              :: a(:,:,:)

      complex(kind=c_float), intent(in)  :: row(:)
      integer(kind=ik)                           :: i, noff, nl

      call obj%timer%start("unpack_row_&
      &complex&
      &_cpu" // &
      &"_single" &
         )

      do i=1,stripe_count
         nl = merge(stripe_width, last_stripe_width, i<stripe_count)
         noff = (i-1)*stripe_width
         a(1:nl,n,i) = row(noff+1:noff+nl)

      enddo

      call obj%timer%stop("unpack_row_&
      &complex&
      &_cpu" // &
      &"_single" &
         )

   end subroutine

end module
