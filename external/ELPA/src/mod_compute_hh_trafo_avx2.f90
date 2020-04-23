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

module compute_hh_trafo
   use elpa_mpi
   implicit none

   public compute_hh_trafo_real_double

   public compute_hh_trafo_complex_double

   public compute_hh_trafo_real_single

   public compute_hh_trafo_complex_single
contains

   !real double precision

   subroutine compute_hh_trafo_&
   &real&
   &_&
   &double &
      (obj, useGPU, wantDebug, a, a_dev, stripe_width, a_dim2, stripe_count, max_threads, &
      a_off, nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
      hh_tau_dev, kernel_flops, kernel_time, n_times, off, ncols, istripe, &
      last_stripe_width, &
      kernel)

      use precision
      use elpa_abstract_impl
      use, intrinsic :: iso_c_binding

      use single_hh_trafo_real

      use cuda_c_kernel
      use cuda_functions

      use elpa_generated_fortran_interfaces

      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      logical, intent(in)                        :: useGPU, wantDebug
      real(kind=c_double), intent(inout)         :: kernel_time  ! MPI_WTIME always needs double
      integer(kind=lik)                          :: kernel_flops
      integer(kind=ik), intent(in)               :: nbw, max_blk_size
      real(kind=c_double)                 :: bcast_buffer(nbw,max_blk_size)
      integer(kind=ik), intent(in)               :: a_off

      integer(kind=ik), intent(in)               :: stripe_width,a_dim2,stripe_count

      integer(kind=ik), intent(in)               :: max_threads
      integer(kind=ik), intent(in)               :: last_stripe_width
!         real(kind=c_double)                :: a(stripe_width,a_dim2,stripe_count)
      real(kind=c_double), pointer        :: a(:,:,:)

      integer(kind=ik), intent(in)               :: kernel

      integer(kind=c_intptr_t)                   :: a_dev
      integer(kind=c_intptr_t)                   :: bcast_buffer_dev
      integer(kind=c_intptr_t)                   :: hh_tau_dev
      integer(kind=c_intptr_t)                   :: dev_offset, dev_offset_1, dev_offset_2

      ! Private variables in OMP regions (my_thread) should better be in the argument list!
      integer(kind=ik)                           :: off, ncols, istripe
      integer(kind=ik)                           :: j, nl, jj, jjj, n_times
      real(kind=c_double)                 :: w(nbw,6)
      real(kind=c_double)                        :: ttt ! MPI_WTIME always needs double

      integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
      &double&
      &_&
      &real

      j = -99

      if (wantDebug) call obj%timer%start("compute_hh_trafo_&
      &real&
      &" // &
      &"_double" &
         )

      ttt = mpi_wtime()

      nl = merge(stripe_width, last_stripe_width, istripe<stripe_count)

! implementation of avx block 2 real case

      do j = ncols, 2, -2
         w(:,1) = bcast_buffer(1:nbw,j+off)
         w(:,2) = bcast_buffer(1:nbw,j+off-1)
         call double_hh_trafo_&
         &real&
         &_avx_avx2_2hv_&
         &double&
         & (c_loc(a(1,j+off+a_off-1,istripe)), w, nbw, nl, stripe_width, nbw)
      enddo

      if (j==1) call single_hh_trafo_&
      &real&
      &_cpu_&
      &double&
      & (a(1:stripe_width,1+off+a_off:1+off+a_off+nbw-1,istripe), bcast_buffer(1:nbw,off+1), nbw, nl,&
         stripe_width)

      kernel_flops = kernel_flops + 4*int(nl,lik)*int(ncols,lik)*int(nbw,lik)
      kernel_time = kernel_time + mpi_wtime()-ttt
      n_times = n_times + 1

      if (wantDebug) call obj%timer%stop("compute_hh_trafo_&
      &real&
      &" // &
      &"_double" &
         )

   end subroutine

   ! real single precision

   subroutine compute_hh_trafo_&
   &real&
   &_&
   &single &
      (obj, useGPU, wantDebug, a, a_dev, stripe_width, a_dim2, stripe_count, max_threads, &
      a_off, nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
      hh_tau_dev, kernel_flops, kernel_time, n_times, off, ncols, istripe, &
      last_stripe_width, &
      kernel)

      use precision
      use elpa_abstract_impl
      use, intrinsic :: iso_c_binding

      use single_hh_trafo_real

      use cuda_c_kernel
      use cuda_functions

      use elpa_generated_fortran_interfaces

      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      logical, intent(in)                        :: useGPU, wantDebug
      real(kind=c_double), intent(inout)         :: kernel_time  ! MPI_WTIME always needs double
      integer(kind=lik)                          :: kernel_flops
      integer(kind=ik), intent(in)               :: nbw, max_blk_size
      real(kind=c_float)                 :: bcast_buffer(nbw,max_blk_size)
      integer(kind=ik), intent(in)               :: a_off

      integer(kind=ik), intent(in)               :: stripe_width,a_dim2,stripe_count

      integer(kind=ik), intent(in)               :: max_threads
      integer(kind=ik), intent(in)               :: last_stripe_width
!         real(kind=c_float)                :: a(stripe_width,a_dim2,stripe_count)
      real(kind=c_float), pointer        :: a(:,:,:)

      integer(kind=ik), intent(in)               :: kernel

      integer(kind=c_intptr_t)                   :: a_dev
      integer(kind=c_intptr_t)                   :: bcast_buffer_dev
      integer(kind=c_intptr_t)                   :: hh_tau_dev
      integer(kind=c_intptr_t)                   :: dev_offset, dev_offset_1, dev_offset_2

      ! Private variables in OMP regions (my_thread) should better be in the argument list!
      integer(kind=ik)                           :: off, ncols, istripe
      integer(kind=ik)                           :: j, nl, jj, jjj, n_times
      real(kind=c_float)                 :: w(nbw,6)
      real(kind=c_double)                        :: ttt ! MPI_WTIME always needs double

      integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
      &single&
      &_&
      &real

      j = -99

      if (wantDebug) call obj%timer%start("compute_hh_trafo_&
      &real&
      &" // &
      &"_single" &
         )

      ttt = mpi_wtime()

      nl = merge(stripe_width, last_stripe_width, istripe<stripe_count)

! implementation of avx block 2 real case

      do j = ncols, 2, -2
         w(:,1) = bcast_buffer(1:nbw,j+off)
         w(:,2) = bcast_buffer(1:nbw,j+off-1)
         call double_hh_trafo_&
         &real&
         &_avx_avx2_2hv_&
         &single&
         & (c_loc(a(1,j+off+a_off-1,istripe)), w, nbw, nl, stripe_width, nbw)
      enddo

      if (j==1) call single_hh_trafo_&
      &real&
      &_cpu_&
      &single&
      & (a(1:stripe_width,1+off+a_off:1+off+a_off+nbw-1,istripe), bcast_buffer(1:nbw,off+1), nbw, nl,&
         stripe_width)

      kernel_flops = kernel_flops + 4*int(nl,lik)*int(ncols,lik)*int(nbw,lik)
      kernel_time = kernel_time + mpi_wtime()-ttt
      n_times = n_times + 1

      if (wantDebug) call obj%timer%stop("compute_hh_trafo_&
      &real&
      &" // &
      &"_single" &
         )

   end subroutine

   !complex double precision

   subroutine compute_hh_trafo_&
   &complex&
   &_&
   &double &
      (obj, useGPU, wantDebug, a, a_dev, stripe_width, a_dim2, stripe_count, max_threads, &
      a_off, nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
      hh_tau_dev, kernel_flops, kernel_time, n_times, off, ncols, istripe, &
      last_stripe_width, &
      kernel)

      use precision
      use elpa_abstract_impl
      use, intrinsic :: iso_c_binding

      use cuda_c_kernel
      use cuda_functions

      use elpa_generated_fortran_interfaces

      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      logical, intent(in)                        :: useGPU, wantDebug
      real(kind=c_double), intent(inout)         :: kernel_time  ! MPI_WTIME always needs double
      integer(kind=lik)                          :: kernel_flops
      integer(kind=ik), intent(in)               :: nbw, max_blk_size
      complex(kind=c_double)              :: bcast_buffer(nbw,max_blk_size)
      integer(kind=ik), intent(in)               :: a_off

      integer(kind=ik), intent(in)               :: stripe_width,a_dim2,stripe_count

      integer(kind=ik), intent(in)               :: max_threads
      integer(kind=ik), intent(in)               :: last_stripe_width
!          complex(kind=c_double)            :: a(stripe_width,a_dim2,stripe_count)
      complex(kind=c_double),pointer     :: a(:,:,:)

      integer(kind=ik), intent(in)               :: kernel

      integer(kind=c_intptr_t)                   :: a_dev
      integer(kind=c_intptr_t)                   :: bcast_buffer_dev
      integer(kind=c_intptr_t)                   :: hh_tau_dev
      integer(kind=c_intptr_t)                   :: dev_offset, dev_offset_1, dev_offset_2

      ! Private variables in OMP regions (my_thread) should better be in the argument list!
      integer(kind=ik)                           :: off, ncols, istripe
      integer(kind=ik)                           :: j, nl, jj, jjj, n_times
      complex(kind=c_double)              :: w(nbw,2)
      real(kind=c_double)                        :: ttt ! MPI_WTIME always needs double

      integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
      &double&
      &_&
      &complex

      j = -99

      if (wantDebug) call obj%timer%start("compute_hh_trafo_&
      &complex&
      &" // &
      &"_double" &
         )

      ttt = mpi_wtime()

      nl = merge(stripe_width, last_stripe_width, istripe<stripe_count)

! avx block1 complex kernel
      do j = ncols, 1, -1
         call single_hh_trafo_&
         &complex&
         &_avx_avx2_1hv_&
         &double&
         & (c_loc(a(1,j+off+a_off,istripe)), bcast_buffer(1,j+off),nbw,nl,stripe_width)
      enddo

      kernel_flops = kernel_flops + 4*int(nl,lik)*int(ncols,lik)*int(nbw,lik)
      kernel_time = kernel_time + mpi_wtime()-ttt
      n_times = n_times + 1

      if (wantDebug) call obj%timer%stop("compute_hh_trafo_&
      &complex&
      &" // &
      &"_double" &
         )

   end subroutine

   ! complex single precision

   subroutine compute_hh_trafo_&
   &complex&
   &_&
   &single &
      (obj, useGPU, wantDebug, a, a_dev, stripe_width, a_dim2, stripe_count, max_threads, &
      a_off, nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
      hh_tau_dev, kernel_flops, kernel_time, n_times, off, ncols, istripe, &
      last_stripe_width, &
      kernel)

      use precision
      use elpa_abstract_impl
      use, intrinsic :: iso_c_binding

      use cuda_c_kernel
      use cuda_functions

      use elpa_generated_fortran_interfaces

      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      logical, intent(in)                        :: useGPU, wantDebug
      real(kind=c_double), intent(inout)         :: kernel_time  ! MPI_WTIME always needs double
      integer(kind=lik)                          :: kernel_flops
      integer(kind=ik), intent(in)               :: nbw, max_blk_size
      complex(kind=c_float)              :: bcast_buffer(nbw,max_blk_size)
      integer(kind=ik), intent(in)               :: a_off

      integer(kind=ik), intent(in)               :: stripe_width,a_dim2,stripe_count

      integer(kind=ik), intent(in)               :: max_threads
      integer(kind=ik), intent(in)               :: last_stripe_width
!          complex(kind=c_float)            :: a(stripe_width,a_dim2,stripe_count)
      complex(kind=c_float),pointer     :: a(:,:,:)

      integer(kind=ik), intent(in)               :: kernel

      integer(kind=c_intptr_t)                   :: a_dev
      integer(kind=c_intptr_t)                   :: bcast_buffer_dev
      integer(kind=c_intptr_t)                   :: hh_tau_dev
      integer(kind=c_intptr_t)                   :: dev_offset, dev_offset_1, dev_offset_2

      ! Private variables in OMP regions (my_thread) should better be in the argument list!
      integer(kind=ik)                           :: off, ncols, istripe
      integer(kind=ik)                           :: j, nl, jj, jjj, n_times
      complex(kind=c_float)              :: w(nbw,2)
      real(kind=c_double)                        :: ttt ! MPI_WTIME always needs double

      integer(kind=c_intptr_t), parameter        :: size_of_datatype = size_of_&
      &single&
      &_&
      &complex

      j = -99

      if (wantDebug) call obj%timer%start("compute_hh_trafo_&
      &complex&
      &" // &
      &"_single" &
         )

      ttt = mpi_wtime()

      nl = merge(stripe_width, last_stripe_width, istripe<stripe_count)

! avx block1 complex kernel
      do j = ncols, 1, -1
         call single_hh_trafo_&
         &complex&
         &_avx_avx2_1hv_&
         &single&
         & (c_loc(a(1,j+off+a_off,istripe)), bcast_buffer(1,j+off),nbw,nl,stripe_width)
      enddo

      kernel_flops = kernel_flops + 4*int(nl,lik)*int(ncols,lik)*int(nbw,lik)
      kernel_time = kernel_time + mpi_wtime()-ttt
      n_times = n_times + 1

      if (wantDebug) call obj%timer%stop("compute_hh_trafo_&
      &complex&
      &" // &
      &"_single" &
         )

   end subroutine

end module
