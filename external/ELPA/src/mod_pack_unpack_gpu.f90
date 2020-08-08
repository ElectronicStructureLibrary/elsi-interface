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

module pack_unpack_gpu
   implicit none

   private

   public pack_row_group_real_gpu_double, &
      unpack_and_prepare_row_group_real_gpu_double, &
      extract_hh_tau_real_gpu_double

   public pack_row_group_complex_gpu_double, &
      unpack_and_prepare_row_group_complex_gpu_double, &
      extract_hh_tau_complex_gpu_double

   public pack_row_group_real_gpu_single, &
      unpack_and_prepare_row_group_real_gpu_single, &
      extract_hh_tau_real_gpu_single

   public pack_row_group_complex_gpu_single, &
      unpack_and_prepare_row_group_complex_gpu_single, &
      extract_hh_tau_complex_gpu_single
contains

   !real double precision

   ! Pack a filled row group (i.e. an array of consecutive rows)

   subroutine pack_row_group_&
   &real&
   &_gpu_&
   &double &
      (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev, &
      rows, n_offset, row_count)
      use cuda_c_kernel
      use cuda_functions
      use precision
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t)     :: row_group_dev, a_dev

      integer(kind=ik), intent(in) :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(in) :: n_offset, row_count
      real(kind=c_double)   :: rows(:,:)
      integer(kind=ik)             :: max_idx
      logical                      :: successCUDA

      ! Use many blocks for higher GPU occupancy
      max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

      ! Use one kernel call to pack the entire row group

!      call my_pack_kernel<<<grid_size, stripe_width>>>(n_offset, max_idx, stripe_width, a_dim2, stripe_count, a_dev, row_group_dev)

      call launch_my_pack_gpu_kernel_&
      &real&
      &_&
      &double &
         (row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev, a_dev, row_group_dev)

      ! Issue one single transfer call for all rows (device to host)
!        rows(:, 1 : row_count) = row_group_dev(:, 1 : row_count)

      successCUDA =  cuda_memcpy(int(loc(rows(:, 1: row_count)),kind=c_intptr_t), row_group_dev , row_count * l_nev * size_of_&
      &double&
      &_&
      &real&
      & , cudaMemcpyDeviceToHost)
      if (.not.(successCUDA)) then
         print *,"pack_row_group_&
         &real&
         &_gpu_&
         &double&
         &: error in cudaMemcpy"
         stop 1
      endif

   end subroutine

   ! Unpack a filled row group (i.e. an array of consecutive rows)
   subroutine unpack_row_group_&
   &real&
   &_gpu_&
   &double &
      (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, &
      a_dim2, l_nev, rows, n_offset, row_count)
      use cuda_c_kernel
      use precision
      use, intrinsic :: iso_c_binding
      use cuda_functions
      implicit none
      integer(kind=c_intptr_t)                     :: row_group_dev, a_dev
      integer(kind=ik), intent(in)                 :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(in)                 :: n_offset, row_count
      real(kind=c_double), intent(in) :: rows(:, :)

      integer(kind=ik)                             :: max_idx
      logical                                      :: successCUDA

      ! Use many blocks for higher GPU occupancy
      max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

      ! Issue one single transfer call for all rows (host to device)
!      row_group_dev(:, 1 : row_count) = rows(:, 1 : row_count)

      successCUDA =  cuda_memcpy( row_group_dev , int(loc(rows(1, 1)),kind=c_intptr_t),row_count * l_nev * &
         size_of_&
      &double&
      &_&
      &real&
      &, cudaMemcpyHostToDevice)
      if (.not.(successCUDA)) then
         print *,"unpack_row_group_&
         &real&
         &_gpu_&
         &double&
         &: error in cudaMemcpy"
         stop 1
      endif

      ! Use one kernel call to pack the entire row group
      !        call my_unpack_kernel<<<grid_size, stripe_width>>>(n_offset, max_idx, stripe_width, a_dim2, stripe_count, row_group_dev, a_dev)

      call launch_my_unpack_gpu_kernel_&
      &real&
      &_&
      &double &
         ( row_count, n_offset, max_idx,stripe_width,a_dim2, stripe_count, l_nev, &
         row_group_dev,a_dev)

   end subroutine

   ! This subroutine must be called before queuing the next row for unpacking; it ensures that an unpacking of the current row group
   ! occurs when the queue is full or when the next row belongs to another group
   subroutine unpack_and_prepare_row_group_&
   &real&
   &_gpu_&
   &double &
      (row_group, row_group_dev, a_dev, stripe_count, stripe_width, &
      last_stripe_width, a_dim2, l_nev, row_group_size, nblk,      &
      unpack_idx, next_unpack_idx, force)

      use, intrinsic :: iso_c_binding
      use precision
      implicit none
      real(kind=c_double)      :: row_group(:,:)
      integer(kind=c_intptr_t)        :: row_group_dev, a_dev
      integer(kind=ik), intent(in)    :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(inout) :: row_group_size
      integer(kind=ik), intent(in)    :: nblk
      integer(kind=ik), intent(inout) :: unpack_idx
      integer(kind=ik), intent(in)    :: next_unpack_idx
      logical, intent(in)             :: force

      if (row_group_size == 0) then
         ! Nothing to flush, just prepare for the upcoming row
         row_group_size = 1
      else
         if (force .or. (row_group_size == nblk) .or. (unpack_idx + 1 /= next_unpack_idx)) then
            ! A flush and a reset must be performed
            call unpack_row_group_&
            &real&
            &_gpu_&
            &double&
               (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, &
               a_dim2, l_nev, row_group(:, :), unpack_idx - row_group_size, row_group_size)
            row_group_size = 1
         else
            ! Just prepare for the upcoming row
            row_group_size = row_group_size + 1
         endif
      endif
      ! Always update the index for the upcoming row
      unpack_idx = next_unpack_idx
   end subroutine

   ! The host wrapper for extracting "tau" from the HH reflectors (see the kernel below)
   subroutine extract_hh_tau_&
   &real&
   &_gpu_&
   &double&
   & (bcast_buffer_dev, hh_tau_dev, nbw, n, is_zero)
      use cuda_c_kernel
      use precision
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t) :: bcast_buffer_dev, hh_tau_dev
      integer(kind=ik), value  :: nbw, n
      logical, value           :: is_zero
      integer(kind=ik)         :: val_is_zero
      if (is_zero) then
         val_is_zero = 1
      else
         val_is_zero = 0
      endif

      call launch_extract_hh_tau_gpu_kernel_&
      &real&
      &_&
      &double&
      & (bcast_buffer_dev, hh_tau_dev, nbw, n, val_is_zero)
   end subroutine

   ! real single precision

   ! Pack a filled row group (i.e. an array of consecutive rows)

   subroutine pack_row_group_&
   &real&
   &_gpu_&
   &single &
      (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev, &
      rows, n_offset, row_count)
      use cuda_c_kernel
      use cuda_functions
      use precision
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t)     :: row_group_dev, a_dev

      integer(kind=ik), intent(in) :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(in) :: n_offset, row_count
      real(kind=c_float)   :: rows(:,:)
      integer(kind=ik)             :: max_idx
      logical                      :: successCUDA

      ! Use many blocks for higher GPU occupancy
      max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

      ! Use one kernel call to pack the entire row group

!      call my_pack_kernel<<<grid_size, stripe_width>>>(n_offset, max_idx, stripe_width, a_dim2, stripe_count, a_dev, row_group_dev)

      call launch_my_pack_gpu_kernel_&
      &real&
      &_&
      &single &
         (row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev, a_dev, row_group_dev)

      ! Issue one single transfer call for all rows (device to host)
!        rows(:, 1 : row_count) = row_group_dev(:, 1 : row_count)

      successCUDA =  cuda_memcpy(int(loc(rows(:, 1: row_count)),kind=c_intptr_t), row_group_dev , row_count * l_nev * size_of_&
      &single&
      &_&
      &real&
      & , cudaMemcpyDeviceToHost)
      if (.not.(successCUDA)) then
         print *,"pack_row_group_&
         &real&
         &_gpu_&
         &single&
         &: error in cudaMemcpy"
         stop 1
      endif

   end subroutine

   ! Unpack a filled row group (i.e. an array of consecutive rows)
   subroutine unpack_row_group_&
   &real&
   &_gpu_&
   &single &
      (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, &
      a_dim2, l_nev, rows, n_offset, row_count)
      use cuda_c_kernel
      use precision
      use, intrinsic :: iso_c_binding
      use cuda_functions
      implicit none
      integer(kind=c_intptr_t)                     :: row_group_dev, a_dev
      integer(kind=ik), intent(in)                 :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(in)                 :: n_offset, row_count
      real(kind=c_float), intent(in) :: rows(:, :)

      integer(kind=ik)                             :: max_idx
      logical                                      :: successCUDA

      ! Use many blocks for higher GPU occupancy
      max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

      ! Issue one single transfer call for all rows (host to device)
!      row_group_dev(:, 1 : row_count) = rows(:, 1 : row_count)

      successCUDA =  cuda_memcpy( row_group_dev , int(loc(rows(1, 1)),kind=c_intptr_t),row_count * l_nev * &
         size_of_&
      &single&
      &_&
      &real&
      &, cudaMemcpyHostToDevice)
      if (.not.(successCUDA)) then
         print *,"unpack_row_group_&
         &real&
         &_gpu_&
         &single&
         &: error in cudaMemcpy"
         stop 1
      endif

      ! Use one kernel call to pack the entire row group
      !        call my_unpack_kernel<<<grid_size, stripe_width>>>(n_offset, max_idx, stripe_width, a_dim2, stripe_count, row_group_dev, a_dev)

      call launch_my_unpack_gpu_kernel_&
      &real&
      &_&
      &single &
         ( row_count, n_offset, max_idx,stripe_width,a_dim2, stripe_count, l_nev, &
         row_group_dev,a_dev)

   end subroutine

   ! This subroutine must be called before queuing the next row for unpacking; it ensures that an unpacking of the current row group
   ! occurs when the queue is full or when the next row belongs to another group
   subroutine unpack_and_prepare_row_group_&
   &real&
   &_gpu_&
   &single &
      (row_group, row_group_dev, a_dev, stripe_count, stripe_width, &
      last_stripe_width, a_dim2, l_nev, row_group_size, nblk,      &
      unpack_idx, next_unpack_idx, force)

      use, intrinsic :: iso_c_binding
      use precision
      implicit none
      real(kind=c_float)      :: row_group(:,:)
      integer(kind=c_intptr_t)        :: row_group_dev, a_dev
      integer(kind=ik), intent(in)    :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(inout) :: row_group_size
      integer(kind=ik), intent(in)    :: nblk
      integer(kind=ik), intent(inout) :: unpack_idx
      integer(kind=ik), intent(in)    :: next_unpack_idx
      logical, intent(in)             :: force

      if (row_group_size == 0) then
         ! Nothing to flush, just prepare for the upcoming row
         row_group_size = 1
      else
         if (force .or. (row_group_size == nblk) .or. (unpack_idx + 1 /= next_unpack_idx)) then
            ! A flush and a reset must be performed
            call unpack_row_group_&
            &real&
            &_gpu_&
            &single&
               (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, &
               a_dim2, l_nev, row_group(:, :), unpack_idx - row_group_size, row_group_size)
            row_group_size = 1
         else
            ! Just prepare for the upcoming row
            row_group_size = row_group_size + 1
         endif
      endif
      ! Always update the index for the upcoming row
      unpack_idx = next_unpack_idx
   end subroutine

   ! The host wrapper for extracting "tau" from the HH reflectors (see the kernel below)
   subroutine extract_hh_tau_&
   &real&
   &_gpu_&
   &single&
   & (bcast_buffer_dev, hh_tau_dev, nbw, n, is_zero)
      use cuda_c_kernel
      use precision
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t) :: bcast_buffer_dev, hh_tau_dev
      integer(kind=ik), value  :: nbw, n
      logical, value           :: is_zero
      integer(kind=ik)         :: val_is_zero
      if (is_zero) then
         val_is_zero = 1
      else
         val_is_zero = 0
      endif

      call launch_extract_hh_tau_gpu_kernel_&
      &real&
      &_&
      &single&
      & (bcast_buffer_dev, hh_tau_dev, nbw, n, val_is_zero)
   end subroutine

   !complex double precision

   ! Pack a filled row group (i.e. an array of consecutive rows)

   subroutine pack_row_group_&
   &complex&
   &_gpu_&
   &double &
      (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev, &
      rows, n_offset, row_count)
      use cuda_c_kernel
      use cuda_functions
      use precision
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t)     :: row_group_dev, a_dev

      integer(kind=ik), intent(in) :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(in) :: n_offset, row_count
      complex(kind=c_double) :: rows(:,:)
      integer(kind=ik)             :: max_idx
      logical                      :: successCUDA

      ! Use many blocks for higher GPU occupancy
      max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

      ! Use one kernel call to pack the entire row group

!      call my_pack_kernel<<<grid_size, stripe_width>>>(n_offset, max_idx, stripe_width, a_dim2, stripe_count, a_dev, row_group_dev)

      call launch_my_pack_gpu_kernel_&
      &complex&
      &_&
      &double &
         (row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev, a_dev, row_group_dev)

      ! Issue one single transfer call for all rows (device to host)
!        rows(:, 1 : row_count) = row_group_dev(:, 1 : row_count)

      successCUDA =  cuda_memcpy(int(loc(rows(:, 1: row_count)),kind=c_intptr_t), row_group_dev , row_count * l_nev * size_of_&
      &double&
      &_&
      &complex&
      & , cudaMemcpyDeviceToHost)
      if (.not.(successCUDA)) then
         print *,"pack_row_group_&
         &complex&
         &_gpu_&
         &double&
         &: error in cudaMemcpy"
         stop 1
      endif

   end subroutine

   ! Unpack a filled row group (i.e. an array of consecutive rows)
   subroutine unpack_row_group_&
   &complex&
   &_gpu_&
   &double &
      (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, &
      a_dim2, l_nev, rows, n_offset, row_count)
      use cuda_c_kernel
      use precision
      use, intrinsic :: iso_c_binding
      use cuda_functions
      implicit none
      integer(kind=c_intptr_t)                     :: row_group_dev, a_dev
      integer(kind=ik), intent(in)                 :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(in)                 :: n_offset, row_count
      complex(kind=c_double), intent(in) :: rows(:, :)

      integer(kind=ik)                             :: max_idx
      logical                                      :: successCUDA

      ! Use many blocks for higher GPU occupancy
      max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

      ! Issue one single transfer call for all rows (host to device)
!      row_group_dev(:, 1 : row_count) = rows(:, 1 : row_count)

      successCUDA =  cuda_memcpy( row_group_dev , int(loc(rows(1, 1)),kind=c_intptr_t),row_count * l_nev * &
         size_of_&
      &double&
      &_&
      &complex&
      &, cudaMemcpyHostToDevice)
      if (.not.(successCUDA)) then
         print *,"unpack_row_group_&
         &complex&
         &_gpu_&
         &double&
         &: error in cudaMemcpy"
         stop 1
      endif

      ! Use one kernel call to pack the entire row group
      !        call my_unpack_kernel<<<grid_size, stripe_width>>>(n_offset, max_idx, stripe_width, a_dim2, stripe_count, row_group_dev, a_dev)

      call launch_my_unpack_gpu_kernel_&
      &complex&
      &_&
      &double &
         ( row_count, n_offset, max_idx,stripe_width,a_dim2, stripe_count, l_nev, &
         row_group_dev,a_dev)

   end subroutine

   ! This subroutine must be called before queuing the next row for unpacking; it ensures that an unpacking of the current row group
   ! occurs when the queue is full or when the next row belongs to another group
   subroutine unpack_and_prepare_row_group_&
   &complex&
   &_gpu_&
   &double &
      (row_group, row_group_dev, a_dev, stripe_count, stripe_width, &
      last_stripe_width, a_dim2, l_nev, row_group_size, nblk,      &
      unpack_idx, next_unpack_idx, force)

      use, intrinsic :: iso_c_binding
      use precision
      implicit none
      complex(kind=c_double)   :: row_group(:,:)
      integer(kind=c_intptr_t)        :: row_group_dev, a_dev
      integer(kind=ik), intent(in)    :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(inout) :: row_group_size
      integer(kind=ik), intent(in)    :: nblk
      integer(kind=ik), intent(inout) :: unpack_idx
      integer(kind=ik), intent(in)    :: next_unpack_idx
      logical, intent(in)             :: force

      if (row_group_size == 0) then
         ! Nothing to flush, just prepare for the upcoming row
         row_group_size = 1
      else
         if (force .or. (row_group_size == nblk) .or. (unpack_idx + 1 /= next_unpack_idx)) then
            ! A flush and a reset must be performed
            call unpack_row_group_&
            &complex&
            &_gpu_&
            &double&
               (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, &
               a_dim2, l_nev, row_group(:, :), unpack_idx - row_group_size, row_group_size)
            row_group_size = 1
         else
            ! Just prepare for the upcoming row
            row_group_size = row_group_size + 1
         endif
      endif
      ! Always update the index for the upcoming row
      unpack_idx = next_unpack_idx
   end subroutine

   ! The host wrapper for extracting "tau" from the HH reflectors (see the kernel below)
   subroutine extract_hh_tau_&
   &complex&
   &_gpu_&
   &double&
   & (bcast_buffer_dev, hh_tau_dev, nbw, n, is_zero)
      use cuda_c_kernel
      use precision
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t) :: bcast_buffer_dev, hh_tau_dev
      integer(kind=ik), value  :: nbw, n
      logical, value           :: is_zero
      integer(kind=ik)         :: val_is_zero
      if (is_zero) then
         val_is_zero = 1
      else
         val_is_zero = 0
      endif

      call launch_extract_hh_tau_gpu_kernel_&
      &complex&
      &_&
      &double&
      & (bcast_buffer_dev, hh_tau_dev, nbw, n, val_is_zero)
   end subroutine

   ! complex single precision

   ! Pack a filled row group (i.e. an array of consecutive rows)

   subroutine pack_row_group_&
   &complex&
   &_gpu_&
   &single &
      (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev, &
      rows, n_offset, row_count)
      use cuda_c_kernel
      use cuda_functions
      use precision
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t)     :: row_group_dev, a_dev

      integer(kind=ik), intent(in) :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(in) :: n_offset, row_count
      complex(kind=c_float) :: rows(:,:)
      integer(kind=ik)             :: max_idx
      logical                      :: successCUDA

      ! Use many blocks for higher GPU occupancy
      max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

      ! Use one kernel call to pack the entire row group

!      call my_pack_kernel<<<grid_size, stripe_width>>>(n_offset, max_idx, stripe_width, a_dim2, stripe_count, a_dev, row_group_dev)

      call launch_my_pack_gpu_kernel_&
      &complex&
      &_&
      &single &
         (row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev, a_dev, row_group_dev)

      ! Issue one single transfer call for all rows (device to host)
!        rows(:, 1 : row_count) = row_group_dev(:, 1 : row_count)

      successCUDA =  cuda_memcpy(int(loc(rows(:, 1: row_count)),kind=c_intptr_t), row_group_dev , row_count * l_nev * size_of_&
      &single&
      &_&
      &complex&
      & , cudaMemcpyDeviceToHost)
      if (.not.(successCUDA)) then
         print *,"pack_row_group_&
         &complex&
         &_gpu_&
         &single&
         &: error in cudaMemcpy"
         stop 1
      endif

   end subroutine

   ! Unpack a filled row group (i.e. an array of consecutive rows)
   subroutine unpack_row_group_&
   &complex&
   &_gpu_&
   &single &
      (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, &
      a_dim2, l_nev, rows, n_offset, row_count)
      use cuda_c_kernel
      use precision
      use, intrinsic :: iso_c_binding
      use cuda_functions
      implicit none
      integer(kind=c_intptr_t)                     :: row_group_dev, a_dev
      integer(kind=ik), intent(in)                 :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(in)                 :: n_offset, row_count
      complex(kind=c_float), intent(in) :: rows(:, :)

      integer(kind=ik)                             :: max_idx
      logical                                      :: successCUDA

      ! Use many blocks for higher GPU occupancy
      max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

      ! Issue one single transfer call for all rows (host to device)
!      row_group_dev(:, 1 : row_count) = rows(:, 1 : row_count)

      successCUDA =  cuda_memcpy( row_group_dev , int(loc(rows(1, 1)),kind=c_intptr_t),row_count * l_nev * &
         size_of_&
      &single&
      &_&
      &complex&
      &, cudaMemcpyHostToDevice)
      if (.not.(successCUDA)) then
         print *,"unpack_row_group_&
         &complex&
         &_gpu_&
         &single&
         &: error in cudaMemcpy"
         stop 1
      endif

      ! Use one kernel call to pack the entire row group
      !        call my_unpack_kernel<<<grid_size, stripe_width>>>(n_offset, max_idx, stripe_width, a_dim2, stripe_count, row_group_dev, a_dev)

      call launch_my_unpack_gpu_kernel_&
      &complex&
      &_&
      &single &
         ( row_count, n_offset, max_idx,stripe_width,a_dim2, stripe_count, l_nev, &
         row_group_dev,a_dev)

   end subroutine

   ! This subroutine must be called before queuing the next row for unpacking; it ensures that an unpacking of the current row group
   ! occurs when the queue is full or when the next row belongs to another group
   subroutine unpack_and_prepare_row_group_&
   &complex&
   &_gpu_&
   &single &
      (row_group, row_group_dev, a_dev, stripe_count, stripe_width, &
      last_stripe_width, a_dim2, l_nev, row_group_size, nblk,      &
      unpack_idx, next_unpack_idx, force)

      use, intrinsic :: iso_c_binding
      use precision
      implicit none
      complex(kind=c_float)   :: row_group(:,:)
      integer(kind=c_intptr_t)        :: row_group_dev, a_dev
      integer(kind=ik), intent(in)    :: stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev
      integer(kind=ik), intent(inout) :: row_group_size
      integer(kind=ik), intent(in)    :: nblk
      integer(kind=ik), intent(inout) :: unpack_idx
      integer(kind=ik), intent(in)    :: next_unpack_idx
      logical, intent(in)             :: force

      if (row_group_size == 0) then
         ! Nothing to flush, just prepare for the upcoming row
         row_group_size = 1
      else
         if (force .or. (row_group_size == nblk) .or. (unpack_idx + 1 /= next_unpack_idx)) then
            ! A flush and a reset must be performed
            call unpack_row_group_&
            &complex&
            &_gpu_&
            &single&
               (row_group_dev, a_dev, stripe_count, stripe_width, last_stripe_width, &
               a_dim2, l_nev, row_group(:, :), unpack_idx - row_group_size, row_group_size)
            row_group_size = 1
         else
            ! Just prepare for the upcoming row
            row_group_size = row_group_size + 1
         endif
      endif
      ! Always update the index for the upcoming row
      unpack_idx = next_unpack_idx
   end subroutine

   ! The host wrapper for extracting "tau" from the HH reflectors (see the kernel below)
   subroutine extract_hh_tau_&
   &complex&
   &_gpu_&
   &single&
   & (bcast_buffer_dev, hh_tau_dev, nbw, n, is_zero)
      use cuda_c_kernel
      use precision
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t) :: bcast_buffer_dev, hh_tau_dev
      integer(kind=ik), value  :: nbw, n
      logical, value           :: is_zero
      integer(kind=ik)         :: val_is_zero
      if (is_zero) then
         val_is_zero = 1
      else
         val_is_zero = 0
      endif

      call launch_extract_hh_tau_gpu_kernel_&
      &complex&
      &_&
      &single&
      & (bcast_buffer_dev, hh_tau_dev, nbw, n, val_is_zero)
   end subroutine

end module
