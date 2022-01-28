










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
last_stripe_width, kernel)

  use precision
  use elpa_abstract_impl
  use, intrinsic :: iso_c_binding
  use single_hh_trafo_real


!#if defined(WITH_REAL_GENERIC_SIMPLE_BLOCK6_KERNEL) && !(defined(1))
!  use real_generic_simple_block6_kernel !, only : double_hh_trafo_generic_simple
!#endif






  !use cuda_c_kernel
  !use cuda_functions
  !use hip_functions
  use gpu_c_kernel
  use elpa_gpu

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
!  real(kind=c_double)                :: a(stripe_width,a_dim2,stripe_count)
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

  if (wantDebug) then
  endif

  ! intel missing
  if (kernel .eq. ELPA_2STAGE_REAL_NVIDIA_GPU .or. &
      kernel .eq. ELPA_2STAGE_REAL_AMD_GPU) then
    ! ncols - indicates the number of HH reflectors to apply; at least 1 must be available
    if (ncols < 1) then
      if (wantDebug) then
        !print *, "Returning early from compute_hh_trafo"
      endif
      return
    endif
  endif

  if (wantDebug) call obj%timer%start("compute_hh_trafo_&
  &real&
  &" // &
  &"_double" &
  )


    ttt = mpi_wtime()


  nl = merge(stripe_width, last_stripe_width, istripe<stripe_count)

! GPU kernel real
  if (kernel .eq. ELPA_2STAGE_REAL_NVIDIA_GPU .or. &
      kernel .eq. ELPA_2STAGE_REAL_AMD_GPU) then
    if (wantDebug) then
      call obj%timer%start("compute_hh_trafo: GPU")
    endif

    dev_offset = ((a_off+off)*stripe_width+(istripe-1)*stripe_width*a_dim2)*size_of_datatype

    dev_offset_1 = off*nbw*size_of_datatype

    dev_offset_2 = off*size_of_datatype

      call launch_compute_hh_trafo_gpu_kernel_&
           &real&
           &_&
           &double&
           &(a_dev + dev_offset, bcast_buffer_dev + dev_offset_1, &
           hh_tau_dev + dev_offset_2, nl, nbw,stripe_width, ncols)


    if (wantDebug) then
      call obj%timer%stop("compute_hh_trafo: GPU")
    endif

  else ! not CUDA kernel

    if (wantDebug) then
      call obj%timer%start("compute_hh_trafo: CPU")
    endif


      !FORTRAN CODE / X86 INRINISIC CODE / BG ASSEMBLER USING 2 HOUSEHOLDER VECTORS
      ! generic kernel real case



        ! generic simple real kernel



        ! sse assembly kernel real case



        ! no sse, vsx, sparc64 sve block1 real kernel






      !no avx block1 real kernel


      ! no avx512 block1 real kernel
      ! no sve512 block1 real kernel


      ! implementation of sparc64 block 2 real case


      ! implementation of neon_arch64 block 2 real case

      ! implementation of neon_arch64 block 2 real case



      ! implementation of vsx block 2 real case


      ! implementation of sse block 2 real case







      ! implementation of avx block 2 real case






      ! implementation of avx512 block 2 real case


        do j = ncols, 2, -2
          w(:,1) = bcast_buffer(1:nbw,j+off)
          w(:,2) = bcast_buffer(1:nbw,j+off-1)
          call double_hh_trafo_&
          &real&
          &_avx512_2hv_&
          &double&
          & (c_loc(a(1,j+off+a_off-1,istripe)), w, nbw, nl, stripe_width, nbw)
        enddo



! implementation of sve512 block 2 real case












      if (j==1) call single_hh_trafo_&
      &real&
      &_cpu_&
      &double&
      & (a(1:stripe_width,1+off+a_off:1+off+a_off+nbw-1,istripe), bcast_buffer(1:nbw,off+1), nbw, nl,&
               stripe_width)



    ! generic simple block4 real kernel



    !real generic simple block6 kernel



    ! sparc64 block 4 real kernel




    ! neon_arch64 block 4 real kernel



    ! sve128 block 4 real kernel



    ! vsx block4 real kernel



    ! sse block4 real kernel




    ! avx block4 real kernel

    ! avx2 block4 real kernel

    ! sve256 block4 real kernel



    ! avx512 block4 real kernel



! sve512 block4 real kernel





    !sparc64 block6 real kernel


    !neon_arch64 block6 real kernel



    !sve128 block6 real kernel


    !vsx block6 real kernel


    !sse block6 real kernel



    ! avx block6 real kernel


    ! avx2 block6 real kernel


    ! sve256 block6 real kernel




    ! avx512 block6 kernel


! sve512 block6 kernel



    if (wantDebug) then
      call obj%timer%stop("compute_hh_trafo: CPU")
    endif
  endif ! GPU_KERNEL

    kernel_flops = kernel_flops + 4*int(nl,8)*int(ncols,8)*int(nbw,8)
    kernel_time = kernel_time + mpi_wtime()-ttt
    n_times = n_times + 1

  if (wantDebug) call obj%timer%stop("compute_hh_trafo_&
  &real&
  &" // &
  &"_double" &
  )

end subroutine

! vim: syntax=fortran

 ! real single precision



















subroutine compute_hh_trafo_&
&real&
&_&
&single &
(obj, useGPU, wantDebug, a, a_dev, stripe_width, a_dim2, stripe_count, max_threads, &
a_off, nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
hh_tau_dev, kernel_flops, kernel_time, n_times, off, ncols, istripe, &
last_stripe_width, kernel)

  use precision
  use elpa_abstract_impl
  use, intrinsic :: iso_c_binding
  use single_hh_trafo_real


!#if defined(WITH_REAL_GENERIC_SIMPLE_BLOCK6_KERNEL) && !(defined(1))
!  use real_generic_simple_block6_kernel !, only : double_hh_trafo_generic_simple
!#endif






  !use cuda_c_kernel
  !use cuda_functions
  !use hip_functions
  use gpu_c_kernel
  use elpa_gpu

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
!  real(kind=c_float)                :: a(stripe_width,a_dim2,stripe_count)
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

  if (wantDebug) then
  endif

  ! intel missing
  if (kernel .eq. ELPA_2STAGE_REAL_NVIDIA_GPU .or. &
      kernel .eq. ELPA_2STAGE_REAL_AMD_GPU) then
    ! ncols - indicates the number of HH reflectors to apply; at least 1 must be available
    if (ncols < 1) then
      if (wantDebug) then
        !print *, "Returning early from compute_hh_trafo"
      endif
      return
    endif
  endif

  if (wantDebug) call obj%timer%start("compute_hh_trafo_&
  &real&
  &" // &
  &"_single" &
  )


    ttt = mpi_wtime()


  nl = merge(stripe_width, last_stripe_width, istripe<stripe_count)

! GPU kernel real
  if (kernel .eq. ELPA_2STAGE_REAL_NVIDIA_GPU .or. &
      kernel .eq. ELPA_2STAGE_REAL_AMD_GPU) then
    if (wantDebug) then
      call obj%timer%start("compute_hh_trafo: GPU")
    endif

    dev_offset = ((a_off+off)*stripe_width+(istripe-1)*stripe_width*a_dim2)*size_of_datatype

    dev_offset_1 = off*nbw*size_of_datatype

    dev_offset_2 = off*size_of_datatype

      call launch_compute_hh_trafo_gpu_kernel_&
           &real&
           &_&
           &single&
           &(a_dev + dev_offset, bcast_buffer_dev + dev_offset_1, &
           hh_tau_dev + dev_offset_2, nl, nbw,stripe_width, ncols)


    if (wantDebug) then
      call obj%timer%stop("compute_hh_trafo: GPU")
    endif

  else ! not CUDA kernel

    if (wantDebug) then
      call obj%timer%start("compute_hh_trafo: CPU")
    endif


      !FORTRAN CODE / X86 INRINISIC CODE / BG ASSEMBLER USING 2 HOUSEHOLDER VECTORS
      ! generic kernel real case



        ! generic simple real kernel



        ! sse assembly kernel real case



        ! no sse, vsx, sparc64 sve block1 real kernel






      !no avx block1 real kernel


      ! no avx512 block1 real kernel
      ! no sve512 block1 real kernel


      ! implementation of sparc64 block 2 real case


      ! implementation of neon_arch64 block 2 real case

      ! implementation of neon_arch64 block 2 real case



      ! implementation of vsx block 2 real case


      ! implementation of sse block 2 real case







      ! implementation of avx block 2 real case






      ! implementation of avx512 block 2 real case


        do j = ncols, 2, -2
          w(:,1) = bcast_buffer(1:nbw,j+off)
          w(:,2) = bcast_buffer(1:nbw,j+off-1)
          call double_hh_trafo_&
          &real&
          &_avx512_2hv_&
          &single&
          & (c_loc(a(1,j+off+a_off-1,istripe)), w, nbw, nl, stripe_width, nbw)
        enddo



! implementation of sve512 block 2 real case












      if (j==1) call single_hh_trafo_&
      &real&
      &_cpu_&
      &single&
      & (a(1:stripe_width,1+off+a_off:1+off+a_off+nbw-1,istripe), bcast_buffer(1:nbw,off+1), nbw, nl,&
               stripe_width)



    ! generic simple block4 real kernel



    !real generic simple block6 kernel



    ! sparc64 block 4 real kernel




    ! neon_arch64 block 4 real kernel



    ! sve128 block 4 real kernel



    ! vsx block4 real kernel



    ! sse block4 real kernel




    ! avx block4 real kernel

    ! avx2 block4 real kernel

    ! sve256 block4 real kernel



    ! avx512 block4 real kernel



! sve512 block4 real kernel





    !sparc64 block6 real kernel


    !neon_arch64 block6 real kernel



    !sve128 block6 real kernel


    !vsx block6 real kernel


    !sse block6 real kernel



    ! avx block6 real kernel


    ! avx2 block6 real kernel


    ! sve256 block6 real kernel




    ! avx512 block6 kernel


! sve512 block6 kernel



    if (wantDebug) then
      call obj%timer%stop("compute_hh_trafo: CPU")
    endif
  endif ! GPU_KERNEL

    kernel_flops = kernel_flops + 4*int(nl,8)*int(ncols,8)*int(nbw,8)
    kernel_time = kernel_time + mpi_wtime()-ttt
    n_times = n_times + 1

  if (wantDebug) call obj%timer%stop("compute_hh_trafo_&
  &real&
  &" // &
  &"_single" &
  )

end subroutine

! vim: syntax=fortran

  !complex double precision



















subroutine compute_hh_trafo_&
&complex&
&_&
&double &
(obj, useGPU, wantDebug, a, a_dev, stripe_width, a_dim2, stripe_count, max_threads, &
a_off, nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
hh_tau_dev, kernel_flops, kernel_time, n_times, off, ncols, istripe, &
last_stripe_width, kernel)

  use precision
  use elpa_abstract_impl
  use, intrinsic :: iso_c_binding




  !use cuda_c_kernel
  !use cuda_functions
  !use hip_functions
  use gpu_c_kernel
  use elpa_gpu

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
!  complex(kind=c_double)            :: a(stripe_width,a_dim2,stripe_count)
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

  if (wantDebug) then
  endif

  ! intel missing
  if (kernel .eq. ELPA_2STAGE_COMPLEX_NVIDIA_GPU .or. &
      kernel .eq. ELPA_2STAGE_COMPLEX_AMD_GPU) then
    ! ncols - indicates the number of HH reflectors to apply; at least 1 must be available
    if (ncols < 1) then
      if (wantDebug) then
        !print *, "Returning early from compute_hh_trafo"
      endif
      return
    endif
  endif

  if (wantDebug) call obj%timer%start("compute_hh_trafo_&
  &complex&
  &" // &
  &"_double" &
  )


    ttt = mpi_wtime()


  nl = merge(stripe_width, last_stripe_width, istripe<stripe_count)

! GPU kernel complex
  if (kernel .eq. ELPA_2STAGE_COMPLEX_NVIDIA_GPU .or. &
      kernel .eq. ELPA_2STAGE_COMPLEX_AMD_GPU) then
    if (wantDebug) then
      call obj%timer%start("compute_hh_trafo: GPU")
    endif

    dev_offset = ((a_off+off)*stripe_width+(istripe-1)*stripe_width*a_dim2)*size_of_datatype

    dev_offset_1 = off*nbw*size_of_datatype

    dev_offset_2 = off*size_of_datatype

      call launch_compute_hh_trafo_gpu_kernel_&
           &complex&
           &_&
           &double&
           &(a_dev + dev_offset, bcast_buffer_dev + dev_offset_1, &
           hh_tau_dev + dev_offset_2, nl, nbw,stripe_width, ncols)


    if (wantDebug) then
      call obj%timer%stop("compute_hh_trafo: GPU")
    endif

  else ! not CUDA kernel

    if (wantDebug) then
      call obj%timer%start("compute_hh_trafo: CPU")
    endif

      !FORTRAN CODE / X86 INRINISIC CODE / BG ASSEMBLER USING 2 HOUSEHOLDER VECTORS

      ! generic kernel complex case



        ! generic simple complex case





        ! sse assembly kernel complex case



        ! sparc64 block1 complex kernel




      ! vsx block1 complex kernel




      ! sse block1 complex kernel

      ! neon_arch64 block1 complex kernel

      ! sve128 block1 complex kernel




      ! avx block1 complex kernel






      ! avx512 block1 complex kernel

        ttt = mpi_wtime()
        do j = ncols, 1, -1
          call single_hh_trafo_&
          &complex&
          &_avx512_1hv_&
          &double&
          & (c_loc(a(1,j+off+a_off,istripe)), bcast_buffer(1,j+off),nbw,nl,stripe_width)
        enddo


      ! sve512 block1 complex kernel








      ! implementation of sparc64 block 2 complex case



      ! implementation of vsx block 2 complex case


      ! implementation of sse block 2 complex case


      ! implementation of neon_arch64 block 2 complex case


      ! implementation of sve128 block 2 complex case





      ! implementation of avx block 2 complex case

      ! implementation of avx2 block 2 complex case

      ! implementation of sve256 block 2 complex case




! implementation of avx512 block 2 complex case

! implementation of vse512 block 2 complex case




      ! complex bgp/bgq kernel implemented













    !no sse block4 complex kernel


    !no avx block4 complex kernel


    !no avx512 block4 complex kernel
    !no sve512 block4 complex kernel








    ! no sse block6 complex kernel


    !no avx block6 complex kernel


    !no avx512 block6 complex kernel
    !no sve512 block6 complex kernel

    if (wantDebug) then
      call obj%timer%stop("compute_hh_trafo: CPU")
    endif
  endif ! GPU_KERNEL

    kernel_flops = kernel_flops + 4*int(nl,8)*int(ncols,8)*int(nbw,8)
    kernel_time = kernel_time + mpi_wtime()-ttt
    n_times = n_times + 1

  if (wantDebug) call obj%timer%stop("compute_hh_trafo_&
  &complex&
  &" // &
  &"_double" &
  )

end subroutine

! vim: syntax=fortran

 ! complex single precision


















subroutine compute_hh_trafo_&
&complex&
&_&
&single &
(obj, useGPU, wantDebug, a, a_dev, stripe_width, a_dim2, stripe_count, max_threads, &
a_off, nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
hh_tau_dev, kernel_flops, kernel_time, n_times, off, ncols, istripe, &
last_stripe_width, kernel)

  use precision
  use elpa_abstract_impl
  use, intrinsic :: iso_c_binding




  !use cuda_c_kernel
  !use cuda_functions
  !use hip_functions
  use gpu_c_kernel
  use elpa_gpu

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
!  complex(kind=c_float)            :: a(stripe_width,a_dim2,stripe_count)
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

  if (wantDebug) then
  endif

  ! intel missing
  if (kernel .eq. ELPA_2STAGE_COMPLEX_NVIDIA_GPU .or. &
      kernel .eq. ELPA_2STAGE_COMPLEX_AMD_GPU) then
    ! ncols - indicates the number of HH reflectors to apply; at least 1 must be available
    if (ncols < 1) then
      if (wantDebug) then
        !print *, "Returning early from compute_hh_trafo"
      endif
      return
    endif
  endif

  if (wantDebug) call obj%timer%start("compute_hh_trafo_&
  &complex&
  &" // &
  &"_single" &
  )


    ttt = mpi_wtime()


  nl = merge(stripe_width, last_stripe_width, istripe<stripe_count)

! GPU kernel complex
  if (kernel .eq. ELPA_2STAGE_COMPLEX_NVIDIA_GPU .or. &
      kernel .eq. ELPA_2STAGE_COMPLEX_AMD_GPU) then
    if (wantDebug) then
      call obj%timer%start("compute_hh_trafo: GPU")
    endif

    dev_offset = ((a_off+off)*stripe_width+(istripe-1)*stripe_width*a_dim2)*size_of_datatype

    dev_offset_1 = off*nbw*size_of_datatype

    dev_offset_2 = off*size_of_datatype

      call launch_compute_hh_trafo_gpu_kernel_&
           &complex&
           &_&
           &single&
           &(a_dev + dev_offset, bcast_buffer_dev + dev_offset_1, &
           hh_tau_dev + dev_offset_2, nl, nbw,stripe_width, ncols)


    if (wantDebug) then
      call obj%timer%stop("compute_hh_trafo: GPU")
    endif

  else ! not CUDA kernel

    if (wantDebug) then
      call obj%timer%start("compute_hh_trafo: CPU")
    endif

      !FORTRAN CODE / X86 INRINISIC CODE / BG ASSEMBLER USING 2 HOUSEHOLDER VECTORS

      ! generic kernel complex case



        ! generic simple complex case





        ! sse assembly kernel complex case



        ! sparc64 block1 complex kernel




      ! vsx block1 complex kernel




      ! sse block1 complex kernel

      ! neon_arch64 block1 complex kernel

      ! sve128 block1 complex kernel




      ! avx block1 complex kernel






      ! avx512 block1 complex kernel

        ttt = mpi_wtime()
        do j = ncols, 1, -1
          call single_hh_trafo_&
          &complex&
          &_avx512_1hv_&
          &single&
          & (c_loc(a(1,j+off+a_off,istripe)), bcast_buffer(1,j+off),nbw,nl,stripe_width)
        enddo


      ! sve512 block1 complex kernel








      ! implementation of sparc64 block 2 complex case



      ! implementation of vsx block 2 complex case


      ! implementation of sse block 2 complex case


      ! implementation of neon_arch64 block 2 complex case


      ! implementation of sve128 block 2 complex case





      ! implementation of avx block 2 complex case

      ! implementation of avx2 block 2 complex case

      ! implementation of sve256 block 2 complex case




! implementation of avx512 block 2 complex case

! implementation of vse512 block 2 complex case




      ! complex bgp/bgq kernel implemented













    !no sse block4 complex kernel


    !no avx block4 complex kernel


    !no avx512 block4 complex kernel
    !no sve512 block4 complex kernel








    ! no sse block6 complex kernel


    !no avx block6 complex kernel


    !no avx512 block6 complex kernel
    !no sve512 block6 complex kernel

    if (wantDebug) then
      call obj%timer%stop("compute_hh_trafo: CPU")
    endif
  endif ! GPU_KERNEL

    kernel_flops = kernel_flops + 4*int(nl,8)*int(ncols,8)*int(nbw,8)
    kernel_time = kernel_time + mpi_wtime()-ttt
    n_times = n_times + 1

  if (wantDebug) call obj%timer%stop("compute_hh_trafo_&
  &complex&
  &" // &
  &"_single" &
  )

end subroutine

! vim: syntax=fortran

!
!  !complex double precision
!#define COMPLEXCASE 1
!#define DOUBLE_PRECISION 1
!#include "../general/precision_macros.h"
!#include "compute_hh_trafo_complex_gpu.F90"
!#undef COMPLEXCASE
!#undef DOUBLE_PRECISION
!
! ! complex single precision
!#if defined(1)
!#define COMPLEXCASE 1
!#define SINGLE_PRECISION 1
!#include "../general/precision_macros.h"
!#include "compute_hh_trafo_complex_gpu.F90"
!#undef COMPLEXCASE
!#undef SINGLE_PRECISION
!#endif
!
end module
