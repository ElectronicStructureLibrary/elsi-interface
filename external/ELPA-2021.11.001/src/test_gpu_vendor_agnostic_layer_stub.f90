












module test_gpu
  !use precision_for_tests
  use precision_for_tests
  use iso_c_binding
!#if TEST_INTEL_GPU == 1
!  use mkl_offload
!#endif
  integer(kind=c_int), parameter :: nvidia_gpu = 1
  integer(kind=c_int), parameter :: amd_gpu = 2
  integer(kind=c_int), parameter :: intel_gpu = 3
  integer(kind=c_int), parameter :: no_gpu = -1
  integer(kind=c_int)            :: use_gpu_vendor
  integer(kind=c_int)            :: gpuHostRegisterDefault    
  integer(kind=c_int)            :: gpuMemcpyHostToDevice    
  integer(kind=c_int)            :: gpuMemcpyDeviceToHost   
  integer(kind=c_int)            :: gpuMemcpyDeviceToDevice
  integer(kind=c_int)            :: gpuHostRegisterMapped
  integer(kind=c_int)            :: gpuHostRegisterPortable

  integer(kind=c_intptr_t), parameter :: size_of_double_real    = 8_rk8
  integer(kind=c_intptr_t), parameter :: size_of_single_real    = 4_rk4

  integer(kind=c_intptr_t), parameter :: size_of_double_complex = 16_ck8
  integer(kind=c_intptr_t), parameter :: size_of_single_complex = 8_ck4

  interface gpu_memcpy
    module procedure gpu_memcpy_intptr
    module procedure gpu_memcpy_cptr
    module procedure gpu_memcpy_mixed
  end interface

  interface gpu_malloc
    module procedure gpu_malloc_intptr
    module procedure gpu_malloc_cptr
  end interface

  interface gpu_free
    module procedure gpu_free_intptr
    module procedure gpu_free_cptr
  end interface

  contains
    function gpu_vendor(set_vendor) result(vendor)
      use precision_for_tests
      implicit none
      integer(kind=c_int)             :: vendor
      integer(kind=c_int), intent(in) :: set_vendor
      ! default
      vendor = no_gpu
      if (set_vendor == nvidia_gpu) then
        vendor = nvidia_gpu
      endif
      if (set_vendor == amd_gpu) then
        vendor = amd_gpu
      endif
!#if TEST_INTEL_GPU == 1
!      vendor = intel_gpu
!#endif
      use_gpu_vendor = vendor
      return
    end function

    subroutine set_gpu_parameters
      implicit none



    end subroutine

    function gpu_malloc_intptr(array, elements) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T)             :: array
      integer(kind=c_intptr_t), intent(in) :: elements
      logical                              :: success


    end function

    function gpu_malloc_cptr(array, elements) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)                          :: array
      integer(kind=c_intptr_t), intent(in) :: elements
      logical                              :: success
      success = .false.



    end function

    function gpu_memcpy_intptr(dst, src, size, dir) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_t)              :: dst
      integer(kind=C_intptr_t)              :: src
      integer(kind=c_intptr_t), intent(in)  :: size
      integer(kind=C_INT), intent(in)       :: dir
      logical :: success



    end function
    
    function gpu_memcpy_cptr(dst, src, size, dir) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)                           :: dst
      type(c_ptr)                           :: src
      integer(kind=c_intptr_t), intent(in)  :: size
      integer(kind=C_INT), intent(in)       :: dir
      logical :: success



    end function
    
    function gpu_memcpy_mixed(dst, src, size, dir) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)                           :: dst
      integer(kind=C_intptr_t)              :: src
      integer(kind=c_intptr_t), intent(in)  :: size
      integer(kind=C_INT), intent(in)       :: dir
      logical :: success



    end function

    function gpu_free_intptr(a) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_intptr_t)                :: a

      logical :: success



    end function

    function gpu_free_cptr(a) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)                :: a

      logical :: success



    end function

end module 
