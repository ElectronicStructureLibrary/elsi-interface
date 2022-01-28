












module elpa_gpu
  use precision
  use iso_c_binding
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

  interface gpu_memcpy2d
    module procedure gpu_memcpy2d_intptr
    module procedure gpu_memcpy2d_cptr
  end interface

  interface gpublas_dcopy
    module procedure gpublas_dcopy_intptr
    module procedure gpublas_dcopy_cptr
  end interface

  interface gpublas_scopy
    module procedure gpublas_scopy_intptr
    module procedure gpublas_scopy_cptr
  end interface

  interface gpublas_zcopy
    module procedure gpublas_zcopy_intptr
    module procedure gpublas_zcopy_cptr
  end interface

  interface gpublas_ccopy
    module procedure gpublas_ccopy_intptr
    module procedure gpublas_ccopy_cptr
  end interface

  interface gpublas_dtrmm
    module procedure gpublas_dtrmm_intptr
    module procedure gpublas_dtrmm_cptr
  end interface
  
  interface gpublas_strmm
    module procedure gpublas_strmm_intptr
    module procedure gpublas_strmm_cptr
  end interface

  interface gpublas_ztrmm
    module procedure gpublas_ztrmm_intptr
    module procedure gpublas_ztrmm_cptr
  end interface

  interface gpublas_ctrmm
    module procedure gpublas_ctrmm_intptr
    module procedure gpublas_ctrmm_cptr
  end interface

  interface gpublas_dtrsm
    module procedure gpublas_dtrsm_intptr
    module procedure gpublas_dtrsm_cptr
  end interface
  
  interface gpublas_strsm
    module procedure gpublas_strsm_intptr
    module procedure gpublas_strsm_cptr
  end interface

  interface gpublas_ztrsm
    module procedure gpublas_ztrsm_intptr
    module procedure gpublas_ztrsm_cptr
  end interface

  interface gpublas_ctrsm
    module procedure gpublas_ctrsm_intptr
    module procedure gpublas_ctrsm_cptr
  end interface


  contains
    function gpu_vendor() result(vendor)
      use precision
      implicit none
      integer(kind=c_int) :: vendor
      ! default
      vendor = no_gpu
      use_gpu_vendor = vendor
      return
    end function

    subroutine set_gpu_parameters
      use cuda_functions
      use hip_functions
      implicit none

      if (use_gpu_vendor == nvidia_gpu) then
        cudaMemcpyHostToDevice   = cuda_memcpyHostToDevice()
        gpuMemcpyHostToDevice    = cudaMemcpyHostToDevice
        cudaMemcpyDeviceToHost   = cuda_memcpyDeviceToHost()
        gpuMemcpyDeviceToHost    = cudaMemcpyDeviceToHost
        cudaMemcpyDeviceToDevice = cuda_memcpyDeviceToDevice()
        gpuMemcpyDeviceToDevice  = cudaMemcpyDeviceToDevice
        cudaHostRegisterPortable = cuda_hostRegisterPortable()
        gpuHostRegisterPortable  = cudaHostRegisterPortable
        cudaHostRegisterMapped   = cuda_hostRegisterMapped()
        gpuHostRegisterMapped    = cudaHostRegisterMapped
        cudaHostRegisterDefault  = cuda_hostRegisterDefault()
        gpuHostRegisterDefault   = cudaHostRegisterDefault
      endif

      if (use_gpu_vendor == amd_gpu) then
        hipMemcpyHostToDevice   = hip_memcpyHostToDevice()
        gpuMemcpyHostToDevice   = hipMemcpyHostToDevice
        hipMemcpyDeviceToHost   = hip_memcpyDeviceToHost()
        gpuMemcpyDeviceToHost   = hipMemcpyDeviceToHost
        hipMemcpyDeviceToDevice = hip_memcpyDeviceToDevice()
        gpuMemcpyDeviceToDevice = hipMemcpyDeviceToDevice
        hipHostRegisterPortable = hip_hostRegisterPortable()
        gpuHostRegisterPortable = hipHostRegisterPortable
        hipHostRegisterMapped   = hip_hostRegisterMapped()
        gpuHostRegisterMapped   = hipHostRegisterMapped
        hipHostRegisterDefault  = hip_hostRegisterDefault()
        gpuHostRegisterDefault  = hipHostRegisterDefault
      endif

    end subroutine

    function gpu_devicesynchronize() result(success)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions
      implicit none
      logical                              :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_devicesynchronize()
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_devicesynchronize()
      endif
    end function


    function gpu_malloc_host(array, elements) result(success)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions
      implicit none
      type(c_ptr)                          :: array
      integer(kind=c_intptr_t), intent(in) :: elements
      logical                              :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_malloc_host(array, elements)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_host_malloc(array, elements)
      endif

    end function

    function gpu_malloc(array, elements) result(success)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions
      implicit none
      integer(kind=C_intptr_T)             :: array
      integer(kind=c_intptr_t), intent(in) :: elements
      logical                              :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_malloc(array, elements)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_malloc(array, elements)
      endif

    end function

    function gpu_host_register(array, elements, flag) result(success)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions
      implicit none
      integer(kind=C_intptr_t)              :: array
      integer(kind=c_intptr_t), intent(in)  :: elements
      integer(kind=C_INT), intent(in)       :: flag
      logical :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_host_register(array, elements, flag)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_host_register(array, elements, flag)
      endif

    end function
    
    function gpu_memcpy_intptr(dst, src, size, dir) result(success)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions
      implicit none
      integer(kind=C_intptr_t)              :: dst
      integer(kind=C_intptr_t)              :: src
      integer(kind=c_intptr_t), intent(in)  :: size
      integer(kind=C_INT), intent(in)       :: dir
      logical :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_memcpy_intptr(dst, src, size, dir)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_memcpy_intptr(dst, src, size, dir)
      endif
    
    end function
    
    function gpu_memcpy_cptr(dst, src, size, dir) result(success)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions
      implicit none
      type(c_ptr)                           :: dst
      type(c_ptr)                           :: src
      integer(kind=c_intptr_t), intent(in)  :: size
      integer(kind=C_INT), intent(in)       :: dir
      logical :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_memcpy_cptr(dst, src, size, dir)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_memcpy_cptr(dst, src, size, dir)
      endif
    
    end function
    
    function gpu_memcpy_mixed(dst, src, size, dir) result(success)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions
      implicit none
      type(c_ptr)                           :: dst
      integer(kind=C_intptr_t)              :: src
      integer(kind=c_intptr_t), intent(in)  :: size
      integer(kind=C_INT), intent(in)       :: dir
      logical :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_memcpy_mixed(dst, src, size, dir)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_memcpy_mixed(dst, src, size, dir)
      endif
    
    end function


    function gpu_memset(a, val, size) result(success)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions
      implicit none
      integer(kind=c_intptr_t)                :: a
      integer(kind=ik)                        :: val
      integer(kind=c_intptr_t), intent(in)      :: size
      integer(kind=C_INT)                     :: istat

      logical :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_memset(a, val, size)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_memset(a, val, size)
      endif

    end function

    function gpu_free(a) result(success)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions
      implicit none
      integer(kind=c_intptr_t)                :: a

      logical :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_free(a)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_free(a)
      endif

    end function

    function gpu_free_host(a) result(success)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions
      implicit none
      type(c_ptr), value          :: a

      logical :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_free_host(a)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_host_free(a)
      endif

    end function

    function gpu_host_unregister(a) result(success)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions
      implicit none
      integer(kind=c_intptr_t)                :: a

      logical :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_host_unregister(a)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_host_unregister(a)
      endif

    end function


    function gpu_memcpy2d_intptr(dst, dpitch, src, spitch, width, height , dir) result(success)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none

      integer(kind=C_intptr_T)           :: dst
      integer(kind=c_intptr_t), intent(in) :: dpitch
      integer(kind=C_intptr_T)           :: src
      integer(kind=c_intptr_t), intent(in) :: spitch
      integer(kind=c_intptr_t), intent(in) :: width
      integer(kind=c_intptr_t), intent(in) :: height
      integer(kind=C_INT), intent(in)    :: dir
      logical                            :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_memcpy2d_intptr(dst, dpitch, src, spitch, width, height , dir)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_memcpy2d_intptr(dst, dpitch, src, spitch, width, height , dir)
      endif
    end function

    function gpu_memcpy2d_cptr(dst, dpitch, src, spitch, width, height , dir) result(success)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none

      type(c_ptr)         :: dst
      integer(kind=c_intptr_t), intent(in) :: dpitch
      type(c_ptr)           :: src
      integer(kind=c_intptr_t), intent(in) :: spitch
      integer(kind=c_intptr_t), intent(in) :: width
      integer(kind=c_intptr_t), intent(in) :: height
      integer(kind=C_INT), intent(in)    :: dir
      logical                            :: success

      if (use_gpu_vendor == nvidia_gpu) then
        success = cuda_memcpy2d_cptr(dst, dpitch, src, spitch, width, height , dir)
      endif

      if (use_gpu_vendor == amd_gpu) then
        success = hip_memcpy2d_cptr(dst, dpitch, src, spitch, width, height , dir)
      endif
    end function

    subroutine gpusolver_dtrtri(uplo, diag, n, a, lda, info)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: uplo, diag
      integer(kind=C_INT64_T)         :: n, lda
      integer(kind=c_intptr_t)        :: a
      integer(kind=c_int)             :: info

      if (use_gpu_vendor == nvidia_gpu) then
        call cusolver_dtrtri(uplo, diag, n, a, lda, info)
      endif
      if (use_gpu_vendor == amd_gpu) then
        !call hipsolver_dtrtri(uplo, diag, n, a, lda, info)
      endif
    end subroutine

    subroutine gpusolver_strtri(uplo, diag, n, a, lda, info)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: uplo, diag
      integer(kind=C_INT64_T)         :: n, lda
      integer(kind=c_intptr_t)        :: a
      integer(kind=c_int)             :: info

      if (use_gpu_vendor == nvidia_gpu) then
        call cusolver_strtri(uplo, diag, n, a, lda, info)
      endif
      if (use_gpu_vendor == amd_gpu) then
        !call hipsolver_strtri(uplo, diag, n, a, lda, info)
      endif
    end subroutine

    subroutine gpusolver_ztrtri(uplo, diag, n, a, lda, info)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: uplo, diag
      integer(kind=C_INT64_T)         :: n, lda
      integer(kind=c_intptr_t)        :: a
      integer(kind=c_int)             :: info

      if (use_gpu_vendor == nvidia_gpu) then
        call cusolver_ztrtri(uplo, diag, n, a, lda, info)
      endif
      if (use_gpu_vendor == amd_gpu) then
        !call hipsolver_ztrtri(uplo, diag, n, a, lda, info)
      endif
    end subroutine

    subroutine gpusolver_ctrtri(uplo, diag, n, a, lda, info)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: uplo, diag
      integer(kind=C_INT64_T)         :: n, lda
      integer(kind=c_intptr_t)        :: a
      integer(kind=c_int)             :: info

      if (use_gpu_vendor == nvidia_gpu) then
        call cusolver_ctrtri(uplo, diag, n, a, lda, info)
      endif
      if (use_gpu_vendor == amd_gpu) then
        !call hipsolver_ctrtri(uplo, diag, n, a, lda, info)
      endif
    end subroutine

    subroutine gpusolver_dpotrf(uplo, n, a, lda, info)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: uplo
      integer(kind=C_INT)             :: n, lda
      integer(kind=c_intptr_t)        :: a
      integer(kind=c_int)             :: info

      if (use_gpu_vendor == nvidia_gpu) then
        call cusolver_dpotrf(uplo, n, a, lda, info)
      endif
      if (use_gpu_vendor == amd_gpu) then
        !call cusolver_dpotrf(uplo, n, a, lda, info)
      endif
    end subroutine

    subroutine gpusolver_spotrf(uplo, n, a, lda, info)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: uplo
      integer(kind=C_INT)             :: n, lda
      integer(kind=c_intptr_t)        :: a
      integer(kind=c_int)             :: info

      if (use_gpu_vendor == nvidia_gpu) then
        call cusolver_spotrf(uplo, n, a, lda, info)
      endif
      if (use_gpu_vendor == amd_gpu) then
        !call cusolver_spotrf(uplo, n, a, lda, info)
      endif
    end subroutine

    subroutine gpusolver_zpotrf(uplo, n, a, lda, info)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: uplo
      integer(kind=C_INT)             :: n, lda
      integer(kind=c_intptr_t)        :: a
      integer(kind=c_int)             :: info

      if (use_gpu_vendor == nvidia_gpu) then
        call cusolver_zpotrf(uplo, n, a, lda, info)
      endif
      if (use_gpu_vendor == amd_gpu) then
        !call cusolver_zpotrf(uplo, n, a, lda, info)
      endif
    end subroutine

    subroutine gpusolver_cpotrf(uplo, n, a, lda, info)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: uplo
      integer(kind=C_INT)             :: n, lda
      integer(kind=c_intptr_t)        :: a
      integer(kind=c_int)             :: info

      if (use_gpu_vendor == nvidia_gpu) then
        call cusolver_cpotrf(uplo, n, a, lda, info)
      endif
      if (use_gpu_vendor == amd_gpu) then
        !call cusolver_cpotrf(uplo, n, a, lda, info)
      endif
    end subroutine

    subroutine gpublas_dgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: cta
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,incx,incy
      real(kind=C_DOUBLE)             :: alpha,beta
      integer(kind=C_intptr_T)        :: a, x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_dgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_dgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      endif
    end subroutine

    subroutine gpublas_sgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: cta
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,incx,incy
      real(kind=C_FLOAT)              :: alpha,beta
      integer(kind=C_intptr_T)        :: a, x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_sgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_sgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      endif

    end subroutine

    subroutine gpublas_zgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: cta
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,incx,incy
      complex(kind=C_DOUBLE_COMPLEX)             :: alpha,beta
      integer(kind=C_intptr_T)        :: a, x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_zgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_zgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      endif

    end subroutine

    subroutine gpublas_cgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: cta
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,incx,incy
      complex(kind=C_FLOAT_COMPLEX)              :: alpha,beta
      integer(kind=C_intptr_T)        :: a, x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_cgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_cgemv(cta, m, n, alpha, a, lda, x, incx, beta, y, incy)
      endif

    end subroutine

    subroutine gpublas_dgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: cta, ctb
      integer(kind=C_INT)             :: m,n,k
      integer(kind=C_INT), intent(in) :: lda,ldb,ldc
      real(kind=C_DOUBLE)             :: alpha,beta
      integer(kind=C_intptr_T)        :: a, b, c

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_dgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_dgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)
      endif

    end subroutine 


    subroutine gpublas_sgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: cta, ctb
      integer(kind=C_INT)             :: m,n,k
      integer(kind=C_INT), intent(in) :: lda,ldb,ldc
      real(kind=C_FLOAT)              :: alpha,beta
      integer(kind=C_intptr_T)        :: a, b, c

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_sgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_sgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)
      endif


    end subroutine

    subroutine gpublas_zgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: cta, ctb
      integer(kind=C_INT)             :: m,n,k
      integer(kind=C_INT), intent(in) :: lda,ldb,ldc
      complex(kind=C_DOUBLE_COMPLEX)          :: alpha,beta
      integer(kind=C_intptr_T)        :: a, b, c

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_zgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_zgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)
      endif

    end subroutine

    subroutine gpublas_cgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: cta, ctb
      integer(kind=C_INT)             :: m,n,k
      integer(kind=C_INT), intent(in) :: lda,ldb,ldc
      complex(kind=C_FLOAT_COMPLEX)           :: alpha,beta
      integer(kind=C_intptr_T)        :: a, b, c

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_cgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_cgemm(cta, ctb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)
      endif
    end subroutine

    subroutine gpublas_dcopy_intptr(n, x, incx, y, incy)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      integer(kind=C_INT)             :: n
      integer(kind=C_INT), intent(in) :: incx, incy
      integer(kind=C_intptr_T)        :: x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_dcopy_intptr(n, x, incx, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_dcopy_intptr(n, x, incx, y, incy)
      endif

    end subroutine

    subroutine gpublas_dcopy_cptr(n, x, incx, y, incy)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      integer(kind=C_INT)             :: n
      integer(kind=C_INT), intent(in) :: incx, incy
      type(c_ptr)        :: x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_dcopy_cptr(n, x, incx, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_dcopy_cptr(n, x, incx, y, incy)
      endif

    end subroutine

    subroutine gpublas_scopy_intptr(n, x, incx, y, incy)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      integer(kind=C_INT)             :: n
      integer(kind=C_INT), intent(in) :: incx, incy
      integer(kind=C_intptr_T)        :: x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_scopy_intptr(n, x, incx, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_dcopy_intptr(n, x, incx, y, incy)
      endif

    end subroutine

    subroutine gpublas_scopy_cptr(n, x, incx, y, incy)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      integer(kind=C_INT)             :: n
      integer(kind=C_INT), intent(in) :: incx, incy
      type(c_ptr)        :: x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_scopy_cptr(n, x, incx, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_scopy_cptr(n, x, incx, y, incy)
      endif

    end subroutine


    subroutine gpublas_zcopy_intptr(n, x, incx, y, incy)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      integer(kind=C_INT)             :: n
      integer(kind=C_INT), intent(in) :: incx, incy
      integer(kind=C_intptr_T)        :: x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_zcopy_intptr(n, x, incx, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_zcopy_intptr(n, x, incx, y, incy)
      endif

    end subroutine

    subroutine gpublas_zcopy_cptr(n, x, incx, y, incy)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      integer(kind=C_INT)             :: n
      integer(kind=C_INT), intent(in) :: incx, incy
      type(c_ptr)        :: x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_zcopy_cptr(n, x, incx, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_dcopy_cptr(n, x, incx, y, incy)
      endif

    end subroutine

    subroutine gpublas_ccopy_intptr(n, x, incx, y, incy)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      integer(kind=C_INT)             :: n
      integer(kind=C_INT), intent(in) :: incx, incy
      integer(kind=C_intptr_T)        :: x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_ccopy_intptr(n, x, incx, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_ccopy_intptr(n, x, incx, y, incy)
      endif

    end subroutine

    subroutine gpublas_ccopy_cptr(n, x, incx, y, incy)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      integer(kind=C_INT)             :: n
      integer(kind=C_INT), intent(in) :: incx, incy
      type(c_ptr)        :: x, y

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_ccopy_cptr(n, x, incx, y, incy)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_ccopy_cptr(n, x, incx, y, incy)
      endif

    end subroutine


    subroutine gpublas_dtrmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      real(kind=C_DOUBLE)             :: alpha
      integer(kind=C_intptr_T)        :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_dtrmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_dtrmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

    end subroutine


    subroutine gpublas_dtrmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      real(kind=C_DOUBLE)             :: alpha
      type(c_ptr)                     :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_dtrmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_dtrmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

    end subroutine


    subroutine gpublas_strmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      real(kind=C_FLOAT)              :: alpha
      integer(kind=C_intptr_T)        :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_strmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_strmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

    end subroutine


    subroutine gpublas_strmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      real(kind=C_FLOAT)              :: alpha
      type(c_ptr)                     :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_strmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_strmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

    end subroutine


    subroutine gpublas_ztrmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      complex(kind=C_DOUBLE_COMPLEX)          :: alpha
      integer(kind=C_intptr_T)        :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_ztrmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_ztrmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif
    end subroutine



    subroutine gpublas_ztrmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      complex(kind=C_DOUBLE_COMPLEX)          :: alpha
      type(c_ptr)                     :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_ztrmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_ztrmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif
    end subroutine


    subroutine gpublas_ctrmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      complex(kind=C_FLOAT_COMPLEX)   :: alpha
      integer(kind=C_intptr_T)        :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_ctrmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_ctrmm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif
    end subroutine



    subroutine gpublas_ctrmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      complex(kind=C_FLOAT_COMPLEX)   :: alpha
      type(c_ptr)                     :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_ctrmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_ctrmm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif
    end subroutine


    subroutine gpublas_dtrsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      real(kind=C_DOUBLE)             :: alpha
      integer(kind=C_intptr_T)        :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_dtrsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_dtrsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

    end subroutine


    subroutine gpublas_dtrsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      real(kind=C_DOUBLE)             :: alpha
      type(c_ptr)                     :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_dtrsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_dtrsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

    end subroutine


    subroutine gpublas_strsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      real(kind=C_FLOAT)              :: alpha
      integer(kind=C_intptr_T)        :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_strsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_strsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

    end subroutine


    subroutine gpublas_strsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      real(kind=C_FLOAT)              :: alpha
      type(c_ptr)                     :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_strsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_strsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

    end subroutine


    subroutine gpublas_ztrsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      complex(kind=C_DOUBLE_COMPLEX)          :: alpha
      integer(kind=C_intptr_T)        :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_ztrsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_ztrsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif
    end subroutine



    subroutine gpublas_ztrsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      complex(kind=C_DOUBLE_COMPLEX)          :: alpha
      type(c_ptr)                     :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_ztrsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_ztrsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif
    end subroutine


    subroutine gpublas_ctrsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      complex(kind=C_FLOAT_COMPLEX)   :: alpha
      integer(kind=C_intptr_T)        :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_ctrsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_ctrsm_intptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif
    end subroutine



    subroutine gpublas_ctrsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)

      use, intrinsic :: iso_c_binding
      use cuda_functions
      use hip_functions

      implicit none
      character(1,C_CHAR),value       :: side, uplo, trans, diag
      integer(kind=C_INT)             :: m,n
      integer(kind=C_INT), intent(in) :: lda,ldb
      complex(kind=C_FLOAT_COMPLEX)   :: alpha
      type(c_ptr)                     :: a, b

      if (use_gpu_vendor == nvidia_gpu) then
        call cublas_ctrsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif

      if (use_gpu_vendor == amd_gpu) then
        call rocblas_ctrsm_cptr(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
      endif
    end subroutine


end module

