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

module cuda_routines

  use precision
  implicit none

  public

  ! TODO: take these values from the definition header files of CUDA !!
  integer(kind=ik), parameter :: cudaMemcpyHostToDevice  = 1
  integer(kind=ik), parameter :: cudaMemcpyDeviceToHost  = 2
  integer(kind=ik), parameter :: cudaHostRegisterPortable  = 3
  integer(kind=ik), parameter :: cudaHostRegisterMapped  = 4
  integer(kind=ik), parameter :: cudaMemcpyDeviceToDevice = 5

  interface
    function cuda_setdevice_c(n) result(istat) &
             bind(C, name="cudaSetDeviceFromC")

      use iso_c_binding
      implicit none
      integer(C_INT), value    :: n
      integer(C_INT)           :: istat
    end function cuda_setdevice_c
  end interface

  interface
    function cuda_getdevicecount_c(n) result(istat) &
             bind(C, name="cudaGetDeviceCountFromC")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(out) :: n
      integer(C_INT)              :: istat
    end function cuda_getdevicecount_c
  end interface

   interface
     function cuda_malloc_c(a, width_height) result(istat) &
              bind(C, name="cudaMallocFromC")

       use iso_c_binding
       implicit none

       integer(C_SIZE_T)                    :: a
       integer(C_SIZE_T), intent(in), value :: width_height
       integer(C_INT)                       :: istat

     end function cuda_malloc_c
   end interface

  interface
    function cuda_free_c(a) result(istat) &
             bind(C, name="cudaFreeFromC")

      use iso_c_binding

      implicit none
      integer(C_SIZE_T), value  :: a
      integer(C_INT)            :: istat

    end function cuda_free_c
  end interface

  interface
    function cuda_memcpy_c(dst, src, size, dir) result(istat) &
             bind(C, name="cudaMemcpyFromC")

      use iso_c_binding

      implicit none
      integer(C_SIZE_T), value              :: dst
      integer(C_SIZE_T), value              :: src
      integer(C_SIZE_T), intent(in), value  :: size
      integer(C_INT), intent(in), value     :: dir
      integer(C_INT)                        :: istat

    end function cuda_memcpy_c
  end interface

  interface
    function cuda_memcpy2d_c(dst, dpitch, src, spitch, width, height , dir) result(istat) &
             bind(C, name="cudaMemcpy2dFromC")

      use iso_c_binding

      implicit none

      integer(C_SIZE_T), value              :: dst
      integer(C_SIZE_T), intent(in), value  :: dpitch
      integer(C_SIZE_T), value              :: src
      integer(C_SIZE_T), intent(in), value  :: spitch
      integer(C_SIZE_T), intent(in), value  :: width
      integer(C_SIZE_T), intent(in), value  :: height
      integer(C_INT), intent(in), value     :: dir
      integer(C_INT)                        :: istat

    end function cuda_memcpy2d_c
  end interface

  interface

 function cuda_memset_c(a, val, size) result(istat) &
          bind(C, name="cudaMemsetFromC")

   use iso_c_binding

   implicit none

   integer(C_SIZE_T), value              :: a
   !integer(C_INT)                       :: val
   integer(C_INT)                        :: val
   integer(C_SIZE_T), intent(in), value  :: size
   integer(C_INT)                        :: istat

 end function cuda_memset_c
 end interface

contains
  function cuda_setdevice(n) result(success)
    use iso_c_binding
    use precision
    implicit none

    integer(kind=ik), intent(in)  :: n
    logical                       :: success

#ifdef WITH_GPU_VERSION
    success = cuda_setdevice_c(int(n,kind=c_int)) /= 0
#else
    success = .true.
#endif
  end function cuda_setdevice

  function cuda_getdevicecount(n) result(success)
    use iso_c_binding
    use precision
    implicit none

    integer(kind=ik), intent(out) :: n
    logical                       :: success

#ifdef WITH_GPU_VERSION
    success = cuda_getdevicecount_c(n) /=0
#else
    success = .true.
    n     = 0
#endif
  end function cuda_getdevicecount

  function cuda_malloc(a, width_height) result(success)
    use iso_c_binding
    implicit none

    integer(C_SIZE_T)             :: a
    integer(C_SIZE_T), intent(in) :: width_height
    logical                       :: success

#ifdef WTIH_GPU_VERSION
    success = cuda_malloc_c(a,width_height) /= 0
#else
    success = .false.
#endif
  end function

  function cuda_memcpy(dst, src, size, dir) result(success)

    use iso_c_binding

    implicit none
    integer(C_SIZE_T)             :: dst
    integer(C_SIZE_T)             :: src
    integer(C_SIZE_T), intent(in) :: size
    integer(C_INT), intent(in)    :: dir
    logical                       :: success

#ifdef WITH_GPU_VERSION
    success = cuda_memcpy_c(dst, src, size, dir) /=0
#else
    success = .false.
#endif

  end function cuda_memcpy

  function cuda_free(a) result(success)

    use iso_c_binding

    implicit none
    integer(C_SIZE_T) :: a
    logical                       :: success

#ifdef WITH_GPU_VERSION
    success = cuda_free_c(a) /= 0
#else
    success = .false.
#endif

  end function cuda_free

  function cuda_memcpy2d(dst, dpitch, src, spitch, width, height , dir) result(success)

    use iso_c_binding

    implicit none

    integer(C_SIZE_T)             :: dst
    integer(C_SIZE_T), intent(in) :: dpitch
    integer(C_SIZE_T)             :: src
    integer(C_SIZE_T), intent(in) :: spitch
    integer(C_SIZE_T), intent(in) :: width
    integer(C_SIZE_T), intent(in) :: height
    integer(C_INT), intent(in)    :: dir
    integer(C_INT)                :: istat
    logical                       :: success

#ifdef WITH_GPU_VERSION
    success = cuda_memcpy2d_c(dst, dpitch, src, spitch, width, height , dir) /= 0
#else
    success = .false.
#endif

  end function cuda_memcpy2d

end module cuda_routines
