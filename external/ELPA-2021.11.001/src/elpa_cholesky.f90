










!    Copyright 2021, A. Marek
!
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

module elpa_cholesky
  use, intrinsic :: iso_c_binding
  use precision
  implicit none

  public

  public :: elpa_cholesky_real_double_impl       !< Cholesky factorization of a double-precision real matrix

  public :: elpa_cholesky_complex_double_impl    !< Cholesky factorization of a double-precision complex matrix

  public :: elpa_cholesky_real_single_impl       !< Cholesky factorization of a single-precision real matrix


  public :: elpa_cholesky_complex_single_impl    !< Cholesky factorization of a single-precision complex matrix

  contains

















!> \brief  elpa_cholesky_real_double_impl: Cholesky factorization of a double-precision real symmetric matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
   function elpa_cholesky_real_double_impl (obj, a) result(success)
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
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
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




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)

  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_cuda, only : copy_double_tmp1_tmp2, &
                              copy_double_a_tmp1
  use cholesky_cuda
  implicit none
!    Copyright 2011, A. Marek
!
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
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
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
!
  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: rck = C_DOUBLE
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)              :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  real(kind=rck)      :: a(obj%local_nrows,*)
  integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)              :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)              :: n, nc, i, info
  integer(kind=BLAS_KIND)       :: infoBLAS
  integer(kind=ik)              :: lcs, lce, lrs, lre
  integer(kind=ik)              :: tile_size, l_rows_tile, l_cols_tile

  real(kind=rck), allocatable    :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                       :: wantDebug
  logical                       :: success
  integer(kind=ik)              :: istat, debug, error
  character(200)                :: errorMessage
  integer(kind=ik)              :: nrThreads, limitThreads
  character(20)                 :: gpuString
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_intptr_t)      :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  type(c_ptr)                      :: tmp1_mpi_dev
  real(kind=rck), pointer :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)      :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                      :: tmatc_mpi_dev
  real(kind=rck), pointer :: tmatc_mpi_fortran_ptr(:,:)
  integer(kind=c_int)              :: gpu_cholesky

  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &double&
                                                            &_&
                                                            &real

  gpu_cholesky = 0
  ! GPU settings
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"ELPA_CHOLESKY: Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif

    call obj%get("gpu_cholesky",gpu_cholesky, error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for gpu_cholesky. Aborting..."
      stop
    endif

  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
  endif

  if (gpu_cholesky .eq. 1) then
    useGPU = (gpu == 1)
  else
    useGPU = .false.
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("elpa_cholesky_&
  &real&
  &_&
  &double&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("mpi_comm_parent", mpi_comm_all, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Error getting option for mpi_comm_all. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_rows",mpi_comm_rows,error )
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for mpi_comm_rows. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_cols",mpi_comm_cols,error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for mpi_comm_cols. Aborting..."
    stop
  endif

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    stop
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), myidMPI, mpierr)


  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  myid    = int(myidMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")
  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then
    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 257,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 260,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 263,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 266,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 269,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 272,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 275,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 278,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: a_dev", 281,  successGPU)

  endif
  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 285,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 288,  istat,  errorMessage)

  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 294,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 297,  istat,  errorMessage)

  tmatr = 0
  tmatc = 0

  if (useGPU) then
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_cholesky 1: memcpy a-> a_dev", 305,  successGPU)
  endif

  do n = 1, na, nblk
    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 340,  successGPU)

          call DPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 347,  successGPU)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false.
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))
      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          call DPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false.
            return
          endif

        endif
      endif ! useGPU
      exit ! Loop
    endif ! (n+nblk > na) 

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 412,  successGPU)

          call DPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 419,  successGPU)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false.
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          call DPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf 2: ",info
            success = .false.
            return
          endif ! useGPU
        endif
 
        if (useGPU) then
          call copy_double_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk)
        else ! useGPU
          nc = 0
          do i=1,nblk
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), tmp1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_cholesky: tmp1_dev to tmp1", 479,  successGPU)

      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      call obj%timer%start("mpi_communication")

      call MPI_Bcast(tmp1, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
                    MPI_REAL8,         &
                    int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")
!#else
!      tmp1_mpi_dev = transfer(tmp1_dev, tmp1_mpi_dev)
!      ! and associate a fortran pointer
!      call c_f_pointer(tmp1_mpi_dev, tmp1_mpi_fortran_ptr, [nblk,nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("cholesky: device_synchronize", 503,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!
!      call MPI_Bcast(tmp1_mpi_fortran_ptr, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
!#if 1 == 1
!                    MPI_REAL8,         &
!#endif
!#if COMPLEXCASE == 1
!                    MPI_COMPLEX_PRECISION,      &
!#endif
!                    int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!
!      call obj%timer%stop("mpi_cuda_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmp1_dev, int(loc(tmp1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_cholesky: tmp1 to tmp1_dev", 524,  successGPU)

      endif
!#endif


      if (useGPU) then
        call copy_double_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_DTRSM('L', 'U', 'T', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows)
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("blas")
        if (l_cols-l_colx+1>0) &
        call DTRSM('L', 'U', 'T', 'N', int(nblk,kind=BLAS_KIND),  &
                            int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
                            int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
        call obj%timer%stop("blas")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        if (l_cols-l_colx+1>0) &
           call copy_double_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)
      endif
    else ! useGPU
      do i=1,nblk
        if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = a(l_row1+i-1,l_colx:l_cols)
      enddo
    endif ! useGPU

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      if (l_cols-l_colx+1 > 0) then
        num = l_cols*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 586,  successGPU)
      endif
    endif
!#endif


!#ifndef WITH_CUDA_AWARE_MPI
    do i=1,nblk
      call obj%timer%start("mpi_communication")
      if (l_cols-l_colx+1>0) &
      call MPI_Bcast(tmatc(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), MPI_REAL8, &
                     int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")
    enddo
!#else
!    tmatc_mpi_dev = transfer(tmatc_dev, tmatc_mpi_dev)
!    ! and associate a fortran pointer
!    call c_f_pointer(tmatc_mpi_dev, tmatc_mpi_fortran_ptr, [l_cols,nblk])
!    
!    if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!    successGPU = gpu_devicesynchronize()
!    call check_memcpy_GPU_f("cholesky: device_synchronize", 610,  successGPU)
!    if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!
!    do i=1,nblk
!      call obj%timer%start("mpi_cuda_communication")
!      if (l_cols-l_colx+1>0) &
!      call MPI_Bcast(tmatc_mpi_fortran_ptr(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), &
!                     MPI_REAL8, &
!                     int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)
!
!      call obj%timer%stop("mpi_cuda_communication")
!    enddo
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      !if (l_cols-l_colx+1 > 0) then
        num = l_cols*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmatc_dev, int(loc(tmatc),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_cholesky: tmatc to tmatc_dev", 632,  successGPU)
      !endif
    endif
!#endif

    if (useGPU) then
      ! can optimize memcpy here with previous ones

      ! a gpu version of elpa_transpose_vectors is needed

!#if !defined(1) || (defined(1) && defined(WITH_CUDA_AWARE_MPI))

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 657,  successGPU)
    endif

    call elpa_transpose_vectors_&
    &real&
    &_&
    &double &
    (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
    tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
    n, na, nblk, nblk, nrThreads, .false.)

    if (useGPU) then
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 672,  successGPU)
    endif


    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_DGEMM('N', 'T', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows)
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call DGEMM('N', 'T', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                            int(nblk,kind=BLAS_KIND), -ONE,  &
                            tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
                            int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
                            ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU

  enddo ! n = 1, na, nblk

  if (useGPU) then
    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 713,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 716,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 719,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 722,  successGPU)
  endif

  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 726,  istat,  errorMessage)


  if (useGPU) then
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 732,  successGPU)
  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  !if (useGPU) then
  !else ! useGPU
    do i=1,na
      if (my_pcol==pcol(i, nblk, np_cols)) then
        ! column i is on local processor
        l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
        l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
        a(l_row1:l_rows,l_col1) = 0
      endif
    enddo
  !endif ! useGPU

  if (useGPU) then
    ! copy back
    !successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
    !                   matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    !call check_memcpy_GPU_f("elpa_cholesky: memcpy a-> d_dev", 752,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_cholesky: a_dev", 755,  successGPU)
  endif

  ! restore original OpenMP settings
  call obj%timer%stop("elpa_cholesky_&
  &real&
  &_&
  &double&
  &"//gpuString)

    end function elpa_cholesky_real_double_impl




















!> \brief  elpa_cholesky_real_single_impl: Cholesky factorization of a double-precision real symmetric matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
   function elpa_cholesky_real_single_impl(obj, a) result(success)
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
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
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




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)

  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_cuda, only : copy_float_tmp1_tmp2, &
                              copy_float_a_tmp1
  use cholesky_cuda
  implicit none
!    Copyright 2011, A. Marek
!
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
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
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
!
  integer, parameter :: rk = C_FLOAT
  integer, parameter :: rck = C_FLOAT
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)              :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  real(kind=rck)      :: a(obj%local_nrows,*)
  integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)              :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)              :: n, nc, i, info
  integer(kind=BLAS_KIND)       :: infoBLAS
  integer(kind=ik)              :: lcs, lce, lrs, lre
  integer(kind=ik)              :: tile_size, l_rows_tile, l_cols_tile

  real(kind=rck), allocatable    :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                       :: wantDebug
  logical                       :: success
  integer(kind=ik)              :: istat, debug, error
  character(200)                :: errorMessage
  integer(kind=ik)              :: nrThreads, limitThreads
  character(20)                 :: gpuString
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_intptr_t)      :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  type(c_ptr)                      :: tmp1_mpi_dev
  real(kind=rck), pointer :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)      :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                      :: tmatc_mpi_dev
  real(kind=rck), pointer :: tmatc_mpi_fortran_ptr(:,:)
  integer(kind=c_int)              :: gpu_cholesky

  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &single&
                                                            &_&
                                                            &real

  gpu_cholesky = 0
  ! GPU settings
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"ELPA_CHOLESKY: Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif

    call obj%get("gpu_cholesky",gpu_cholesky, error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for gpu_cholesky. Aborting..."
      stop
    endif

  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
  endif

  if (gpu_cholesky .eq. 1) then
    useGPU = (gpu == 1)
  else
    useGPU = .false.
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("elpa_cholesky_&
  &real&
  &_&
  &single&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("mpi_comm_parent", mpi_comm_all, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Error getting option for mpi_comm_all. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_rows",mpi_comm_rows,error )
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for mpi_comm_rows. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_cols",mpi_comm_cols,error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for mpi_comm_cols. Aborting..."
    stop
  endif

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    stop
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), myidMPI, mpierr)


  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  myid    = int(myidMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")
  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then
    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 257,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 260,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 263,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 266,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 269,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 272,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 275,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 278,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: a_dev", 281,  successGPU)

  endif
  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 285,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 288,  istat,  errorMessage)

  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 294,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 297,  istat,  errorMessage)

  tmatr = 0
  tmatc = 0

  if (useGPU) then
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_cholesky 1: memcpy a-> a_dev", 305,  successGPU)
  endif

  do n = 1, na, nblk
    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 340,  successGPU)

          call SPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 347,  successGPU)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false.
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))
      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          call SPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false.
            return
          endif

        endif
      endif ! useGPU
      exit ! Loop
    endif ! (n+nblk > na) 

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 412,  successGPU)

          call SPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 419,  successGPU)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf: ",info
            success = .false.
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          call SPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &real&

            &: Error in dpotrf 2: ",info
            success = .false.
            return
          endif ! useGPU
        endif
 
        if (useGPU) then
          call copy_float_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk)
        else ! useGPU
          nc = 0
          do i=1,nblk
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), tmp1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_cholesky: tmp1_dev to tmp1", 479,  successGPU)

      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      call obj%timer%start("mpi_communication")

      call MPI_Bcast(tmp1, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
                    MPI_REAL4,         &
                    int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")
!#else
!      tmp1_mpi_dev = transfer(tmp1_dev, tmp1_mpi_dev)
!      ! and associate a fortran pointer
!      call c_f_pointer(tmp1_mpi_dev, tmp1_mpi_fortran_ptr, [nblk,nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("cholesky: device_synchronize", 503,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!
!      call MPI_Bcast(tmp1_mpi_fortran_ptr, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
!#if 1 == 1
!                    MPI_REAL4,         &
!#endif
!#if COMPLEXCASE == 1
!                    MPI_COMPLEX_PRECISION,      &
!#endif
!                    int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!
!      call obj%timer%stop("mpi_cuda_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmp1_dev, int(loc(tmp1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_cholesky: tmp1 to tmp1_dev", 524,  successGPU)

      endif
!#endif


      if (useGPU) then
        call copy_float_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_STRSM('L', 'U', 'T', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows)
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("blas")
        if (l_cols-l_colx+1>0) &
        call STRSM('L', 'U', 'T', 'N', int(nblk,kind=BLAS_KIND),  &
                            int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
                            int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
        call obj%timer%stop("blas")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        if (l_cols-l_colx+1>0) &
           call copy_float_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)
      endif
    else ! useGPU
      do i=1,nblk
        if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = a(l_row1+i-1,l_colx:l_cols)
      enddo
    endif ! useGPU

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      if (l_cols-l_colx+1 > 0) then
        num = l_cols*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 586,  successGPU)
      endif
    endif
!#endif


!#ifndef WITH_CUDA_AWARE_MPI
    do i=1,nblk
      call obj%timer%start("mpi_communication")
      if (l_cols-l_colx+1>0) &
      call MPI_Bcast(tmatc(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), MPI_REAL4, &
                     int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")
    enddo
!#else
!    tmatc_mpi_dev = transfer(tmatc_dev, tmatc_mpi_dev)
!    ! and associate a fortran pointer
!    call c_f_pointer(tmatc_mpi_dev, tmatc_mpi_fortran_ptr, [l_cols,nblk])
!    
!    if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!    successGPU = gpu_devicesynchronize()
!    call check_memcpy_GPU_f("cholesky: device_synchronize", 610,  successGPU)
!    if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!
!    do i=1,nblk
!      call obj%timer%start("mpi_cuda_communication")
!      if (l_cols-l_colx+1>0) &
!      call MPI_Bcast(tmatc_mpi_fortran_ptr(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), &
!                     MPI_REAL4, &
!                     int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)
!
!      call obj%timer%stop("mpi_cuda_communication")
!    enddo
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      !if (l_cols-l_colx+1 > 0) then
        num = l_cols*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmatc_dev, int(loc(tmatc),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_cholesky: tmatc to tmatc_dev", 632,  successGPU)
      !endif
    endif
!#endif

    if (useGPU) then
      ! can optimize memcpy here with previous ones

      ! a gpu version of elpa_transpose_vectors is needed

!#if !defined(1) || (defined(1) && defined(WITH_CUDA_AWARE_MPI))

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 657,  successGPU)
    endif

    call elpa_transpose_vectors_&
    &real&
    &_&
    &single &
    (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
    tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
    n, na, nblk, nblk, nrThreads, .false.)

    if (useGPU) then
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 672,  successGPU)
    endif


    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_SGEMM('N', 'T', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows)
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call SGEMM('N', 'T', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                            int(nblk,kind=BLAS_KIND), -ONE,  &
                            tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
                            int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
                            ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU

  enddo ! n = 1, na, nblk

  if (useGPU) then
    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 713,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 716,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 719,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 722,  successGPU)
  endif

  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 726,  istat,  errorMessage)


  if (useGPU) then
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 732,  successGPU)
  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  !if (useGPU) then
  !else ! useGPU
    do i=1,na
      if (my_pcol==pcol(i, nblk, np_cols)) then
        ! column i is on local processor
        l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
        l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
        a(l_row1:l_rows,l_col1) = 0
      endif
    enddo
  !endif ! useGPU

  if (useGPU) then
    ! copy back
    !successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
    !                   matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    !call check_memcpy_GPU_f("elpa_cholesky: memcpy a-> d_dev", 752,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_cholesky: a_dev", 755,  successGPU)
  endif

  ! restore original OpenMP settings
  call obj%timer%stop("elpa_cholesky_&
  &real&
  &_&
  &single&
  &"//gpuString)

    end function elpa_cholesky_real_single_impl




















!> \brief  elpa_cholesky_complex_double_impl: Cholesky factorization of a double-precision complex hermitian matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
    function elpa_cholesky_complex_double_impl(obj, a) result(success)

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
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
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




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)

  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_cuda, only : copy_double_complex_tmp1_tmp2, &
                              copy_double_complex_a_tmp1
  use cholesky_cuda
  implicit none
!    Copyright 2011, A. Marek
!
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
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
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
!

  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: ck = C_DOUBLE_COMPLEX
  integer, parameter :: rck = C_DOUBLE_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)              :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  complex(kind=rck)      :: a(obj%local_nrows,*)
  integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)              :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)              :: n, nc, i, info
  integer(kind=BLAS_KIND)       :: infoBLAS
  integer(kind=ik)              :: lcs, lce, lrs, lre
  integer(kind=ik)              :: tile_size, l_rows_tile, l_cols_tile

  complex(kind=rck), allocatable    :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                       :: wantDebug
  logical                       :: success
  integer(kind=ik)              :: istat, debug, error
  character(200)                :: errorMessage
  integer(kind=ik)              :: nrThreads, limitThreads
  character(20)                 :: gpuString
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_intptr_t)      :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  type(c_ptr)                      :: tmp1_mpi_dev
  complex(kind=rck), pointer :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)      :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                      :: tmatc_mpi_dev
  complex(kind=rck), pointer :: tmatc_mpi_fortran_ptr(:,:)
  integer(kind=c_int)              :: gpu_cholesky

  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &double&
                                                            &_&
                                                            &complex

  gpu_cholesky = 0
  ! GPU settings
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"ELPA_CHOLESKY: Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif

    call obj%get("gpu_cholesky",gpu_cholesky, error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for gpu_cholesky. Aborting..."
      stop
    endif

  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
  endif

  if (gpu_cholesky .eq. 1) then
    useGPU = (gpu == 1)
  else
    useGPU = .false.
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("elpa_cholesky_&
  &complex&
  &_&
  &double&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("mpi_comm_parent", mpi_comm_all, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Error getting option for mpi_comm_all. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_rows",mpi_comm_rows,error )
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for mpi_comm_rows. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_cols",mpi_comm_cols,error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for mpi_comm_cols. Aborting..."
    stop
  endif

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    stop
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), myidMPI, mpierr)


  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  myid    = int(myidMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")
  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then
    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 257,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 260,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 263,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 266,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 269,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 272,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 275,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 278,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: a_dev", 281,  successGPU)

  endif
  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 285,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 288,  istat,  errorMessage)

  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 294,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 297,  istat,  errorMessage)

  tmatr = 0
  tmatc = 0

  if (useGPU) then
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_cholesky 1: memcpy a-> a_dev", 305,  successGPU)
  endif

  do n = 1, na, nblk
    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 340,  successGPU)

          call ZPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 347,  successGPU)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false.
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))
      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          call ZPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false.
            return
          endif

        endif
      endif ! useGPU
      exit ! Loop
    endif ! (n+nblk > na) 

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 412,  successGPU)

          call ZPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 419,  successGPU)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false.
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          call ZPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf 2: ",info

            success = .false.
            return
          endif ! useGPU
        endif
 
        if (useGPU) then
          call copy_double_complex_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk)
        else ! useGPU
          nc = 0
          do i=1,nblk
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), tmp1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_cholesky: tmp1_dev to tmp1", 479,  successGPU)

      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      call obj%timer%start("mpi_communication")

      call MPI_Bcast(tmp1, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
                    MPI_DOUBLE_COMPLEX,      &
                    int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")
!#else
!      tmp1_mpi_dev = transfer(tmp1_dev, tmp1_mpi_dev)
!      ! and associate a fortran pointer
!      call c_f_pointer(tmp1_mpi_dev, tmp1_mpi_fortran_ptr, [nblk,nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("cholesky: device_synchronize", 503,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!
!      call MPI_Bcast(tmp1_mpi_fortran_ptr, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
!#if REALCASE == 1
!                    MPI_REAL8,         &
!#endif
!#if 1 == 1
!                    MPI_DOUBLE_COMPLEX,      &
!#endif
!                    int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!
!      call obj%timer%stop("mpi_cuda_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmp1_dev, int(loc(tmp1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_cholesky: tmp1 to tmp1_dev", 524,  successGPU)

      endif
!#endif


      if (useGPU) then
        call copy_double_complex_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_ZTRSM('L', 'U', 'C', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows)
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("blas")
        if (l_cols-l_colx+1>0) &
        call ZTRSM('L', 'U', 'C', 'N', int(nblk,kind=BLAS_KIND),  &
                            int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
                            int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
        call obj%timer%stop("blas")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        if (l_cols-l_colx+1>0) &
           call copy_double_complex_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)
      endif
    else ! useGPU
      do i=1,nblk
        if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = conjg(a(l_row1+i-1,l_colx:l_cols))
      enddo
    endif ! useGPU

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      if (l_cols-l_colx+1 > 0) then
        num = l_cols*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 586,  successGPU)
      endif
    endif
!#endif


!#ifndef WITH_CUDA_AWARE_MPI
    do i=1,nblk
      call obj%timer%start("mpi_communication")
      if (l_cols-l_colx+1>0) &
      call MPI_Bcast(tmatc(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                     int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")
    enddo
!#else
!    tmatc_mpi_dev = transfer(tmatc_dev, tmatc_mpi_dev)
!    ! and associate a fortran pointer
!    call c_f_pointer(tmatc_mpi_dev, tmatc_mpi_fortran_ptr, [l_cols,nblk])
!    
!    if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!    successGPU = gpu_devicesynchronize()
!    call check_memcpy_GPU_f("cholesky: device_synchronize", 610,  successGPU)
!    if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!
!    do i=1,nblk
!      call obj%timer%start("mpi_cuda_communication")
!      if (l_cols-l_colx+1>0) &
!      call MPI_Bcast(tmatc_mpi_fortran_ptr(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), &
!                     MPI_DOUBLE_COMPLEX, &
!                     int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)
!
!      call obj%timer%stop("mpi_cuda_communication")
!    enddo
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      !if (l_cols-l_colx+1 > 0) then
        num = l_cols*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmatc_dev, int(loc(tmatc),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_cholesky: tmatc to tmatc_dev", 632,  successGPU)
      !endif
    endif
!#endif

    if (useGPU) then
      ! can optimize memcpy here with previous ones

      ! a gpu version of elpa_transpose_vectors is needed

!#if !defined(1) || (defined(1) && defined(WITH_CUDA_AWARE_MPI))

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 657,  successGPU)
    endif

    call elpa_transpose_vectors_&
    &complex&
    &_&
    &double &
    (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
    tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
    n, na, nblk, nblk, nrThreads, .false.)

    if (useGPU) then
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 672,  successGPU)
    endif


    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_ZGEMM('N', 'C', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows)
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call ZGEMM('N', 'C', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                            int(nblk,kind=BLAS_KIND), -ONE,  &
                            tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
                            int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
                            ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU

  enddo ! n = 1, na, nblk

  if (useGPU) then
    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 713,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 716,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 719,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 722,  successGPU)
  endif

  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 726,  istat,  errorMessage)


  if (useGPU) then
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 732,  successGPU)
  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  !if (useGPU) then
  !else ! useGPU
    do i=1,na
      if (my_pcol==pcol(i, nblk, np_cols)) then
        ! column i is on local processor
        l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
        l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
        a(l_row1:l_rows,l_col1) = 0
      endif
    enddo
  !endif ! useGPU

  if (useGPU) then
    ! copy back
    !successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
    !                   matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    !call check_memcpy_GPU_f("elpa_cholesky: memcpy a-> d_dev", 752,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_cholesky: a_dev", 755,  successGPU)
  endif

  ! restore original OpenMP settings
  call obj%timer%stop("elpa_cholesky_&
  &complex&
  &_&
  &double&
  &"//gpuString)

    end function elpa_cholesky_complex_double_impl


















!> \brief  elpa_cholesky_complex_single_impl: Cholesky factorization of a single-precision complex hermitian matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
    function elpa_cholesky_complex_single_impl(obj, a) result(success)

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
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
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




!cannot use "../src/cholesky/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)

  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu
  use mod_check_for_gpu
  use invert_trm_cuda, only : copy_float_complex_tmp1_tmp2, &
                              copy_float_complex_a_tmp1
  use cholesky_cuda
  implicit none
!    Copyright 2011, A. Marek
!
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
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
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
!

  integer, parameter :: rk = C_FLOAT
  integer, parameter :: ck = C_FLOAT_COMPLEX
  integer, parameter :: rck = C_FLOAT_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)              :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  complex(kind=rck)      :: a(obj%local_nrows,*)
  integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)              :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)              :: n, nc, i, info
  integer(kind=BLAS_KIND)       :: infoBLAS
  integer(kind=ik)              :: lcs, lce, lrs, lre
  integer(kind=ik)              :: tile_size, l_rows_tile, l_cols_tile

  complex(kind=rck), allocatable    :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
  logical                       :: wantDebug
  logical                       :: success
  integer(kind=ik)              :: istat, debug, error
  character(200)                :: errorMessage
  integer(kind=ik)              :: nrThreads, limitThreads
  character(20)                 :: gpuString
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_intptr_t)      :: tmp1_dev, tmatc_dev, tmatr_dev, a_dev, tmp2_dev
  type(c_ptr)                      :: tmp1_mpi_dev
  complex(kind=rck), pointer :: tmp1_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)      :: a_off, tmatc_off, tmatr_off
  type(c_ptr)                      :: tmatc_mpi_dev
  complex(kind=rck), pointer :: tmatc_mpi_fortran_ptr(:,:)
  integer(kind=c_int)              :: gpu_cholesky

  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &single&
                                                            &_&
                                                            &complex

  gpu_cholesky = 0
  ! GPU settings
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"ELPA_CHOLESKY: Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif

    call obj%get("gpu_cholesky",gpu_cholesky, error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for gpu_cholesky. Aborting..."
      stop
    endif

  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_CHOLESKY: Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
  endif

  if (gpu_cholesky .eq. 1) then
    useGPU = (gpu == 1)
  else
    useGPU = .false.
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("elpa_cholesky_&
  &complex&
  &_&
  &single&
  &"//gpuString)

  nrThreads=1

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("mpi_comm_parent", mpi_comm_all, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Error getting option for mpi_comm_all. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_rows",mpi_comm_rows,error )
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for mpi_comm_rows. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_cols",mpi_comm_cols,error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for mpi_comm_cols. Aborting..."
    stop
  endif

  call obj%get("debug",debug,error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_CHOLESKY: Problem getting option for debug settings. Aborting..."
    stop
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .false.
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), myidMPI, mpierr)


  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  myid    = int(myidMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")
  success = .true.

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"ELPA_CHOLESKY: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU


  if (useGPU) then
    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp1_dev", 257,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp1_dev", 260,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmp2_dev", 263,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmp2_dev", 266,  successGPU)

    successGPU = gpu_malloc(tmatc_dev, l_cols*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatc_dev", 269,  successGPU)

    successGPU = gpu_memset(tmatc_dev, 0, l_cols*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatc_dev", 272,  successGPU)

    successGPU = gpu_malloc(tmatr_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: tmatr_dev", 275,  successGPU)

    successGPU = gpu_memset(tmatr_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_cholesky: memset tmatr_dev", 278,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_cholesky: a_dev", 281,  successGPU)

  endif
  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp1", 285,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmp2", 288,  istat,  errorMessage)

  tmp1 = 0
  tmp2 = 0

  allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatr", 294,  istat,  errorMessage)

  allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_cholesky: tmatc", 297,  istat,  errorMessage)

  tmatr = 0
  tmatc = 0

  if (useGPU) then
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_cholesky 1: memcpy a-> a_dev", 305,  successGPU)
  endif

  do n = 1, na, nblk
    ! Calculate first local row and column of the still remaining matrix
    ! on the local processor

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

    if (n+nblk > na) then

      ! This is the last step, just do a Cholesky-Factorization
      ! of the remaining block

      if (useGPU) then
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 340,  successGPU)

          call CPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 347,  successGPU)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false.
            return
          endif ! info
        endif ! (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols))
      else ! useGPU
        if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
          call obj%timer%start("blas")

          call CPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                             int(matrixRows,kind=BLAS_KIND), infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false.
            return
          endif

        endif
      endif ! useGPU
      exit ! Loop
    endif ! (n+nblk > na) 

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then

        if (useGPU) then
          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 412,  successGPU)

          call CPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t), &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("elpa_cholesky: memcpy a_dev-> a", 419,  successGPU)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf: ",info
            success = .false.
            return
          endif ! info
        else ! useGPU
          ! The process owning the upper left remaining block does the
          ! Cholesky-Factorization of this block
          call obj%timer%start("blas")

          call CPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                               int(matrixRows,kind=BLAS_KIND) , infoBLAS )
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")

          if (info/=0) then
            if (wantDebug) write(error_unit,*) "elpa_cholesky_&
            &complex&

            &: Error in zpotrf 2: ",info

            success = .false.
            return
          endif ! useGPU
        endif
 
        if (useGPU) then
          call copy_float_complex_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nblk)
        else ! useGPU
          nc = 0
          do i=1,nblk
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! (my_pcol==pcol(n, nblk, np_cols))

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), tmp1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_cholesky: tmp1_dev to tmp1", 479,  successGPU)

      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      call obj%timer%start("mpi_communication")

      call MPI_Bcast(tmp1, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
                    MPI_COMPLEX,      &
                    int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")
!#else
!      tmp1_mpi_dev = transfer(tmp1_dev, tmp1_mpi_dev)
!      ! and associate a fortran pointer
!      call c_f_pointer(tmp1_mpi_dev, tmp1_mpi_fortran_ptr, [nblk,nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("cholesky: device_synchronize", 503,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!
!      call MPI_Bcast(tmp1_mpi_fortran_ptr, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
!#if REALCASE == 1
!                    MPI_REAL4,         &
!#endif
!#if 1 == 1
!                    MPI_COMPLEX,      &
!#endif
!                    int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!
!      call obj%timer%stop("mpi_cuda_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmp1_dev, int(loc(tmp1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_cholesky: tmp1 to tmp1_dev", 524,  successGPU)

      endif
!#endif


      if (useGPU) then
        call copy_float_complex_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nblk)
      else ! useGPU
        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1-1 + (l_colx-1)*matrixRows) * size_of_datatype
          call gpublas_CTRSM('L', 'U', 'C', 'N', nblk, l_cols-l_colx+1, ONE, &
                            tmp2_dev, nblk, a_dev+a_off, matrixRows)
        endif
        call obj%timer%stop("gpublas")

      else ! useGPU

        call obj%timer%start("blas")
        if (l_cols-l_colx+1>0) &
        call CTRSM('L', 'U', 'C', 'N', int(nblk,kind=BLAS_KIND),  &
                            int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
                            int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
        call obj%timer%stop("blas")
      endif ! useGPU
    endif ! (my_prow==prow(n, nblk, np_rows))


    if (useGPU) then
      if (my_prow==prow(n, nblk, np_rows)) then
        ! if l_cols-l_colx+1 == 0 kernel launch with 0 blocks => raises error
        if (l_cols-l_colx+1>0) &
           call copy_float_complex_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)
      endif
    else ! useGPU
      do i=1,nblk
        if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = conjg(a(l_row1+i-1,l_colx:l_cols))
      enddo
    endif ! useGPU

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      if (l_cols-l_colx+1 > 0) then
        num = l_cols*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmatc),kind=c_intptr_t), tmatc_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_cholesky: tmatc_dev to tmatc", 586,  successGPU)
      endif
    endif
!#endif


!#ifndef WITH_CUDA_AWARE_MPI
    do i=1,nblk
      call obj%timer%start("mpi_communication")
      if (l_cols-l_colx+1>0) &
      call MPI_Bcast(tmatc(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), MPI_COMPLEX, &
                     int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")
    enddo
!#else
!    tmatc_mpi_dev = transfer(tmatc_dev, tmatc_mpi_dev)
!    ! and associate a fortran pointer
!    call c_f_pointer(tmatc_mpi_dev, tmatc_mpi_fortran_ptr, [l_cols,nblk])
!    
!    if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!    successGPU = gpu_devicesynchronize()
!    call check_memcpy_GPU_f("cholesky: device_synchronize", 610,  successGPU)
!    if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!
!    do i=1,nblk
!      call obj%timer%start("mpi_cuda_communication")
!      if (l_cols-l_colx+1>0) &
!      call MPI_Bcast(tmatc_mpi_fortran_ptr(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), &
!                     MPI_COMPLEX, &
!                     int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)
!
!      call obj%timer%stop("mpi_cuda_communication")
!    enddo
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      !if (l_cols-l_colx+1 > 0) then
        num = l_cols*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmatc_dev, int(loc(tmatc),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_cholesky: tmatc to tmatc_dev", 632,  successGPU)
      !endif
    endif
!#endif

    if (useGPU) then
      ! can optimize memcpy here with previous ones

      ! a gpu version of elpa_transpose_vectors is needed

!#if !defined(1) || (defined(1) && defined(WITH_CUDA_AWARE_MPI))

      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(int(loc(tmatr),kind=c_intptr_t), tmatr_dev, num, &
                              gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("elpa_cholesky: tmatr_dev to tmatr", 657,  successGPU)
    endif

    call elpa_transpose_vectors_&
    &complex&
    &_&
    &single &
    (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
    tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
    n, na, nblk, nblk, nrThreads, .false.)

    if (useGPU) then
      num = l_rows*nblk*size_of_datatype
      successGPU = gpu_memcpy(tmatr_dev, int(loc(tmatr),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("elpa_cholesky: tmat to tmatr_dev", 672,  successGPU)
    endif


    if (useGPU) then
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce  < lcs .or. lre < lrs) cycle
        call obj%timer%start("gpublas")
        tmatr_off = (lrs-1 + (1-1)*l_rows) * size_of_datatype
        tmatc_off = (lcs-1 + (1-1)*l_cols) * size_of_datatype
        a_off = (lrs-1 + (lcs-1)*matrixRows) * size_of_datatype
        call gpublas_CGEMM('N', 'C', lre-lrs+1, lce-lcs+1, nblk, &
                            -ONE, tmatr_dev+tmatr_off, l_rows, tmatc_dev+tmatc_off, l_cols, ONE, &
                            a_dev+a_off, matrixRows)
        call obj%timer%stop("gpublas")
      enddo
    else !useGPU
      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle
        call obj%timer%start("blas")
        call CGEMM('N', 'C', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                            int(nblk,kind=BLAS_KIND), -ONE,  &
                            tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
                            int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
                            ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      enddo
    endif ! useGPU

  enddo ! n = 1, na, nblk

  if (useGPU) then
    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 713,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmp1_dev", 716,  successGPU)

    successGPU = gpu_free(tmatc_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatc_dev", 719,  successGPU)

    successGPU = gpu_free(tmatr_dev)
    call check_dealloc_GPU_f("elpa_cholesky: tmatr_dev", 722,  successGPU)
  endif

  deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 726,  istat,  errorMessage)


  if (useGPU) then
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                     matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_cholesky: memcpy 2 a-> a_dev", 732,  successGPU)
  endif
  ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

  !if (useGPU) then
  !else ! useGPU
    do i=1,na
      if (my_pcol==pcol(i, nblk, np_cols)) then
        ! column i is on local processor
        l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
        l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
        a(l_row1:l_rows,l_col1) = 0
      endif
    enddo
  !endif ! useGPU

  if (useGPU) then
    ! copy back
    !successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
    !                   matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    !call check_memcpy_GPU_f("elpa_cholesky: memcpy a-> d_dev", 752,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_cholesky: a_dev", 755,  successGPU)
  endif

  ! restore original OpenMP settings
  call obj%timer%stop("elpa_cholesky_&
  &complex&
  &_&
  &single&
  &"//gpuString)

    end function elpa_cholesky_complex_single_impl



end module
