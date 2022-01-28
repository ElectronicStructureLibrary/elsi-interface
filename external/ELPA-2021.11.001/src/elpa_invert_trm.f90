










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
module elpa_invert_trm
  use, intrinsic :: iso_c_binding
  use precision
  implicit none

  public
  public :: elpa_invert_trm_real_double_impl    !< Invert double-precision real triangular matrix
  public :: elpa_invert_trm_complex_double_impl  !< Invert double-precision complex triangular matrix

  public :: elpa_invert_trm_real_single_impl     !< Invert single-precision real triangular matrix

  public :: elpa_invert_trm_complex_single_impl  !< Invert single-precision complex triangular matrix
  contains

















!> \brief  elpa_invert_trm_real_double: Inverts a double-precision real upper triangular matrix
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
    function elpa_invert_trm_real_double_impl(obj, a) result(success)
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
! ELPA1 -- Faster replacements for ScaLAPACK symmetric eigenvalue routines
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".




!cannot use "../src/invert_trm/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)



















  use precision
  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use elpa_abstract_impl
  use elpa_gpu
  use mod_check_for_gpu
  use elpa_blas_interfaces
  use invert_trm_cuda

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
  integer(kind=ik)             :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  integer(kind=ik)             :: mpi_comm_all
  real(kind=rck)      :: a(obj%local_nrows,*)
  integer :: ii, jj

  integer(kind=ik)             :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)       :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)             :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)             :: n, nc, i, info, ns, nb
  integer(kind=BLAS_KIND)      :: infoBLAS
  real(kind=rck), allocatable   :: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)
  logical                      :: wantDebug
  logical                      :: success
  integer(kind=ik)             :: istat, debug, error
  character(200)               :: errorMessage
  character(20)                 :: gpuString
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=c_intptr_t)      :: tmat1_dev, tmat2_dev, a_dev, tmp1_dev, tmp2_dev, zero_dev
  type(c_ptr)                      :: tmp1_mpi_dev
  real(kind=rck), pointer :: tmp1_mpi_fortran_ptr(:)
  type(c_ptr)                      :: tmat1_mpi_dev, tmat2_mpi_dev
  real(kind=rck), pointer :: tmat1_mpi_fortran_ptr(:,:), tmat2_mpi_fortran_ptr(:,:)

  type(c_ptr)                   :: tmp2_mpi_dev, a_mpi_dev
  integer(kind=c_intptr_t)      :: a_off, tmat2_off, tmp1_off, tmp2_off
   real(kind=rck), pointer :: a_mpi_deviceptr(:,:), initializer_ptr(:) !DEB
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_int)           :: gpu_invert_trm
  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &double&
                                                            &_&
                                                            &real


  ! GPU settings
  gpu_invert_trm = 0
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"ELPA_INVERT_TRM: Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif
    call obj%get("gpu_invert_trm",gpu_invert_trm,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for gpu_cholesky. Aborting..."
      stop
    endif

  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
  endif

  if (gpu_invert_trm .eq. 1) then
    useGPU = (gpu == 1)
  else
    useGPU = .false.
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("elpa_invert_trm_&
  &real&
  &_&
  &double&
  &"//gpuString)

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("mpi_comm_parent", mpi_comm_all, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_all. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_rows", mpi_comm_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_rows. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_cols", mpi_comm_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_cols. Aborting..."
    stop
  endif

  call obj%get("debug", debug, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for debug. Aborting..."
    stop
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .true.
  endif
  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), myidMPI, mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)
  myid    = int(myidMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")


  success = .true.

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
     call obj%timer%start("check_for_gpu")
     call obj%set("use_gpu_id", myid, error)
    if (check_for_gpu(obj, myid, numGPU, .TRUE.)) then
       ! set the neccessary parameters       
      call set_gpu_parameters()
    else
      print *,"ELPA_INVERT_TRM: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU

  if (useGPU) then
    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmp1_dev", 237,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmp1_dev", 240,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmp2_dev", 243,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmp2_dev", 246,  successGPU)

    successGPU = gpu_malloc(tmat1_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmat1_dev", 249,  successGPU)

    successGPU = gpu_memset(tmat1_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmat1_dev", 252,  successGPU)

    successGPU = gpu_malloc(tmat2_dev, nblk*l_cols*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmat1_dev", 255,  successGPU)

    successGPU = gpu_memset(tmat2_dev, 0, nblk*l_cols*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmat2_dev", 258,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: a_dev", 261,  successGPU)

    successGPU = gpu_malloc(zero_dev, 1*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: zero_dev", 264,  successGPU)

    successGPU = gpu_memset(zero_dev, 0, 1*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset zero_dev", 267,  successGPU)
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmp1", 272,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmp2", 275,  istat,  errorMessage)

  tmp1 = 0
  tmp2 = 0

  allocate(tmat1(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat1", 281,  istat,  errorMessage)

  allocate(tmat2(nblk,l_cols), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat2", 284,  istat,  errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat2", 285,  istat,  errorMessage)

  tmat1 = 0
  tmat2 = 0

  if (useGPU) then
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t),  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_invert_trm: memcpy a-> d_dev", 293,  successGPU)
  endif


  ns = ((na-1)/nblk)*nblk + 1

  do n = ns,1,-nblk

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    nb = nblk
    if (na-n+1 < nblk) nb = na-n+1

    l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then
        if (useGPU) then

         
          ! still have to use cpu blas -> a generic GPU implementation would be needed

          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev, &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("invert_trm: memcpy a_dev -> a", 334,  successGPU)

          call DTRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                             infoBLAS)
          info = int(infoBLAS,kind=ik)

          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t),  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("invert_trm: memcpy a -> a_dev", 342,  successGPU)
          call obj%timer%stop("blas")

        else ! useGPU
          call obj%timer%start("blas")

          call DTRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                             infoBLAS)
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")
        endif ! useGPU
        if (info/=0) then
          if (wantDebug) write(error_unit,*) "elpa_invert_trm_&
            &real&

          &: Error in DTRTRI"

          success = .false.
          call obj%timer%stop("elpa_invert_trm_&
          &real&
          &_&
          &double&
          &"//gpuString)
          return
        endif

        if (useGPU) then
          call copy_double_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb)
        else ! useGPU
          nc = 0
          do i=1,nb
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! my_pcol==pcol(n, nblk, np_cols)

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), tmp1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmp1_dev to tmp1", 391,  successGPU)

      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      call obj%timer%start("mpi_communication")
      call MPI_Bcast(tmp1, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_REAL8,       &
                     int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
      call obj%timer%stop("mpi_communication")
!#else
!      tmp1_mpi_dev = transfer(tmp1_dev, tmp1_mpi_dev) 
!      ! and associate a fortran pointer
!      call c_f_pointer(tmp1_mpi_dev, tmp1_mpi_fortran_ptr, [nblk*nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 407,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!
!
!
!      if (wantDebug) call obj%timer%start("cuda_mpi_communication")
!      call MPI_Bcast(tmp1_mpi_fortran_ptr, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_REAL8,       &
!                     int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!      if (wantDebug) call obj%timer%stop("cuda_mpi_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if ((useGPU)) then  
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmp1_dev, int(loc(tmp1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmp1 to tmp1_dev", 423,  successGPU)

      endif
!#endif
      
      if (useGPU) then
        call copy_double_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nb)
      else ! useGPU
        nc = 0
        do i=1,nb
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1 -1 + (l_colx-1)*matrixRows) * size_of_datatype

          call gpublas_DTRMM('L', 'U', 'N', 'N', nb, l_cols-l_colx+1, ONE, tmp2_dev, &
                                      nblk, a_dev+a_off, matrixRows)
        
        endif
        call obj%timer%stop("gpublas")

        if (l_colx <= l_cols) then
          call copy_double_a_tmat2 (a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, & 
                                       l_row1, nb)
        endif

        if (my_pcol==pcol(n, nblk, np_cols)) then
           ! tmp2 has the lower left triangle 0
          call copy_double_tmp2_tmat2 (tmp2_dev, tmat2_dev, nblk, l_col1, nb) 
        endif
      else ! useGPU
        call obj%timer%start("blas")
        if (l_cols-l_colx+1>0) &
        call DTRMM('L', 'U', 'N', 'N', int(nb,kind=BLAS_KIND), int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, &
                              tmp2, int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
        if (l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
        if (my_pcol==pcol(n, nblk, np_cols)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0
      endif ! useGPU

    endif ! (my_prow==prow(n, nblk, np_rows)

    if (l_row1>1) then
      if (my_pcol==pcol(n, nblk, np_cols)) then
        if (useGPU) then
          call copy_double_a_tmat1 (a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1, zero_dev)
        else
          tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
          a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
        endif
      endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = l_rows*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmat1),kind=c_intptr_t), tmat1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat1_dev to tmat1", 487,  successGPU)
      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      do i=1,nb
        call obj%timer%start("mpi_communication")
        call MPI_Bcast(tmat1(1,i), int(l_row1-1,kind=MPI_KIND), MPI_REAL8, &
                       int(pcol(n, nblk, np_cols),kind=MPI_KIND), & 
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)

        call obj%timer%stop("mpi_communication")
      enddo
!#else
!      tmat1_mpi_dev = transfer(tmat1_dev, tmat1_mpi_dev)
!      ! and associate a fortran pointer
!      call c_f_pointer(tmat1_mpi_dev, tmat1_mpi_fortran_ptr, [l_rows,nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 508,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!      do i=1,nb
!        call MPI_Bcast(tmat1_mpi_fortran_ptr(1,i), int(l_row1-1,kind=MPI_KIND), MPI_REAL8, &
!                       int(pcol(n, nblk, np_cols),kind=MPI_KIND), & 
!                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!
!      enddo
!      call obj%timer%stop("mpi_cuda_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        ! cuda aware MPI here
        num = l_rows*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmat1_dev, int(loc(tmat1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat1 to tmat1_dev", 528,  successGPU)

      endif
!#endif
    endif ! (l_row1>1)

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      
      if (l_cols-l_col1+1 > 0) then
        num = nblk*l_cols*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmat2),kind=c_intptr_t), tmat2_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat2_dev to tmat2", 543,  successGPU)
      endif
    endif

    call obj%timer%start("mpi_communication")
    if (l_cols-l_col1+1 > 0) &
    call MPI_Bcast(tmat2(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), MPI_REAL8, &
                   int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

    call obj%timer%stop("mpi_communication")

    if (useGPU) then
      if (l_cols-l_col1+1 > 0) then
        num = nblk*l_cols*size_of_datatype
        successGPU = gpu_memcpy(tmat2_dev, int(loc(tmat2),kind=c_intptr_t), num, &
                                gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat2 to tmat2_dev", 559,  successGPU)
      endif
    endif
!#else
!      tmat2_mpi_dev = transfer(tmat2_dev, tmat2_mpi_dev)     
!      call c_f_pointer(tmat2_mpi_dev, tmat2_mpi_fortran_ptr, [nblk,l_cols])
!      
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 568,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!    if (l_cols-l_col1+1 > 0) &
!        call MPI_Bcast(tmat2_mpi_fortran_ptr(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), & 
!                       MPI_REAL8, int(prow(n, nblk, np_rows),kind=MPI_KIND), & 
!                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
!      call obj%timer%stop("mpi_cuda_communication")
!
!#endif

    if (useGPU) then
      call obj%timer%start("gpublas")

      tmat2_off = (1 - 1 + (l_col1-1) * nblk) * size_of_datatype      
      a_off = (1 - 1 + (l_col1-1) * matrixRows) * size_of_datatype
      if (l_row1>1 .and. l_cols-l_col1+1>0) &
        call gpublas_DGEMM('N', 'N', l_row1-1, l_cols-l_col1+1, nb, -ONE, &
                                    tmat1_dev, l_rows, tmat2_dev + tmat2_off, &
                                    nblk, ONE, a_dev+a_off, matrixRows)

      call obj%timer%stop("gpublas")

    else ! useGPU
      call obj%timer%start("blas")
      if (l_row1>1 .and. l_cols-l_col1+1>0) &
        call DGEMM('N', 'N', int(l_row1-1,kind=BLAS_KIND), int(l_cols-l_col1+1,kind=BLAS_KIND), &
                            int(nb,kind=BLAS_KIND), -ONE, &
                             tmat1, int(ubound(tmat1,dim=1),kind=BLAS_KIND), tmat2(1,l_col1), &
                             int(ubound(tmat2,dim=1),kind=BLAS_KIND), ONE, &
                              a(1,l_col1), int(matrixRows,kind=BLAS_KIND) )

      call obj%timer%stop("blas")
    endif ! useGPU
  enddo

  if (useGPU) then
  ! copy results back
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_invert_trm: memcpy a-> d_dev", 609,  successGPU)

  endif

  if (useGPU) then
    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmp1_dev", 615,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmp2_dev", 618,  successGPU)

    successGPU = gpu_free(tmat1_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmat1_dev", 621,  successGPU)

    successGPU = gpu_free(tmat2_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmat2_dev", 624,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: a_dev", 627,  successGPU)

    successGPU = gpu_free(zero_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: zero_dev", 630,  successGPU)

    !successGPU = gpu_host_unregister(int(loc(b),kind=c_intptr_t))
    !call check_host_unregister_GPU_f("elpa_multiply_a_b: b", 633,  successGPU)
  endif ! useGPU

  deallocate(tmp1, tmp2, tmat1, tmat2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_invert_trm: tmp1, tmp2, tmat1, tmat2", 637,  istat,  errorMessage)

  call obj%timer%stop("elpa_invert_trm_&
  &real&
  &_&
  &double&
  &"//gpuString)
     end function elpa_invert_trm_real_double_impl




















!> \brief  elpa_invert_trm_real_single_impl: Inverts a single-precision real upper triangular matrix
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

    function elpa_invert_trm_real_single_impl(obj, a) result(success)
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
! ELPA1 -- Faster replacements for ScaLAPACK symmetric eigenvalue routines
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".




!cannot use "../src/invert_trm/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)





















  use precision
  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use elpa_abstract_impl
  use elpa_gpu
  use mod_check_for_gpu
  use elpa_blas_interfaces
  use invert_trm_cuda

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
  integer(kind=ik)             :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  integer(kind=ik)             :: mpi_comm_all
  real(kind=rck)      :: a(obj%local_nrows,*)
  integer :: ii, jj

  integer(kind=ik)             :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)       :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)             :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)             :: n, nc, i, info, ns, nb
  integer(kind=BLAS_KIND)      :: infoBLAS
  real(kind=rck), allocatable   :: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)
  logical                      :: wantDebug
  logical                      :: success
  integer(kind=ik)             :: istat, debug, error
  character(200)               :: errorMessage
  character(20)                 :: gpuString
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=c_intptr_t)      :: tmat1_dev, tmat2_dev, a_dev, tmp1_dev, tmp2_dev, zero_dev
  type(c_ptr)                      :: tmp1_mpi_dev
  real(kind=rck), pointer :: tmp1_mpi_fortran_ptr(:)
  type(c_ptr)                      :: tmat1_mpi_dev, tmat2_mpi_dev
  real(kind=rck), pointer :: tmat1_mpi_fortran_ptr(:,:), tmat2_mpi_fortran_ptr(:,:)

  type(c_ptr)                   :: tmp2_mpi_dev, a_mpi_dev
  integer(kind=c_intptr_t)      :: a_off, tmat2_off, tmp1_off, tmp2_off
   real(kind=rck), pointer :: a_mpi_deviceptr(:,:), initializer_ptr(:) !DEB
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_int)           :: gpu_invert_trm
  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &single&
                                                            &_&
                                                            &real


  ! GPU settings
  gpu_invert_trm = 0
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"ELPA_INVERT_TRM: Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif
    call obj%get("gpu_invert_trm",gpu_invert_trm,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for gpu_cholesky. Aborting..."
      stop
    endif

  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
  endif

  if (gpu_invert_trm .eq. 1) then
    useGPU = (gpu == 1)
  else
    useGPU = .false.
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("elpa_invert_trm_&
  &real&
  &_&
  &single&
  &"//gpuString)

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("mpi_comm_parent", mpi_comm_all, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_all. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_rows", mpi_comm_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_rows. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_cols", mpi_comm_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_cols. Aborting..."
    stop
  endif

  call obj%get("debug", debug, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for debug. Aborting..."
    stop
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .true.
  endif
  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), myidMPI, mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)
  myid    = int(myidMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")


  success = .true.

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
     call obj%timer%start("check_for_gpu")
     call obj%set("use_gpu_id", myid, error)
    if (check_for_gpu(obj, myid, numGPU, .TRUE.)) then
       ! set the neccessary parameters       
      call set_gpu_parameters()
    else
      print *,"ELPA_INVERT_TRM: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU

  if (useGPU) then
    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmp1_dev", 237,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmp1_dev", 240,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmp2_dev", 243,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmp2_dev", 246,  successGPU)

    successGPU = gpu_malloc(tmat1_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmat1_dev", 249,  successGPU)

    successGPU = gpu_memset(tmat1_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmat1_dev", 252,  successGPU)

    successGPU = gpu_malloc(tmat2_dev, nblk*l_cols*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmat1_dev", 255,  successGPU)

    successGPU = gpu_memset(tmat2_dev, 0, nblk*l_cols*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmat2_dev", 258,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: a_dev", 261,  successGPU)

    successGPU = gpu_malloc(zero_dev, 1*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: zero_dev", 264,  successGPU)

    successGPU = gpu_memset(zero_dev, 0, 1*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset zero_dev", 267,  successGPU)
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmp1", 272,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmp2", 275,  istat,  errorMessage)

  tmp1 = 0
  tmp2 = 0

  allocate(tmat1(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat1", 281,  istat,  errorMessage)

  allocate(tmat2(nblk,l_cols), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat2", 284,  istat,  errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat2", 285,  istat,  errorMessage)

  tmat1 = 0
  tmat2 = 0

  if (useGPU) then
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t),  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_invert_trm: memcpy a-> d_dev", 293,  successGPU)
  endif


  ns = ((na-1)/nblk)*nblk + 1

  do n = ns,1,-nblk

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    nb = nblk
    if (na-n+1 < nblk) nb = na-n+1

    l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then
        if (useGPU) then

         
          ! still have to use cpu blas -> a generic GPU implementation would be needed

          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev, &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("invert_trm: memcpy a_dev -> a", 334,  successGPU)

          call STRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                             infoBLAS)
          info = int(infoBLAS,kind=ik)

          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t),  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("invert_trm: memcpy a -> a_dev", 342,  successGPU)
          call obj%timer%stop("blas")

        else ! useGPU
          call obj%timer%start("blas")

          call STRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                             infoBLAS)
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")
        endif ! useGPU
        if (info/=0) then
          if (wantDebug) write(error_unit,*) "elpa_invert_trm_&
            &real&

          &: Error in DTRTRI"

          success = .false.
          call obj%timer%stop("elpa_invert_trm_&
          &real&
          &_&
          &single&
          &"//gpuString)
          return
        endif

        if (useGPU) then
          call copy_float_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb)
        else ! useGPU
          nc = 0
          do i=1,nb
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! my_pcol==pcol(n, nblk, np_cols)

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), tmp1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmp1_dev to tmp1", 391,  successGPU)

      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      call obj%timer%start("mpi_communication")
      call MPI_Bcast(tmp1, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_REAL4,       &
                     int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
      call obj%timer%stop("mpi_communication")
!#else
!      tmp1_mpi_dev = transfer(tmp1_dev, tmp1_mpi_dev) 
!      ! and associate a fortran pointer
!      call c_f_pointer(tmp1_mpi_dev, tmp1_mpi_fortran_ptr, [nblk*nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 407,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!
!
!
!      if (wantDebug) call obj%timer%start("cuda_mpi_communication")
!      call MPI_Bcast(tmp1_mpi_fortran_ptr, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_REAL4,       &
!                     int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!      if (wantDebug) call obj%timer%stop("cuda_mpi_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if ((useGPU)) then  
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmp1_dev, int(loc(tmp1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmp1 to tmp1_dev", 423,  successGPU)

      endif
!#endif
      
      if (useGPU) then
        call copy_float_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nb)
      else ! useGPU
        nc = 0
        do i=1,nb
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1 -1 + (l_colx-1)*matrixRows) * size_of_datatype

          call gpublas_STRMM('L', 'U', 'N', 'N', nb, l_cols-l_colx+1, ONE, tmp2_dev, &
                                      nblk, a_dev+a_off, matrixRows)
        
        endif
        call obj%timer%stop("gpublas")

        if (l_colx <= l_cols) then
          call copy_float_a_tmat2 (a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, & 
                                       l_row1, nb)
        endif

        if (my_pcol==pcol(n, nblk, np_cols)) then
           ! tmp2 has the lower left triangle 0
          call copy_float_tmp2_tmat2 (tmp2_dev, tmat2_dev, nblk, l_col1, nb) 
        endif
      else ! useGPU
        call obj%timer%start("blas")
        if (l_cols-l_colx+1>0) &
        call STRMM('L', 'U', 'N', 'N', int(nb,kind=BLAS_KIND), int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, &
                              tmp2, int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
        if (l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
        if (my_pcol==pcol(n, nblk, np_cols)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0
      endif ! useGPU

    endif ! (my_prow==prow(n, nblk, np_rows)

    if (l_row1>1) then
      if (my_pcol==pcol(n, nblk, np_cols)) then
        if (useGPU) then
          call copy_float_a_tmat1 (a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1, zero_dev)
        else
          tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
          a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
        endif
      endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = l_rows*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmat1),kind=c_intptr_t), tmat1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat1_dev to tmat1", 487,  successGPU)
      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      do i=1,nb
        call obj%timer%start("mpi_communication")
        call MPI_Bcast(tmat1(1,i), int(l_row1-1,kind=MPI_KIND), MPI_REAL4, &
                       int(pcol(n, nblk, np_cols),kind=MPI_KIND), & 
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)

        call obj%timer%stop("mpi_communication")
      enddo
!#else
!      tmat1_mpi_dev = transfer(tmat1_dev, tmat1_mpi_dev)
!      ! and associate a fortran pointer
!      call c_f_pointer(tmat1_mpi_dev, tmat1_mpi_fortran_ptr, [l_rows,nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 508,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!      do i=1,nb
!        call MPI_Bcast(tmat1_mpi_fortran_ptr(1,i), int(l_row1-1,kind=MPI_KIND), MPI_REAL4, &
!                       int(pcol(n, nblk, np_cols),kind=MPI_KIND), & 
!                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!
!      enddo
!      call obj%timer%stop("mpi_cuda_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        ! cuda aware MPI here
        num = l_rows*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmat1_dev, int(loc(tmat1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat1 to tmat1_dev", 528,  successGPU)

      endif
!#endif
    endif ! (l_row1>1)

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      
      if (l_cols-l_col1+1 > 0) then
        num = nblk*l_cols*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmat2),kind=c_intptr_t), tmat2_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat2_dev to tmat2", 543,  successGPU)
      endif
    endif

    call obj%timer%start("mpi_communication")
    if (l_cols-l_col1+1 > 0) &
    call MPI_Bcast(tmat2(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), MPI_REAL4, &
                   int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

    call obj%timer%stop("mpi_communication")

    if (useGPU) then
      if (l_cols-l_col1+1 > 0) then
        num = nblk*l_cols*size_of_datatype
        successGPU = gpu_memcpy(tmat2_dev, int(loc(tmat2),kind=c_intptr_t), num, &
                                gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat2 to tmat2_dev", 559,  successGPU)
      endif
    endif
!#else
!      tmat2_mpi_dev = transfer(tmat2_dev, tmat2_mpi_dev)     
!      call c_f_pointer(tmat2_mpi_dev, tmat2_mpi_fortran_ptr, [nblk,l_cols])
!      
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 568,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!    if (l_cols-l_col1+1 > 0) &
!        call MPI_Bcast(tmat2_mpi_fortran_ptr(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), & 
!                       MPI_REAL4, int(prow(n, nblk, np_rows),kind=MPI_KIND), & 
!                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
!      call obj%timer%stop("mpi_cuda_communication")
!
!#endif

    if (useGPU) then
      call obj%timer%start("gpublas")

      tmat2_off = (1 - 1 + (l_col1-1) * nblk) * size_of_datatype      
      a_off = (1 - 1 + (l_col1-1) * matrixRows) * size_of_datatype
      if (l_row1>1 .and. l_cols-l_col1+1>0) &
        call gpublas_SGEMM('N', 'N', l_row1-1, l_cols-l_col1+1, nb, -ONE, &
                                    tmat1_dev, l_rows, tmat2_dev + tmat2_off, &
                                    nblk, ONE, a_dev+a_off, matrixRows)

      call obj%timer%stop("gpublas")

    else ! useGPU
      call obj%timer%start("blas")
      if (l_row1>1 .and. l_cols-l_col1+1>0) &
        call SGEMM('N', 'N', int(l_row1-1,kind=BLAS_KIND), int(l_cols-l_col1+1,kind=BLAS_KIND), &
                            int(nb,kind=BLAS_KIND), -ONE, &
                             tmat1, int(ubound(tmat1,dim=1),kind=BLAS_KIND), tmat2(1,l_col1), &
                             int(ubound(tmat2,dim=1),kind=BLAS_KIND), ONE, &
                              a(1,l_col1), int(matrixRows,kind=BLAS_KIND) )

      call obj%timer%stop("blas")
    endif ! useGPU
  enddo

  if (useGPU) then
  ! copy results back
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_invert_trm: memcpy a-> d_dev", 609,  successGPU)

  endif

  if (useGPU) then
    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmp1_dev", 615,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmp2_dev", 618,  successGPU)

    successGPU = gpu_free(tmat1_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmat1_dev", 621,  successGPU)

    successGPU = gpu_free(tmat2_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmat2_dev", 624,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: a_dev", 627,  successGPU)

    successGPU = gpu_free(zero_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: zero_dev", 630,  successGPU)

    !successGPU = gpu_host_unregister(int(loc(b),kind=c_intptr_t))
    !call check_host_unregister_GPU_f("elpa_multiply_a_b: b", 633,  successGPU)
  endif ! useGPU

  deallocate(tmp1, tmp2, tmat1, tmat2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_invert_trm: tmp1, tmp2, tmat1, tmat2", 637,  istat,  errorMessage)

  call obj%timer%stop("elpa_invert_trm_&
  &real&
  &_&
  &single&
  &"//gpuString)
    end function elpa_invert_trm_real_single_impl




















!> \brief  elpa_invert_trm_complex_double_impl: Inverts a double-precision complex upper triangular matrix
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
     function elpa_invert_trm_complex_double_impl(obj, a) result(success)
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
! ELPA1 -- Faster replacements for ScaLAPACK symmetric eigenvalue routines
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".




!cannot use "../src/invert_trm/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)





















  use precision
  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use elpa_abstract_impl
  use elpa_gpu
  use mod_check_for_gpu
  use elpa_blas_interfaces
  use invert_trm_cuda

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
  integer(kind=ik)             :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  integer(kind=ik)             :: mpi_comm_all
  complex(kind=rck)      :: a(obj%local_nrows,*)
  integer :: ii, jj

  integer(kind=ik)             :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)       :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)             :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)             :: n, nc, i, info, ns, nb
  integer(kind=BLAS_KIND)      :: infoBLAS
  complex(kind=rck), allocatable   :: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)
  logical                      :: wantDebug
  logical                      :: success
  integer(kind=ik)             :: istat, debug, error
  character(200)               :: errorMessage
  character(20)                 :: gpuString
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=c_intptr_t)      :: tmat1_dev, tmat2_dev, a_dev, tmp1_dev, tmp2_dev, zero_dev
  type(c_ptr)                      :: tmp1_mpi_dev
  complex(kind=rck), pointer :: tmp1_mpi_fortran_ptr(:)
  type(c_ptr)                      :: tmat1_mpi_dev, tmat2_mpi_dev
  complex(kind=rck), pointer :: tmat1_mpi_fortran_ptr(:,:), tmat2_mpi_fortran_ptr(:,:)

  type(c_ptr)                   :: tmp2_mpi_dev, a_mpi_dev
  integer(kind=c_intptr_t)      :: a_off, tmat2_off, tmp1_off, tmp2_off
   complex(kind=rck), pointer :: a_mpi_deviceptr(:,:), initializer_ptr(:) !DEB
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_int)           :: gpu_invert_trm
  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &double&
                                                            &_&
                                                            &complex


  ! GPU settings
  gpu_invert_trm = 0
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"ELPA_INVERT_TRM: Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif
    call obj%get("gpu_invert_trm",gpu_invert_trm,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for gpu_cholesky. Aborting..."
      stop
    endif

  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
  endif

  if (gpu_invert_trm .eq. 1) then
    useGPU = (gpu == 1)
  else
    useGPU = .false.
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("elpa_invert_trm_&
  &complex&
  &_&
  &double&
  &"//gpuString)

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("mpi_comm_parent", mpi_comm_all, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_all. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_rows", mpi_comm_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_rows. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_cols", mpi_comm_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_cols. Aborting..."
    stop
  endif

  call obj%get("debug", debug, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for debug. Aborting..."
    stop
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .true.
  endif
  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), myidMPI, mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)
  myid    = int(myidMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")


  success = .true.

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
     call obj%timer%start("check_for_gpu")
     call obj%set("use_gpu_id", myid, error)
    if (check_for_gpu(obj, myid, numGPU, .TRUE.)) then
       ! set the neccessary parameters       
      call set_gpu_parameters()
    else
      print *,"ELPA_INVERT_TRM: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU

  if (useGPU) then
    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmp1_dev", 237,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmp1_dev", 240,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmp2_dev", 243,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmp2_dev", 246,  successGPU)

    successGPU = gpu_malloc(tmat1_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmat1_dev", 249,  successGPU)

    successGPU = gpu_memset(tmat1_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmat1_dev", 252,  successGPU)

    successGPU = gpu_malloc(tmat2_dev, nblk*l_cols*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmat1_dev", 255,  successGPU)

    successGPU = gpu_memset(tmat2_dev, 0, nblk*l_cols*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmat2_dev", 258,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: a_dev", 261,  successGPU)

    successGPU = gpu_malloc(zero_dev, 1*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: zero_dev", 264,  successGPU)

    successGPU = gpu_memset(zero_dev, 0, 1*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset zero_dev", 267,  successGPU)
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmp1", 272,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmp2", 275,  istat,  errorMessage)

  tmp1 = 0
  tmp2 = 0

  allocate(tmat1(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat1", 281,  istat,  errorMessage)

  allocate(tmat2(nblk,l_cols), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat2", 284,  istat,  errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat2", 285,  istat,  errorMessage)

  tmat1 = 0
  tmat2 = 0

  if (useGPU) then
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t),  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_invert_trm: memcpy a-> d_dev", 293,  successGPU)
  endif


  ns = ((na-1)/nblk)*nblk + 1

  do n = ns,1,-nblk

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    nb = nblk
    if (na-n+1 < nblk) nb = na-n+1

    l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then
        if (useGPU) then

         
          ! still have to use cpu blas -> a generic GPU implementation would be needed

          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev, &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("invert_trm: memcpy a_dev -> a", 334,  successGPU)

          call ZTRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                             infoBLAS)
          info = int(infoBLAS,kind=ik)

          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t),  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("invert_trm: memcpy a -> a_dev", 342,  successGPU)
          call obj%timer%stop("blas")

        else ! useGPU
          call obj%timer%start("blas")

          call ZTRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                             infoBLAS)
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")
        endif ! useGPU
        if (info/=0) then
          if (wantDebug) write(error_unit,*) "elpa_invert_trm_&
            &complex&

          &: Error in ZTRTRI"

          success = .false.
          call obj%timer%stop("elpa_invert_trm_&
          &complex&
          &_&
          &double&
          &"//gpuString)
          return
        endif

        if (useGPU) then
          call copy_double_complex_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb)
        else ! useGPU
          nc = 0
          do i=1,nb
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! my_pcol==pcol(n, nblk, np_cols)

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), tmp1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmp1_dev to tmp1", 391,  successGPU)

      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      call obj%timer%start("mpi_communication")
      call MPI_Bcast(tmp1, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,       &
                     int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
      call obj%timer%stop("mpi_communication")
!#else
!      tmp1_mpi_dev = transfer(tmp1_dev, tmp1_mpi_dev) 
!      ! and associate a fortran pointer
!      call c_f_pointer(tmp1_mpi_dev, tmp1_mpi_fortran_ptr, [nblk*nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 407,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!
!
!
!      if (wantDebug) call obj%timer%start("cuda_mpi_communication")
!      call MPI_Bcast(tmp1_mpi_fortran_ptr, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,       &
!                     int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!      if (wantDebug) call obj%timer%stop("cuda_mpi_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if ((useGPU)) then  
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmp1_dev, int(loc(tmp1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmp1 to tmp1_dev", 423,  successGPU)

      endif
!#endif
      
      if (useGPU) then
        call copy_double_complex_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nb)
      else ! useGPU
        nc = 0
        do i=1,nb
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1 -1 + (l_colx-1)*matrixRows) * size_of_datatype

          call gpublas_ZTRMM('L', 'U', 'N', 'N', nb, l_cols-l_colx+1, ONE, tmp2_dev, &
                                      nblk, a_dev+a_off, matrixRows)
        
        endif
        call obj%timer%stop("gpublas")

        if (l_colx <= l_cols) then
          call copy_double_complex_a_tmat2 (a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, & 
                                       l_row1, nb)
        endif

        if (my_pcol==pcol(n, nblk, np_cols)) then
           ! tmp2 has the lower left triangle 0
          call copy_double_complex_tmp2_tmat2 (tmp2_dev, tmat2_dev, nblk, l_col1, nb) 
        endif
      else ! useGPU
        call obj%timer%start("blas")
        if (l_cols-l_colx+1>0) &
        call ZTRMM('L', 'U', 'N', 'N', int(nb,kind=BLAS_KIND), int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, &
                              tmp2, int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
        if (l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
        if (my_pcol==pcol(n, nblk, np_cols)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0
      endif ! useGPU

    endif ! (my_prow==prow(n, nblk, np_rows)

    if (l_row1>1) then
      if (my_pcol==pcol(n, nblk, np_cols)) then
        if (useGPU) then
          call copy_double_complex_a_tmat1 (a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1, zero_dev)
        else
          tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
          a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
        endif
      endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = l_rows*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmat1),kind=c_intptr_t), tmat1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat1_dev to tmat1", 487,  successGPU)
      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      do i=1,nb
        call obj%timer%start("mpi_communication")
        call MPI_Bcast(tmat1(1,i), int(l_row1-1,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                       int(pcol(n, nblk, np_cols),kind=MPI_KIND), & 
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)

        call obj%timer%stop("mpi_communication")
      enddo
!#else
!      tmat1_mpi_dev = transfer(tmat1_dev, tmat1_mpi_dev)
!      ! and associate a fortran pointer
!      call c_f_pointer(tmat1_mpi_dev, tmat1_mpi_fortran_ptr, [l_rows,nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 508,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!      do i=1,nb
!        call MPI_Bcast(tmat1_mpi_fortran_ptr(1,i), int(l_row1-1,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
!                       int(pcol(n, nblk, np_cols),kind=MPI_KIND), & 
!                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!
!      enddo
!      call obj%timer%stop("mpi_cuda_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        ! cuda aware MPI here
        num = l_rows*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmat1_dev, int(loc(tmat1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat1 to tmat1_dev", 528,  successGPU)

      endif
!#endif
    endif ! (l_row1>1)

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      
      if (l_cols-l_col1+1 > 0) then
        num = nblk*l_cols*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmat2),kind=c_intptr_t), tmat2_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat2_dev to tmat2", 543,  successGPU)
      endif
    endif

    call obj%timer%start("mpi_communication")
    if (l_cols-l_col1+1 > 0) &
    call MPI_Bcast(tmat2(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                   int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

    call obj%timer%stop("mpi_communication")

    if (useGPU) then
      if (l_cols-l_col1+1 > 0) then
        num = nblk*l_cols*size_of_datatype
        successGPU = gpu_memcpy(tmat2_dev, int(loc(tmat2),kind=c_intptr_t), num, &
                                gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat2 to tmat2_dev", 559,  successGPU)
      endif
    endif
!#else
!      tmat2_mpi_dev = transfer(tmat2_dev, tmat2_mpi_dev)     
!      call c_f_pointer(tmat2_mpi_dev, tmat2_mpi_fortran_ptr, [nblk,l_cols])
!      
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 568,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!    if (l_cols-l_col1+1 > 0) &
!        call MPI_Bcast(tmat2_mpi_fortran_ptr(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), & 
!                       MPI_DOUBLE_COMPLEX, int(prow(n, nblk, np_rows),kind=MPI_KIND), & 
!                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
!      call obj%timer%stop("mpi_cuda_communication")
!
!#endif

    if (useGPU) then
      call obj%timer%start("gpublas")

      tmat2_off = (1 - 1 + (l_col1-1) * nblk) * size_of_datatype      
      a_off = (1 - 1 + (l_col1-1) * matrixRows) * size_of_datatype
      if (l_row1>1 .and. l_cols-l_col1+1>0) &
        call gpublas_ZGEMM('N', 'N', l_row1-1, l_cols-l_col1+1, nb, -ONE, &
                                    tmat1_dev, l_rows, tmat2_dev + tmat2_off, &
                                    nblk, ONE, a_dev+a_off, matrixRows)

      call obj%timer%stop("gpublas")

    else ! useGPU
      call obj%timer%start("blas")
      if (l_row1>1 .and. l_cols-l_col1+1>0) &
        call ZGEMM('N', 'N', int(l_row1-1,kind=BLAS_KIND), int(l_cols-l_col1+1,kind=BLAS_KIND), &
                            int(nb,kind=BLAS_KIND), -ONE, &
                             tmat1, int(ubound(tmat1,dim=1),kind=BLAS_KIND), tmat2(1,l_col1), &
                             int(ubound(tmat2,dim=1),kind=BLAS_KIND), ONE, &
                              a(1,l_col1), int(matrixRows,kind=BLAS_KIND) )

      call obj%timer%stop("blas")
    endif ! useGPU
  enddo

  if (useGPU) then
  ! copy results back
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_invert_trm: memcpy a-> d_dev", 609,  successGPU)

  endif

  if (useGPU) then
    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmp1_dev", 615,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmp2_dev", 618,  successGPU)

    successGPU = gpu_free(tmat1_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmat1_dev", 621,  successGPU)

    successGPU = gpu_free(tmat2_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmat2_dev", 624,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: a_dev", 627,  successGPU)

    successGPU = gpu_free(zero_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: zero_dev", 630,  successGPU)

    !successGPU = gpu_host_unregister(int(loc(b),kind=c_intptr_t))
    !call check_host_unregister_GPU_f("elpa_multiply_a_b: b", 633,  successGPU)
  endif ! useGPU

  deallocate(tmp1, tmp2, tmat1, tmat2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_invert_trm: tmp1, tmp2, tmat1, tmat2", 637,  istat,  errorMessage)

  call obj%timer%stop("elpa_invert_trm_&
  &complex&
  &_&
  &double&
  &"//gpuString)
    end function elpa_invert_trm_complex_double_impl



















!> \brief  elpa_invert_trm_complex_single_impl: Inverts a single-precision complex upper triangular matrix
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
    function elpa_invert_trm_complex_single_impl(obj, a) result(success)
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
! ELPA1 -- Faster replacements for ScaLAPACK symmetric eigenvalue routines
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".




!cannot use "../src/invert_trm/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




















  use precision
  use elpa1_compute
  use elpa_utilities
  use elpa_mpi
  use elpa_abstract_impl
  use elpa_gpu
  use mod_check_for_gpu
  use elpa_blas_interfaces
  use invert_trm_cuda

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
  integer(kind=ik)             :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  integer(kind=ik)             :: mpi_comm_all
  complex(kind=rck)      :: a(obj%local_nrows,*)
  integer :: ii, jj

  integer(kind=ik)             :: my_prow, my_pcol, np_rows, np_cols, myid
  integer(kind=MPI_KIND)       :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, myidMPI
  integer(kind=ik)             :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
  integer(kind=ik)             :: n, nc, i, info, ns, nb
  integer(kind=BLAS_KIND)      :: infoBLAS
  complex(kind=rck), allocatable   :: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)
  logical                      :: wantDebug
  logical                      :: success
  integer(kind=ik)             :: istat, debug, error
  character(200)               :: errorMessage
  character(20)                 :: gpuString
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=c_intptr_t)      :: tmat1_dev, tmat2_dev, a_dev, tmp1_dev, tmp2_dev, zero_dev
  type(c_ptr)                      :: tmp1_mpi_dev
  complex(kind=rck), pointer :: tmp1_mpi_fortran_ptr(:)
  type(c_ptr)                      :: tmat1_mpi_dev, tmat2_mpi_dev
  complex(kind=rck), pointer :: tmat1_mpi_fortran_ptr(:,:), tmat2_mpi_fortran_ptr(:,:)

  type(c_ptr)                   :: tmp2_mpi_dev, a_mpi_dev
  integer(kind=c_intptr_t)      :: a_off, tmat2_off, tmp1_off, tmp2_off
   complex(kind=rck), pointer :: a_mpi_deviceptr(:,:), initializer_ptr(:) !DEB
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_int)           :: gpu_invert_trm
  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &single&
                                                            &_&
                                                            &complex


  ! GPU settings
  gpu_invert_trm = 0
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"ELPA_INVERT_TRM: Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif
    call obj%get("gpu_invert_trm",gpu_invert_trm,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for gpu_cholesky. Aborting..."
      stop
    endif

  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"ELPA_INVERT_TRM: Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
  endif

  if (gpu_invert_trm .eq. 1) then
    useGPU = (gpu == 1)
  else
    useGPU = .false.
  endif

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("elpa_invert_trm_&
  &complex&
  &_&
  &single&
  &"//gpuString)

  na         = obj%na
  matrixRows = obj%local_nrows
  nblk       = obj%nblk
  matrixCols = obj%local_ncols

  call obj%get("mpi_comm_parent", mpi_comm_all, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_all. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_rows", mpi_comm_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_rows. Aborting..."
    stop
  endif
  call obj%get("mpi_comm_cols", mpi_comm_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for mpi_comm_cols. Aborting..."
    stop
  endif

  call obj%get("debug", debug, error)
  if (error .ne. ELPA_OK) then
    print *,"ELPA_INVERT_TRM: Error getting option for debug. Aborting..."
    stop
  endif
  if (debug == 1) then
    wantDebug = .true.
  else
    wantDebug = .true.
  endif
  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), myidMPI, mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)
  myid    = int(myidMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")


  success = .true.

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

  if (useGPU) then
     call obj%timer%start("check_for_gpu")
     call obj%set("use_gpu_id", myid, error)
    if (check_for_gpu(obj, myid, numGPU, .TRUE.)) then
       ! set the neccessary parameters       
      call set_gpu_parameters()
    else
      print *,"ELPA_INVERT_TRM: GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")
  else ! useGPU
  endif ! useGPU

  if (useGPU) then
    successGPU = gpu_malloc(tmp1_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmp1_dev", 237,  successGPU)

    successGPU = gpu_memset(tmp1_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmp1_dev", 240,  successGPU)

    successGPU = gpu_malloc(tmp2_dev, nblk*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmp2_dev", 243,  successGPU)

    successGPU = gpu_memset(tmp2_dev, 0, nblk*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmp2_dev", 246,  successGPU)

    successGPU = gpu_malloc(tmat1_dev, l_rows*nblk*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmat1_dev", 249,  successGPU)

    successGPU = gpu_memset(tmat1_dev, 0, l_rows*nblk*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmat1_dev", 252,  successGPU)

    successGPU = gpu_malloc(tmat2_dev, nblk*l_cols*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: tmat1_dev", 255,  successGPU)

    successGPU = gpu_memset(tmat2_dev, 0, nblk*l_cols*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset tmat2_dev", 258,  successGPU)

    successGPU = gpu_malloc(a_dev, matrixRows*matrixCols*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: a_dev", 261,  successGPU)

    successGPU = gpu_malloc(zero_dev, 1*size_of_datatype)
    call check_alloc_GPU_f("elpa_invert_trm: zero_dev", 264,  successGPU)

    successGPU = gpu_memset(zero_dev, 0, 1*size_of_datatype)
    call check_memcpy_GPU_f("elpa_invert_trm: memset zero_dev", 267,  successGPU)
  endif ! useGPU


  allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmp1", 272,  istat,  errorMessage)

  allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmp2", 275,  istat,  errorMessage)

  tmp1 = 0
  tmp2 = 0

  allocate(tmat1(l_rows,nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat1", 281,  istat,  errorMessage)

  allocate(tmat2(nblk,l_cols), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat2", 284,  istat,  errorMessage)
  call check_allocate_f("elpa_invert_trm: tmat2", 285,  istat,  errorMessage)

  tmat1 = 0
  tmat2 = 0

  if (useGPU) then
    successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t),  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_invert_trm: memcpy a-> d_dev", 293,  successGPU)
  endif


  ns = ((na-1)/nblk)*nblk + 1

  do n = ns,1,-nblk

    l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
    l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

    nb = nblk
    if (na-n+1 < nblk) nb = na-n+1

    l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
    l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)

    if (my_prow==prow(n, nblk, np_rows)) then

      if (my_pcol==pcol(n, nblk, np_cols)) then
        if (useGPU) then

         
          ! still have to use cpu blas -> a generic GPU implementation would be needed

          call obj%timer%start("blas")
          successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev, &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("invert_trm: memcpy a_dev -> a", 334,  successGPU)

          call CTRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                             infoBLAS)
          info = int(infoBLAS,kind=ik)

          successGPU = gpu_memcpy(a_dev, int(loc(a(1,1)),kind=c_intptr_t),  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("invert_trm: memcpy a -> a_dev", 342,  successGPU)
          call obj%timer%stop("blas")

        else ! useGPU
          call obj%timer%start("blas")

          call CTRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                             infoBLAS)
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")
        endif ! useGPU
        if (info/=0) then
          if (wantDebug) write(error_unit,*) "elpa_invert_trm_&
            &complex&

          &: Error in ZTRTRI"

          success = .false.
          call obj%timer%stop("elpa_invert_trm_&
          &complex&
          &_&
          &single&
          &"//gpuString)
          return
        endif

        if (useGPU) then
          call copy_float_complex_a_tmp1 (a_dev, tmp1_dev, l_row1, l_col1, matrixRows, nb)
        else ! useGPU
          nc = 0
          do i=1,nb
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif ! useGPU
      endif ! my_pcol==pcol(n, nblk, np_cols)

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), tmp1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmp1_dev to tmp1", 391,  successGPU)

      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      call obj%timer%start("mpi_communication")
      call MPI_Bcast(tmp1, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_COMPLEX,       &
                     int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
      call obj%timer%stop("mpi_communication")
!#else
!      tmp1_mpi_dev = transfer(tmp1_dev, tmp1_mpi_dev) 
!      ! and associate a fortran pointer
!      call c_f_pointer(tmp1_mpi_dev, tmp1_mpi_fortran_ptr, [nblk*nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 407,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!
!
!
!      if (wantDebug) call obj%timer%start("cuda_mpi_communication")
!      call MPI_Bcast(tmp1_mpi_fortran_ptr, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_COMPLEX,       &
!                     int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!      if (wantDebug) call obj%timer%stop("cuda_mpi_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if ((useGPU)) then  
        num = nblk*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmp1_dev, int(loc(tmp1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmp1 to tmp1_dev", 423,  successGPU)

      endif
!#endif
      
      if (useGPU) then
        call copy_float_complex_tmp1_tmp2 (tmp1_dev, tmp2_dev, nblk, nb)
      else ! useGPU
        nc = 0
        do i=1,nb
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo
      endif ! useGPU

      if (useGPU) then
        call obj%timer%start("gpublas")
        if (l_cols-l_colx+1 > 0) then
          a_off = (l_row1 -1 + (l_colx-1)*matrixRows) * size_of_datatype

          call gpublas_CTRMM('L', 'U', 'N', 'N', nb, l_cols-l_colx+1, ONE, tmp2_dev, &
                                      nblk, a_dev+a_off, matrixRows)
        
        endif
        call obj%timer%stop("gpublas")

        if (l_colx <= l_cols) then
          call copy_float_complex_a_tmat2 (a_dev, tmat2_dev, nblk, matrixRows, l_cols, l_colx, & 
                                       l_row1, nb)
        endif

        if (my_pcol==pcol(n, nblk, np_cols)) then
           ! tmp2 has the lower left triangle 0
          call copy_float_complex_tmp2_tmat2 (tmp2_dev, tmat2_dev, nblk, l_col1, nb) 
        endif
      else ! useGPU
        call obj%timer%start("blas")
        if (l_cols-l_colx+1>0) &
        call CTRMM('L', 'U', 'N', 'N', int(nb,kind=BLAS_KIND), int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, &
                              tmp2, int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND))
        call obj%timer%stop("blas")
        if (l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
        if (my_pcol==pcol(n, nblk, np_cols)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0
      endif ! useGPU

    endif ! (my_prow==prow(n, nblk, np_rows)

    if (l_row1>1) then
      if (my_pcol==pcol(n, nblk, np_cols)) then
        if (useGPU) then
          call copy_float_complex_a_tmat1 (a_dev, tmat1_dev, l_rows, matrixRows, nb, l_row1, l_col1, zero_dev)
        else
          tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
          a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
        endif
      endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        num = l_rows*nblk*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmat1),kind=c_intptr_t), tmat1_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat1_dev to tmat1", 487,  successGPU)
      endif
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      do i=1,nb
        call obj%timer%start("mpi_communication")
        call MPI_Bcast(tmat1(1,i), int(l_row1-1,kind=MPI_KIND), MPI_COMPLEX, &
                       int(pcol(n, nblk, np_cols),kind=MPI_KIND), & 
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)

        call obj%timer%stop("mpi_communication")
      enddo
!#else
!      tmat1_mpi_dev = transfer(tmat1_dev, tmat1_mpi_dev)
!      ! and associate a fortran pointer
!      call c_f_pointer(tmat1_mpi_dev, tmat1_mpi_fortran_ptr, [l_rows,nblk])
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 508,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!      do i=1,nb
!        call MPI_Bcast(tmat1_mpi_fortran_ptr(1,i), int(l_row1-1,kind=MPI_KIND), MPI_COMPLEX, &
!                       int(pcol(n, nblk, np_cols),kind=MPI_KIND), & 
!                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
!
!      enddo
!      call obj%timer%stop("mpi_cuda_communication")
!#endif

!#ifndef WITH_CUDA_AWARE_MPI
      if (useGPU) then
        ! cuda aware MPI here
        num = l_rows*nblk*size_of_datatype
        successGPU = gpu_memcpy(tmat1_dev, int(loc(tmat1),kind=c_intptr_t), num, &
                              gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat1 to tmat1_dev", 528,  successGPU)

      endif
!#endif
    endif ! (l_row1>1)

!#ifndef WITH_CUDA_AWARE_MPI
    if (useGPU) then
      
      if (l_cols-l_col1+1 > 0) then
        num = nblk*l_cols*size_of_datatype
        successGPU = gpu_memcpy(int(loc(tmat2),kind=c_intptr_t), tmat2_dev, num, &
                              gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat2_dev to tmat2", 543,  successGPU)
      endif
    endif

    call obj%timer%start("mpi_communication")
    if (l_cols-l_col1+1 > 0) &
    call MPI_Bcast(tmat2(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), MPI_COMPLEX, &
                   int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

    call obj%timer%stop("mpi_communication")

    if (useGPU) then
      if (l_cols-l_col1+1 > 0) then
        num = nblk*l_cols*size_of_datatype
        successGPU = gpu_memcpy(tmat2_dev, int(loc(tmat2),kind=c_intptr_t), num, &
                                gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("elpa_invert_trm: tmat2 to tmat2_dev", 559,  successGPU)
      endif
    endif
!#else
!      tmat2_mpi_dev = transfer(tmat2_dev, tmat2_mpi_dev)     
!      call c_f_pointer(tmat2_mpi_dev, tmat2_mpi_fortran_ptr, [nblk,l_cols])
!      
!      if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
!      successGPU = gpu_devicesynchronize()
!      call check_memcpy_GPU_f("invert_trm: device_synchronize", 568,  successGPU)
!      if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
!      call obj%timer%start("mpi_cuda_communication")
!    if (l_cols-l_col1+1 > 0) &
!        call MPI_Bcast(tmat2_mpi_fortran_ptr(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), & 
!                       MPI_COMPLEX, int(prow(n, nblk, np_rows),kind=MPI_KIND), & 
!                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
!      call obj%timer%stop("mpi_cuda_communication")
!
!#endif

    if (useGPU) then
      call obj%timer%start("gpublas")

      tmat2_off = (1 - 1 + (l_col1-1) * nblk) * size_of_datatype      
      a_off = (1 - 1 + (l_col1-1) * matrixRows) * size_of_datatype
      if (l_row1>1 .and. l_cols-l_col1+1>0) &
        call gpublas_CGEMM('N', 'N', l_row1-1, l_cols-l_col1+1, nb, -ONE, &
                                    tmat1_dev, l_rows, tmat2_dev + tmat2_off, &
                                    nblk, ONE, a_dev+a_off, matrixRows)

      call obj%timer%stop("gpublas")

    else ! useGPU
      call obj%timer%start("blas")
      if (l_row1>1 .and. l_cols-l_col1+1>0) &
        call CGEMM('N', 'N', int(l_row1-1,kind=BLAS_KIND), int(l_cols-l_col1+1,kind=BLAS_KIND), &
                            int(nb,kind=BLAS_KIND), -ONE, &
                             tmat1, int(ubound(tmat1,dim=1),kind=BLAS_KIND), tmat2(1,l_col1), &
                             int(ubound(tmat2,dim=1),kind=BLAS_KIND), ONE, &
                              a(1,l_col1), int(matrixRows,kind=BLAS_KIND) )

      call obj%timer%stop("blas")
    endif ! useGPU
  enddo

  if (useGPU) then
  ! copy results back
    successGPU = gpu_memcpy(int(loc(a(1,1)),kind=c_intptr_t), a_dev,  &
                       matrixRows*matrixCols* size_of_datatype, gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("elpa_invert_trm: memcpy a-> d_dev", 609,  successGPU)

  endif

  if (useGPU) then
    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmp1_dev", 615,  successGPU)

    successGPU = gpu_free(tmp2_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmp2_dev", 618,  successGPU)

    successGPU = gpu_free(tmat1_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmat1_dev", 621,  successGPU)

    successGPU = gpu_free(tmat2_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: tmat2_dev", 624,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: a_dev", 627,  successGPU)

    successGPU = gpu_free(zero_dev)
    call check_dealloc_GPU_f("elpa_invert_trm: zero_dev", 630,  successGPU)

    !successGPU = gpu_host_unregister(int(loc(b),kind=c_intptr_t))
    !call check_host_unregister_GPU_f("elpa_multiply_a_b: b", 633,  successGPU)
  endif ! useGPU

  deallocate(tmp1, tmp2, tmat1, tmat2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_invert_trm: tmp1, tmp2, tmat1, tmat2", 637,  istat,  errorMessage)

  call obj%timer%stop("elpa_invert_trm_&
  &complex&
  &_&
  &single&
  &"//gpuString)
    end function elpa_invert_trm_complex_single_impl

end module
